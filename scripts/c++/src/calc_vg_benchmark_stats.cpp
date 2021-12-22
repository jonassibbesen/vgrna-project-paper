
/*
calc_vg_benchmark_stats
Calculates benchmarking statistics for vg simulated reference mapped reads.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <assert.h>

#include "SeqLib/RefGenome.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/GenomicRegionCollection.h"

#include "utils.hpp"

using namespace SeqLib;


unordered_map<string, pair<string, BamRecord> > parseTranscriptAlignments(const string & transcript_bam_file) {

    unordered_map<string, pair<string, BamRecord> > transcript_alignments;

    BamReader bam_reader;
    bam_reader.Open(transcript_bam_file);
    assert(bam_reader.IsOpen());

    BamRecord bam_record;

    while (bam_reader.GetNextRecord(bam_record)) { 

        assert(bam_record.GetCigar().NumQueryConsumed() == bam_record.Length());
        assert(transcript_alignments.emplace(bam_record.Qname(), make_pair(bam_record.ChrName(bam_reader.Header()), bam_record)).second);
    }

    bam_reader.Close();

    return transcript_alignments;
}

unordered_map<string, pair<string, uint32_t> > parseReadsTranscriptInfo(const string & read_transcript_file) {

    unordered_map<string, pair<string, uint32_t> > read_transcript_info;

    ifstream read_istream(read_transcript_file);
    assert(read_istream.is_open());

    string line;

    while (read_istream.good()) {

        getline(read_istream, line);

        if (line.empty()) {

            continue;
        }

        auto line_split = splitString(line, '\t');
        assert(line_split.size() == 4);

        if (line_split.at(2) == "offset") {

            continue;
        }

        assert(read_transcript_info.emplace(line_split.front(), make_pair(line_split.at(1), stoi(line_split.at(2)))).second);
    }

    read_istream.close();

    return read_transcript_info;
}

void writeEmptyEvaluation(const BamRecord & bam_record, const BamReader & bam_reader, const bool debug_output, bool count_variants, unordered_map<string, uint32_t> * benchmark_stats) {
    
    stringstream benchmark_stats_ss;
    benchmark_stats_ss << '0';                              // TruthAlignmentLength
    benchmark_stats_ss << '\t' << bam_record.MappedFlag();  // IsMapped
    benchmark_stats_ss << '\t' << bam_record.MapQuality();  // MapQ
    benchmark_stats_ss << '\t' << getAllelicMapQ(bam_record);             // AllelicMapQ
    benchmark_stats_ss << '\t' << bam_record.Length();      // Length
    benchmark_stats_ss << '\t' << '0';                      // SoftClipLength
    benchmark_stats_ss << '\t' << '0';                      // Overlap
    
    if (count_variants) {
        benchmark_stats_ss << "\t0\t0\t0\t0";               // SubstitutionBP / IndelBP
    }

    if (debug_output) {

        cout << bam_record.Qname();
        cout << '\t' << bam_record.ChrName(bam_reader.Header()) << ':';
        cout << '\t' << ':';
        cout << '\t' << benchmark_stats_ss.str();
        cout << endl;

    } else {         

        auto benchmark_stats_it = benchmark_stats->emplace(benchmark_stats_ss.str(), 0);
        benchmark_stats_it.first->second++;  
    }
}

int main(int argc, char* argv[]) {

    if (!(argc == 7 || argc == 8)) {

        cerr << "Usage: calc_vg_benchmark_stats <read_bam> <transcript_bam> <read_transcript_file> <min_base_quality> <vcf1,vcf2,...> <sample> (<enable_debug_output>) > statistics.txt" << endl;
        return 1;
    }
    
    printScriptHeader(argc, argv);

    BamReader bam_reader;
    bam_reader.Open(argv[1]);
    assert(bam_reader.IsOpen());
    
    auto transcript_alignments = parseTranscriptAlignments(argv[2]);
    cerr << "Number of transcript alignments: " << transcript_alignments.size() << endl;

    auto read_transcript_info = parseReadsTranscriptInfo(argv[3]);
    cerr << "Number of reads: " << read_transcript_info.size() << "\n" << endl;
    
    const uint32_t min_base_quality = stoi(argv[4]);

    // read in VCFs
    string sample_name = argv[6];
    vector<string> vcf_filenames = splitString(argv[5], ',');
    unordered_map<string, int> contig_to_vcf;
    vector<tuple<htsFile*, bcf_hdr_t*, tbx_t*, int>> vcfs = initializeVCFs(vcf_filenames, sample_name, contig_to_vcf);

    const bool debug_output = (argc == 8);

    stringstream base_header; 
    base_header << "TruthAlignmentLength";
    base_header << "\t" << "IsMapped";
    base_header << "\t" << "MapQ";
    base_header << "\t" << "AllelicMapQ";
    base_header << "\t" << "Length";
    base_header << "\t" << "SoftClipLength";
    base_header << "\t" << "Overlap";

    if (!vcf_filenames.empty()) {
        base_header << "\t" << "SubstitutionBP1";
        base_header << "\t"  << "IndelBP1";
        base_header << "\t" << "SubstitutionBP2";
        base_header << "\t"  << "IndelBP2";
    }
    
    if (debug_output) {

        cout << "Name" << "\t" << "Alignment" << "\t" << "TruthAlignment" << "\t" << base_header.str() << endl;
    } 

    BamRecord bam_record;

    unordered_map<string, uint32_t> benchmark_stats;

    uint32_t num_reads = 0;
    double sum_overlap = 0;

    while (bam_reader.GetNextRecord(bam_record)) { 

        if (bam_record.SecondaryFlag()) {

            continue;
        }
        
        num_reads++;
        
        int32_t trimmed_start = 0;
        int32_t trimmed_end = 0;

        bam_record.QualityTrimmedSequence(min_base_quality, trimmed_start, trimmed_end);
        assert(trimmed_end <= bam_record.Length());

        if (trimmed_end < 0) {

            writeEmptyEvaluation(bam_record, bam_reader, debug_output, !vcf_filenames.empty(), &benchmark_stats);
            sum_overlap += 1;

            continue;
        }

        assert(trimmed_end > trimmed_start);

        uint32_t trimmed_length = trimmed_end - trimmed_start;
        assert(trimmed_length <= bam_record.Length());

        string read_name = bam_record.Qname();

        auto read_transcript_info_it = read_transcript_info.find(read_name);
        
        if (read_transcript_info_it == read_transcript_info.end()) {

            if (bam_record.FirstFlag()) {

                read_name += "_1";
            
            } else {

                read_name += "_2";            
            }

            read_transcript_info_it = read_transcript_info.find(read_name);
        }

        assert(read_transcript_info_it != read_transcript_info.end());

        auto read_transcript_id = read_transcript_info_it->second.first;

        auto transcript_alignments_it = transcript_alignments.find(read_transcript_id);
        assert(transcript_alignments_it != transcript_alignments.end());

        auto read_transcript_pos = read_transcript_info_it->second.second;

        if (transcript_alignments_it->second.second.ReverseFlag()) {

            read_transcript_pos = transcript_alignments_it->second.second.Length() - read_transcript_pos - bam_record.Length();
        }

        auto transcript_read_cigar = trimCigar(transcript_alignments_it->second.second.GetCigar(), read_transcript_pos + trimmed_start, trimmed_length);

        if (transcript_read_cigar.first.NumReferenceConsumed() == 0) {

            writeEmptyEvaluation(bam_record, bam_reader, debug_output, !vcf_filenames.empty(), &benchmark_stats);
            sum_overlap += 1;

            continue;
        }

        uint32_t soft_clip_length = 0;
        double overlap = 0;

        string read_genomic_regions_str = "";
        auto transcript_cigar_genomic_regions = cigarToGenomicRegions(transcript_read_cigar.first, 0, transcript_alignments_it->second.second.Position() + transcript_read_cigar.second);

        if (bam_record.MappedFlag()) {

            assert(trimmed_start + trimmed_length <= bam_record.GetCigar().NumQueryConsumed());
            auto trimmed_cigar = trimCigar(bam_record.GetCigar(), trimmed_start, trimmed_length);     

            assert(trimmed_cigar.first.NumQueryConsumed() >= transcript_read_cigar.first.NumQueryConsumed());
            soft_clip_length = cigarTypeLength(trimmed_cigar.first, 'S');

            auto read_cigar_genomic_regions = cigarToGenomicRegions(trimmed_cigar.first, 0, bam_record.Position() + trimmed_cigar.second);

            if (debug_output) {

                read_genomic_regions_str = genomicRegionsToString(read_cigar_genomic_regions);
            }

            if (bam_record.ChrName(bam_reader.Header()) == transcript_alignments_it->second.first) {

                read_cigar_genomic_regions.CreateTreeMap();
                transcript_cigar_genomic_regions.CreateTreeMap();

                auto cigar_genomic_regions_intersection = transcript_cigar_genomic_regions.Intersection(read_cigar_genomic_regions, true);
                overlap = cigar_genomic_regions_intersection.TotalWidth() / static_cast<double>(transcript_cigar_genomic_regions.TotalWidth());
            }
        }

        stringstream benchmark_stats_ss;

        benchmark_stats_ss << transcript_cigar_genomic_regions.TotalWidth();
        benchmark_stats_ss << '\t' << bam_record.MappedFlag();
        benchmark_stats_ss << '\t' << bam_record.MapQuality();
        benchmark_stats_ss << '\t' << getAllelicMapQ(bam_record);
        benchmark_stats_ss << '\t' << trimmed_length;
        benchmark_stats_ss << '\t' << soft_clip_length;
        benchmark_stats_ss << '\t' << overlap;
        
        
        if (!vcf_filenames.empty()) {
            
            int32_t subs_bp_1 = 0;
            int32_t indel_bp_1 = 0;
            int32_t subs_bp_2 = 0;
            int32_t indel_bp_2 = 0;
            
            
            string contig = bam_record.ChrName(bam_reader.Header());

            // we tolerate the contig being because it can happen when chrY vcfs are left out
            // for XX samples
            if (contig_to_vcf.count(contig)) {
                htsFile* vcf;
                bcf_hdr_t* header;
                tbx_t* tabix_index;
                int samp_idx;
                tie(vcf, header, tabix_index, samp_idx) = vcfs.at(contig_to_vcf.at(contig));
                
                if (samp_idx < 0) {
                    cerr << "error: truth alignment for " << bam_record.Qname() << " is to contig " << contig << " in VCF file " << vcf_filenames[contig_to_vcf.at(contig)] << " that does not contain sample " << sample_name << endl;
                    return 1;
                }
                
                tie(subs_bp_1, indel_bp_1, subs_bp_2, indel_bp_2) = countIndelsAndSubs(transcript_alignments_it->second.first,
                                                                                       transcript_cigar_genomic_regions,
                                                                                       vcf, header, tabix_index, samp_idx);
            }
            benchmark_stats_ss << '\t' << subs_bp_1;
            benchmark_stats_ss << '\t' << indel_bp_1;
            benchmark_stats_ss << '\t' << subs_bp_2;
            benchmark_stats_ss << '\t' << indel_bp_2;
        }
        
        if (debug_output) {

            cout << bam_record.Qname();
            cout << '\t' << bam_record.ChrName(bam_reader.Header()) << ':' << read_genomic_regions_str;
            cout << '\t' << transcript_alignments_it->second.first << ':' << genomicRegionsToString(transcript_cigar_genomic_regions);
            cout << '\t' << benchmark_stats_ss.str();
            cout << '\n';

        } else {         

            auto benchmark_stats_it = benchmark_stats.emplace(benchmark_stats_ss.str(), 0);
            benchmark_stats_it.first->second++;  
        }

        sum_overlap += overlap;

//        if (bam_record.Qname() == "seed_7640106_fragment_662692_2") {
//
//            cerr << endl;
//            cerr << read_transcript_pos << endl;
//            cerr << bam_record.ReverseFlag() << endl;
//            cerr << transcript_alignments_it->second.second.ReverseFlag() << endl;
//            cerr << trimmed_start << endl;
//            cerr << trimmed_length << endl;
//            cerr << bam_record.GetCigar() << endl;
//            cerr << trimCigar(bam_record.GetCigar(), trimmed_start, trimmed_length).first << endl;
//            cerr << trimCigar(bam_record.GetCigar(), trimmed_start, trimmed_length).second << endl;
//            cerr << bam_record.Qname() << endl;
//            cerr << bam_record.Position() << endl;
//            cerr << "\t" << bam_record.ChrName(bam_reader.Header()) << ":" << read_genomic_regions_str;
//            cerr << "\t" << transcript_alignments_it->second.first << ":" << genomicRegionsToString(transcript_cigar_genomic_regions);
//            cerr << "\t" << benchmark_stats_ss.str();
//            cerr << endl;
//        }

        if (num_reads % 10000000 == 0) {

            cerr << "Number of analysed reads: " << num_reads << endl;
        }
    }

    bam_reader.Close();

    if (!debug_output) {

        cout << "Count" << "\t" << base_header.str() << endl;

        for (auto & stats: benchmark_stats) {

            cout << stats.second << "\t" << stats.first << endl;
        }
    }

    cerr << "\nTotal number of analysed reads: " << num_reads << endl;
    cerr << "Average overlap: " << sum_overlap/num_reads << endl;
    
    for (auto vcf_rec : vcfs) {
        tbx_destroy(get<2>(vcf_rec));
        bcf_hdr_destroy(get<1>(vcf_rec));
        hts_close(get<0>(vcf_rec));
    }

	return 0;
}
