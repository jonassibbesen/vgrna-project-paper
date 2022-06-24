
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

struct ReadTranscriptInfo {

    const string transcript_id;
    const uint32_t position;
    const bool is_reverse;

    ReadTranscriptInfo(const string & transcript_id_in, const uint32_t position_in, const bool is_reverse_in) : transcript_id(transcript_id_in), position(position_in), is_reverse(is_reverse_in) {}
};

unordered_map<string, ReadTranscriptInfo> parseReadsTranscriptInfo(const string & read_transcript_file) {

    unordered_map<string, ReadTranscriptInfo> read_transcript_info;

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

        assert(read_transcript_info.emplace(line_split.front(), ReadTranscriptInfo(line_split.at(1), stoi(line_split.at(2)), stoi(line_split.at(3)))).second);
    }

    read_istream.close();

    return read_transcript_info;
}

stringstream createEmptyStats(const BamRecord & bam_record, const BamReader & bam_reader, const bool debug_output) {
    
    stringstream benchmark_stats_ss;

    if (debug_output) {

        benchmark_stats_ss << bam_record.Qname();
        benchmark_stats_ss << '\t' << ':';
        benchmark_stats_ss << '\t' << bam_record.ChrName(bam_reader.Header()) << ':';
        benchmark_stats_ss << '\t';
    }

    benchmark_stats_ss << '0';                                  // TruthAlignmentLength
    benchmark_stats_ss << '\t' << bam_record.MappedFlag();      // IsMapped
    benchmark_stats_ss << '\t' << bam_record.MapQuality();      // MapQ
    benchmark_stats_ss << '\t' << getAllelicMapQ(bam_record);   // AllelicMapQ
    benchmark_stats_ss << '\t' << getGroupMapQ(bam_record);     // GroupMapQ
    benchmark_stats_ss << '\t' << bam_record.Length();          // Length
    benchmark_stats_ss << '\t' << '0';                          // SoftClipLength
    benchmark_stats_ss << '\t' << '0';                          // Overlap
    benchmark_stats_ss << '\t' << '0';                          // NonAnnoSJ
    benchmark_stats_ss << '\t' << '0';                          // ComplexFrac
    benchmark_stats_ss << "\t0\t0\t0\t0";                       // SubstitutionBP / IndelBP

    return benchmark_stats_ss;
}

void addStats(unordered_map<string, pair<uint32_t, double> > * benchmark_stats, const BamRecord & bam_record, const double overlap, const stringstream & benchmark_stats_ss, string * prev_read_name, uint32_t * prev_mapq, double * prev_overlap, string * prev_output_string, uint32_t * cur_primary_mapq, double * cur_primary_overlap) {

    string read_name = bam_record.Qname();

    if (isFirstRead(bam_record)) {

        read_name += "_1";
            
    } else {

        read_name += "_2";            
    }

    if (*prev_output_string == "") {

        *prev_read_name = read_name;
        *prev_mapq = bam_record.MapQuality();
        *prev_overlap = overlap;
        *prev_output_string = benchmark_stats_ss.str();

    } else if (read_name == *prev_read_name) {

        if (overlap > *prev_overlap || (doubleCompare(overlap, *prev_overlap) && bam_record.MapQuality() > *prev_mapq)) {

            *prev_read_name = read_name;
            *prev_mapq = bam_record.MapQuality();
            *prev_overlap = overlap;
            *prev_output_string = benchmark_stats_ss.str();
        } 

    } else {

        stringstream primary_ss;
        primary_ss << "\t" << *cur_primary_mapq;
        primary_ss << "\t" << *cur_primary_overlap;

        *prev_output_string += primary_ss.str();

        auto benchmark_stats_it = benchmark_stats->emplace(*prev_output_string, pair<uint32_t, double>(0, 0));

        benchmark_stats_it.first->second.first++;  
        benchmark_stats_it.first->second.second += *prev_overlap;  

        *prev_read_name = read_name;
        *prev_mapq = bam_record.MapQuality();
        *prev_overlap = overlap;
        *prev_output_string = benchmark_stats_ss.str();
    }  

    if (!bam_record.SecondaryFlag()) {

        *cur_primary_mapq = bam_record.MapQuality();
        *cur_primary_overlap = overlap;
    }
}

int main(int argc, char* argv[]) {

    if (!(argc == 9 || argc == 10)) {

        cerr << "Usage: calc_vg_benchmark_stats <read_bam> <transcript_bam> <read_transcript_file> <min_base_quality> <transcript_gff> <complex_bed> <vcf1,vcf2,...> <sample> (<enable_debug_output>) > statistics.txt" << endl;
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

    GRC splice_junctions;

    if (string(argv[5]) != "null") {

        splice_junctions = parseSpliceJunctions(argv[5], bam_reader.Header());
    }

    GRC complex_regions;

    if (string(argv[6]) != "null") {

        complex_regions = parseRegionsBed(argv[6], bam_reader.Header());
    }

    cerr << "Total number of splice junction: " << splice_junctions.size() << endl;
    cerr << "Total length of complex regions: " << complex_regions.TotalWidth() << "\n" << endl;

    vector<string> vcf_filenames;
    unordered_map<string, int> contig_to_vcf;
    vector<tuple<htsFile*, bcf_hdr_t*, tbx_t*, int> > vcfs;

    string sample_name = argv[8];

    if (string(argv[7]) != "null") {

        // read in VCFs
        vcf_filenames = splitString(argv[7], ',');
        vcfs = initializeVCFs(vcf_filenames, sample_name, contig_to_vcf);
    }

    const bool debug_output = (argc == 10);

    stringstream base_header; 
    base_header << "TruthAlignmentLength";
    base_header << "\t" << "IsMapped";
    base_header << "\t" << "MapQ";
    base_header << "\t" << "AllelicMapQ";
    base_header << "\t" << "GroupMapQ";
    base_header << "\t" << "Length";
    base_header << "\t" << "SoftClipLength";
    base_header << "\t" << "Overlap";
    base_header << "\t" << "NonAnnoSJ";
    base_header << "\t" << "ComplexFrac";
    base_header << "\t" << "SubstitutionBP1";
    base_header << "\t" << "IndelBP1";
    base_header << "\t" << "SubstitutionBP2";
    base_header << "\t" << "IndelBP2";
    base_header << "\t" << "PrimaryMapq";
    base_header << "\t" << "PrimaryOverlap";

    if (debug_output) {

        cout << "Name" << "\t" << "TruthAlignment" << "\t" << "Alignment" << "\t" << base_header.str() << endl;
    } 

    BamRecord bam_record;

    unordered_map<string, pair<uint32_t, double> > benchmark_stats;

    string prev_read_name = "";

    uint32_t prev_mapq = 0;    
    double prev_overlap = 0;

    string prev_output_string = "";

    uint32_t cur_primary_mapq = 0;
    double cur_primary_overlap = 0;

    uint32_t num_reads = 0;

    while (bam_reader.GetNextRecord(bam_record)) { 

        num_reads++;

        int32_t trimmed_start = 0;
        int32_t trimmed_end = 0;

        bam_record.QualityTrimmedSequence(min_base_quality, trimmed_start, trimmed_end);
        assert(trimmed_end <= bam_record.Length());

        if (trimmed_end < 0) {

            auto benchmark_stats_ss = createEmptyStats(bam_record, bam_reader, debug_output);

            if (debug_output) {

                cout << benchmark_stats_ss.str() << "\t" << bam_record.MapQuality() << "\t" << "0" << endl;
            
            } else {

                addStats(&benchmark_stats, bam_record, 0, benchmark_stats_ss, &prev_read_name, &prev_mapq, &prev_overlap, &prev_output_string, &cur_primary_mapq, &cur_primary_overlap);
            }

            continue;
        }

        assert(trimmed_end > trimmed_start);

        uint32_t trimmed_length = trimmed_end - trimmed_start;
        assert(trimmed_start + trimmed_length <= bam_record.Length());

        string read_name = bam_record.Qname();

        auto read_transcript_info_it = read_transcript_info.find(read_name);
        
        if (read_transcript_info_it == read_transcript_info.end()) {

            if (isFirstRead(bam_record)) {

                read_name += "_1";
            
            } else {

                read_name += "_2";            
            }

            read_transcript_info_it = read_transcript_info.find(read_name);
        }

        assert(read_transcript_info_it != read_transcript_info.end());

        auto read_transcript_id = read_transcript_info_it->second.transcript_id;

        auto transcript_alignments_it = transcript_alignments.find(read_transcript_id);
        assert(transcript_alignments_it != transcript_alignments.end());

        auto read_transcript_pos = read_transcript_info_it->second.position;

        auto read_transcript_trimmed_start = trimmed_start;
        auto read_transcript_trimmed_length = trimmed_length;

        if (transcript_alignments_it->second.second.ReverseFlag()) {

            read_transcript_pos = transcript_alignments_it->second.second.Length() - read_transcript_pos - bam_record.Length();

            if (bam_record.ReverseFlag() == read_transcript_info_it->second.is_reverse) {

                read_transcript_trimmed_start = bam_record.Length() - read_transcript_trimmed_start - read_transcript_trimmed_length;
                assert(read_transcript_trimmed_start >= 0);
            }

        } else {

            if (bam_record.ReverseFlag() != read_transcript_info_it->second.is_reverse) {

                read_transcript_trimmed_start = bam_record.Length() - read_transcript_trimmed_start - read_transcript_trimmed_length;
                assert(read_transcript_trimmed_start >= 0);
            }
        }

        assert(read_transcript_trimmed_start + read_transcript_trimmed_length <= bam_record.Length());

        auto transcript_read_cigar = trimCigar(transcript_alignments_it->second.second.GetCigar(), read_transcript_pos + read_transcript_trimmed_start, read_transcript_trimmed_length);

        if (transcript_read_cigar.first.NumReferenceConsumed() == 0) {

            auto benchmark_stats_ss = createEmptyStats(bam_record, bam_reader, debug_output);

            if (debug_output) {

                cout << benchmark_stats_ss.str() << "\t" << bam_record.MapQuality() << "\t" << "0" << endl;
            
            } else {

                addStats(&benchmark_stats, bam_record, 0, benchmark_stats_ss, &prev_read_name, &prev_mapq, &prev_overlap, &prev_output_string, &cur_primary_mapq, &cur_primary_overlap);
            }

            continue;
        }

        auto transcript_cigar_genomic_regions = cigarToGenomicRegions(transcript_read_cigar.first, bam_reader.Header().Name2ID(transcript_alignments_it->second.first), transcript_alignments_it->second.second.Position() + transcript_read_cigar.second);
        transcript_cigar_genomic_regions.CreateTreeMap();

        uint32_t soft_clip_length = 0;
        double overlap = 0;

        string read_genomic_regions_str = "";

        if (bam_record.MappedFlag()) {

            assert(trimmed_start + trimmed_length <= bam_record.GetCigar().NumQueryConsumed());
            auto trimmed_cigar = trimCigar(bam_record.GetCigar(), trimmed_start, trimmed_length);     

            assert(trimmed_cigar.first.NumQueryConsumed() >= transcript_read_cigar.first.NumQueryConsumed());
            soft_clip_length = cigarTypeLength(trimmed_cigar.first, 'S');

            auto read_cigar_genomic_regions = cigarToGenomicRegions(trimmed_cigar.first, bam_record.ChrID(), bam_record.Position() + trimmed_cigar.second);
            read_cigar_genomic_regions.CreateTreeMap();

            if (debug_output) {

                read_genomic_regions_str = genomicRegionsToString(read_cigar_genomic_regions);
            }

            auto intersection = transcript_cigar_genomic_regions.Intersection(read_cigar_genomic_regions, true);
            intersection.MergeOverlappingIntervals();

            overlap = intersection.TotalWidth() / static_cast<double>(transcript_cigar_genomic_regions.TotalWidth());
        }

        stringstream benchmark_stats_ss;

        benchmark_stats_ss << transcript_cigar_genomic_regions.TotalWidth();
        benchmark_stats_ss << '\t' << bam_record.MappedFlag();
        benchmark_stats_ss << '\t' << bam_record.MapQuality();
        benchmark_stats_ss << '\t' << getAllelicMapQ(bam_record);
        benchmark_stats_ss << '\t' << getGroupMapQ(bam_record);
        benchmark_stats_ss << '\t' << trimmed_length;
        benchmark_stats_ss << '\t' << soft_clip_length;
        benchmark_stats_ss << '\t' << overlap;

        uint32_t non_anno_splice_junctions = 0;
        uint32_t transcript_read_cigar_genomic_pos = transcript_alignments_it->second.second.Position() + transcript_read_cigar.second;

        for (auto & field: transcript_read_cigar.first) {

            assert(field.Type() != 'H');
            assert(field.Length() > 0);

            if (!field.ConsumesQuery()) { 

                assert(field.Type() == 'D' || field.Type() == 'N');
                auto del_region = GenomicRegion(bam_reader.Header().Name2ID(transcript_alignments_it->second.first), transcript_read_cigar_genomic_pos, transcript_read_cigar_genomic_pos + field.Length() - 1);

                auto sj_indexes = splice_junctions.FindOverlappedIntervals(del_region, true);
                uint32_t best_sj_dist = std::numeric_limits<uint32_t>::max();

                for (auto & sj_idx: sj_indexes) {

                    auto sj = splice_junctions.at(sj_idx);
                    best_sj_dist = min(best_sj_dist, static_cast<uint32_t>(max(sj.DistanceBetweenStarts(del_region), sj.DistanceBetweenEnds(del_region))));
                }

                if (field.Length() >= 20 && best_sj_dist > 5) {

                    non_anno_splice_junctions++;
                }
            }
            
            if (field.ConsumesReference()) { 

                transcript_read_cigar_genomic_pos += field.Length();
            }
        }

        benchmark_stats_ss << '\t' << non_anno_splice_junctions;

        auto complex_regions_intersection = complex_regions.Intersection(transcript_cigar_genomic_regions, true);
        complex_regions_intersection.MergeOverlappingIntervals();

        double complex_frac = complex_regions_intersection.TotalWidth() / static_cast<double>(transcript_cigar_genomic_regions.TotalWidth());;

        benchmark_stats_ss << '\t' << complex_frac;
                            
        int32_t subs_bp_1 = 0;
        int32_t indel_bp_1 = 0;
        int32_t subs_bp_2 = 0;
        int32_t indel_bp_2 = 0;

        if (!vcf_filenames.empty()) {
        
            // we tolerate the contig being because it can happen when chrY vcfs are left out
            // for XX samples
            if (contig_to_vcf.count(transcript_alignments_it->second.first)) {
                htsFile* vcf;
                bcf_hdr_t* header;
                tbx_t* tabix_index;
                int samp_idx;
                tie(vcf, header, tabix_index, samp_idx) = vcfs.at(contig_to_vcf.at(transcript_alignments_it->second.first));
                
                if (samp_idx < 0) {
                    cerr << "error: truth alignment for " << bam_record.Qname() << " is to contig " << transcript_alignments_it->second.first << " in VCF file " << vcf_filenames[contig_to_vcf.at(transcript_alignments_it->second.first)] << " that does not contain sample " << sample_name << endl;
                    return 1;
                }
                
                tie(subs_bp_1, indel_bp_1, subs_bp_2, indel_bp_2) = countIndelsAndSubs(transcript_alignments_it->second.first,
                                                                                       transcript_cigar_genomic_regions,
                                                                                       vcf, header, tabix_index, samp_idx);
            }
        }

        benchmark_stats_ss << '\t' << subs_bp_1;
        benchmark_stats_ss << '\t' << indel_bp_1;
        benchmark_stats_ss << '\t' << subs_bp_2;
        benchmark_stats_ss << '\t' << indel_bp_2;

        if (debug_output) {

            cout << bam_record.Qname();
            cout << '\t' << transcript_alignments_it->second.first << ':' << genomicRegionsToString(transcript_cigar_genomic_regions);
            cout << '\t' << bam_record.ChrName(bam_reader.Header()) << ':' << read_genomic_regions_str;
            cout << '\t' << benchmark_stats_ss.str();
            cout << '\t' << bam_record.MapQuality();
            cout << '\t' << overlap;
            cout << endl;

        } else {   

            addStats(&benchmark_stats, bam_record, overlap, benchmark_stats_ss, &prev_read_name, &prev_mapq, &prev_overlap, &prev_output_string, &cur_primary_mapq, &cur_primary_overlap);  
        }

        if (num_reads % 10000000 == 0) {

            cerr << "Number of analysed reads: " << num_reads << endl;
        }
    }

    if (!debug_output) {

        stringstream benchmark_stats_ss;
        prev_read_name = "";

        addStats(&benchmark_stats, bam_record, 0, benchmark_stats_ss, &prev_read_name, &prev_mapq, &prev_overlap, &prev_output_string, &cur_primary_mapq, &cur_primary_overlap);  
    }

    bam_reader.Close();

    for (auto vcf_rec : vcfs) {
        tbx_destroy(get<2>(vcf_rec));
        bcf_hdr_destroy(get<1>(vcf_rec));
        hts_close(get<0>(vcf_rec));
    }

    if (!debug_output) {

        cout << "Count" << "\t" << base_header.str() << endl;

        uint32_t num_unique_reads = 0;
        double sum_overlap = 0;

        for (auto & stats: benchmark_stats) {

            cout << stats.second.first << "\t" << stats.first << endl;

            num_unique_reads += stats.second.first;
            sum_overlap += stats.second.second;

        }

        cerr << "\nNumber of unique reads: " << num_unique_reads << endl;
        cerr << "Average overlap: " << sum_overlap/num_unique_reads << endl;
    }

    cerr << "\nTotal number of analysed reads: " << num_reads << endl;

	return 0;
}
