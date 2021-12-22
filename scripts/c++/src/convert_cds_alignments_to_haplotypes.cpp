
/*
convert_cds_alignments_to_haplotypes
Add intron reference seqeunces to cds allele alignments and
write the resulting genomic haplotype sequences as fasta. 
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>
#include <assert.h>

#include "SeqLib/FastqReader.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/GenomicRegionCollection.h"

#include "utils.hpp"

using namespace SeqLib;

unordered_map<int32_t, string> parseGenome(const string & genome_file, const BamHeader & bam_header) {

    unordered_map<int32_t, string> genome;

    FastqReader fasta_reader(genome_file);    
    UnalignedSequence chrom;

    while (fasta_reader.GetNextSequence(chrom)) {
            
        assert(genome.emplace(bam_header.Name2ID(chrom.Name), chrom.Seq).second);
    }

    return genome;
}

GRC parseSpliceJunctions(const string & transcripts_file, const BamHeader & bam_header) {

    GRC splice_junctions;

    ifstream transcripts_istream(transcripts_file);
    assert(transcripts_istream.is_open());

    string line;

    smatch regex_id_match;
    regex regex_id_exp("transcript_id\\s{1}\"?([^\"]*)\"?;?");

    string prev_transcript_id = "";
    pair<uint32_t, uint32_t> prev_exon(0, 0);

    while (transcripts_istream.good()) {

        getline(transcripts_istream, line, '\n');

        if (line.empty() || line.front() == '#') {

            continue;
        }

        auto line_split = splitString(line, '\t');

        if (line_split.at(2) != "exon") {

            continue;
        }

        string transcript_id = "";

        if (std::regex_search(line_split.at(8), regex_id_match, regex_id_exp)) {

            assert(regex_id_match.size() == 2);
            transcript_id = regex_id_match[1];
        }

        pair<uint32_t, uint32_t> exon(stoi(line_split.at(3)), stoi(line_split.at(4)));

        if (prev_transcript_id == transcript_id) {

            auto splice_start = prev_exon.second;
            auto splice_end = exon.first;

            if (line_split.at(6) == "-") {

                if (exon.first < prev_exon.second) {

                    splice_start = exon.second;
                    splice_end = prev_exon.first;
                }
            }

            splice_junctions.add(GenomicRegion(bam_header.Name2ID(line_split.at(0)), splice_start, splice_end - 2));
        }

        prev_transcript_id = transcript_id;
        prev_exon = exon;
    }

    transcripts_istream.close();

    splice_junctions.CoordinateSort();
    splice_junctions.CreateTreeMap();

    return splice_junctions;
}

int main(int argc, char* argv[]) {

    if (argc != 5) {

        cerr << "Usage: convert_cds_alignments_to_haplotypes <cds_bam> <genome_fasta> <transcripts_gtf> <max_sj_dist> > haplotypes.fa" << endl;
        return 1;
    }

    printScriptHeader(argc, argv);

    BamReader bam_reader;
    bam_reader.Open(argv[1]);
    assert(bam_reader.IsOpen());

    auto genome = parseGenome(argv[2], bam_reader.Header());
    cerr << "Number of chromosomes: " << genome.size() << endl;

    auto splice_junctions = parseSpliceJunctions(argv[3], bam_reader.Header());
    cerr << "Number of splice-junctions: " << splice_junctions.size() << "\n" << endl;

    const uint32_t max_sj_dist = stoi(argv[4]); 

    BamRecord bam_record;
    uint32_t num_alignments = 0;

    while (bam_reader.GetNextRecord(bam_record)) { 

        num_alignments++;

        assert(bam_record.MappedFlag());    
        assert(bam_record.GetCigar().NumQueryConsumed() == bam_record.Length());

        string genomic_query_seq = "";

        uint32_t query_pos = 0;
        uint32_t genomic_pos = bam_record.Position();

        for (auto & field: bam_record.GetCigar()) {

            assert(field.Type() != 'H');
            assert(field.Length() > 0);

            if (field.ConsumesQuery()) { 

                genomic_query_seq += bam_record.Sequence().substr(query_pos, field.Length());
                query_pos += field.Length();
            
            } else {

                assert(field.Type() == 'D' || field.Type() == 'N');
                auto del_region = GenomicRegion(bam_record.ChrID(), genomic_pos, genomic_pos + field.Length() - 1);

                auto sj_indexes = splice_junctions.FindOverlappedIntervals(del_region, true);

                uint32_t best_sj_idx = 0;
                uint32_t best_sj_dist = std::numeric_limits<uint32_t>::max();

                for (auto & sj_idx: sj_indexes) {

                    auto sj = splice_junctions.at(sj_idx);
                    uint32_t sj_dist = max(sj.DistanceBetweenStarts(del_region), sj.DistanceBetweenEnds(del_region));

                    if (sj_dist < best_sj_dist) {

                        best_sj_idx = sj_idx;
                        best_sj_dist = sj_dist;
                    }
                }

                if (best_sj_dist <= max_sj_dist) {

                    auto best_sj = splice_junctions.at(best_sj_idx);
                    assert(best_sj.pos1 <= best_sj.pos2);

                    genomic_query_seq += genome.at(best_sj.chr).substr(best_sj.pos1, best_sj.pos2 - best_sj.pos1 + 1);

                    if (field.Length() < 20) {

                        cerr << "Deletion (<20) in " << bam_record.Qname() << " CDS alignment found (length: " << field.Length() << ", distance: " << best_sj_dist << ")" << endl; 
                    }

                } else if (field.Length() >= 20) {

                    cerr << "Intron (≥20) in " << bam_record.Qname() << " CDS alignment not found (length: " << field.Length() << ", distance: " << best_sj_dist << ")" << endl; 
                }
            }
            
            if (field.ConsumesReference()) { 

                genomic_pos += field.Length();
            }
        }

        if (bam_record.ReverseFlag()) {

            rcomplement(genomic_query_seq);
        }

        cout << ">" + bam_record.Qname() << "\n";
        cout << genomic_query_seq << endl;
    }

    bam_reader.Close();

    cerr << "\nNumber of converted CDS alignments: " << num_alignments << endl;

	return 0;
}
