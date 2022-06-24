
/*
calc_read_transcript_overlap_stats
Calculates overlapping statistics between mapped reads and 
regions in bed format.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <assert.h>

#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/GenomicRegionCollection.h"

#include "utils.hpp"

using namespace SeqLib;


int main(int argc, char* argv[]) {

    if (!(argc == 3 || argc == 4)) {

        cerr << "Usage: calc_read_transcript_overlap_stats <read_bam> <region_bed> (<enable_debug_output>) > statistics.txt" << endl;
        return 1;
    }

    printScriptHeader(argc, argv);

    BamReader bam_reader;
    bam_reader.Open(argv[1]);
    assert(bam_reader.IsOpen());

    auto transcript_regions = parseRegionsBed(argv[2], bam_reader.Header());

    cerr << "Total length of transcript regions: " << transcript_regions.TotalWidth() << "\n" << endl;

    bool debug_output = (argc == 4);

    stringstream base_header; 
    base_header << "IsMapped" << "\t" << "MapQ" << "\t" << "AllelicMapQ" << "\t" << "Length" << "\t" << "InsertionLength" << "\t" << "SoftClipLength" << "\t" << "Overlap";

    if (debug_output) {

        cout << "Name" << "\t" << base_header.str() << endl;
    }

    BamRecord bam_record;

    unordered_map<string, uint32_t> overlap_stats;

    uint32_t num_reads = 0;
    double sum_overlap = 0;

    while (bam_reader.GetNextRecord(bam_record)) { 

        if (bam_record.SecondaryFlag()) {

            continue;
        }

        num_reads++;

        uint32_t insertion_length = 0;
        uint32_t soft_clip_length = 0;

        double overlap = 0;

        if (bam_record.MappedFlag()) { 
            
            assert(bam_record.GetCigar().NumQueryConsumed() == bam_record.Length());

            insertion_length = cigarTypeLength(bam_record.GetCigar(), 'I');
            soft_clip_length = cigarTypeLength(bam_record.GetCigar(), 'S');

            auto read_cigar_genomic_regions = cigarToGenomicRegions(bam_record.GetCigar(), bam_record.ChrID(), bam_record.Position());
            read_cigar_genomic_regions.CreateTreeMap();

            auto cigar_genomic_regions_intersection = transcript_regions.Intersection(read_cigar_genomic_regions, true);
            overlap = cigar_genomic_regions_intersection.TotalWidth() / static_cast<double>(bam_record.Length());
        }

        stringstream overlap_stats_ss;
        
        overlap_stats_ss << bam_record.MappedFlag();
        overlap_stats_ss << "\t" << bam_record.MapQuality();
        overlap_stats_ss << "\t" << getAllelicMapQ(bam_record);
        overlap_stats_ss << "\t" << bam_record.Length();
        overlap_stats_ss << "\t" << insertion_length;
        overlap_stats_ss << "\t" << soft_clip_length;
        overlap_stats_ss << "\t" << overlap;   

        if (debug_output) {

            cout << bam_record.Qname();
            cout << "\t" << overlap_stats_ss.str();
            cout << endl;

        } else {

            auto overlap_stats_it = overlap_stats.emplace(overlap_stats_ss.str(), 0);
            overlap_stats_it.first->second++;  
        }

        sum_overlap += overlap;

        if (num_reads % 10000000 == 0) {

            cerr << "Number of analysed reads: " << num_reads << endl;
        }        
    }

    bam_reader.Close();

    if (!debug_output) {

        cout << "Count" << "\t" << base_header.str() << endl;

        for (auto & stats: overlap_stats) {

            cout << stats.second << "\t" << stats.first << endl;
        }
    }

    cerr << "\nTotal number of analysed reads: " << num_reads << endl;
    cerr << "Average overlap: " << sum_overlap/num_reads << endl;

	return 0;
}
