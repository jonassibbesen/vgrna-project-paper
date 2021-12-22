
/*
calc_allele_rsem_read_coverage
Calculates mapped read coverage across alleles in a diploid vcf file. 
The read names should end with the haplotype origin (e.g. "_h1")
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <assert.h>

#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"

#include "utils.hpp"

using namespace SeqLib;

void printAlleleReadCoverage(BamReader * bam_reader, const vector<string> & line_split, const uint32_t up_pos, const uint32_t down_pos, const uint32_t allele_id, const string & allele_type, const int32_t rel_allele_length) {

    BamRecord bam_record;

    unordered_map<pair<uint32_t, uint32_t>, pair<uint32_t, uint32_t>, MapQPairHash> mapq_read_counts;

    GenomicRegion genomic_region_up(line_split.at(0), to_string(up_pos), to_string(up_pos + 1), bam_reader->Header());
    assert(bam_reader->SetRegion(genomic_region_up));

    while (bam_reader->GetNextRecord(bam_record)) { 

        if (bam_record.SecondaryFlag()) {

            continue;
        }

        auto mapq_read_counts_it = mapq_read_counts.emplace(make_pair(bam_record.MapQuality(), getAllelicMapQ(bam_record)), make_pair(0, 0));
        mapq_read_counts_it.first->second.first++;
    }

    GenomicRegion genomic_region_down(line_split.at(0), to_string(down_pos), to_string(down_pos + 1), bam_reader->Header());
    assert(bam_reader->SetRegion(genomic_region_down));

    while (bam_reader->GetNextRecord(bam_record)) { 

        if (bam_record.SecondaryFlag()) {

            continue;
        }

        auto mapq_read_counts_it = mapq_read_counts.emplace(make_pair(bam_record.MapQuality(), getAllelicMapQ(bam_record)), make_pair(0, 0));
        mapq_read_counts_it.first->second.second++;
    }

    for (auto & mapq_count: mapq_read_counts) {

        cout << "1";
        cout << "\t" << mapq_count.first.first;
        cout << "\t" << mapq_count.first.second;
        cout << "\t" << line_split.at(0) << ":" << stoi(line_split.at(1));
        cout << "\t" << allele_id;
        cout << "\t" << allele_type;
        cout << "\t" << rel_allele_length;
        cout << "\t" << mapq_count.second.first;
        cout << "\t" << mapq_count.second.second;
        cout << endl;
    }
}


int main(int argc, char* argv[]) {

    if (argc != 4) {

        cerr << "Usage: calc_allele_read_coverage <read_bam_h1> <read_bam_h2> <variant_vcf> > coverage.txt" << endl;
        return 1;
    }

    printScriptHeader(argc, argv);

    BamReader bam_reader_h1;
    bam_reader_h1.Open(argv[1]);
    assert(bam_reader_h1.IsOpen());


    BamReader bam_reader_h2;
    bam_reader_h2.Open(argv[2]);
    assert(bam_reader_h2.IsOpen());

    ifstream vcf_istream(argv[3]);
    assert(vcf_istream.is_open());

    cout << "Count" << "\t" << "MapQ" << "\t" << "AllelicMapQ" << "\t" << "VariantPosition" << "\t" << "AlleleId" << "\t" << "AlleleType" << "\t" << "RelativeAlleleLength" << "\t" << "UpReadCount" << "\t" << "DownReadCount" << endl;

    string line;
    uint32_t num_variants = 0;

    while (vcf_istream.good()) {

        getline(vcf_istream, line);

        if (line.empty() || line.front() == '#') {

            continue;
        }

        num_variants++;

        auto line_split = splitString(line, '\t');
        assert(line_split.size() == 10);

        auto genotype = parseGenotype(line_split.at(9));
        assert(genotype.size() == 2);

        if (genotype.front() == "." || genotype.back() == ".") {

            continue;
        }

        auto pos = stoi(line_split.at(1)) - 1;

        auto ref_seq_h1 = line_split.at(3);
        auto ref_seq_h2 = line_split.at(3);

        auto allele_seq_h1 = alleleIdxToSequence(genotype.front(), line_split);
        auto allele_seq_h2 = alleleIdxToSequence(genotype.back(), line_split);

        auto up_trim_h1 = trimAlleles(&ref_seq_h1, &allele_seq_h1);
        auto up_trim_h2 = trimAlleles(&ref_seq_h2, &allele_seq_h2);

        auto up_pos = pos + min(up_trim_h1, up_trim_h2) - 1;
        auto down_pos = up_pos + max(ref_seq_h1.size(), ref_seq_h2.size()) + 1;

        printAlleleReadCoverage(&bam_reader_h1, line_split, up_pos, down_pos, 1, getAlleleType(ref_seq_h1, allele_seq_h1), allele_seq_h1.size() - ref_seq_h1.size());
        printAlleleReadCoverage(&bam_reader_h2, line_split, up_pos, down_pos, 2, getAlleleType(ref_seq_h2, allele_seq_h2), allele_seq_h2.size() - ref_seq_h2.size());

        if (num_variants % 10000 == 0) {

            cerr << "Number of analysed variants: " << num_variants << endl;
        }        
    }

    bam_reader_h1.Close();
    bam_reader_h2.Close();

    vcf_istream.close();

    cerr << "\nTotal number of analysed variants: " << num_variants << endl;

	return 0;
}
