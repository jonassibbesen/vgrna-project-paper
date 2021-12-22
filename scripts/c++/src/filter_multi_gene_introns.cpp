
/*
filter_multi_gene_introns
Filter introns overlapping multiple genes.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <assert.h>

#include "SeqLib/GenomicRegion.h"
#include "SeqLib/GenomicRegionCollection.h"

#include "utils.hpp"

using namespace SeqLib;

GRC parseGenes(const string & genes_file, unordered_map<string, uint32_t> * contig_name_to_id_map) {

    GRC gene_intervals;

    ifstream genes_istream(genes_file);
    assert(genes_istream.is_open());

    string line;

    while (genes_istream.good()) {

        getline(genes_istream, line, '\n');

        if (line.empty() || line.front() == '#') {

            continue;
        }

        auto line_split = splitString(line, '\t');

        if (line_split.at(2) == "gene") {

            auto contig_name_to_id_map_it = contig_name_to_id_map->emplace(line_split.at(0), contig_name_to_id_map->size());
            gene_intervals.add(GenomicRegion(contig_name_to_id_map_it.first->second, stoi(line_split.at(3)), stoi(line_split.at(4))));
        }
    }

    genes_istream.close();

    gene_intervals.MergeOverlappingIntervals();
    gene_intervals.CoordinateSort();

    return gene_intervals;
}

int main(int argc, char* argv[]) {

    if (argc != 4) {

        cerr << "Usage: filter_multi_gene_introns <introns_bed_name> <genes_gtf_name> <max_genes_overlap> > introns.bed" << endl;
        return 1;
    }

    printScriptHeader(argc, argv);

    unordered_map<string, uint32_t> contig_name_to_id_map;

    auto gene_intervals = parseGenes(argv[2], &contig_name_to_id_map);
    cerr << "Number of non-overlapping gene intervals: " << gene_intervals.size() << endl;

    uint32_t max_genes_overlap = stoi(argv[3]);

    gene_intervals.CreateTreeMap();

    ifstream introns_istream(argv[1]);
    assert(introns_istream.is_open());

    string line;

    uint32_t num_introns = 0;
    uint32_t num_introns_filt = 0;

    while (introns_istream.good()) {

        getline(introns_istream, line, '\n');

        if (line.empty() || line.front() == '#') {

            continue;
        }

        num_introns++;

        auto line_split = splitString(line, '\t');
        auto contig_name_to_id_map_it = contig_name_to_id_map.find(line_split.at(0));

        if (contig_name_to_id_map_it != contig_name_to_id_map.end()) {

            auto intron_interval = GenomicRegion(contig_name_to_id_map_it->second, stoi(line_split.at(1)) + 1, stoi(line_split.at(2)));
            auto gene_introns_intersection = gene_intervals.FindOverlappedIntervals(intron_interval, true);

            if (gene_introns_intersection.size() > max_genes_overlap) {

                num_introns_filt++;
                continue;
            }
        }

        cout << line << endl;
    }

    introns_istream.close();

    cerr << "\nNumber of introns: " << num_introns << endl;
    cerr << "Number of introns filtered: " << num_introns_filt << endl;

	return 0;
}
