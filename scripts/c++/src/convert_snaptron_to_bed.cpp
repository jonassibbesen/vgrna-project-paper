
/*
convert_snaptron_to_bed
Convert introns and filter introns based on number 
of samples and read coverage.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <assert.h>

#include "utils.hpp"

using namespace SeqLib;


uint32_t updateEndSiteClusters(map<uint32_t, vector<pair<uint32_t, vector<string> > > > * end_site_clusters, vector<pair<uint32_t, vector<string> > > * start_site_cluster, const uint32_t max_splice_edges) {

    sort(start_site_cluster->begin(), start_site_cluster->end());

    uint32_t splice_edges_added = 0;
    auto start_site_cluster_rit = start_site_cluster->rbegin();

    while (start_site_cluster_rit != start_site_cluster->rend()) {

        if (splice_edges_added == max_splice_edges) {

            break;
        }

        auto end_site_clusters_it = end_site_clusters->emplace(stoi(start_site_cluster_rit->second.at(3)), vector<pair<uint32_t, vector<string> > >());
        end_site_clusters_it.first->second.emplace_back(*start_site_cluster_rit);

        ++splice_edges_added;
        ++start_site_cluster_rit;
    }

    auto num_start_sites = start_site_cluster->size();
    assert(num_start_sites >= splice_edges_added);

    start_site_cluster->clear();
    return (num_start_sites - splice_edges_added);
}

uint32_t writeIntrons(map<uint32_t, vector<pair<uint32_t, vector<string> > > > * end_site_clusters, const uint32_t start_pos, const uint32_t max_splice_edges) {

    uint32_t introns_filtered = 0;
    auto end_site_clusters_it = end_site_clusters->begin();

    while (end_site_clusters_it != end_site_clusters->end()) {

        if (end_site_clusters_it->first >= start_pos) {

            break;
        }

        sort(end_site_clusters_it->second.begin(), end_site_clusters_it->second.end());

        uint32_t introns_written = 0;
        auto end_site_cluster_it = end_site_clusters_it->second.rbegin();

        while (end_site_cluster_it != end_site_clusters_it->second.rend()) {

            if (introns_written == max_splice_edges) {

                break;
            }

            auto & intron_line = end_site_cluster_it->second;

            cout << intron_line.at(1) << "\t" << stoi(intron_line.at(2)) - 1 << "\t" << intron_line.at(3) << "\t" << intron_line.front() << "\t" << end_site_cluster_it->first << "\t" << intron_line.at(5) << endl;

            ++introns_written;
            ++end_site_cluster_it;
        }

        auto num_end_sites = end_site_clusters_it->second.size();
        assert(num_end_sites >= introns_written);

        introns_filtered += num_end_sites - introns_written;

        auto del_end_site_clusters_it = end_site_clusters_it;
        ++end_site_clusters_it;

        end_site_clusters->erase(del_end_site_clusters_it);
    }

    return introns_filtered;
}

int main(int argc, char* argv[]) {

    if (argc != 5) {

        cerr << "Usage: convert_snaptron_to_bed <introns_tsv_name> <min_num_samples> <min_num_reads_per_sample> <max_splice_edges> > introns.bed" << endl;
        return 1;
    }

    printScriptHeader(argc, argv);

    ifstream introns_istream(argv[1]);
    assert(introns_istream.is_open());

    uint32_t min_num_samples = stoi(argv[2]);
    assert(min_num_samples > 0);

    uint32_t min_num_reads_per_sample = stoi(argv[3]);
    assert(min_num_reads_per_sample > 0);

    uint32_t max_splice_edges = stoi(argv[4]);
    assert(max_splice_edges > 0);

    string line;

    uint32_t num_introns = 0;
    uint32_t num_introns_filt_samples = 0;
    uint32_t num_introns_filt_edges = 0;

    vector<pair<uint32_t, vector<string> > > start_site_cluster;
    map<uint32_t, vector<pair<uint32_t, vector<string> > > > end_site_clusters;

    string prev_chrom = "";
    uint32_t prev_start_pos = 0;

    while (introns_istream.good()) {

        getline(introns_istream, line, '\n');

        if (line.empty() || line.front() == '#') {

            continue;
        }

        num_introns++;

        auto line_split = splitString(line, '\t');
        
        uint32_t start_pos = stoi(line_split.at(2));

        if (prev_chrom == line_split.at(1)) {

            assert(prev_start_pos <= start_pos);
        }

        if (prev_chrom != line_split.at(1) || prev_start_pos < start_pos) {

            num_introns_filt_edges += updateEndSiteClusters(&end_site_clusters, &start_site_cluster, max_splice_edges);
            assert(start_site_cluster.empty());
        }

        if (prev_chrom != line_split.at(1)) {

            num_introns_filt_edges += writeIntrons(&end_site_clusters, std::numeric_limits<uint32_t>::max(), max_splice_edges);
            assert(end_site_clusters.empty());

        } else {

            num_introns_filt_edges += writeIntrons(&end_site_clusters, start_pos, max_splice_edges);
        }

        prev_chrom = line_split.at(1);
        prev_start_pos = stoi(line_split.at(2));

        auto samples_split = splitString(line_split.at(11), ',');

        assert(samples_split.front() == "");
        assert(samples_split.size() > 1);

        if ((samples_split.size() - 1) < min_num_samples) {

            ++num_introns_filt_samples;
            continue;
        }

        uint32_t num_read_samples = 0;

        auto samples_split_it = samples_split.begin();
        assert(samples_split_it != samples_split.end());

        ++samples_split_it;
        assert(samples_split_it != samples_split.end());

        while (samples_split_it != samples_split.end()) {

            auto sample_split = splitString(*samples_split_it, ':');
            
            assert(sample_split.size() == 2);
            assert(stoi(sample_split.back()) >= 1);

            if (stoi(sample_split.back()) >= min_num_reads_per_sample) {

                ++num_read_samples;
            }

            ++samples_split_it;
        }

        if (num_read_samples < min_num_samples) {

            ++num_introns_filt_samples;
            continue;
        }

        start_site_cluster.emplace_back(num_read_samples, line_split);
    }

    num_introns_filt_edges += updateEndSiteClusters(&end_site_clusters, &start_site_cluster, max_splice_edges);
    num_introns_filt_edges += writeIntrons(&end_site_clusters, std::numeric_limits<uint32_t>::max(), max_splice_edges);

    assert(start_site_cluster.empty());
    assert(end_site_clusters.empty());

    introns_istream.close();

    cerr << "\nNumber of introns: " << num_introns << endl;
    cerr << "Number of introns filtered (samples): " << num_introns_filt_samples << endl;
    cerr << "Number of introns filtered (edges): " << num_introns_filt_edges << endl;

	return 0;
}
