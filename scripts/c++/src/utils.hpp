
#ifndef VGRNA_PROJECT_SCRIPTS_UTILS_HPP
#define VGRNA_PROJECT_SCRIPTS_UTILS_HPP

#include <string>
#include <vector>
#include <sstream>
#include <assert.h>
#include <sys/stat.h>
#include <limits>
#include <regex>

#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"
#include "SeqLib/GenomicRegion.h"
#include "SeqLib/GenomicRegionCollection.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"


using namespace std;

// Precision used when comparing double variables.
static const double double_precision = numeric_limits<double>::epsilon() * 100;

template<class T>
ostream & operator<<(ostream & os, const vector<T> & values) {

    auto values_it = values.cbegin();

    if (values_it == values.cend()) {

        return os;
    }

    os << *values_it;
    ++values_it;

    while (values_it != values.cend()) {

        os << " " << *values_it;
        ++values_it;
    }

    return os;
}

struct MapQPairHash {

    size_t operator()(const pair<uint32_t, uint32_t> & mapqs) const {

        return hash<uint32_t>{}(mapqs.first) ^ (hash<uint32_t>{}(mapqs.second) << 1); 
    }
};

// Compare double variables using above precision.
bool doubleCompare(const double a, const double b) {

    assert(isfinite(a));
    assert(isfinite(b));

    return ((a == b) or (abs(a - b) < abs(min(a, b)) * double_precision));
}

vector<string> splitString(const string & str, const char delim) {

    stringstream ss(str);
    vector<string> elems;

    for (string item; getline(ss, item, delim);) {

        elems.push_back(item);
    }

    return elems;
}

void printScriptHeader(int argc, char * argv[]) {

    vector<string> argv_vec;
    argv_vec.assign(argv, argv + argc);

    cerr << GIT_COMMIT << endl;
    cerr << argv_vec << "\n" << endl;
}

pair<SeqLib::Cigar, uint32_t> trimCigar(const SeqLib::Cigar & cigar, const uint32_t trim_start, const uint32_t trim_length) {

    if (trim_start > 0 || trim_length < cigar.NumQueryConsumed()) { 

        SeqLib::Cigar trimmed_cigar;
        uint32_t num_shifted_genomic_bases = 0;

        uint32_t cur_query_length = 0;

        for (auto & field: cigar) {

            if (field.ConsumesQuery()) {

                uint32_t new_field_length = 0;

                if (cur_query_length <= trim_start && trim_start <= cur_query_length + field.Length() - 1) {

                    assert(trimmed_cigar.TotalLength() == 0);

                    new_field_length = min(trim_length, cur_query_length + field.Length() - trim_start);
                    assert(new_field_length > 0);

                    if (field.ConsumesReference()) {

                        num_shifted_genomic_bases += (trim_start - cur_query_length);
                    }

                } else if (trimmed_cigar.size() > 0) {

                    new_field_length = min(static_cast<uint32_t>(trim_length - trimmed_cigar.NumQueryConsumed()), field.Length());
                    assert(new_field_length > 0);
                }

                assert(new_field_length <= field.Length());
                cur_query_length += field.Length();

                if (new_field_length > 0) {

                    trimmed_cigar.add(SeqLib::CigarField(field.Type(), new_field_length));
                
                    if (trimmed_cigar.NumQueryConsumed() == trim_length) {

                        break;
                    }   
                } 

            } else if (trimmed_cigar.size() > 0) {

                trimmed_cigar.add(field); 
            }
            
            if (field.ConsumesReference() && trimmed_cigar.size() == 0) {

                num_shifted_genomic_bases += field.Length();
            }
        }
    
        assert(trimmed_cigar.NumQueryConsumed() <= trim_length); 
        return make_pair(trimmed_cigar, num_shifted_genomic_bases);

    } else {

        return make_pair(cigar, 0);
    }
}

void cigarToGenomicRegions(SeqLib::GRC * cigar_genomic_regions, const SeqLib::Cigar & cigar, const uint32_t chrom_idx, const uint32_t start_pos) {
    
    uint32_t cur_length = 0;

    for (auto & field: cigar) {

        if (field.Type() == 'M' || field.Type() == '=' || field.Type() == 'X') {

            cigar_genomic_regions->add(SeqLib::GenomicRegion(chrom_idx, start_pos + cur_length, start_pos + cur_length + field.Length() - 1));
            cur_length += field.Length();

        } else if (field.Type() == 'D' || field.Type() == 'N') {

            cur_length += field.Length();
        }
    }
}

SeqLib::GRC cigarToGenomicRegions(const SeqLib::Cigar & cigar, const uint32_t chrom_idx, const uint32_t start_pos) {

    SeqLib::GRC cigar_genomic_regions;
    cigarToGenomicRegions(&cigar_genomic_regions, cigar, chrom_idx, start_pos);

    return cigar_genomic_regions;
}

uint32_t cigarTypeLength(const SeqLib::Cigar & cigar, const char type) {

    uint32_t type_length = 0;

    for (auto & field: cigar) {

        if (field.Type() == type) {

            type_length += field.Length();
        }
    }

    return type_length;
}

string genomicRegionsToString(const SeqLib::GRC & genomic_regions) {

    stringstream genomic_regions_ss; 

    bool is_first = true;
    for (auto & region: genomic_regions.AsGenomicRegionVector()) {

        if (!is_first) {

            genomic_regions_ss << ",";
        }

        genomic_regions_ss << region.pos1 + 1 << "-" << region.pos2 + 1;
        is_first = false;
    }

    return genomic_regions_ss.str();
}

vector<tuple<htsFile*, bcf_hdr_t*, tbx_t*, int>> initializeVCFs(const vector<string>& vcf_filenames,
                                                                const string& sample_name,
                                                                unordered_map<string, int>& contig_to_vcf_out) {
    
    vector<tuple<htsFile*, bcf_hdr_t*, tbx_t*, int>> vcfs;
    for (const string& vcf_filename : vcf_filenames) {
        
        // make sure the VCF and tabis exist
        string tabix_filename = vcf_filename + ".tbi";
        struct stat stat_tbi, stat_vcf;
        if (stat(vcf_filename.c_str(), &stat_vcf) != 0) {
            cerr << "VCF file " << vcf_filename << " not found." << endl;
            exit(1);
        }
        if (stat(tabix_filename.c_str(), &stat_tbi) != 0) {
            cerr << "Tabix file " << tabix_filename << " not found. Must tabix index VCF file " << vcf_filename << " before running benchmark." << endl;
            exit(1);
        }
        
        // load them up
        htsFile* vcf = bcf_open(vcf_filename.c_str(), "r");
        if (vcf == nullptr) {
            cerr << "error: could not load VCF file " << vcf_filename << endl;
            exit(1);
        }
        
        bcf_hdr_t* header = bcf_hdr_read(vcf);
        if (header == nullptr) {
            cerr << "error: could not read header for VCF file " << vcf_filename << endl;
            exit(1);
        }
        
        tbx_t* tabix_index = tbx_index_load(tabix_filename.c_str());
        if (tabix_index == nullptr) {
            cerr << "error: could not load tabix file " << tabix_filename << endl;
            exit(1);
        }
        
        
        // find the index of the sample we want
        // TODO: there should be a way to do this using the dictionary in the
        // header, but i can't find it...
        int idx = -1;
        for (int i = 0, n = bcf_hdr_nsamples(header); i < n; ++i) {
            if (header->samples[i] == sample_name) {
                idx = i;
                break;
            }
        }
        vcfs.emplace_back(vcf, header, tabix_index, idx);
        
        // record which contigs occur in this VCF
        int num_seq_names = 0;
        const char** contig_names = tbx_seqnames(tabix_index, &num_seq_names);
        for (int j = 0; j < num_seq_names; ++j) {
            string contig = contig_names[j];
            contig_to_vcf_out[contig] = vcfs.size() - 1;
        }
        free(contig_names);
    }
    
    return vcfs;
}
    

inline tuple<int32_t, int32_t, int32_t, int32_t>
    countIndelsAndSubs(const string& chrom, const SeqLib::GRC & regions, htsFile * vcf,
                       bcf_hdr_t * bcf_header, tbx_t * tabix_index, int sample_idx) {
    
    int32_t indel_bps_1 = 0;
    int32_t subs_bps_1 = 0;
    int32_t indel_bps_2 = 0;
    int32_t subs_bps_2 = 0;
        
    for (const SeqLib::GenomicRegion& region : regions) {
        
        // init iteration variables
        int tid = tbx_name2id(tabix_index, chrom.c_str());
        hts_itr_t * itr = tbx_itr_queryi(tabix_index, tid, region.pos1, region.pos2);
        bcf1_t * bcf_record = bcf_init();
        kstring_t kstr = {0, 0, 0};
        
        // iterate over VCF lines
        while (tbx_itr_next(vcf, tabix_index, itr, &kstr) >= 0) {
            
            vcf_parse(&kstr, bcf_header, bcf_record);
            
            // init a genotype array
            int32_t* genotypes = nullptr;
            int arr_size = 0;
            // and query it
            int num_genotypes = bcf_get_genotypes(bcf_header, bcf_record, &genotypes, &arr_size);
            int ploidy = num_genotypes / bcf_hdr_nsamples(bcf_header);
            assert(ploidy <= 2);
            int allele_1 = -1;
            int allele_2 = -1;
            // look at the genotype of this sample
            for (int i = ploidy * sample_idx, n = ploidy * (sample_idx + 1); i < n; ++i) {
                if (genotypes[i] == bcf_int32_vector_end) {
                    // sample has lower ploidy
                    break;
                }
                if (bcf_gt_is_missing(genotypes[i])) {
                    continue;
                }
                if (i == 0) {
                    allele_1 = bcf_gt_allele(genotypes[i]);
                }
                else {
                    allele_2 = bcf_gt_allele(genotypes[i]);
                }
            }
            free(genotypes);
            
            for (bool do_allele_1 : {true, false}) {
                auto allele = do_allele_1 ? allele_1 : allele_2;
                if (allele < 0) {
                    // there is an allele for this haplotype (might not exist from
                    // missing calls or ploidy 1 chromosomes: X, Y, MT)
                    continue;
                }
                
                if (allele == 0) {
                    // this is the reference allele
                    continue;
                }
                
                auto& subs_bps = do_allele_1 ? subs_bps_1 : subs_bps_2;
                auto& indel_bps = do_allele_1 ? indel_bps_1 : indel_bps_2;
    
                // make sure the lazily-unpacked metadata is unpacked through the alleles
                bcf_unpack(bcf_record, BCF_UN_STR);
                
                int ref_len = strlen(bcf_record->d.allele[0]);
                int var_len = strlen(bcf_record->d.allele[allele]);
                
                if (ref_len == var_len) {
                    // substitution
                    subs_bps += max(ref_len, var_len);
                }
                else if (max(ref_len, var_len) > 1 && min(ref_len, var_len) == 1) {
                    // indel
                    indel_bps += max(ref_len, var_len) - 1;
                }
                // else: complex variant, we ignore it for simplicity
            }
        }
        
        bcf_destroy(bcf_record);
        tbx_itr_destroy(itr);
        if (kstr.s) {
            free(kstr.s);
        }
    }
        
    auto return_val = make_tuple(subs_bps_1, indel_bps_1, subs_bps_2, indel_bps_2);
    return return_val;
}

vector<string> parseGenotype(const string & sample) {

    auto genotype_str = splitString(sample, ':');

    if (genotype_str.front().find('/') != std::string::npos) {

        assert(genotype_str.front().find('|') == std::string::npos);
        return splitString(genotype_str.front(), '/');

    } else {

        return splitString(genotype_str.front(), '|');

    }
}

string alleleIdxToSequence(const string allele_idx_str, const vector<string> & variant) {

    if (allele_idx_str == ".") {

        return "";        
    
    } else {

        auto allele_idx = stoi(allele_idx_str); 
        
        if (allele_idx == 0) {

            return variant.at(3);

        } else {

            return splitString(variant.at(4), ',').at(allele_idx - 1);
        }
    }
}

void rightTrim(string * allele, const uint32_t trim_length) {

    if (trim_length >= allele->size()) {

        *allele = "";
    
    } else {

        *allele = allele->substr(0, allele->size() - trim_length);
    }
}

void leftTrim(string * allele, const uint32_t trim_length) {

    if (trim_length >= allele->size()) {

        *allele = "";
    
    } else {

        *allele = allele->substr(trim_length);
    }
}

uint32_t trimAlleles(string * ref_allele, string * alt_allele) {

    assert(!ref_allele->empty());
    assert(!alt_allele->empty());

    if (*ref_allele == *alt_allele) {

        *ref_allele = "";
        *alt_allele = "";

        return 1;
    }

    uint32_t right_trim_len = 0;

    for (size_t i = 0; i < min(ref_allele->size(), alt_allele->size()); ++i) {

        if (ref_allele->at(ref_allele->size() - i - 1) == alt_allele->at(alt_allele->size() - i - 1)) {

            right_trim_len++;
        
        } else {

            break;
        }
    }

    if (right_trim_len > 0) {

        rightTrim(ref_allele, right_trim_len);
        rightTrim(alt_allele, right_trim_len);
    }

    uint32_t left_trim_len = 0;

    for (size_t i = 0; i < min(ref_allele->size(), alt_allele->size()); ++i) {

        if (ref_allele->at(i) == alt_allele->at(i)) {

            left_trim_len++;
        
        } else {

            break;
        }
    }

    if (left_trim_len > 0) {

        leftTrim(ref_allele, left_trim_len);
        leftTrim(alt_allele, left_trim_len);
    }

    return left_trim_len;
}

string getAlleleType(string ref_allele, string alt_allele) {

    if (ref_allele == alt_allele) {

        return "REF";

    } else if (ref_allele.size() == alt_allele.size() and ref_allele.size() == 1) {

        return "SNV";

    } else if (ref_allele.empty()) {

        assert(!alt_allele.empty());
        return "INS";

    } else if (alt_allele.empty()) {

        assert(!ref_allele.empty());
        return "DEL";        

    } else {

        assert(!ref_allele.empty());
        assert(!alt_allele.empty());
        return "COM";
    }
}

uint32_t getAllelicMapQ(const SeqLib::BamRecord & bam_record) {

    int32_t allelic_mapq = -1;

    if (!bam_record.GetIntTag("AQ", allelic_mapq)) {

        // if allelic mapq isn't annotated, it's assumed to be the overall mapq
        allelic_mapq = bam_record.MapQuality();
    }

    assert(allelic_mapq >= 0);
    return allelic_mapq;
}

uint32_t getGroupMapQ(const SeqLib::BamRecord & bam_record) {

    int32_t group_mapq = -1;

    if (!bam_record.GetIntTag("GM", group_mapq)) {

        // if group mapq isn't annotated, it's assumed to be the overall mapq
        group_mapq = bam_record.MapQuality();
    }

    assert(group_mapq >= 0);
    return group_mapq;
}

SeqLib::GRC parseRegionsBed(const string & region_bed_file, const SeqLib::BamHeader & bam_header) {

    SeqLib::GRC regions;
    assert(regions.ReadBED(region_bed_file, bam_header));

    auto regions_it = regions.begin();

    while (regions_it != regions.end()) {

        regions_it->pos2--;
        ++regions_it;
    }

    regions.MergeOverlappingIntervals();
    regions.CreateTreeMap();

    return regions;
}

SeqLib::GRC parseSpliceJunctions(const string & transcripts_file, const SeqLib::BamHeader & bam_header) {

    SeqLib::GRC splice_junctions;

    ifstream transcripts_istream(transcripts_file);
    assert(transcripts_istream.is_open());

    string line;

    smatch regex_id_match;
    regex regex_id_exp("transcript_id=([^;]*);?");

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

            splice_junctions.add(SeqLib::GenomicRegion(bam_header.Name2ID(line_split.at(0)), splice_start, splice_end - 2));
        }

        prev_transcript_id = transcript_id;
        prev_exon = exon;
    }

    transcripts_istream.close();

    splice_junctions.CoordinateSort();
    splice_junctions.CreateTreeMap();

    return splice_junctions;
}

bool isFirstRead(const SeqLib::BamRecord & bam_record) {

    if (bam_record.PairedFlag()) {

        return bam_record.FirstFlag();
    }

    auto read_name = bam_record.Qname();

    if (read_name.size() > 2) {

        if (read_name.substr(read_name.size() - 2, 1) == "/" || read_name.substr(read_name.size() - 2, 1) == "_") {

            if (read_name.substr(read_name.size() - 1, 1) == "1") {

                return true;
            
            } else {

                assert(read_name.substr(read_name.size() - 1, 1) == "2");
                return false;
            }
        } 
    } 

    return true;
}

#endif
