
/*
convert_rsem_sim_info
Converts 
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


vector<pair<string, uint32_t> > parseTranscriptInfo(const string & rsem_expression_file) {

    vector<pair<string, uint32_t> > transcripts;

    ifstream rsem_istream(rsem_expression_file);
    assert(rsem_istream.is_open());

    string id;
    string length;

    while (rsem_istream.good()) {

        getline(rsem_istream, id, '\t');

        if (id.empty() || id == "transcript_id") {

            rsem_istream.ignore(numeric_limits<streamsize>::max(), '\n');
            continue;
        }

        getline(rsem_istream, length, '\t');
        getline(rsem_istream, length, '\t');

        transcripts.emplace_back(id, stoi(length));
        rsem_istream.ignore(numeric_limits<streamsize>::max(), '\n');
    }

    rsem_istream.close();

    return transcripts;
}

int main(int argc, char* argv[]) {

    if (argc != 3) {

        cerr << "Usage: convert_rsem_sim_info <read_bam> <rsem_expression_file> > info.txt" << endl;
        return 1;
    }

    printScriptHeader(argc, argv);

    BamReader bam_reader;
    bam_reader.Open(argv[1]);
    assert(bam_reader.IsOpen());

    auto transcripts = parseTranscriptInfo(argv[2]);
    cerr << "Number of transcripts: " << transcripts.size() << "\n" << endl;

    BamRecord bam_record;

    unordered_set<string> read_info;
    uint32_t num_reads = 0;

    while (bam_reader.GetNextRecord(bam_record)) { 

        if (bam_record.SecondaryFlag()) {

            continue;
        }

        num_reads++;

        auto read_name_split = splitString(bam_record.Qname(), '_');
        assert(read_name_split.size() == 5);

        auto read_transcript_id = transcripts.at(stoi(read_name_split.at(2)) - 1).first;
        uint32_t read_transcript_pos = stoi(read_name_split.at(3));

        auto read_name_end_split = splitString(read_name_split.at(4), '/');
        assert(read_name_end_split.size() <= 2);

        bool is_first = isFirstRead(bam_record); 

        if (!is_first) {

            read_transcript_pos += stoi(read_name_end_split.front()) - bam_record.Length();
        }

        if (read_name_split.at(1) == "1") {

            read_transcript_pos = transcripts.at(stoi(read_name_split.at(2)) - 1).second - read_transcript_pos - bam_record.Length();
        }

        string read_name = bam_record.Qname();

        if (is_first) {

            read_name += "_1";
        
        } else {

            read_name += "_2";            
        }

        string reverse = "0";

        if ((read_name_split.at(1) == "0" && !is_first) || (read_name_split.at(1) == "1" && is_first)) {

            reverse = "1";
        }  

        stringstream read_info_ss;
        read_info_ss << read_name << "\t" << read_transcript_id << "\t" << read_transcript_pos << "\t" << reverse;  

        read_info.emplace(read_info_ss.str());   

        if (num_reads % 10000000 == 0) {

            cerr << "Number of converted reads: " << num_reads << endl;
        }
    }

    bam_reader.Close();

    cout << "read" << "\t" << "path" << "\t" << "offset" << "\t" << "reverse" << endl;

    for (auto & read: read_info) {

        cout << read << endl;
    }

    cerr << "\nNumber of converted unique reads: " << read_info.size() << endl;
    cerr << "Total number of converted reads: " << num_reads << endl;

	return 0;
}
