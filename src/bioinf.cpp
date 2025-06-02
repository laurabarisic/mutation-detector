#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

string DATA_DIR = "../data/";

// Structure for storing SAM record information
struct SamRecord {
    string qname;
    int flag;
    string rname;
    int64_t pos;
    string cigar;
    string seq;
};

// Structure for storing votes at a position in the reference genome
struct PosVotes {
    int none = 0;
    int deleted = 0;
    int inserted = 0;
    int substituted = 0;
    vector<char> substitutionBases;
    vector<char> insertionBases;
};

// Function for reverse complementing a sequence
string reversee(const string &seq) {
    string reversed_seq(seq.length(), ' ');
    size_t index = 0;
    for (int i = (int)seq.length() - 1; i >= 0; i--) {
        char base = seq[i];
        char new_base;
        switch (base) {
        case 'A':
            new_base = 'T';
            break;
        case 'T':
            new_base = 'A';
            break;
        case 'G':
            new_base = 'C';
            break;
        case 'C':
            new_base = 'G';
            break;
        default:
            new_base = base;
        }
        reversed_seq[index++] = new_base;
    }
    return reversed_seq;
}

// Function to read a FASTA file and return the sequence as a string
string read_fasta(const string &filename) {
    ifstream file(filename);
    if (!file) {
        cout << "Greška pri otvaranju FASTA datoteke." << endl;
        return "";
    }

    string line, sequence;
    while (getline(file, line)) {
        if (!line.empty() && line[0] != '>') {
            sequence += line;
        }
    }

    return sequence;
}

// Function to parse a SAM line and fill the SamRecord structure
// Returns true if the line is successfully parsed, otherwise false
bool parse_sam_line(const string &line, SamRecord &record) {
    istringstream iss(line);
    vector<string> fields;
    string field;

    while (getline(iss, field, '\t')) {
        fields.push_back(field);
    }

    if (fields.size() < 11 || stoi(fields[1]) & 4) {
        return false;
    }

    record.qname = fields[0];
    record.flag = stoi(fields[1]);
    record.rname = fields[2];
    record.pos = stoll(fields[3]);
    record.cigar = fields[5];
    record.seq = fields[9];

    // If this is commented we get better results
    // but it is not correct according to the SAM specification
    if (record.flag & 16)
        record.seq = reversee(record.seq);

    return true;
}

// Function to read a SAM file and return a vector of SamRecord structures
vector<SamRecord> read_sam(const string &filename) {
    ifstream file(filename);
    if (!file) {
        cout << "Greška pri otvaranju SAM datoteke." << endl;
        return {};
    }

    vector<SamRecord> records;
    string line;

    while (getline(file, line)) {
        if (!line.empty() && line[0] == '@')
            continue;

        SamRecord record;
        if (parse_sam_line(line, record)) {
            records.push_back(record);
        }
    }

    return records;
}

// Function to perform voting on the mutations
// at each position in the reference genome
void voting(unordered_map<int64_t, PosVotes> &dict,
            unordered_map<int64_t, pair<string, string>> &final_dict) {
    string max_votes;
    string max_base;
    for (const auto &[pos, votes] : dict) {
        max_votes = "none";
        max_base = "-";

        int total_votes =
            votes.none + votes.deleted + votes.inserted + votes.substituted;

        if (votes.none >= votes.substituted && votes.none >= votes.inserted &&
            votes.none >= votes.deleted && votes.none > 2 &&
            votes.none >= ceil(0.4 * total_votes)) {
            continue;
        } else if (votes.substituted >= votes.inserted &&
                   votes.substituted >= votes.deleted &&
                   votes.substituted >= votes.none && votes.substituted > 2 &&
                   votes.substituted >= ceil(0.4 * total_votes)) {
            max_votes = "X";
            unordered_map<char, int> freq;
            char max_elem = '\0';
            int max_count = 0;
            for (char c : votes.substitutionBases) {
                int count = ++freq[c];
                if (count > max_count) {
                    max_count = count;
                    max_elem = c;
                }
            }
            max_base = string(1, max_elem);
        } else if (votes.inserted >= votes.substituted &&
                   votes.inserted >= votes.deleted &&
                   votes.inserted >= votes.none && votes.inserted > 2 &&
                   votes.inserted >= ceil(0.4 * total_votes)) {
            max_votes = "I";
            unordered_map<char, int> freq;
            char max_elem = '\0';
            int max_count = 0;
            for (char c : votes.insertionBases) {
                int count = ++freq[c];
                if (count > max_count) {
                    max_count = count;
                    max_elem = c;
                }
            }
            max_base = string(1, max_elem);
        } else if (votes.deleted >= votes.substituted &&
                   votes.deleted >= votes.inserted &&
                   votes.deleted >= votes.none && votes.deleted > 2 &&
                   votes.deleted >= ceil(0.4 * total_votes)) {
            max_votes = "D";
            max_base = "-";
        }

        final_dict[pos] = {max_votes, max_base};
    }
}

// Function to process and identify mutations
void mutations(const vector<SamRecord> &sam_records,
               unordered_map<int64_t, PosVotes> &dict,
               const string &fasta_sequence,
               unordered_map<int64_t, pair<string, string>> &final_dict) {
    ofstream matchingFile(DATA_DIR + "matching.txt");
    if (!matchingFile) {
        cerr << "Greška pri otvaranju datoteke matching.txt" << endl;
        return;
    }
    for (const SamRecord &record : sam_records) {

        int64_t refPos = record.pos - 1;
        int64_t readPos = 0;
        int i = 0;

        while (i < record.cigar.length()) {
            int64_t length = 0;
            while (i < record.cigar.length() && isdigit(record.cigar[i])) {
                length *= 10;
                length += record.cigar[i] - '0';
                i++;
            }

            if (i >= record.cigar.length())
                break;
            char op = record.cigar[i];
            i++;

            // Check mutation operation and process accordingly
            if (op == 'M') {
                for (int j = 0; j < length; ++j) {
                    if (refPos >= fasta_sequence.size() ||
                        readPos >= record.seq.size())
                        break;

                    char refBase = fasta_sequence[refPos];
                    char readBase = record.seq[readPos];

                    if (refBase == readBase) {
                        matchingFile << "refpos " << refPos << " - REF_BASE "
                                     << refBase << " | " << "readpos "
                                     << readPos + 1 << " - READ_BASE "
                                     << readBase << " [MATCH]" << endl;
                        dict[refPos].none++;
                    } else {
                        matchingFile << "refpos " << refPos << " - REF_BASE "
                                     << refBase << " | " << "readpos "
                                     << readPos + 1 << " - READ_BASE "
                                     << readBase << " [MISS]" << endl;
                        dict[refPos].substituted++;
                        dict[refPos].substitutionBases.push_back(readBase);
                    }
                    refPos++;
                    readPos++;
                }

            } else if (op == 'I') {
                for (int j = 0; j < length; ++j) {
                    if (readPos + j >= record.seq.size())
                        break;

                    char refBase = fasta_sequence[refPos];
                    char readBase = record.seq[readPos + j];

                    matchingFile << "refpos " << refPos << " - REF_BASE "
                                 << refBase << " | " << "readpos "
                                 << (readPos + j + 1) << " - READ_BASE "
                                 << readBase << " [INSERT]" << endl;

                    dict[refPos].insertionBases.push_back(readBase);
                    dict[refPos].inserted++;
                }
                readPos += length;
            } else if (op == 'D') {
                for (int j = 0; j < length; ++j) {
                    if (refPos >= fasta_sequence.size() ||
                        readPos >= record.seq.size())
                        break;

                    char refBase = fasta_sequence[refPos];
                    char readBase = record.seq[readPos];

                    matchingFile << "refpos " << refPos << " - REF_BASE "
                                 << refBase << " | " << "readpos "
                                 << readPos + 1 << " - READ_BASE " << readBase
                                 << " [DELETE]" << endl;
                    dict[refPos].deleted++;
                    refPos++;
                }
            } else if (op == 'S') {
                readPos += length;
            }
        }
    }
    matchingFile.close();

    ofstream votingFile(DATA_DIR + "voting.txt");
    // Writing votes to the file
    votingFile << "\n--- Glasovi po pozicijama u referentnom genomu ---\n";
    for (const auto &[pos, votes] : dict) {
        votingFile << "Pozicija: " << pos << "\n";
        votingFile << "  - Podudaranja (none): " << votes.none << "\n";
        votingFile << "  - Supstitucije: " << votes.substituted << " ";
        if (!votes.substitutionBases.empty()) {
            votingFile << "[";
            for (char base : votes.substitutionBases)
                votingFile << base << ",";
            votingFile << "]";
        }
        votingFile << "\n";
        votingFile << "  - Brisanja (deleted): " << votes.deleted << "\n";
        votingFile << "  - Umetanja (inserted): " << votes.inserted << " ";
        if (!votes.insertionBases.empty()) {
            votingFile << "[";
            for (char base : votes.insertionBases)
                votingFile << base << ",";
            votingFile << "]";
        }
        votingFile << "\n\n";
    }
    // Closing the voting file
    votingFile.close();

    voting(dict, final_dict);
}

int main() {
    auto start_time = std::chrono::high_resolution_clock::now();

    string fasta_path = DATA_DIR + "lambda.fasta";
    string sam_path = DATA_DIR + "lambda.sam";

    string fasta_sequence = read_fasta(fasta_path);
    cout << "FASTA duljina: " << fasta_sequence.size() << " znakova\n";

    vector<SamRecord> sam_records = read_sam(sam_path);
    cout << "\nUkupno očitanih SAM zapisa: " << sam_records.size() << "\n";

    for (const auto &r : sam_records) {
        cout << "QNAME: " << r.qname << " | FLAG: " << r.flag
             << " | RNAME: " << r.rname << " | POS: " << r.pos
             << " | CIGAR: " << r.cigar << " | SEQ: " << r.seq << "\n";
    }

    unordered_map<int64_t, PosVotes> dict;
    unordered_map<int64_t, pair<string, string>> final_dict;
    mutations(sam_records, dict, fasta_sequence, final_dict);

    // Writing to CSV file
    ofstream outfile(DATA_DIR + "mutations.csv");
    if (!outfile) {
        cerr << "Greška pri otvaranju datoteke za pisanje mutacija." << endl;
        return 1;
    }

    vector<pair<int64_t, pair<string, string>>> sorted_mutations(
        final_dict.begin(), final_dict.end());
    sort(sorted_mutations.begin(), sorted_mutations.end());

    for (const auto &[pos, result] : sorted_mutations) {
        if (result.first != "none") {
            outfile << result.first << "," << pos << "," << result.second
                    << "\n";
        }
    }
    outfile.close();

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;

    cout << "Vrijeme izvođenja: " << elapsed_seconds.count() << " sekundi\n";
    return 0;
}