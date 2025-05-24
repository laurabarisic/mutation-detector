#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

string DATA_DIR = "../data/";

struct SamRecord {
    string qname;
    int flag;
    string rname;
    int64_t pos;
    string cigar;
    string seq;
};

struct PosVotes {
    int none = 0;
    int deleted = 0;
    int inserted = 0;
    int substituted = 0;
    vector<char> substitutionBases;
    vector<char> insertionBases;
};

string reverse(const string& seq)
{
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

string read_fasta(const string& filename)
{
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

bool parse_sam_line(const string& line, SamRecord& record)
{
    istringstream iss(line);
    vector<string> fields;
    string field;

    while (getline(iss, field, '\t')) {
        fields.push_back(field);
    }

    if (fields.size() < 11 || stoi(fields[1]) == 4) {
        return false;
    }

    record.qname = fields[0];
    record.flag = stoi(fields[1]);
    record.rname = fields[2];
    record.pos = stoll(fields[3]);
    record.cigar = fields[5];
    record.seq = fields[9];

    if (record.flag == 16)
        record.seq = reverse(record.seq);

    return true;
}

vector<SamRecord> read_sam(const string& filename)
{
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

void voting(unordered_map<int64_t, PosVotes>& dict,
    unordered_map<int64_t, pair<string, string>>& final_dict)
{
    string max_votes;
    string max_base;
    for (const auto& [pos, votes] : dict) {
        max_votes = "none";
        max_base = "-";

        if (votes.none >= (votes.deleted + votes.inserted + votes.substituted)) {
            continue;
        } else if (votes.deleted >= (votes.none + votes.inserted + votes.substituted)) {
            max_votes = "D";
            max_base = "-";
        } else if (votes.inserted >= (votes.none + votes.deleted + votes.substituted)) {
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
        } else if (votes.substituted >= (votes.none + votes.inserted + votes.deleted)) {
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
        }

        final_dict[pos] = { max_votes, max_base };
    }
}

void mutations(const vector<SamRecord>& sam_records,
    unordered_map<int64_t, PosVotes>& dict,
    const string& fasta_sequence,
    unordered_map<int64_t, pair<string, string>>& final_dict)
{
    ofstream matchingFile("matching.txt");
    if (!matchingFile) {
        cerr << "Greška pri otvaranju datoteke matching.txt" << endl;
        return;
    }
    for (const SamRecord& record : sam_records) {

        int64_t refPos = record.pos - 1;
        int64_t readPos = 0;
        int i = 0;

        while (i < record.cigar.length()) {
            string num = "";
            // int64_t length = 0;
            // while (i < record.cigar.length() && isdigit(record.cigar[i])) {
            //     length *= 10;
            //     length += record.cigar[i] - '0';
            //     i++;
            // }

            while (i < record.cigar.length() && isdigit(record.cigar[i])) {
                num += record.cigar[i];
                i++;
            }

            if (i >= record.cigar.length())
                break;
            char op = record.cigar[i];
            int64_t length = stoll(num);
            i++;

            if (op == 'M') {
                for (int j = 0; j < length; ++j) {
                    if (refPos >= fasta_sequence.size() || readPos >= record.seq.size())
                        break;

                    char refBase = fasta_sequence[refPos];
                    char readBase = record.seq[readPos];

                    // cout << refPos << " - REF_BASE " << refBase << " | " << readPos << " - READ_BASE " << readBase << endl;
                    matchingFile << refPos << " - REF_BASE " << refBase << " | " << readPos << " - READ_BASE " << readBase << endl;

                    if (refBase == readBase) {
                        dict[refPos].none++;
                    } else {
                        dict[refPos].substituted++;
                        dict[refPos].substitutionBases.push_back(readBase);
                    }
                    refPos++;
                    readPos++;
                }

            } else if (op == 'I') {
                for (int j = 0; j < length; ++j) {
                    if (readPos >= record.seq.size())
                        break;

                    dict[refPos].inserted++;
                    dict[refPos].insertionBases.push_back(record.seq[readPos]);
                    readPos++;
                }
            } else if (op == 'D') {
                for (int j = 0; j < length; ++j) {
                    if (refPos >= fasta_sequence.size() || readPos >= record.seq.size())
                        break;

                    dict[refPos].deleted++;
                    refPos++;
                }
            } else if (op == 'S') {
                readPos += length;
            }
        }
    }
    matchingFile.close();

    ofstream votingFile("voting.txt");
    // Zapisivanje glasova u datoteku
    votingFile << "\n--- Glasovi po pozicijama u referentnom genomu ---\n";
    for (const auto& [pos, votes] : dict) {
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

    // Zatvaranje datoteke
    votingFile.close();

    voting(dict, final_dict);
}

int main()
{
    string fasta_path = DATA_DIR + "lambda.fasta";
    string sam_path = DATA_DIR + "lambda.sam";

    string fasta_sequence = read_fasta(fasta_path);
    cout << "FASTA duljina: " << fasta_sequence.size() << " znakova\n";

    vector<SamRecord> sam_records = read_sam(sam_path);
    cout << "\nUkupno očitanih SAM zapisa: " << sam_records.size() << "\n";

    for (const auto& r : sam_records) {
        cout << "QNAME: " << r.qname << " | FLAG: " << r.flag
             << " | RNAME: " << r.rname << " | POS: " << r.pos
             << " | CIGAR: " << r.cigar << " | SEQ: " << r.seq << "\n";
    }

    unordered_map<int64_t, PosVotes> dict;
    unordered_map<int64_t, pair<string, string>> final_dict;
    mutations(sam_records, dict, fasta_sequence, final_dict);

       // zapis u CSV datoteku
    ofstream outfile("mutations.csv");
    if (!outfile) {
        cerr << "Greška pri otvaranju datoteke za pisanje mutacija." << endl;
        return 1;
    }

    vector<pair<int64_t, pair<string, string>>> sorted_mutations(final_dict.begin(),
        final_dict.end());
    sort(sorted_mutations.begin(), sorted_mutations.end());

    for (const auto& [pos, result] : sorted_mutations) {
        if (result.first != "none") {
            outfile << result.first << "," << pos << "," << result.second << "\n";
        }
    }
    outfile.close();

    return 0;
}
