#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

struct SamRecord {
    string qname;
    int flag;
    string rname;
    long pos;
    string cigar;
    string seq;
};

struct PosVotes {
    int none = 0;
    int deleted = 0;
    int inserted = 0;
    int substituted = 0;
};

string reverse(const string& seq){
    string reversed_seq = "";
    for (int i = seq.length()-1; i >= 0; i--){
        char base = seq[i];
        char new_base;

        if (base == 'A') new_base = 'T';
        else if (base == 'T') new_base = 'A';
        else if (base == 'G') new_base = 'C';
        else if (base == 'C') new_base = 'G';
        else new_base = base; // za nepoznate baze (npr. N)

        reversed_seq += new_base;
    } 
    cout << "original: " << seq << ", reversed: " << reversed_seq << endl;
    return reversed_seq;
}


// čitanje FASTA datoteke (spajanje svih linija bez zaglavlja u jedan string)
string read_fasta(const string& filename) {
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

// parsiranje linije iz SAM datoteke
bool parse_sam_line(const string& line, SamRecord& record) {
    istringstream iss(line);
    vector<string> fields;
    string field;

    while (getline(iss, field, '\t')) {
        fields.push_back(field);
    }

    if (fields.size() < 11 || stoi(fields[1]) == 4) return false;

    record.qname = fields[0];
    record.flag = stoi(fields[1]);
    record.rname = fields[2];
    record.pos = stol(fields[3]);
    record.cigar = fields[5];
    record.seq = fields[9];

    if (record.flag == 16) record.seq  = reverse(record.seq);

    return true;
}

// čitanje SAM datoteke i spremanje zapisa
vector<SamRecord> read_sam(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Greška pri otvaranju SAM datoteke." << endl;
        return {};
    }

    vector<SamRecord> records;
    string line;

    while (getline(file, line)) {
        if (!line.empty() && line[0] == '@') continue;

        SamRecord record;
        if (parse_sam_line(line, record)) {
            records.push_back(record);
        }
    }

    return records;
}

int main() {
    string fasta_path = "/home/laura/mutation-detector/data/lambda.fasta";
    string sam_path = "/home/laura/mutation-detector/data/lambda.sam";

    // učitanje FASTA
    string fasta_sequence = read_fasta(fasta_path);
    cout << "FASTA duljina: " << fasta_sequence.size() << " znakova\n";
    cout << fasta_sequence << "\n";

    // učitanje SAM
    vector<SamRecord> sam_records = read_sam(sam_path);
    cout << "\nUkupno očitanih SAM zapisa: " << sam_records.size() << "\n";

    // ispisivanje zapisa iz SAM datoteke
    for (size_t i = 0; i < sam_records.size(); ++i) {
        const SamRecord& r = sam_records[i];
        cout << "QNAME: " << r.qname << " | FLAG: " << r.flag << " | RNAME: " << r.rname
             << " | POS: " << r.pos << " | CIGAR: " << r.cigar << " | SEQ: " << r.seq << "\n";
    }

    string test_seq = "ACGTTGCA";
    string reversed = reverse(test_seq);

    
    return 0;
}
