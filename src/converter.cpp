#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

string DATA_DIR = "../data/";

// Determines mutation type based on REF and ALT values from VCF
// Returns:
//   "X" - substitution (same length, single base change)
//   "I" - insertion (ALT longer than REF)
//   "D" - deletion (REF longer than ALT)
string tip_mutacije(const string& ref, const string& alt_raw) {
    string alt = alt_raw;

    // If multiple ALT values (comma-separated), only consider the first
    size_t comma_pos = alt.find(',');
    if (comma_pos != string::npos) {
        alt = alt.substr(0, comma_pos);
    }

    // Determine mutation type by comparing lengths of REF and ALT
    if (ref.length() == 1 && alt.length() == 1) {
        return "X";
    } else if (ref.length() < alt.length()) {
        return "I"; 
    } else if (ref.length() > alt.length()) {
        return "D"; 
    } else {
        return "X"; 
    }
}

int main() {
    // Open input VCF file and output CSV file
    ifstream vcf(DATA_DIR + "freebayes.vcf");
    ofstream csv_out(DATA_DIR + "freebayes_mutations.csv");

    if (!vcf.is_open() || !csv_out.is_open()) {
        cerr << "Error opening files!" << endl;
        return 1;
    }

    // Write CSV header
    csv_out << "Position,Type,REF,ALT\n";

    string line;
    while (getline(vcf, line)) {
        // Skip empty lines and VCF metadata (lines starting with '#')
        if (line.empty() || line[0] == '#') continue;

        stringstream ss(line);
        string column;
        string columns[5];
        int i = 0;

        // Extract the first 5 tab-separated columns from the VCF line
        while (getline(ss, column, '\t') && i < 5) {
            columns[i++] = column;
        }

        if (i < 5) continue; 

        string position = columns[1];
        string ref = columns[3];
        string alt = columns[4];
        string mut_type = tip_mutacije(ref, alt);

        // Write mutation data to CSV
        csv_out << position << "," << mut_type << "," << ref << "," << alt << "\n";
    }

    // Close files
    vcf.close();
    csv_out.close();

    return 0;
}
