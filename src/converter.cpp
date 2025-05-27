#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

string DATA_DIR = "../data/";

string tip_mutacije(const string& ref, const string& alt_raw) {
    string alt = alt_raw;
    size_t comma_pos = alt.find(',');
    if (comma_pos != string::npos) {
        alt = alt.substr(0, comma_pos);
    }

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
    ifstream vcf(DATA_DIR + "freebayes.vcf");
    ofstream csv_out(DATA_DIR + "freebayes_mutations.csv");

    if (!vcf.is_open() || !csv_out.is_open()) {
        cerr << "Greška pri otvaranju datoteka!" << endl;
        return 1;
    }

    csv_out << "Position,Type,REF,ALT\n";

    string line;
    while (getline(vcf, line)) {
        if (line.empty() || line[0] == '#') continue;

        stringstream ss(line);
        string kolona;
        string kolone[5];
        int i = 0;

        while (getline(ss, kolona, '\t') && i < 5) {
            kolone[i++] = kolona;
        }

        if (i < 5) continue;

        string poz = kolone[1];
        string ref = kolone[3];
        string alt = kolone[4];
        string mut_tip = tip_mutacije(ref, alt);

        csv_out << poz << "," << mut_tip << "," << ref << "," << alt << "\n";
    }

    vcf.close();
    csv_out.close();

    return 0;
}
