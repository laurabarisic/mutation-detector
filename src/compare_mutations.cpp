#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>

using namespace std;

string DATA_DIR = "../data/";

// Učitavanje custom formata CSV: tip, pozicija, baza
unordered_set<string> loadCustomCSV(const string &filename) {
    unordered_set<string> variants;
    ifstream file(filename);
    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        string type, pos, base;
        getline(ss, type, ',');
        getline(ss, pos, ',');
        getline(ss, base, ',');
        if (!type.empty() && !pos.empty())
            variants.insert(type + "|" + pos);
    }
    return variants;
}

// Učitavanje freebayes formata: pozicija, tip, ref, alt
unordered_set<string> loadFreebayesCSV(const string &filename) {
    unordered_set<string> variants;
    ifstream file(filename);
    string line;
    getline(file, line); // preskoči header
    while (getline(file, line)) {
        stringstream ss(line);
        string pos, type, ref, alt;
        getline(ss, pos, ',');
        getline(ss, type, ',');
        getline(ss, ref, ',');
        getline(ss, alt, ',');
        if (!type.empty() && !pos.empty())
            variants.insert(type + "|" + pos);
    }
    return variants;
}

unordered_set<string> getIntersection(const unordered_set<string> &a,
                                      const unordered_set<string> &b) {
    unordered_set<string> result;
    for (const auto &v : a)
        if (b.count(v))
            result.insert(v);
    return result;
}

unordered_set<string> getIntersectionThree(const unordered_set<string> &a,
                                           const unordered_set<string> &b,
                                           const unordered_set<string> &c) {
    unordered_set<string> result;
    for (const auto &v : a)
        if (b.count(v) && c.count(v))
            result.insert(v);
    return result;
}

unordered_set<string> getUnique(const unordered_set<string> &a,
                                const unordered_set<string> &b,
                                const unordered_set<string> &c) {
    unordered_set<string> result;
    for (const auto &v : a)
        if (!b.count(v) && !c.count(v))
            result.insert(v);
    return result;
}

unordered_set<string> getUnion(const unordered_set<string> &a,
                               const unordered_set<string> &b,
                               const unordered_set<string> &c) {
    unordered_set<string> result = a;
    for (const auto &v : b)
        result.insert(v);
    for (const auto &v : c)
        result.insert(v);
    return result;
}

int main() {
    auto lambda_mutated = loadCustomCSV(DATA_DIR + "lambda_mutated.csv");
    auto lambda_mutations = loadCustomCSV(DATA_DIR + "lambda_mutations.csv");
    auto lambda_freebayes = loadFreebayesCSV(DATA_DIR + "lambda_freebayes_mutations_filtered.csv");

    auto ecoli_mutated = loadCustomCSV(DATA_DIR + "ecoli_mutated.csv");
    auto ecoli_mutations = loadCustomCSV(DATA_DIR + "ecoli_mutations.csv");
    auto ecoli_freebayes = loadFreebayesCSV(DATA_DIR + "ecoli_freebayes_filtered_mutations.csv");

    cout << "LAMBDA PODACI\n";
    cout << "Ukupan broj varijanti:\n";
    cout << "lambda_mutated: " << lambda_mutated.size() << endl;
    cout << "lambda_mutations: " << lambda_mutations.size() << endl;
    cout << "lambda_freebayes: " << lambda_freebayes.size() << endl << endl;

    cout << "Zajednicke varijante:\n";
    cout << "lambda_mutated & lambda_mutations: "
         << getIntersection(lambda_mutated, lambda_mutations).size() << endl;
    cout << "lambda_mutated & lambda_freebayes: "
         << getIntersection(lambda_mutated, lambda_freebayes).size() << endl;
    cout << "lambda_mutations & lambda_freebayes: "
         << getIntersection(lambda_mutations, lambda_freebayes).size() << endl
         << endl;

    cout << "Jedinstvene varijante:\n";
    cout << "Samo u lambda_mutated: "
         << getUnique(lambda_mutated, lambda_mutations, lambda_freebayes).size()
         << endl;
    cout << "Samo u lambda_mutations: "
         << getUnique(lambda_mutations, lambda_mutated, lambda_freebayes).size()
         << endl;
    cout << "Samo u lambda_freebayes: "
         << getUnique(lambda_freebayes, lambda_mutated, lambda_mutations).size()
         << endl;

    cout << "U svim lambda datotekama: "
         << getIntersectionThree(lambda_freebayes, lambda_mutated,
                                 lambda_mutations)
                .size()
         << endl;
    cout << "Ukupan broj mutacija svih lambda datoteka: "
         << getUnion(lambda_freebayes, lambda_mutated, lambda_mutations).size()
         << endl;

    cout << "\nECOLI PODACI\n";
    cout << "Ukupan broj varijanti:\n";
    cout << "ecoli_mutated: " << ecoli_mutated.size() << endl;
    cout << "ecoli_mutations: " << ecoli_mutations.size() << endl;
    cout << "ecoli_freebayes: " << ecoli_freebayes.size() << endl << endl;

    cout << "Zajednicke varijante:\n";
    cout << "ecoli_mutated & ecoli_mutations: "
         << getIntersection(ecoli_mutated, ecoli_mutations).size() << endl;
    cout << "ecoli_mutated & ecoli_freebayes: "
         << getIntersection(ecoli_mutated, ecoli_freebayes).size() << endl;
    cout << "ecoli_mutations & ecoli_freebayes: "
         << getIntersection(ecoli_mutations, ecoli_freebayes).size() << endl
         << endl;

    cout << "Jedinstvene varijante:\n";
    cout << "Samo u ecoli_mutated: "
         << getUnique(ecoli_mutated, ecoli_mutations, ecoli_freebayes).size()
         << endl;
    cout << "Samo u ecoli_mutations: "
         << getUnique(ecoli_mutations, ecoli_mutated, ecoli_freebayes).size()
         << endl;
    cout << "Samo u ecoli_freebayes: "
         << getUnique(ecoli_freebayes, ecoli_mutated, ecoli_mutations).size()
         << endl;
    cout << "U svim ecoli datotekama: "
         << getIntersectionThree(ecoli_freebayes, ecoli_mutated,
                                 ecoli_mutations)
                .size()
         << endl;
    cout << "Ukupan broj mutacija svih ecoli datoteka: "
         << getUnion(ecoli_freebayes, ecoli_mutated, ecoli_mutations).size()
         << endl;

    return 0;
}
