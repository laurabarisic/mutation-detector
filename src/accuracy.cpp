#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

string DATA_PATH = "../data/";

// Represents a single mutation
struct Mutation {
    string type; // "X", "I", "D"
    int position;
    string newValue;
};

// Aliases for better readability
using MutationList = vector<Mutation>;
using MutationMap =
    unordered_map<int, pair<string, string>>; // position -> (type, newValue)

// Parse a CSV line into a Mutation object
Mutation parseLine(const string &line) {
    stringstream ss(line);
    string type, posStr, val;
    getline(ss, type, ',');
    getline(ss, posStr, ',');
    getline(ss, val, ',');
    return {type, stoi(posStr), val};
}

// Load mutations from a CSV file (ignores header)
MutationList loadFromFile(const string &filename) {
    MutationList list;
    ifstream file(filename);
    if (!file) {
        cerr << "Error: cannot open file " << filename << endl;
        return list;
    }

    string line;
    getline(file, line); // Skip header
    while (getline(file, line)) {
        list.push_back(parseLine(line));
    }

    return list;
}

// Create a hash map for fast mutation lookup by position
MutationMap indexByPosition(const MutationList &list) {
    MutationMap map;
    for (const auto &m : list) {
        map[m.position] = {m.type, m.newValue};
    }
    return map;
}

// Compare predicted mutations against reference mutations
double evaluate(const MutationList &predicted, const MutationList &reference) {

    const int scorePos = 1;
    const int scoreType = 1;
    const int scoreValue = 1;

    int total = 0;
    int earned = 0;

    MutationMap predictedMap = indexByPosition(predicted);

    for (const auto &ref : reference) {
        total += scorePos + scoreType + scoreValue;

        auto it = predictedMap.find(ref.position);
        if (it != predictedMap.end()) {
            const auto &[pType, pVal] = it->second;
            earned += scorePos;

            if (pType == ref.type) {
                earned += scoreType;
                if (pVal == ref.newValue) {
                    earned += scoreValue;
                }
            }
        }
    }

    double accuracy = (total > 0) ? (earned * 100.0 / total) : 0.0;

    return accuracy;
}

int main() {
    const string predFile = DATA_PATH + "mutations.csv";
    const string refFile = DATA_PATH + "lambda_mutated.csv";

    MutationList predicted = loadFromFile(predFile);
    MutationList reference = loadFromFile(refFile);

    double accuracy = evaluate(predicted, reference);

    cout << fixed << "Mutation detection accuracy: " << accuracy << "%" << endl;
    return 0;
}
