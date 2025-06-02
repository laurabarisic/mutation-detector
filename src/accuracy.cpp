#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <iomanip>

using namespace std;

string DATA_DIR = "../data/";

struct GenomicMutation {
    string mutationType;
    int genomePosition;
    string newBase;
};

using MutationList = vector<GenomicMutation>;
using MutationLookup = unordered_map<int, pair<string, string>>;

GenomicMutation parseMutationRow(const string& csvLine) {
    stringstream lineStream(csvLine);
    string type, positionStr, base;

    getline(lineStream, type, ',');
    getline(lineStream, positionStr, ',');
    getline(lineStream, base, ',');

    return {type, stoi(positionStr), base};
}

MutationList loadMutationsFromFile(const string& filename) {
    MutationList mutations;
    ifstream inputFile(filename);

    if (!inputFile) {
        cerr << "Greška pri otvaranju datoteke: " << filename << endl;
        return mutations;
    }

    string line;
    getline(inputFile, line); 

    while (getline(inputFile, line)) {
        mutations.push_back(parseMutationRow(line));
    }

    return mutations;
}

MutationLookup buildMutationIndex(const MutationList& mutations) {
    MutationLookup index;
    for (const auto& mutation : mutations) {
        index[mutation.genomePosition] = {mutation.mutationType, mutation.newBase};
    }
    return index;
}

double evaluateMutationAccuracy(const MutationList& predicted, const MutationList& reference, const string& logPath) {
    ofstream logFile(logPath);
    if (!logFile) {
        cerr << "Ne mogu otvoriti log datoteku: " << logPath << endl;
        return 0.0;
    }

    const int positionScore = 1;
    const int typeScore = 1;
    const int baseScore = 1;

    int earnedPoints = 0;
    int totalPossiblePoints = 0;

    MutationLookup predictionIndex = buildMutationIndex(predicted);

    logFile << "Rezultati evaluacije\n====================\n";

    for (const auto& refMutation : reference) {
        totalPossiblePoints += positionScore + typeScore + baseScore;
        logFile << "Provjera referentne mutacije: "
                << refMutation.mutationType << ","
                << refMutation.genomePosition << ","
                << refMutation.newBase << "\n";

        auto found = predictionIndex.find(refMutation.genomePosition);
        if (found != predictionIndex.end()) {
            const auto& [predType, predBase] = found->second;
            logFile << "  Pronađena predikcija: "
                    << predType << "," << refMutation.genomePosition << "," << predBase << "\n";

            earnedPoints += positionScore;

            if (predType == refMutation.mutationType) {
                earnedPoints += typeScore;
                if (predBase == refMutation.newBase) {
                    earnedPoints += baseScore;
                    logFile << "    Sve odgovara: +" << (positionScore + typeScore + baseScore) << " bodova\n";
                } else {
                    logFile << "    Tip odgovara, ali baza ne: +" << (positionScore + typeScore) << " bodova\n";
                }
            } else {
                logFile << "    Pozicija odgovara, ali tip ne: +" << positionScore << " bodova\n";
            }
        } else {
            logFile << "  Nema predikcije za tu poziciju\n";
        }
    }

    double accuracyPercentage = (totalPossiblePoints > 0)
        ? (static_cast<double>(earnedPoints) / totalPossiblePoints) * 100.0
        : 0.0;

    logFile << "====================\n";
    logFile << "Osvojeni bodovi: " << earnedPoints << "\n";
    logFile << "Maksimalni bodovi: " << totalPossiblePoints << "\n";
    logFile << fixed << setprecision(2)
            << "Točnost: " << accuracyPercentage << "%\n";

    return accuracyPercentage;
}

int main() {
    const string userFile = DATA_DIR + "mutations.csv";
    const string referenceFile = DATA_DIR + "lambda_mutated.csv";
    const string evaluationLog = DATA_DIR + "evaluation_log.txt";

    MutationList userMutations = loadMutationsFromFile(userFile);
    MutationList referenceMutations = loadMutationsFromFile(referenceFile);

    double finalAccuracy = evaluateMutationAccuracy(userMutations, referenceMutations, evaluationLog);

    cout << fixed << setprecision(2)
         << "Točnost detekcije mutacija: " << finalAccuracy << "%" << endl;

    return 0;
}
