#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <chrono>
#include <sys/resource.h>
#include <algorithm>
#include <cctype>
#include <tuple>
#include <iomanip>

static std::string trim(const std::string &s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

static long getMemoryUsageKB() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss;
}

class DNAAligner 
{
public:
    DNAAligner();

    std::pair<std::string,std::string> fetchSequences(const std::string &path) const;

    std::tuple<double,long,std::vector<std::string>> alignSequences(
        const std::string &seq1,
        const std::string &seq2
    ) const;

    void documentResults(
        const std::string &outPath,
        const std::vector<std::string> &outcomes,
        double duration,
        long memKB
    ) const;

private:
    int indelPenalty;
    std::unordered_map<char,std::unordered_map<char,long>> matchScores;

    std::string amplifyDNA(const std::string &seq, const std::vector<int> &spots) const;
};

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " input.txt output.txt\n";
        return 1;
    }
    DNAAligner aligner;
    auto seqs = aligner.fetchSequences(argv[1]);

    double duration;
    long memKB;
    std::vector<std::string> outcomes;
    std::tie(duration, memKB, outcomes)
        = aligner.alignSequences(seqs.first, seqs.second);

    aligner.documentResults(argv[2], outcomes, duration, memKB);
    return 0;
}

DNAAligner::DNAAligner()
  : indelPenalty(30)
{
    matchScores = {
        {'A', {{'A',0},{'C',110},{'G',48},{'T',94}}},
        {'C', {{'A',110},{'C',0},{'G',118},{'T',48}}},
        {'G', {{'A',48},{'C',118},{'G',0},{'T',110}}},
        {'T', {{'A',94},{'C',48},{'G',110},{'T',0}}}
    };
}

std::pair<std::string,std::string>
DNAAligner::fetchSequences(const std::string &path) const 
{
    std::ifstream infile(path);
    if (!infile) {
        std::cerr << "Error opening input file\n";
        std::exit(1);
    }
    std::string line, initialX, initialY;
    std::getline(infile, line);
    initialX = trim(line);

    std::vector<int> spotsX;
    while (std::getline(infile, line)) {
        std::string t = trim(line);
        if (t.empty() || !std::all_of(t.begin(), t.end(), ::isdigit)) {
            initialY = t;
            break;
        }
        spotsX.push_back(std::stoi(t));
    }

    std::vector<int> spotsY;
    while (std::getline(infile, line)) {
        std::string t = trim(line);
        if (t.empty()) continue;
        spotsY.push_back(std::stoi(t));
    }

    auto seqX = amplifyDNA(initialX, spotsX);
    auto seqY = amplifyDNA(initialY, spotsY);
    std::transform(seqX.begin(), seqX.end(), seqX.begin(), ::toupper);
    std::transform(seqY.begin(), seqY.end(), seqY.begin(), ::toupper);
    return {seqX, seqY};
}

std::tuple<double,long,std::vector<std::string>>
DNAAligner::alignSequences(
    const std::string &seq1,
    const std::string &seq2
) const {
    auto t0 = std::chrono::high_resolution_clock::now();
    int n = seq1.size()+1, m = seq2.size()+1;
    std::vector<std::vector<long>> mat(n, std::vector<long>(m));
    for (int i = 0; i < n; ++i) mat[i][0] = i*indelPenalty;
    for (int j = 1; j < m; ++j) mat[0][j] = j*indelPenalty;

    for (int i = 1; i < n; ++i) {
        for (int j = 1; j < m; ++j) {
            long match = mat[i-1][j-1]
                + matchScores.at(seq1[i-1]).at(seq2[j-1]);
            long del   = mat[i-1][j]   + indelPenalty;
            long ins   = mat[i][j-1]   + indelPenalty;
            mat[i][j] = std::min({match, del, ins});
        }
    }

    int i = n-1, j = m-1;
    std::string alignX, alignY;
    while (i>0 && j>0) {
        if (mat[i-1][j-1]
            + matchScores.at(seq1[i-1]).at(seq2[j-1])
            == mat[i][j])
        {
            alignX.push_back(seq1[i-1]);
            alignY.push_back(seq2[j-1]);
            --i; --j;
        }
        else if (mat[i][j-1] + indelPenalty == mat[i][j]) {
            alignX.push_back('_');
            alignY.push_back(seq2[j-1]);
            --j;
        } else {
            alignX.push_back(seq1[i-1]);
            alignY.push_back('_');
            --i;
        }
    }
    while (i>0) {
        alignX.push_back(seq1[i-1]);
        alignY.push_back('_');
        --i;
    }
    while (j>0) {
        alignX.push_back('_');
        alignY.push_back(seq2[j-1]);
        --j;
    }
    std::reverse(alignX.begin(), alignX.end());
    std::reverse(alignY.begin(), alignY.end());

    auto t1 = std::chrono::high_resolution_clock::now();
    double duration =
      std::chrono::duration<double,std::milli>(t1 - t0).count();
    long memKB = getMemoryUsageKB();

    return {duration, memKB,
            { std::to_string(mat[n-1][m-1]), alignX, alignY }};
}

void DNAAligner::documentResults(
    const std::string &outPath,
    const std::vector<std::string> &outcomes,
    double duration,
    long memKB
) const {
    std::ofstream outfile(outPath);
    if (!outfile) {
        std::cerr << "Error opening output file\n";
        std::exit(1);
    }
    for (auto &line : outcomes) outfile << line << "\n";
    outfile << std::fixed << std::setprecision(6) << duration << "\n";
    outfile << memKB;
}

std::string DNAAligner::amplifyDNA(
    const std::string &seq,
    const std::vector<int> &spots
) const {
    std::string res = seq;
    for (int p : spots) {
        res = res.substr(0,p+1) + res + res.substr(p+1);
    }
    return res;
}
