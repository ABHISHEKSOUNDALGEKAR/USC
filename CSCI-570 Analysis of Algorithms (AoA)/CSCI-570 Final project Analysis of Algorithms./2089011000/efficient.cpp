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

class SequenceProcessor {
public:
    SequenceProcessor();

    std::pair<std::string,std::string> readData(const std::string &path) const;

    std::tuple<double,long,long,std::string,std::string> refinedAlignment(
        const std::string &s1,
        const std::string &s2
    ) const;

    void saveResults(
        const std::string &outPath,
        long cost,
        const std::string &aX,
        const std::string &aY,
        double runtime,
        long memKB
    ) const;

private:
    int gapPenalty;
    std::unordered_map<char,std::unordered_map<char,long>> matchScores;

    std::string extendSequence(const std::string &seq, const std::vector<int> &ins) const;

    std::vector<long> alignmentStrategy(
        const std::string &s1,
        const std::string &s2
    ) const;

    std::vector<std::string> standardAlignment(
        const std::string &s1,
        const std::string &s2
    ) const;

    std::vector<std::string> divideAndConquer(
        const std::string &s1,
        const std::string &s2
    ) const;
};

int main(int argc, char *argv[]) 
{
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " input.txt output.txt\n";
        return 1;
    }
    SequenceProcessor proc;
    auto seqs = proc.readData(argv[1]);
    double runtime;
    long memKB, cost;
    std::string ax, ay;
    std::tie(runtime, memKB, cost, ax, ay)
        = proc.refinedAlignment(seqs.first, seqs.second);
    proc.saveResults(argv[2], cost, ax, ay, runtime, memKB);
    return 0;
}

// ─────────────────────────────────────────────────────────────────────────────
// Definitions of SequenceProcessor members *after* main()
// ─────────────────────────────────────────────────────────────────────────────

SequenceProcessor::SequenceProcessor()
  : gapPenalty(30)
{
    matchScores = {
        {'A', {{'A',0},{'C',110},{'G',48},{'T',94}}},
        {'C', {{'A',110},{'C',0},{'G',118},{'T',48}}},
        {'G', {{'A',48},{'C',118},{'G',0},{'T',110}}},
        {'T', {{'A',94},{'C',48},{'G',110},{'T',0}}}
    };
}

std::pair<std::string,std::string>
SequenceProcessor::readData(const std::string &path) const {
    std::ifstream infile(path);
    if (!infile) {
        std::cerr << "Error opening input file\n";
        std::exit(1);
    }
    std::string line, seqX, seqY;
    std::getline(infile, line);
    seqX = trim(line);

    std::vector<int> ptsX;
    while (std::getline(infile, line)) {
        std::string t = trim(line);
        if (t.empty() || !std::all_of(t.begin(), t.end(), ::isdigit)) {
            seqY = t;
            break;
        }
        ptsX.push_back(std::stoi(t));
    }

    std::vector<int> ptsY;
    while (std::getline(infile, line)) {
        std::string t = trim(line);
        if (t.empty()) continue;
        ptsY.push_back(std::stoi(t));
    }

    auto eX = extendSequence(seqX, ptsX);
    auto eY = extendSequence(seqY, ptsY);
    std::transform(eX.begin(), eX.end(), eX.begin(), ::toupper);
    std::transform(eY.begin(), eY.end(), eY.begin(), ::toupper);
    return {eX, eY};
}

std::tuple<double,long,long,std::string,std::string>
SequenceProcessor::refinedAlignment(
    const std::string &s1,
    const std::string &s2
) const {
    auto t0 = std::chrono::high_resolution_clock::now();
    auto result = divideAndConquer(s1, s2);
    auto t1 = std::chrono::high_resolution_clock::now();
    double runtime = std::chrono::duration<double,std::milli>(t1 - t0).count();
    long memKB = getMemoryUsageKB();
    long cost = std::stol(result[0]);
    return {runtime, memKB, cost, result[1], result[2]};
}

void SequenceProcessor::saveResults(
    const std::string &outPath,
    long cost,
    const std::string &aX,
    const std::string &aY,
    double runtime,
    long memKB
) const {
    std::ofstream outfile(outPath);
    if (!outfile) {
        std::cerr << "Error opening output file\n";
        std::exit(1);
    }
    outfile << cost << "\n";
    outfile << aX   << "\n";
    outfile << aY   << "\n";
    outfile << std::fixed << std::setprecision(6) << runtime << "\n";
    outfile << memKB;
}

std::string SequenceProcessor::extendSequence(
    const std::string &seq,
    const std::vector<int> &ins
) const {
    std::string res = seq;
    for (int p : ins) {
        res = res.substr(0,p+1) + res + res.substr(p+1);
    }
    return res;
}

std::vector<long> SequenceProcessor::alignmentStrategy(
    const std::string &s1,
    const std::string &s2
) const {
    int n = s1.size()+1, m = s2.size()+1;
    std::vector<long> prev(m);
    for (int j = 0; j < m; ++j) prev[j] = j*gapPenalty;
    for (int i = 1; i < n; ++i) {
        std::vector<long> cur(m);
        cur[0] = i*gapPenalty;
        for (int j = 1; j < m; ++j) {
            long match = prev[j-1] + matchScores.at(s1[i-1]).at(s2[j-1]);
            long del   = prev[j]   + gapPenalty;
            long ins   = cur[j-1]  + gapPenalty;
            cur[j] = std::min({match, del, ins});
        }
        prev.swap(cur);
    }
    return prev;
}

std::vector<std::string> SequenceProcessor::standardAlignment(
    const std::string &s1,
    const std::string &s2
) const {
    int n = s1.size()+1, m = s2.size()+1;
    std::vector<std::vector<long>> mat(n, std::vector<long>(m));
    for (int i=0; i<n; ++i) mat[i][0] = i*gapPenalty;
    for (int j=1; j<m; ++j) mat[0][j] = j*gapPenalty;
    for (int i=1; i<n; ++i) {
        for (int j=1; j<m; ++j) {
            long match = mat[i-1][j-1] + matchScores.at(s1[i-1]).at(s2[j-1]);
            long del   = mat[i-1][j]   + gapPenalty;
            long ins   = mat[i][j-1]   + gapPenalty;
            mat[i][j] = std::min({match, del, ins});
        }
    }
    int i=n-1, j=m-1;
    std::string ax, ay;
    while (i>0 && j>0) {
        if (mat[i-1][j-1] + matchScores.at(s1[i-1]).at(s2[j-1]) == mat[i][j]) {
            ax.push_back(s1[i-1]);
            ay.push_back(s2[j-1]);
            --i; --j;
        } else if (mat[i][j-1] + gapPenalty == mat[i][j]) {
            ax.push_back('_');
            ay.push_back(s2[j-1]);
            --j;
        } else {
            ax.push_back(s1[i-1]);
            ay.push_back('_');
            --i;
        }
    }
    while (i>0) { ax.push_back(s1[i-1]); ay.push_back('_'); --i; }
    while (j>0) { ax.push_back('_');    ay.push_back(s2[j-1]); --j; }
    std::reverse(ax.begin(), ax.end());
    std::reverse(ay.begin(), ay.end());
    return { std::to_string(mat[n-1][m-1]), ax, ay };
}

std::vector<std::string> SequenceProcessor::divideAndConquer(
    const std::string &s1,
    const std::string &s2
) const {
    int n = s1.size(), m = s2.size();
    if (n < 2 || m < 2) return standardAlignment(s1, s2);

    auto leftCost = alignmentStrategy(s1.substr(0, n/2), s2);
    std::string rs1 = s1.substr(n/2); std::reverse(rs1.begin(), rs1.end());
    std::string rs2 = s2;            std::reverse(rs2.begin(), rs2.end());
    auto rightCost = alignmentStrategy(rs1, rs2);
    std::reverse(rightCost.begin(), rightCost.end());

    int cut = 0; long best = leftCost[0] + rightCost[0];
    for (int i = 1; i <= m; ++i) {
        long c = leftCost[i] + rightCost[i];
        if (c < best) { best = c; cut = i; }
    }

    auto L = divideAndConquer(s1.substr(0, n/2), s2.substr(0, cut));
    auto R = divideAndConquer(s1.substr(n/2),    s2.substr(cut));
    long cL = std::stol(L[0]), cR = std::stol(R[0]);
    return { std::to_string(cL + cR), L[1] + R[1], L[2] + R[2] };
}
