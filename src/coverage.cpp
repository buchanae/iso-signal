#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "api/BamAux.h"

#include "Alignment.h"
#include "coverage.h"

using GFF::Feature;
using std::endl;
using std::vector;
using std::string;

void diff_signals(vector<double>& a, vector<double>& b, vector<double>& diff)
{
    // TODO some check for whether a.size() == b.size()

    for (int i = 0; i < a.size(); i++)
        diff.push_back(fabs(a.at(i) - b.at(i)));
}

void comp_coverage(string& id, vector<int>& coverage1, vector<int>& coverage2) {
    vector<double> signal1, signal2, diff;
    int WIN = 20, SLIDE = 10, total_above = 0;
    //int WIN = 1, SLIDE = 1, total_above = 0;

    gen_signal(coverage1, signal1, WIN, SLIDE);

    gen_signal(coverage2, signal2, WIN, SLIDE);

    for (int i = 0; i < signal1.size(); i++) {
        double sig_diff = fabs(signal1.at(i) - signal2.at(i));
        diff.push_back(sig_diff);
        if (sig_diff > 0.1) total_above++;
    }

    string file = "plots/" + id + ".plot.R";
    //string file = "test.plot.R";

    std::ofstream out(file.c_str(), std::ios::out | std::ios::trunc);

    out << "x1 <- c(";
    for (int i = 0; i < signal1.size(); i++) {
        out << signal1.at(i);
        if (i < signal1.size() - 1) out << ",";
    }
    out << ")" << endl << endl;

    out << "x2 <- c(";
    for (int i = 0; i < signal2.size(); i++) {
        out << signal2.at(i);
        if (i < signal2.size() - 1) out << ",";
    }
    out << ")" << endl << endl;

    out << "x3 <- c(";
    for (int i = 0; i < diff.size(); i++) {
        out << diff.at(i);
        if (i < diff.size() - 1) out << ",";
    }
    out << ")" << endl << endl;

    out << "plot(x1, type = 'l')" << endl;
    out << "lines(x2)" << endl;
    out << "lines(x3)" << endl;

    out.close();
}

void gen_signal (vector<int>& coverage, vector<double>& signal, int window, int slide) {
    int loops = static_cast<int>((double)(coverage.size() - window + 1) / (double)slide + 0.5);
    double total = 0;

    for (int loop = 1; loop <= loops; loop++) {
        int start = (loop  - 1) * slide + 1;
        int end = start + window - 1;
        for (int index = start; index <= end; index++) total += (double)coverage.at(index - 1);
    }

    for (int loop = 1; loop <= loops; loop++) {
        int start = (loop - 1) * slide + 1;
        int end = start + window - 1;
        double nt = 0;
        for (int index = start; index <= end; index++) nt += (double)coverage.at(index - 1);
        signal.push_back(nt / total * 100); // Show expression as a percentage of total expression
        //signal.push_back(nt); // Show expression as a raw total value
    }
}

void gen_coverage(vector<int>& coverage, vector<Alignment>& alignments, Feature& feature) {
    for (int i = 0; i < feature.end - feature.start + 1; i++)
        coverage.push_back(0);

    int total_alignments = alignments.size();
    for (int i = 0; i < total_alignments; i++) {
        Alignment al = alignments.at(i);

        add_coverage(feature, coverage, al.CigarData, al.Position + 1);
    }
}

void add_coverage(Feature& feature, vector<int>& coverage, vector<BamTools::CigarOp>& CigarData, int pos) {
    int start, end;

    start = pos;
    end = start + CigarData.at(0).Length - 1;
    for (int index = start; index <= end; index++)
        if (index >= feature.start && index <= feature.end) coverage.at(index - feature.start)++;

    if (CigarData.size() == 3) {
        int offset = CigarData.at(0).Length + CigarData.at(1).Length;
        start = pos + offset;
        end = start + CigarData.at(2).Length - 1;
        for (int index = start; index <= end; index++)
            if (index >= feature.start && index <= feature.end) coverage.at(index - feature.start)++;
    }
}
