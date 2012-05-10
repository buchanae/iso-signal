#ifndef _ISOSIGNAL_COVERAGE_H
#define _ISOSIGNAL_COVERAGE_H

#include <string>
#include <vector>

#include "api/BamAux.h"
#include "Feature.h"

#include "Alignment.h"

using GFF::Feature;
using std::vector;
using std::string;

void gen_coverage(vector<int>&, vector<Alignment>&, Feature&);

void add_coverage(Feature&, vector<int>&, vector<BamTools::CigarOp>&, int);

void comp_coverage(string&, vector<int>&, vector<int>&);

void gen_signal(vector<int>&, vector<double>&, int, int);

#endif
