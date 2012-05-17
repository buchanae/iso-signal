#ifndef _ISOSIGNAL_COVERAGE_H
#define _ISOSIGNAL_COVERAGE_H

#include <map>
#include <string>
#include <vector>

#include "api/BamAux.h"
#include "Feature.h"

#include "Alignment.h"

using GFF::Feature;
using std::vector;
using std::string;

class Coverage
{
    public:
        std::map<std::string, vector<int> > coverages;

        // be careful with get/set/increment, they will throw out of range errors
        // if you haven't initialized the references
        int get(string ref_name, int pos);
        void set(string ref_name, int pos, int value);
        void increment(string ref_name, int pos, int value = 1);

        void setMinReferenceLength(string name, int length);

        void add(Alignment& alignment);
        void add(string ref_name, int start, int length);
};

void loadCoverage(std::istream& coverage_stream, Coverage& coverage);

void formatGMBCoverage(Coverage& coverage, std::ostream& coverage_stream);

void formatGMBCoverage(Coverage& coverage, std::string& output);

#endif
