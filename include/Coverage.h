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

        void add(Alignment& alignment);
        void add(string ref_name, int start, int length);

        void setMinReferenceLength(string name, int length);
};

#endif
