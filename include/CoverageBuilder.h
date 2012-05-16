#ifndef ISOSIGNAL_COVERAGEBUILDER_H
#define ISOSIGNAL_COVERAGEBUILDER_H

#include <istream>
#include <map>

#include "BamReader.h"
#include "Coverage.h"
#include "JunctionFilter.h"

using std::istream;

class CoverageBuilder
{
    JunctionFilter filter;
    std::map<std::string, Coverage*> coverages;
    
    public:
        void addJunctionsFromStack(istream& stack_stream);
        void addJunctionsFromGFF(istream& gff_stream);
        void addCoverageFromBam(BamReader& reader);
        std::string coverageString(void);
};
#endif
