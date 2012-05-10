#ifndef _ISOSIGNAL_JUNCTIONINDEX_H
#define _ISOSIGNAL_JUNCTIONINDEX_H

#include <set>
#include <vector>

#include "Feature.h"
#include "Index.h"

using GFF::Feature;

class JunctionIndex : public GFF::IndexBase
{
    struct JunctionComparison
    {
        bool operator() (const Feature& a, const Feature& b) const
        {
            return (a.seqid < b.seqid)
                   || (a.seqid == b.seqid && a.start < b.start)
                   || (a.seqid == b.seqid && a.start == b.start
                       && a.end < b.end);
        }
    };

    std::set<Feature, JunctionComparison> junctions;

    public:
        void overlappingFeature(Feature&, std::vector<Feature>&);
        bool contains(Feature&);

    private:
        void index(Feature&);
};

#endif
