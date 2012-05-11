#include <algorithm>
#include <iterator>
#include <set>
#include <vector>

// TODO make this #include "gff/Feature";
#include "Feature.h"
#include "JunctionIndex.h"

using GFF::Feature;

// TODO fix GFF lib to make this const Feature& f
void JunctionIndex::index(Feature& f)
{
    junctions.insert(f);
}

void JunctionIndex::overlappingFeature(Feature& target,
                                       std::vector<Feature>& overlaps)
{
    Feature lower;
    lower.seqid = target.seqid;
    lower.start = target.start;
    lower.end = target.start;

    Feature upper;
    upper.seqid = target.seqid;
    upper.start = target.end;
    upper.end = target.end;

    typedef std::set<Feature>::iterator iter;
    iter it = junctions.lower_bound(lower);
    iter end = junctions.upper_bound(upper);
    for (; it != end; ++it)
    {
        if ((*it).end < target.end) overlaps.push_back(*it);
    }
}

bool JunctionIndex::contains(Feature& query)
{
    return junctions.find(query) != junctions.end();
}

void JunctionIndex::unique(std::vector<Feature>& juncs)
{
    std::copy(junctions.begin(), junctions.end(), std::back_inserter(juncs));
}
