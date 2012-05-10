#include <vector>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "Feature.h"

#include "JunctionIndex.h"

using GFF::Feature;
using testing::ElementsAre;
using testing::WhenSorted;

TEST(JunctionIndexTest, overlappingFeature)
{
    JunctionIndex index;

    // a and b are valid junctions and will be returned from the query
    Feature a;
    a.seqid = "one";
    a.source = "A";
    a.start = 20;
    a.end = 40;
    index.add(a);

    Feature b;
    b.start = 60;
    b.end = 80;
    b.seqid = "one";
    b.source = "B";
    index.add(b);

    // The query range will fall in the middle of this junction,
    // so it won't be returned.
    // i.e. valid junctions must fall entirely within the query range
    Feature c;
    c.seqid = "one";
    c.source = "C";
    c.start = 100;
    c.end = 150;
    index.add(c);

    // duplicate of Feature a.
    // JunctionIndex only returns unique junctions,
    // so this won't be returned as a duplicate
    Feature d;
    d.seqid = "one";
    d.start = 20;
    d.end = 40;
    index.add(d);

    // will fall in a valid range, but has a different seqid than the query
    // so won't be returned as overlapping query feature
    Feature e;
    e.seqid = "two";
    e.start = 20;
    e.end = 40;
    index.add(e);

    Feature query;
    query.seqid = "one";
    query.start = 10;
    query.end = 110;

    std::vector<Feature> ret;
    index.overlappingFeature(query, ret);
    
    std::vector<std::string> sources;
    for (std::vector<Feature>::iterator it = ret.begin(); it != ret.end(); ++it)
    {
        sources.push_back((*it).source);
    }

    EXPECT_THAT(sources, WhenSorted(ElementsAre("A", "B")));

    ret.clear();
}
