#include <string>
#include <vector>

#include "api/BamAux.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "Alignment.h"
#include "Feature.h"
#include "JunctionFilter.h"

using GFF::Feature;
using testing::ElementsAre;

TEST(JunctionFilterTest, regular_alignment)
{
    // Regular alignments (i.e. not splice junction alignments) will always pass
    // the junction filter.  The junction filter only returns false for splice junction
    // alignments that aren't supported by splice junction data.
    JunctionFilter filter;

    Alignment al;
    al.RefName = "foo";

    CigarOp op;
    op.Type = 'M';
    op.Length = 51;
    al.CigarData.push_back(op);

    EXPECT_TRUE(filter(al));
}

TEST(JunctionFilterTest, supported_alignment_passes)
{
    JunctionFilter filter;
    Feature junction;
    junction.seqid = "foo";

    junction.start = 10;
    junction.end = 20;
    filter.junction_index.add(junction);

    Alignment al;
    al.RefName = "foo";
    CigarOp op;

    op.Type = 'M';
    op.Length = 10;
    al.CigarData.push_back(op);

    op.Type = 'N';
    op.Length = 9;
    al.CigarData.push_back(op);

    op.Type = 'M';
    op.Length = 10;
    al.CigarData.push_back(op);

    al.position(1);
    EXPECT_TRUE(filter(al));
}

TEST(JunctionFilterTest, unsupported_alignment_fails)
{
    JunctionFilter filter;
    Feature junction;
    junction.seqid = "foo";

    junction.start = 10;
    junction.end = 20;
    filter.junction_index.add(junction);

    Alignment al;
    al.RefName = "foo";
    CigarOp op;

    op.Type = 'M';
    op.Length = 10;
    al.CigarData.push_back(op);

    op.Type = 'N';
    op.Length = 9;
    al.CigarData.push_back(op);

    op.Type = 'M';
    op.Length = 10;
    al.CigarData.push_back(op);

    al.position(30);
    EXPECT_FALSE(filter(al));
}
