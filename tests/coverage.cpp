#include <vector>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "Alignment.h"
#include "coverage.h"

using GFF::Feature;
using testing::DoubleEq;
using testing::ElementsAre;

TEST(CoverageTest, diff_signals)
{
    vector<double> a, b, diff;

    a.push_back(5.5);
    b.push_back(4.3);

    a.push_back(7.5);
    b.push_back(4.3);

    diff_signals(a, b, diff);

    EXPECT_THAT(diff, ElementsAre(DoubleEq(1.2), DoubleEq(3.2)));
}

TEST(CoverageTest, Coverage_init)
{
    Coverage c(10);
    EXPECT_THAT(c.data, ElementsAre(0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
}

TEST(CoverageTest, Coverage_add_position)
{
    Coverage c(10);
    c.add(3, 7);

    EXPECT_THAT(c.data, ElementsAre(0, 0, 1, 1, 1, 1, 1, 0, 0, 0));
}

TEST(CoverageTest, Coverage_add_alignment)
{
    Alignment a;
    a.position(3);

    CigarOp op;

    a.CigarData.clear();  
    op.Type = 'M';
    op.Length = 2;
    a.CigarData.push_back(op);

    op.Type = 'N';
    op.Length = 3;
    a.CigarData.push_back(op);

    op.Type = 'M';
    op.Length = 2;
    a.CigarData.push_back(op);

    Coverage c(10);
    c.add(a);

    EXPECT_THAT(c.data, ElementsAre(0, 0, 1, 1, 0, 0, 0, 1, 1, 0));
}
