#include <vector>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

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
