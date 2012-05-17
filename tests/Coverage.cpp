#include <sstream>
#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "Alignment.h"
#include "Coverage.h"

using GFF::Feature;
using testing::DoubleEq;
using testing::ElementsAre;

/*
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
*/

TEST(CoverageTest, Coverage_setMinReferenceLength)
{
    Coverage c;

    EXPECT_EQ(0, c.coverages.size());
    c.setMinReferenceLength("foo", 10);

    EXPECT_EQ(1, c.coverages.size());
    EXPECT_THAT(c.coverages.find("foo")->second,
                ElementsAre(0, 0, 0, 0, 0, 0, 0, 0, 0, 0));

    c.setMinReferenceLength("foo", 5);
    EXPECT_THAT(c.coverages.find("foo")->second,
                ElementsAre(0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
}

TEST(CoverageTest, Coverage_add)
{
    Coverage c;
    c.add("foo", 2, 3);

    EXPECT_THAT(c.coverages.find("foo")->second,
                ElementsAre(0, 1, 1, 1));

    c.add("foo", 2, 2);

    EXPECT_THAT(c.coverages.find("foo")->second,
                ElementsAre(0, 2, 2, 1));

    c.add("foo", 6, 2);

    EXPECT_THAT(c.coverages.find("foo")->second,
                ElementsAre(0, 2, 2, 1, 0, 1, 1));

    c.setMinReferenceLength("foo", 10);

    EXPECT_THAT(c.coverages.find("foo")->second,
                ElementsAre(0, 2, 2, 1, 0, 1, 1, 0, 0, 0));
}

TEST(CoverageTest, Coverage_add_alignment)
{
    Alignment a;
    a.RefName = "foo";
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

    Coverage c;
    c.add(a);

    EXPECT_THAT(c.coverages.find("foo")->second, ElementsAre(0, 0, 1, 1, 0, 0, 0, 1, 1));
}

TEST(HelpersTest, loadCoverage)
{
    Coverage c;

    std::stringstream coverage_str("bar\t6\n1\n1\n1\n0\n0\n0\nfoo\t5\n0\n0\n1\n1\n0\n");
    loadCoverage(coverage_str, c);

    EXPECT_THAT(c.coverages.find("bar")->second, ElementsAre(1, 1, 1, 0, 0, 0));
    EXPECT_THAT(c.coverages.find("foo")->second, ElementsAre(0, 0, 1, 1, 0));
}

TEST(HelpersTest, formatGMBCoverage)
{
    // Note that references are output in sorted order on reference name

    Coverage c;
    c.setMinReferenceLength("foo", 5);
    c.setMinReferenceLength("bar", 6);
    c.add("foo", 3, 2);
    c.add("bar", 1, 3);

    std::string expected = "bar\t6\n1\n1\n1\n0\n0\n0\nfoo\t5\n0\n0\n1\n1\n0\n";

    std::stringstream out;
    formatGMBCoverage(c, out);

    EXPECT_EQ(expected, out.str());

    std::string out_string;
    formatGMBCoverage(c, out_string);

    EXPECT_EQ(expected, out_string);
}
