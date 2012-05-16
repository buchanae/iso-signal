#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "Coverage.h"
#include "helpers.h"

using GFF::Feature;
using testing::ElementsAre;

TEST(HelpersTest, getGenesAndTranscriptFromGFF)
{
    std::ifstream gff_input_stream("dummies/features.gff");
    EXPECT_TRUE(gff_input_stream.is_open());

    std::vector<Feature> all;
    std::vector<Feature> genes;
    std::vector<Feature> transcripts;

    getGenesAndTranscriptsFromGFF(gff_input_stream, all, genes, transcripts);

    EXPECT_EQ(18, all.size());
    EXPECT_EQ(3, genes.size());
    EXPECT_EQ(9, transcripts.size());
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
