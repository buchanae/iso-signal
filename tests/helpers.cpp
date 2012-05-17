#include <iostream>
#include <vector>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

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
