#include <iostream>
#include <sstream>
#include <string>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "Feature.h"
#include "StackReader.h"

using GFF::Feature;
using std::string;

TEST(StackReaderTest, getNextFeature)
{
    std::stringstream data;
    data << "@Chr\tAT-CG\t10\t20\t11" << std::endl;
    data << "ATCG\t1" << std::endl;
    data << "ATCG\t1" << std::endl;
    data << "@Chr2\tAT-CG\t30\t40\t11" << std::endl;

    Feature f;

    EXPECT_TRUE(StackReader::getNextFeature(data, f));
    EXPECT_EQ("Chr", f.seqid);
    EXPECT_EQ("splice_junction", f.type);
    EXPECT_EQ(10, f.start);
    EXPECT_EQ(20, f.end);

    EXPECT_TRUE(StackReader::getNextFeature(data, f));
    EXPECT_EQ("Chr2", f.seqid);
    EXPECT_EQ("splice_junction", f.type);
    EXPECT_EQ(30, f.start);
    EXPECT_EQ(40, f.end);

    EXPECT_FALSE(StackReader::getNextFeature(data, f));
}
