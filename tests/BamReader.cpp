#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "Alignment.h"
#include "BamReader.h"

using GFF::Feature;
using testing::ElementsAre;
using testing::WhenSorted;

TEST(BamReaderTest, set_RefName)
{
    Alignment a;
    BamReader r;
    r.Open("dummies/test.bam");

    r.GetNextAlignment(a);
    EXPECT_EQ("one", a.RefName);

    r.GetNextAlignment(a);
    EXPECT_EQ("two", a.RefName);
}
