#include <vector>

#include "api/BamAlignment.h"
#include "api/BamAux.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "Feature.h"

#include "Alignment.h"

using GFF::Feature;
using testing::ElementsAre;
using testing::WhenSorted;

TEST(AlignmentTest, position)
{
    BamTools::BamAlignment a;
    a.Position = 100; // Bamtools uses 0-based positioning

    Alignment b(a);
    EXPECT_EQ(101, b.position()); // we want 1-based positions, for consistency

    b.position(100);
    EXPECT_EQ(100, b.position());
}

TEST(AlignmentTest, getFeature)
{
    Alignment a;
    a.position(100);

    Feature f;
    CigarOp op;

    EXPECT_FALSE(a.getJunction(f));

    a.CigarData.clear();  
    op.Type = 'M';
    op.Length = 10;
    a.CigarData.push_back(op);

    EXPECT_FALSE(a.getJunction(f));

    op.Type = 'N';
    op.Length = 50;
    a.CigarData.push_back(op);

    op.Type = 'M';
    op.Length = 10;
    a.CigarData.push_back(op);

    EXPECT_TRUE(a.getJunction(f));
    EXPECT_EQ(110, f.start);
    EXPECT_EQ(160, f.end);
}
