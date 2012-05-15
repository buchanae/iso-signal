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

TEST(AlignmentTest, has_RefName)
{
    Alignment a;
    EXPECT_EQ("", a.RefName);
}

TEST(AlignmentTest, position)
{
    BamTools::BamAlignment a;
    a.Position = 100; // Bamtools uses 0-based positioning

    Alignment b(a);
    EXPECT_EQ(101, b.position()); // we want 1-based positions, for consistency

    b.position(100);
    EXPECT_EQ(100, b.position());
}

TEST(AlignmentTest, getJunction)
{
    /*
    First exon is 10 bases.
    Splice junction gap is 20 bases.
    Second exon is 10 bases.

    Splice junction start (A) is the last base of the first exon.
    Splice junction end (B) is the first base of the second exon.

             A                    B
    xxxxxxxxxx--------------------xxxxxxxxxx
    ^        ^                    ^
    1        10                   31
    */
    Alignment a;
    a.RefName = "ref";
    a.position(1);

    Feature f;
    CigarOp op;

    EXPECT_FALSE(a.getJunction(f));

    a.CigarData.clear();  
    op.Type = 'M';
    op.Length = 10;
    a.CigarData.push_back(op);

    EXPECT_FALSE(a.getJunction(f));

    op.Type = 'N';
    op.Length = 20;
    a.CigarData.push_back(op);

    op.Type = 'M';
    op.Length = 10;
    a.CigarData.push_back(op);

    EXPECT_TRUE(a.getJunction(f));
    EXPECT_EQ("ref", f.seqid);
    EXPECT_EQ(10, f.start);
    EXPECT_EQ(31, f.end);
}

TEST(AlignmentTest, getJunction_real_world)
{
    /*
    A real-world example from TAIR10, AT1G01100.1

    GFF definition:
    ...other exons...
    Chr1  TAIR10  exon  50419 50631 . - . Parent=AT1G01100.1
    Chr1  TAIR10  exon  50883 50963 . - . Parent=AT1G01100.1
    ...other exons...

    notice the junction between 50631 and 50883

    Valid alignment, SAM format:
    ID  163 Chr1  50621 255 11M251N40M  = 50942 121 ATCG
    */

    Alignment a;
    a.position(50621);

    Feature f;
    CigarOp op;

    a.CigarData.clear();  
    op.Type = 'M';
    op.Length = 11;
    a.CigarData.push_back(op);

    op.Type = 'N';
    op.Length = 251;
    a.CigarData.push_back(op);

    op.Type = 'M';
    op.Length = 40;
    a.CigarData.push_back(op);

    EXPECT_TRUE(a.getJunction(f));
    EXPECT_EQ(50631, f.start);
    EXPECT_EQ(50883, f.end);
}
