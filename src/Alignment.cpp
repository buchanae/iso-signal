#include <vector>

#include "api/BamAux.h"
#include "Feature.h"

#include "Alignment.h"

Alignment::Alignment(void) : BamTools::BamAlignment() {}
Alignment::Alignment(BamTools::BamAlignment& other) : BamTools::BamAlignment(other) {}

int Alignment::position(void)
{
    return Position + 1;
}

void Alignment::position(int pos)
{
    Position = pos - 1;
}

bool Alignment::getJunction(GFF::Feature& junction)
{
    // this all assumes there is only one gap, i.e. only one 'N' CigarOp

    int len = 0;
    int gap_len = 0;

    for (std::vector<CigarOp>::iterator it = CigarData.begin(); 
         it != CigarData.end(); ++it)
    {
        CigarOp op = *it;
        if (op.Type == 'N')
        {
            gap_len = op.Length;
            break;
        }

        len += op.Length;
    }

    if (gap_len == 0) return false;

    junction.start = position() + len - 1;
    junction.end = junction.start + gap_len + 1;

    return true;
}
