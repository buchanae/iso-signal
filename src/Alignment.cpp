#include <vector>

#include "api/BamAux.h"
#include "Feature.h"

#include "Alignment.h"

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

    junction.start = Position + len;
    junction.end = junction.start + gap_len;

    return true;
}
