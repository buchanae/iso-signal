#ifndef _ISOSIGNAL_ALIGNMENT_H
#define _ISOSIGNAL_ALIGNMENT_H

#include "api/BamAlignment.h"

#include "Feature.h"

using BamTools::CigarOp;

class Alignment : public BamTools::BamAlignment
{
    public:
        bool getJunction(GFF::Feature&);
};
#endif
