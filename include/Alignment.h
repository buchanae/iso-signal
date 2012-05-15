#ifndef _ISOSIGNAL_ALIGNMENT_H
#define _ISOSIGNAL_ALIGNMENT_H

#include <string>

#include "api/BamAlignment.h"

#include "Feature.h"

using BamTools::CigarOp;

class Alignment : public BamTools::BamAlignment
{
    public:
        std::string RefName;

        Alignment(void);
        Alignment(BamTools::BamAlignment& other);
        bool getJunction(GFF::Feature&);
        int position(void);
        void position(int);

    private:
        // Prevent public access to Position, because it's zero-based,
        // which can cause confusion and difficult debugging.
        // Use position() instead.
        BamTools::BamAlignment::Position;
};
#endif
