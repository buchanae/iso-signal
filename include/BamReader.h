#ifndef ISOSIGNAL_BAMREADER_H
#define ISOSIGNAL_BAMREADER_H

#include "api/BamReader.h"

#include "Alignment.h"

class BamReader : public BamTools::BamReader
{
    public:
        bool GetNextAlignment(Alignment& a)
        {
            if (!BamTools::BamReader::GetNextAlignment(a)) return false;

            a.RefName = GetReferenceData().at(a.RefID).RefName;
            return true;
        }
};

#endif
