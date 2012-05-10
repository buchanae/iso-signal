#ifndef _ISOSIGNAL_STACKREADER_H
#define _ISOSIGNAL_STACKREADER_H
#include <iostream>

#include "Feature.h"

namespace StackReader
{
    bool getNextFeature(std::istream&, GFF::Feature&);
}
#endif
