#ifndef ISOSIGNAL_JUNCTION_FILTER_H
#define ISOSIGNAL_JUNCTION_FILTER_H

#include <string>

#include "Alignment.h"
#include "Feature.h"
#include "JunctionIndex.h"

class JunctionFilter
{
    public:
        JunctionIndex junction_index;

        bool operator() (Alignment& al)
        {
            Feature junction;
            if (al.getJunction(junction))
            {
                return junction_index.contains(junction);
            }
            return true;
        }
};

#endif
