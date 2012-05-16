#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "api/BamAux.h"

#include "Alignment.h"
#include "Coverage.h"

using GFF::Feature;
using std::endl;
using std::vector;
using std::string;

void Coverage::add(Alignment& alignment)
{
    int pos = alignment.position();
    for (std::vector<CigarOp>::iterator op = alignment.CigarData.begin(); 
         op != alignment.CigarData.end(); ++op)
    {
        if (op->Type == 'M')
        {
            add(alignment.RefName, pos, (int)op->Length);
        }
        pos += op->Length;
    }
}

void Coverage::add(string ref_name, int start, int length)
{
    setMinReferenceLength(ref_name, start + length - 1);

    for (int i = start - 1; i < start + length - 1; ++i)
    {
        coverages.find(ref_name)->second.at(i) += 1;
    }
}

void Coverage::setMinReferenceLength(string name, int length)
{
    std::map<string, vector<int> >::iterator coverage = coverages.find(name);
    if (coverage == coverages.end())
    {
        coverages.insert(make_pair(name, vector<int>(length, 0)));
    }
    else if (coverage->second.size() < length)
    {
        coverage->second.resize(length, 0);
    }
}
