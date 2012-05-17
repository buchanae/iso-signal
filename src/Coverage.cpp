#include <cstdlib>
#include <iostream>
#include <sstream>
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

int Coverage::get(string ref_name, int pos)
{
    return coverages.find(ref_name)->second.at(pos - 1);
}

void Coverage::set(string ref_name, int pos, int value)
{
    coverages.find(ref_name)->second.at(pos - 1) = value;
}

void Coverage::increment(string ref_name, int pos, int value)
{
    set(ref_name, pos, get(ref_name, pos) + value);
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

void Coverage::add(Alignment& alignment)
{
    int pos = alignment.position();
    for (vector<CigarOp>::iterator op = alignment.CigarData.begin(); 
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

    for (int i = start; i < start + length; ++i)
    {
        increment(ref_name, i);
    }
}

void loadCoverage(std::istream& coverage_stream, Coverage& coverage)
{
    string line;
    string ref_name;
    int i = 1;

    while (std::getline(coverage_stream, line).good())
    {
        size_t pos = line.find("\t");
        if (pos != string::npos)
        {
            i = 1;
            ref_name = line.substr(0, pos);
            int length = atoi(line.substr(pos + 1).c_str());
            coverage.setMinReferenceLength(ref_name, length);
        }
        else if (!line.empty())
        {
            coverage.set(ref_name, i, atoi(line.c_str()));
            ++i;
        }
    }
}

void formatGMBCoverage(Coverage& coverage, std::ostream& coverage_stream)
{
    for (std::map<string, vector<int> >::iterator it = coverage.coverages.begin();
         it != coverage.coverages.end(); ++it)
    {
        coverage_stream << it->first << "\t" << it->second.size() << std::endl;
        for (int i = 0; i < it->second.size(); ++i)
        {
            coverage_stream << it->second.at(i) << std::endl;
        }
    }
}

void formatGMBCoverage(Coverage& coverage, std::string& output)
{
    std::stringstream stream;
    formatGMBCoverage(coverage, stream);
    output = stream.str();
}
