#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
#include <string>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include "Feature.h"
#include "StackReader.h"

using std::string;

typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
typedef tokenizer::iterator token_iter;
typedef boost::char_separator<char> separator;

bool StackReader::getNextFeature(std::istream& input_stream, GFF::Feature& f)
{
    string line;

    while (std::getline(input_stream, line).good())
    {
        if (boost::starts_with(line, "@"))
        {
            std::vector<string> cols;
            tokenizer tokens(line, separator("\t"));
            std::copy(tokens.begin(), tokens.end(), std::back_inserter(cols));

            f.seqid = cols.at(0).substr(1);            
            f.type = "splice_junction";
            f.start = boost::lexical_cast<int>(cols.at(2));
            f.end = boost::lexical_cast<int>(cols.at(3));

            return true;
        }
    }
    return false;
}
