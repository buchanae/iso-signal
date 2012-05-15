#ifndef _ISOSIGNAL_HELPERS_H
#define _ISOSIGNAL_HELPERS_H

#include <vector>
#include <iostream>

#include "api/BamReader.h"

#include "Feature.h"
#include "JunctionIndex.h"

#include "Alignment.h"

using std::vector;

void getGenesAndTranscriptsFromGFF(std::istream& gff_input_stream, 
                                   vector<Feature>& all,
                                   vector<Feature>& genes,
                                   vector<Feature>& transcripts);

#endif
