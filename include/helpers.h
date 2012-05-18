#ifndef _ISOSIGNAL_HELPERS_H
#define _ISOSIGNAL_HELPERS_H

#include <vector>
#include <iostream>

#include "api/BamReader.h"

#include "Feature.h"
#include "JunctionIndex.h"

#include "Alignment.h"

using std::istream;
using std::vector;

void getGenesAndTranscriptsFromGFF(istream& gff_input_stream, 
                                   vector<Feature>& all,
                                   vector<Feature>& genes,
                                   vector<Feature>& transcripts);

void indexJunctionsFromGFF(istream& gff_stream, JunctionIndex& index);

void indexJunctionsFromStack(istream& stack_stream, JunctionIndex& index);

void indexJunctionsFromTranscripts(vector<Feature>& transcripts, JunctionIndex& index);

bool isTranscriptType(Feature& f);

bool isExonType(Feature& f);

#endif
