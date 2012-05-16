#include <istream>
#include <map>
#include <string>
#include <vector>

#include "BamReader.h"
#include "CoverageBuilder.h"
#include "Feature.h"
#include "helpers.h"
#include "StackReader.h"

using std::istream;
using std::string;
using std::vector;

void CoverageBuilder::addJunctionsFromGFF(istream& gff_stream)
{
    vector<Feature> all;
    vector<Feature> genes;
    vector<Feature> transcripts;
    getGenesAndTranscriptsFromGFF(gff_stream, all, genes, transcripts);

    vector<string> exon_types;
    exon_types.push_back("exon");
    exon_types.push_back("pseudogenic_exon");

    vector<Feature> juncs;
    for (vector<Feature>::iterator transcript = transcripts.begin(); 
         transcript != transcripts.end(); ++transcript)
    {
        transcript->spliceJunctions(juncs, exon_types);
    }

    // TODO optimize to add directly to index
    filter.junction_index.add(juncs.begin(), juncs.end());
}

void CoverageBuilder::addJunctionsFromStack(istream& stack_stream)
{
    Feature j;
    while (StackReader::getNextFeature(stack_stream, j))
    {
        filter.junction_index.add(j);
    }
}

void CoverageBuilder::addCoverageFromBam(BamReader& reader)
{
    BamTools::RefVector ref_vec = reader.GetReferenceData();
    for (int i = 0; i < ref_vec.size(); ++i)
    {
        BamTools::RefData data = ref_vec.at(i);
        if (coverages.find(data.RefName) == coverages.end())
        {
            coverages.insert(std::make_pair(data.RefName, new Coverage(data.RefLength)));
        }
    }

    Alignment al;
    while (reader.GetNextAlignment(al) && filter(al))
    {
        coverages[al.RefName]->add(al);
    }
}

string CoverageBuilder::coverageString(void)
{
}
