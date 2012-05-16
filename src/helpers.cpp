#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "Feature.h"
#include "Index.h"
#include "Reader.h"

#include "Coverage.h"
#include "helpers.h"
#include "StackReader.h"

using std::vector;

void getGenesAndTranscriptsFromGFF(std::istream& gff_input_stream, 
                                   vector<Feature>& all,
                                   vector<Feature>& genes,
                                   vector<Feature>& transcripts)
{
    // read all features from GFF reference
    GFF::Reader::readAllAndLinkChildren(gff_input_stream, all);

    // index all features from GFF reference by type
    GFF::TypeIndex types;
    types.add(all.begin(), all.end());

    types.type("gene", genes);
    types.type("pseudogene", genes);
    types.type("transposable_element_gene", genes);

    types.type("mRNA", transcripts);
    types.type("mRNA_TE_gene", transcripts);
    types.type("ncRNA", transcripts);
    types.type("miRNA", transcripts);
    types.type("snoRNA", transcripts);
    types.type("snRNA", transcripts);
    types.type("rRNA", transcripts);
    types.type("tRNA", transcripts);
    types.type("pseudogenic_transcript", transcripts);
}

void indexJunctionsFromGFF(istream& gff_stream, JunctionIndex& index)
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
    index.add(juncs.begin(), juncs.end());
}

void indexJunctionsFromStack(istream& stack_stream, JunctionIndex& index)
{
    Feature j;
    while (StackReader::getNextFeature(stack_stream, j))
    {
        index.add(j);
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
