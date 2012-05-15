#include <iostream>
#include <vector>

#include "Feature.h"
#include "Index.h"
#include "Reader.h"

#include "helpers.h"

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
