#include <iostream>
#include <vector>

#include <boost/tokenizer.hpp>
#include <tclap/CmdLine.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>

#include "Index.h"
#include "Reader.h"
#include "Feature.h"

#include "Alignment.h"
#include "coverage.h"
#include "JunctionIndex.h"
#include "StackReader.h"

#define VERSION "0.1"

using GFF::Feature;
using GFF::TypeIndex;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

void getValidAlignmentsToFeature(Feature&, BamTools::BamReader&, JunctionIndex&,
                                 vector<Alignment>&);

int main (int argc, char* argv[])
{
    string gff_file_path, input_bam1, input_bam2;
    vector<string> stack_file_paths;

    try
    {
        TCLAP::CmdLine cmd("Program description", ' ', VERSION);

        TCLAP::MultiArg<string> inputSTACKS("s", "stack-file", "Stack file", false, "foo.stacks", cmd);
        TCLAP::ValueArg<string> inputGFF("g", "gff-file", "Input GFF file", true, "", "input_file.gff", cmd);
        TCLAP::ValueArg<string> inputBAM2("2", "bam-file2", "Input BAM file 2", true, "", "input_file2.bam", cmd);
        TCLAP::ValueArg<string> inputBAM1("1", "bam-file1", "Input BAM file 1", true, "", "input_file1.bam", cmd);

        cmd.parse(argc, argv);

        gff_file_path = inputGFF.getValue();
        input_bam1 = inputBAM1.getValue();
        input_bam2 = inputBAM2.getValue();
        stack_file_paths = inputSTACKS.getValue();

    } catch (TCLAP::ArgException &e) {
        cerr << "Error: " << e.error() << " " << e.argId() << endl;
    }

    BamTools::BamReader bam_file1, bam_file2;
    if(bam_file1.Open(input_bam1) && bam_file2.Open(input_bam2)) {
        if (!bam_file1.LocateIndex() || !bam_file2.LocateIndex()) {
            cerr << "Error: could not locate bam indexes. Exiting." << endl;
            exit(0);
        }
    } else {
        cerr << "Error opening the bam files. Exiting." << endl;
        exit(0);
    }

    BamTools::BamWriter AlignmentsOut;
    AlignmentsOut.Open("new_alignments.bam", bam_file1.GetHeader(), 
                       bam_file1.GetReferenceData());

    // TODO: ensure all reference data for BAM files are the same

    cerr << "Loading GFF reference" << endl;

    // open GFF reference file
    std::ifstream gff_input_stream(gff_file_path.c_str());
    if (!gff_input_stream.is_open())
    {
        cerr << "Error opening reference GFF file. Exiting." << endl;
        exit(0);
    }

    // read all features from GFF reference
    vector<Feature> all_features;
    GFF::Reader::readAllAndLinkChildren(gff_input_stream, all_features);

    cerr << "Indexing reference features by type" << endl;

    // index all features from GFF reference by type
    TypeIndex types;
    types.add(all_features.begin(), all_features.end());

    JunctionIndex junction_index;

    cerr << "Searching reference for splice junctions" << endl;

    vector<Feature> transcripts;
    types.type("mRNA", transcripts);
    types.type("mRNA_TE_gene", transcripts);
    types.type("ncRNA", transcripts);
    types.type("miRNA", transcripts);
    types.type("snoRNA", transcripts);
    types.type("snRNA", transcripts);
    types.type("rRNA", transcripts);
    types.type("tRNA", transcripts);
    types.type("pseudogenic_transcript", transcripts);

    vector<string> exon_types;
    exon_types.push_back("exon");
    exon_types.push_back("pseudogenic_exon");

    for (vector<Feature>::iterator it = transcripts.begin(); 
         it != transcripts.end(); ++it)
    {
        Feature cur = *it;
        vector<Feature> juncs;
        cur.spliceJunctions(juncs, exon_types);
        // TODO optimize to add directly to index
        junction_index.add(juncs.begin(), juncs.end());
    }

    /*
      TODO could be useful to have something to list all splice junctions,
           so maybe find a place for this to live.
    vector<Feature> unique_juncs;
    junction_index.unique(unique_juncs);
    for (vector<Feature>::iterator iit = unique_juncs.begin();
         iit != unique_juncs.end(); ++iit)
    {
        Feature cur_junc = *iit;
        cout << cur_junc.seqid << "\t" << cur_junc.start << "\t" << cur_junc.end << endl;
    }
    */

    // load splice junctions from stack files
    for (vector<string>::iterator it = stack_file_paths.begin();
         it != stack_file_paths.end(); ++it)
    {
        std::ifstream input_stream((*it).c_str());
        if (!input_stream.is_open())
        {
            cerr << "Error opening stack file: " << *it << endl;
            cerr << "Skipping file." << endl;
        } else {

            int count = 0;

            cerr << "Loading stack file: " << *it << endl;

            Feature j;
            while (StackReader::getNextFeature(input_stream, j))
            {
                count++;
                junction_index.add(j);
            }

            cerr << "Loaded " << count << " splice junctions" << endl;
        }
    }

    vector<Feature> genes;
    types.type("gene", genes);
    types.type("pseudogene", genes);
    types.type("transposable_element_gene", genes);

    int count = 0;

    cerr << "Building coverage for " << genes.size() << " genes" << endl;

    for (vector<Feature>::iterator it = genes.begin(); it != genes.end(); ++it)
    {
        Feature feature = *it;
        vector<Feature> local_juncs;
        JunctionIndex local_juncs_index;

        junction_index.overlappingFeature(feature, local_juncs);
        // TODO optimize this so that overlappingFeature adds directly to the index
        local_juncs_index.add(local_juncs.begin(), local_juncs.end());

        vector<Alignment> alignments1, alignments2;
        getValidAlignmentsToFeature(feature, bam_file1, local_juncs_index, alignments1);
        for (vector<Alignment>::iterator it = alignments1.begin();
             it != alignments1.end(); ++it)
        {
            AlignmentsOut.SaveAlignment(*it);
        }

        getValidAlignmentsToFeature(feature, bam_file2, local_juncs_index, alignments2);

        vector<int> coverage1, coverage2;
        double nt_cov1 = 0, nt_cov2 = 0;

        gen_coverage(coverage1, alignments1, feature);
        gen_coverage(coverage2, alignments2, feature);

        for (int i = 0; i < coverage1.size(); i++) nt_cov1 += coverage1.at(i);
        nt_cov1 /= (double)coverage1.size();

        for (int i = 0; i < coverage2.size(); i++) nt_cov2 += coverage2.at(i);
        nt_cov2 /= (double)coverage2.size();

        if (nt_cov1 >= 5.0 && nt_cov2 >= 5.0){
            string ID = "no ID";
            feature.attributes.get("ID", ID);
            comp_coverage(ID, coverage1, coverage2);
        }

        count++;
    }
    cerr << "Done." << endl;

    AlignmentsOut.Close();
    bam_file1.Close();
    bam_file2.Close();

    return 0;
}

// TODO note that splice junction alignments are only valid if they're in junctions arg
void getValidAlignmentsToFeature(Feature& feature,
                                 BamTools::BamReader& reader,
                                 JunctionIndex& junction_index,
                                 vector<Alignment>& alignments)
{
    int ref_ID = reader.GetReferenceID(feature.seqid);

    BamTools::BamRegion region;
    region = BamTools::BamRegion(ref_ID, feature.start, ref_ID, feature.end);

    reader.SetRegion(region);

    Alignment al;
    while (reader.GetNextAlignment(al))
    {
        string ref_name = reader.GetReferenceData().at(al.RefID).RefName;

        Feature junction;
        junction.seqid = ref_name;
        if (al.getJunction(junction))
        {
            if (junction_index.contains(junction))
            {
                alignments.push_back(al);
            }
        } else {
            alignments.push_back(al);
        }
    }
}
