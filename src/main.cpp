#include <iostream>
#include <vector>

#include <tclap/CmdLine.h>
#include <api/BamReader.h>

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
    string gff_file_path, input_bam1, input_bam2, stack_file_list;
    vector<string> stack_files;

    try
    {
        TCLAP::CmdLine cmd("Program description", ' ', VERSION);

        TCLAP::ValueArg<string> inputSTACKS("s", "stack-files", "List of stack files", false, "", "stacks1,stacks2,stacks3,etc.", cmd);
        TCLAP::ValueArg<string> inputGFF("g", "gff-file", "Input GFF file", true, "", "input_file.gff", cmd);
        TCLAP::ValueArg<string> inputBAM2("2", "bam-file2", "Input BAM file 2", true, "", "input_file2.bam", cmd);
        TCLAP::ValueArg<string> inputBAM1("1", "bam-file1", "Input BAM file 1", true, "", "input_file1.bam", cmd);

        cmd.parse(argc, argv);

        gff_file_path = inputGFF.getValue();
        input_bam1 = inputBAM1.getValue();
        input_bam2 = inputBAM2.getValue();
        stack_file_list = inputSTACKS.getValue();

        // Clean up input file names
        for (int i = 0; i < stack_files.size(); i++) {
            if (stack_files.at(i).size()) {
                bool cont = true;
                string* val = &stack_files.at(i);
                while (val->size() && cont) {
                   if (val->at(0) == ' ') val->erase(0, 1);
                   else cont = false;
                }
                cont = true;
                while (val->size() && cont) {
                   if (val->at(val->size() - 1) == ' ') val->erase(val->size() - 1, 1);
                   else cont = false;
                }
            }
        }
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

    // TODO: ensure all reference data for BAM files are the same

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

    // index all features from GFF reference by type
    TypeIndex types;
    types.add(all_features.begin(), all_features.end());

    JunctionIndex junction_index;

    // load splice junctions from stack files
    for (vector<string>::iterator it = stack_files.begin();
         it != stack_files.end(); ++it)
    {
        std::ifstream input_stream((*it).c_str());
        if (!input_stream.is_open())
        {
            cerr << "Error opening stack file: " << *it << endl;
            cerr << "Skipping file." << endl;
        } else {
            Feature j;
            StackReader::getNextFeature(input_stream, j);
            junction_index.add(j);
        }
    }

    vector<Feature> genes;
    types.type("gene", genes);
    types.type("pseudogene", genes);
    types.type("transposable_element_gene", genes);

    int count = 0;

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

        if (count % 5 == 0) cerr << "\r" << count;
    }
    cerr << "\r" << count << endl << "Done." << endl;

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
