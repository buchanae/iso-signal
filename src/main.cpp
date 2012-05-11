#include <iostream>
#include <vector>

#include <tclap/CmdLine.h>
#include <api/BamReader.h>

#include "Index.h"
#include "Feature.h"

#include "Alignment.h"
#include "coverage.h"
#include "JunctionIndex.h"
#include "StackReader.h"
#include "helpers.h"

#define VERSION "0.1"

using GFF::Feature;
using GFF::TypeIndex;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

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

    // TODO: ensure all reference data for BAM files are the same
    // TODO needed? do one coverage file at a time?

    // open GFF reference file
    std::ifstream gff_input_stream(gff_file_path.c_str());
    if (!gff_input_stream.is_open())
    {
        cerr << "Error opening reference GFF file. Exiting." << endl;
        exit(0);
    }

    vector<Feature> all;
    vector<Feature> genes;
    vector<Feature> transcripts;
    getGenesAndTranscriptsFromGFF(gff_input_stream, all, genes, transcripts);

    cerr << "Searching reference for splice junctions" << endl;

    JunctionIndex junction_index;
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
            Feature j;
            StackReader::getNextFeature(input_stream, j);
            junction_index.add(j);
        }
    }

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
