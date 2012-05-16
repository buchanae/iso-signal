#include <iostream>
#include <map>
#include <vector>

#include <tclap/CmdLine.h>
#include <api/BamReader.h>

#include "Index.h"
#include "Feature.h"

#include "Alignment.h"
#include "coverage.h"
#include "CoverageBuilder.h"
#include "JunctionFilter.h"
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
    string gff_file_path, bam_file_path;
    vector<string> stack_file_paths;

    // TODO allow multiple bam files?
    try
    {
        TCLAP::CmdLine cmd("Program description", ' ', VERSION);

        TCLAP::MultiArg<string> inputSTACKS("s", "stack-file", "Stack file", false, "foo.stacks", cmd);
        TCLAP::ValueArg<string> inputGFF("g", "gff-file", "Input GFF file", true, "", "input_file.gff", cmd);
        TCLAP::ValueArg<string> inputBAM1("b", "bam-file", "Input BAM file", true, "", "input_file.bam", cmd);

        cmd.parse(argc, argv);

        gff_file_path = inputGFF.getValue();
        bam_file_path = inputBAM1.getValue();
        stack_file_paths = inputSTACKS.getValue();

    } catch (TCLAP::ArgException &e) {
        cerr << "Error: " << e.error() << " " << e.argId() << endl;
    }

    BamReader reader;
    if(!reader.Open(bam_file_path))
    {
        cerr << "Error opening the bam file. Exiting." << endl;
        return 0;
    }

    cerr << "Loading the reference GFF" << endl;

    // open GFF reference file
    std::ifstream gff_input_stream(gff_file_path.c_str());
    if (!gff_input_stream.is_open())
    {
        cerr << "Error opening reference GFF file. Exiting." << endl;
        return 0;
    }

    CoverageBuilder builder;
    builder.addJunctionsFromGFF(gff_input_stream);

    cerr << "Loading splice junctions from stack files" << endl;

    // load splice junctions from stack files
    for (vector<string>::iterator it = stack_file_paths.begin();
         it != stack_file_paths.end(); ++it)
    {
        std::ifstream input_stream(it->c_str());
        if (!input_stream.is_open())
        {
            cerr << "Error opening stack file: " << *it << endl;
            cerr << "Skipping file." << endl;
        } else {
            builder.addJunctionsFromStack(input_stream);
        }
    }

    builder.addCoverageFromBam(reader);
    reader.Close();

    cout << builder.coverageString();

    return 0;
}
