#include <iostream>
#include <map>
#include <vector>

#include <tclap/CmdLine.h>

#include "Index.h"
#include "Feature.h"

#include "Alignment.h"
#include "BamReader.h"
#include "Coverage.h"
#include "JunctionFilter.h"
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
    std::ifstream gff_stream(gff_file_path.c_str());
    if (!gff_stream.is_open())
    {
        cerr << "Error opening reference GFF file. Exiting." << endl;
        return 0;
    }

    JunctionFilter filter;
    indexJunctionsFromGFF(gff_stream, filter.junction_index);

    cerr << "Loading splice junctions from stack files" << endl;

    // load splice junctions from stack files
    for (vector<string>::iterator it = stack_file_paths.begin();
         it != stack_file_paths.end(); ++it)
    {
        std::ifstream stack_stream(it->c_str());
        if (!stack_stream.is_open())
        {
            cerr << "Error opening stack file: " << *it << endl;
            cerr << "Skipping file." << endl;
        }
        else
        {
            indexJunctionsFromStack(stack_stream, filter.junction_index);
        }
    }

    Coverage coverage;

    // initialize references
    BamTools::RefVector ref_vec = reader.GetReferenceData();
    for (int i = 0; i < ref_vec.size(); ++i)
    {
        BamTools::RefData data = ref_vec.at(i);
        coverage.setMinReferenceLength(data.RefName, data.RefLength);
    }

    // read and filter alignments, adding to coverages
    Alignment al;
    while (reader.GetNextAlignment(al) && filter(al))
    {
        coverage.add(al);
    }

    reader.Close();

    return 0;
}
