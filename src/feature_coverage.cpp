#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include <tclap/CmdLine.h>

#include "Feature.h"
#include "Reader.h"

#include "Coverage.h"
#include "JunctionIndex.h"
#include "helpers.h"

#define VERSION "0.1"

using GFF::Feature;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

typedef TCLAP::ValueArg<string> string_arg;

int main (int argc, char* argv[])
{
    string gff_file_path, coverage_file_path, output_file_path;

    try
    {
        TCLAP::CmdLine cmd("Program description", ' ', VERSION);

        string_arg coverageFileArg("c", "coverage", "Coverage file", false, "",
                                   "foo.coverage", cmd);

        string_arg gffFileArg("g", "gff", "GFF file", true, "", "reference.gff", cmd);

        string_arg outputFileArg("o", "output", "Output file", true, "", 
                                 "output.coverage", cmd);

        cmd.parse(argc, argv);

        gff_file_path = gffFileArg.getValue();
        coverage_file_path = coverageFileArg.getValue();
        output_file_path = outputFileArg.getValue();

    } catch (TCLAP::ArgException &e) {
        cerr << "Error: " << e.error() << " " << e.argId() << endl;
    }

    std::ostream* output_stream;
    std::ofstream output_file_stream;

    if (output_file_path == "-")
    {
        cerr << "Outputting to standard out." << endl;
        output_stream = &cout;
    }
    else
    {
        output_file_stream.open(output_file_path.c_str(),
                                std::ios::out | std::ios::trunc);

        if (!output_file_stream.is_open())
        {
            cerr << "Error opening output file. Exiting." << endl;
            return 0;
        }
        output_stream = &output_file_stream;
    }

    // open GFF reference file
    std::ifstream coverage_stream(coverage_file_path.c_str());
    if (!coverage_stream.is_open())
    {
        cerr << "Error opening the coverage file. Exiting." << endl;
        return 0;
    }

    cerr << "Loading the coverage file." << endl;

    Coverage coverage;
    loadCoverage(coverage_stream, coverage);

    // open GFF reference file
    std::ifstream gff_stream(gff_file_path.c_str());
    if (!gff_stream.is_open())
    {
        cerr << "Error opening the reference GFF file. Exiting." << endl;
        return 0;
    }

    cerr << "Loading the reference GFF." << endl;

    vector<Feature> transcripts;
    GFF::ParentChildIndex exons_index;

    Feature f;
    while (GFF::Reader::getNextFeature(gff_stream, f))
    {
        if (isTranscriptType(f))
        {
            transcripts.push_back(f);
        }

        if (isExonType(f))
        {
            exons_index.add(f);
        }
    }

    for (vector<Feature>::iterator transcript = transcripts.begin(); 
         transcript != transcripts.end(); ++transcript)
    {
        vector<Feature> junctions;
        vector<Feature> exons;
        exons_index.childrenOf(*transcript, exons);
        spliceJunctions(exons, junctions);

        int transcript_total = 0;
        vector<int> exon_totals(exons.size(), 0);
        vector<int> junction_totals(junctions.size(), 0);

        for (int j = 0; j < exons.size(); ++j)
        {
            Feature exon = exons.at(j);
            for (int i = exon.start; i <= exon.end; ++i)
            {
                int c = coverage.get(exon.seqid, i);
                exon_totals[j] += c;
                transcript_total += c;
            }
        }

        for (int j = 0; j < junctions.size(); ++j)
        {
            Feature junction = junctions.at(j);
            for (int i = junction.start + 1; i <= junction.end - 1; ++i)
            {
                int c = coverage.get(junction.seqid, i);
                junction_totals[j] += c;
                transcript_total += c;
            }
        }

        int transcript_length = transcript->getLength();
        double transcript_coverage = (double) transcript_total / transcript_length;

        for (int j = 0; j < junctions.size(); ++j)
        {
            double percent_of_expected = 0;

            if (junction_totals[j] > 0)
            {
                Feature junction = junctions.at(j);
                double coverage = (double) junction_totals[j] / transcript_length;

                double percent_length_of_transcript = 
                    (double) junction.getLength() / transcript_length;

                double expected = percent_length_of_transcript * transcript_coverage;
                percent_of_expected = coverage / expected;
            }

            if (percent_of_expected > 0)
            {
                string ID;
                transcript->attributes.get("ID", ID);
                cout << ID << "-intron-" << j + 1 << "\t" << percent_of_expected;
                cout << endl;
            }
        }
    }

    if (output_file_stream.is_open())
    {
        output_file_stream.close();
    }

    return 0;
}
