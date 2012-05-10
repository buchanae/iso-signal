#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <stdlib.h>

#include <types.h>
#include <utils.h>

struct Record {
    unsigned int count;
    std::string seq;
    int frag1_length;
    int frag2_length;

    Record(): count(0),
              frag1_length(-1),
              frag2_length(-1) 
    { }
};

struct Stack {
   splice_junction sj;
   unsigned int count;
   std::vector<Record> records;

   Stack(): count(0)
   { }
};

class STACKReader {

    // Public interface
    public:
        //ctor
        STACKReader ();
        STACKReader (const std::vector<std::string>);

        //dtor
        ~STACKReader (void);

        // file handling
        bool open(std::vector<std::string>);
        bool open(void);
        bool close(void);

        // Retrieve features in your GFF file
        void loadStacks(void);

        // Retrieve splice junctions
        void get_splice_junctions(std::vector<splice_junction>&, std::string, int, int);

    // Member variables
    private:
        std::vector<std::string> filenames;
        int file_index;
        std::fstream stack_file;
        std::multimap<splice_junction, Stack> stacks;

    // Private functions
    private:
        void _add_stack(Stack&);
        bool _parse_stack(Stack&);
        bool _parse_splice_junction(splice_junction&, std::string&);
        bool _parse_record(Record&, std::string&);
};
