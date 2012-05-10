#include <STACKReader.h>

STACKReader::STACKReader () :
                          file_index(0){
}

STACKReader::STACKReader (const std::vector<std::string> filenames_p):
                          file_index(0),
                          filenames(filenames_p) {
}

STACKReader::~STACKReader () {
    if (stack_file.is_open()) close();
}

bool STACKReader::open(void) {
    if(stack_file.is_open()) close();

    if (file_index < filenames.size()) {
        if (filenames.at(file_index).size() > 0) {
            stack_file.open(filenames.at(file_index).c_str(), std::ios::in);
            if (stack_file.is_open()) {
                file_index++;
                return true;
            }
        }
    }

    return false;
}

bool STACKReader::open(std::vector<std::string> filenames_p) {
    filenames = filenames_p;
    file_index = 0;
    if (stack_file.is_open()) stack_file.close();

    return open();
}

bool STACKReader::close(void) {
    if (stack_file.is_open()) {
        stack_file.close();
        stack_file.clear();
        return true;
    }
    return false;
}

void STACKReader::loadStacks(void) {
    // Flags to control parsing
    bool more_files = true, cont_parsing = true;

    // If no stack file is prepped for parsing, attempt to do that now 
    if (!stack_file.is_open())
        if (!open()) more_files = false;

    // Summary statistics for all stacks parsed
    int total_stacks = 0, total_counts = 0;
    while (more_files) {
        int stack_count = 0;
        std::cerr << "[ Parsing stack file: " << filenames.at(file_index - 1) << " ]";

        while (cont_parsing) {
            Stack stack;
            if (_parse_stack(stack)) {
                _add_stack(stack);

                stack_count++;
                total_counts += stack.count;

                if (stack_count % 1000 == 0) {
                    std::cerr << "\r" << "[ Parsing file: " << filenames.at(file_index - 1);
                    std::cerr << " (" << stack_count << ") ]";
                }
            } else cont_parsing = false;
        }

        total_stacks += stack_count;

        std::cerr << "\r[ Parsing file: " << filenames.at(file_index - 1);
        std::cerr << " (" << stack_count << ") ]" << std::endl;

        if (!open()) more_files = false;
        else cont_parsing = true;
    }

    std::cerr << "[ Processed: " << total_stacks << " stacks (Total Alignments: " << total_counts << ") ]" << std::endl;

    std::cerr << "[ Found " << stacks.size() << " unique splice junctions. ]" << std::endl;
}

void STACKReader::get_splice_junctions(std::vector<splice_junction>& sjs, std::string ref, int start, int end) {

    std::multimap<splice_junction, Stack>::iterator stack_it;
    std::pair<
        std::multimap<splice_junction, Stack>::iterator, 
        std::multimap<splice_junction, Stack>::iterator
    > range;

    for (int pos = start; pos < end; pos++)
    {
        splice_junction sj = {ref, "", pos, -1};
        range = stacks.equal_range(sj);

        for (stack_it = range.first; stack_it != range.second; stack_it++)
        {
            if (stack_it->second.sj.end <= end) sjs.push_back(stack_it->second.sj);
        }
    }
}

void STACKReader::_add_stack(Stack& stack) {
   std::multimap<splice_junction, Stack>::iterator stack_it;
   std::pair< std::multimap <splice_junction, Stack>::iterator, std::multimap <splice_junction, Stack>::iterator> range;

   range = stacks.equal_range(stack.sj);

   // Add new stack
   if (range.first == range.second) {
       stacks.insert(std::pair <splice_junction, Stack> (stack.sj, stack));

   // Merge into existing stack 
   } else {
       bool merged = false;
       for (stack_it = range.first; stack_it != range.second; stack_it++) {
           if (stack_it->second.sj.end == stack.sj.end) {
               for (int i = 0; i < stack.records.size(); i++ ) {
                   stack_it->second.count += stack.records.at(i).count;
                   stack_it->second.records.push_back(stack.records.at(i));
               }
               merged = true;
           }
       }

       if (!merged) stacks.insert(std::pair <splice_junction, Stack> (stack.sj, stack));
   }
}

bool STACKReader::_parse_stack(Stack& stack) {
    if (stack_file.is_open()) {
        std::string line;
        bool cont = true;
        if (getline(stack_file, line).good()) {
            if(!_parse_splice_junction(stack.sj, line)) {
                std::cerr << "Error: splice junction data improperly formatted. Exiting." << std::endl;
                exit(0);
            }

            // Parse and store all record information denoting reads that support this splice junction
            do {
                std::streampos pos = stack_file.tellg();
                if(getline(stack_file, line).good()) {
                    if (line.size() > 0) {
                        Record record;
                        if (line.at(0) == '@') { // Done parsing this stack
                            // Stop parsing, and rewind back to beginning of record
                            cont = false;
                            stack_file.seekg(pos);  
                        } else if (_parse_record(record, line)) {
                            stack.count += record.count;
                            stack.records.push_back(record);
                        } else {
                            std::cerr << "Error: improper stacker format. Exiting." << std::endl;
                        }
                    } 
                } else cont = false;
            } while (cont);

            return true;
        }
    }

    return false;
}

bool STACKReader::_parse_record(Record& record, std::string& line) {
    std::vector<std::string> data;

    _split_string(data, line, '\t');

    if (data.size() == 4) {
        record.count = atoi(data.at(0).c_str());
        record.seq = data.at(1);
        record.frag1_length = atoi(data.at(2).c_str());
        record.frag2_length = atoi(data.at(3).c_str());
    } else false;

    return true;
}

bool STACKReader::_parse_splice_junction(splice_junction& sj, std::string& line) {
    std::vector<std::string> data;

    _split_string(data, line, '\t');

    if (data.size() == 5) { // Properly formatted stacker junction
        data.at(0).erase(0, 1); // Remove '@' symbol
        sj.ref = data.at(0);
        sj.donor_acceptor = data.at(1);
        sj.start = atoi(data.at(2).c_str());
        sj.end = atoi(data.at(3).c_str());
    } else return false;

    return true;
}
