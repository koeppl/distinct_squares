#include "main.hpp"
#include <gflags/gflags.h>


DEFINE_bool(quiet, false, "Do not output the found squares");
DEFINE_bool(verbose, false, "Output integer values of the characters");




int main(int argc, char* argv[])
{
	google::ParseCommandLineFlags(&argc, &argv, true);

    if(argc == 1 || strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0 || strcmp(argv[1],"-help") == 0) {
//        google::ShowUsageWithFlagsRestrict(argv[0], __FILE__); //shortcut
        std::vector<google::CommandLineFlagInfo> info;
        google::GetAllFlags(&info);

        std::cout << argv[0] << " [options] {file to analyze}" << std::endl;
        std::cout << "You need to provide at least a file whose distinct squares should be computed." << std::endl << std::endl;
        std::cout
            << std::setw(20) << std::setiosflags(std::ios::left) << "Parameter"
            << std::setw(10) << "Type"
            << std::setw(20) << "Default"
            << "Description" << std::endl;
        std::cout << std::endl;
         for(auto it = info.cbegin(); it != info.cend(); ++it) {
             if(it->filename != __FILE__) continue;
                 std::cout
                    << std::setw(20) << std::setiosflags(std::ios::left)<< (std::string("--")+ it->name)
                    << std::setw(10) << it->type
                    << std::setw(20) << (std::string("(") + it->default_value + ")")
                    << it->description << std::endl;
         }
        return 0;
    }

	std::ifstream t(argv[1]);
	
	if(!((bool)t)) {
		std::cerr << "Could not open file " << argv[1] << std::endl;
		return 1;
	}

	std::string text((std::istreambuf_iterator<char>(t)),
			std::istreambuf_iterator<char>());
	
	size_t square_counter = 0;
	if(FLAGS_quiet) {
		compute_distinct_squares(text, [&text,&square_counter] (len_t, len_t ) {
				++square_counter;
				});
	} else {
		compute_distinct_squares(text, [&text,&square_counter] (len_t pos, len_t period) {
				std::cout << "T[" << pos << "," << (pos+period*2-1) << "] = " << text.substr(pos,period) << "," << text.substr(pos+period,period) << " | ";
				++square_counter;
				if(FLAGS_verbose) {
					for(len_t i = pos; i < pos+period*2; ++i) { 
						if(i == pos+period) DVLOG(1) << ","; 
						std::cout << "(" << ((size_t)text[i]) << ")"; 
					}
				}
				std::cout << std::endl;
			});
	}
	return 0;
}
