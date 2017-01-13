#include "main.hpp"
#include <gflags/gflags.h>


DEFINE_bool(quiet, false, "Do not output the found squares");


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
	
	const Vector<Square> s = compute_distinct_squares(text);
	std::cout << "Number of Squares: " << s.size() << std::endl;
	if(!FLAGS_quiet) {
		for(size_t i = 0; i < s.size(); ++i) {
			std::cout << "T[" << s[i].start << "," << (s[i].start+2*s[i].period ) << "]: " << text.substr(s[i].start,2*s[i].period) << std::endl;
		}
	}

	return 0;
}
