#include "naive.hpp"

int main(int argc, char* argv[])
{
    if(argc < 2) {
		std::cout << "Usage: "<< argv[0] << " file" << std::endl;
        return 1;
    }
	std::ifstream t(argv[1]);
	std::string text((std::istreambuf_iterator<char>(t)),
			std::istreambuf_iterator<char>());
	check_square_algo(text);
}

