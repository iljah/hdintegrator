#include "cstdlib"
#include "ctime"
#include "iostream"
#include "string"
#include "sstream"
#include "vector"

int main(int, char**)
{
	std::srand(std::time(0));

	while (true) {
		if (std::rand() % 2 == 0) {
			std::string line;
			std::getline(std::cin, line);
			std::istringstream iss(line);
			double number;
			iss >> number;

			std::vector<double> mins, maxs;
			while (iss >> number) {
				mins.push_back(number);
				if (iss >> number) {
					maxs.push_back(number);
				} else {
					break;
				}
			}
			if (mins.size() == 0) {
				break;
			}
			if (mins.size() != maxs.size()) {
				std::cerr << "Number of minimum and maximum extents differs" << std::endl;
				return EXIT_FAILURE;
			}

			std::cout << "0 0 0" << std::endl;
		} else {
			break;
		}
	}

	return 1;
}
