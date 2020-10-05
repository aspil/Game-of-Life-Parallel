#include <iostream>
#include <iostream>
// #include <cstdlib>
#include <string>
#include <fstream>
#include <cstdlib>
#include <ctime>

// #define FILE_NAME "initial_grid.txt"
#define DEAD '-'
#define ALIVE 'x'

int main(int argc, char **argv){
	if (argc != 3) {
		printf("Usage ./generate <N> <filename>\n");
		return 1;
	}
	std::string filename = "test_grids/" + std::string(argv[2]);
	std::ofstream OutFile(filename);

	srand(time(0));
	int width = std::stoi(argv[1]);
    for(int i = 0; i < width; i++){
        for(int j = 0; j < width; j++){
            if((rand()% 100) < 50)
                OutFile << DEAD;
            else
				OutFile << ALIVE;
		}
		OutFile << '\n';
	}
	std::cout << argv[2] << " was added to ./test_grids/ directory." << std::endl;
	return 0;
}