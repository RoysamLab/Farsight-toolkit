#include <fstream>
#include <iostream>
#include <string>
#include <fregl/fregl_util.h>

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		std::cerr << "Usage: " << argv[0] << " pair_list_file.txt" << std::endl;
		return -1;
	}
	std::ifstream pair_list_file;
	pair_list_file.open(argv[1], std::ifstream::in);

	while (!pair_list_file.eof())
	{
		std::string from_image_name;
		std::string to_image_name;

		pair_list_file >> from_image_name >> to_image_name >> std::ws;
		std::cout << "from_image: " << from_image_name << std::endl;
		std::cout << "to_image: " << to_image_name << std::endl << std::endl;

		fregl_util< unsigned char >::fregl_util_read_image( from_image_name, 0, 0);
	}

	
}
