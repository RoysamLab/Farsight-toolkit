import os, subprocess, sys

if __name__ == '__main__':
	image_list_file=open(sys.argv[1], 'r');
	image_list = montage_tile_file.readlines()

	number_of_lines = len(image_list)
	for x in range(number_of_lines):
                command = 'TIF16BitTo8BitConvertor.exe ' + image_list[x].rstrip('\n')
                print(command)
		os.system(command)
