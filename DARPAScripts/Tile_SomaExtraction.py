import os, subprocess, sys

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print 'Usage: '+sys.argv[0]+' <ImageFileList.txt> <SomaCentroidsFileList.txt> <SomaImageFileList.txt>\n'
        sys.exit(1)
    montage_tile_file=open(sys.argv[1], 'r');
    soma_tile_file=open(sys.argv[3], 'w');
    soma_centroids_file=open(sys.argv[2], 'w');

    montage_tile_lines = montage_tile_file.readlines()

    number_of_lines = len(montage_tile_lines)
    for x in range(number_of_lines):
	print('SomaExtraction.exe ' + montage_tile_lines[x].rstrip('\n') + " Tile_" + montage_tile_lines[x].rstrip('\n') + ".txt ParamFile Soma_" + montage_tile_lines[x].rstrip('\n') + " 5 15000 128")
        os.system('SomaExtraction.exe ' + montage_tile_lines[x].rstrip('\n') + " Tile_" + montage_tile_lines[x].rstrip('\n') + ".txt ParamFile Soma_" + montage_tile_lines[x].rstrip('\n') + " 5 15000 128")
        soma_tile_file.write("Soma_" + montage_tile_lines[x].rstrip('\n') + '\n' )
        soma_centroids_file.write(" Tile_" + montage_tile_lines[x].rstrip('\n') + ".txt\n")
