1. 	Create a directory - put "src" and "bin" in that directory.
3.	Open CMake. "Configure" and "Generate" using the src and bin directories. 
The code needs 2 packages, ITK and libxml2 (http://xmlsoft.org/downloads.html). 
Please make sure that they are there in the set correctly while configuring CMake.
4.	Build the files.
5. 	Run the program with syntax - 
		(for Linux)		./ftl2d_SETracing.exe parafile.txt	 	
		(for Windows)	ftl2d_SETracing.exe parafile.txt	 	

It should run and give the output similar to the output image attached.
