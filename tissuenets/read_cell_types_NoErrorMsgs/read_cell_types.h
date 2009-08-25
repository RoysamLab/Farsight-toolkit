
#include "tinyxml.h"
#include <set> //store histogram data in a multiset
#include <string>
#include <iostream>

using namespace std;

// Implemented by Aytekin Vargun:
// This class is used for determining the list of all different cell types in an XML file.
// NOTE: This version is aimed to be used from the TissueNets python version
// Therefore no error messages will be produced. If there are any errors in the input we simply return nothing
// Please use read_cell_types if error messages are required too


class CellTypes
{
public:
	bool ReadCellTypes(char* fileName);
	void Display();
private:
	set<string> cellTypes;
	set<string>::iterator pos;	
};