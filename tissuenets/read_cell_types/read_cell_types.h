
#include "tinyxml.h"
#include <set> //store histogram data in a multiset
#include <string>
#include <iostream>

using namespace std;

// Aytekin Vargun:
// This class is used for determining the list of all different cell types in an XML file.
//


class CellTypes
{
public:
	bool ReadCellTypes(char* fileName);
	void Display();
private:
	set<string> cellTypes;
	set<string>::iterator pos;	
};