/*****************************************************************************************/
// This program reads a text file and converts it into XML file which
// can be read by the TissueNets program.
// The file should be structured like this:
// The first line contains variable names
// Values start from the second row
// Values are delimited by comma. Variable names ID,X,Y,Z,Class_Membership 
// are required to be given in the first row. Additional variable names and their 
// values can be given
// 
// Here is an example file
// 
// ID, X,Y,Z,Size,Class_Membership
// 2,74,76,8,1869,Neuron
// 3,127,132,23,1534,Neuron// 4,127,132,23,1534,Neuron
// 5,127,132,23,1534,Neuron
// 6,127,132,23,1000,Neuron
/*****************************************************************************************/
// Use this program to read FARSIGHT's text output 
#include <iostream>
#include <iomanip>
#include <fstream>
#include <set>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

vector<string> parse(string textLine) {
	
	vector<string> tokens;
	//tokens.reserve(30);
	vector<string>::iterator i;
	size_t pos = 0;
	string s,textLineWithoutSpace;	
	int c,j=0;

	//Remove the empty spaces first. We cannot use them as variable names
	for (c=0;c<textLine.size();c++) {	
		if (textLine.compare(c,1," ") != 0) {			
			textLineWithoutSpace.append(textLine.substr(c,1));
			j++;
		}
	}

	while( true ) {
		size_t nextPos = textLineWithoutSpace.find(',', pos);		
		if( nextPos == textLineWithoutSpace.npos ) {			
			tokens.push_back( string( textLineWithoutSpace.substr( pos, (textLineWithoutSpace.size() - pos))));
			break;
		};
		tokens.push_back( string( textLineWithoutSpace.substr( pos, nextPos - pos ) ) );
		pos = nextPos + 1;
	}

	/*
	// List the tokens for testing purposes
	for (i=tokens.begin();i!=tokens.end();i++) 
		cout<<*i<<endl;
	*/
	return tokens;
	
}

bool findv(vector<string> varNames, string s) {
	int i;
	for (i=0;i<varNames.size(); i++)
		if (varNames[i] == s) return true;
	
	return false;
}	

bool checkVariableNames(vector<string> s) {

	vector<string> mustVarNames;
	bool flag=true;
	// These variable names and corresponding values have to be read from the input file
	mustVarNames.push_back("ID");
	mustVarNames.push_back("X");
	mustVarNames.push_back("Y");
	mustVarNames.push_back("Z");
	mustVarNames.push_back("Class_Membership");
	if (s.size() < 5) {
		cerr << "There are less then 5 variables"<<endl;
		cerr << "ID, X, Y, Z, and, Class_Membership are required"<<endl;
		return false;
	} else {
		int i;
		for (i=0;i<mustVarNames.size(); i++)
			if (!findv(s,mustVarNames[i])) {
				cerr<<mustVarNames[i]<<" is required"<<endl;
				flag=false;				
			}
		if (!flag) return false;
	}
	return true;
}
		
void writeXML(vector<vector<string>> values, char* fileName) {	
	ofstream outFile;
	outFile.open(fileName);
	int i,j;
    if (outFile.fail()) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }

	string firstLine="<Farsight_Output xml_tb_version=\"3.1\">";	
	//  + "3.1" + ">";
	stringstream ss1, ss2;
	string s1,s2;
	// Find the number of elements
	//The first element contains variable names. Don't count it
	ss1 << values.size() - 1;
	s1 = ss1.str();
	// Find the number of features (=number of variable names)
	ss2<<values[0].size();
	s2=ss2.str();

	string secondLine="  <Nuclei Number_Of_Nuclei=\"" + s1 + "\" Number_Of_Features=\"" + s2 + "\">";
	outFile<< firstLine<<endl;
	outFile<< secondLine<<endl;
	string element;
	for (i=1;i<values.size();i++) {
		element= "\t<Nuclear_Features ";
		// values[i].size would be better but if there is an extra column
		// we will get error. It is safer like this.
		for (j=0;j<values[0].size();j++)
			element= element + values[0][j] + "=\"" + values[i][j]+"\" ";		
		element=element + ">";
		outFile<<element<<endl;
		outFile<<"\t</Nuclear_Features>"<<endl;
	}
	outFile<<"  </Nuclei>"<<endl;
	outFile<<"</Farsight_Output>"<<endl;
	outFile.close();
}


// Use it like this: FindMedian(sizes);
void FindMedian(multiset<double> numbers) {

	int numOfElements,j;
	double num1, num2;
	multiset<double>::iterator i;
	
	if (numbers.size() > 0) {
		numOfElements = numbers.size();		
		j=0;
		for (i=numbers.begin();i!=numbers.end();i++) {
			j++;
			if (j == numOfElements/2)
				num1 = *i;
			else if (j == ((numOfElements/2) + 1)) {
				num2 = *i;
				break;
			}
		}
		if ((numOfElements % 2) == 0) {			
			cout<<"The median is "<< (num1+num2)/2<<endl;
			cout<<"burda"<<endl;
		} else {
			cout<<"The median is "<< num2<<endl;
		}
	} else cout<<"There are no numbers"<<endl;

};

int main(int argc, char* argv[]) {

	multiset<double> sizes;
	multiset<double>::iterator i;	
	vector<string> variableNames;
	vector<string> valueTokens;
	valueTokens.reserve(30);
	vector<vector<string>> values;
	vector<vector<string>>::iterator k;
	vector<string>::iterator v;	
	string value;
	string varNamesLine;
	char XMLfileName[1024]; //The name of the XML file

	if (argc <2) {
			cout<<"Incorrect Usage! Enter the text file's name" <<endl;
			return 0;
	};

	ifstream inFile;
    inFile.open(argv[1]);
    if (!inFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
	
	// The first line of the file should contain variable names
	while(getline(inFile, value)) {		
		valueTokens=parse(value);		
		values.push_back(valueTokens);
		valueTokens.erase(valueTokens.begin(),valueTokens.end());
		//valueTokens.clear();
	}

	// For clarity, copy the first vector element
	// Since the deletion is expensive with vectors we don't delete 
	// the copied element. But we know that values of these variables 
	// start from the second element
	variableNames=values[0];
	//Check the variable names now. There has to be at least 5 variables 
	//
	if (checkVariableNames(variableNames)) {		
		strcat(XMLfileName,argv[1]);
		strcat(XMLfileName,".XML");
		writeXML(values, XMLfileName);
		//for (t=0;t<variableNames.size();t++)
		//	cout<<variableNames[t]<<endl;	
	} else {
	    inFile.close();
		exit(1);
	}
  
    inFile.close();

    return 0;
}
