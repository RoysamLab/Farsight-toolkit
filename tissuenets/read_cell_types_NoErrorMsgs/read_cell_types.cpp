#include "read_cell_types.h" 

bool CellTypes::ReadCellTypes(char* fileName)
{  
  TiXmlDocument doc(fileName);
  doc.LoadFile();
  TiXmlHandle docHandle( &doc );
  TiXmlElement* levelOneElement =
  docHandle.FirstChild("Farsight_Output").Element();
  //docHandle.FirstChild("Trace").FirstChild().Element();
  TiXmlElement *levelTwoElement;
  TiXmlElement *levelThreeElement;
  const char *nodeName;
 
//  while(levelOneElement)
//    {
  if (levelOneElement) {
	nodeName = levelOneElement->Value();

	//Check if this is a histogram_data file
    if (strcmp(nodeName,"Farsight_Output") == 0)
      {
      levelTwoElement = levelOneElement->FirstChildElement();
      //while(levelTwoElement)
	  if (levelTwoElement)
        {
			nodeName = (char*)levelTwoElement->Value();
			if (strcmp(nodeName,"Nuclei") == 0)
				{
					levelThreeElement = levelTwoElement->FirstChildElement();
					while(levelThreeElement) {
					    nodeName = levelThreeElement->Value();
						//Check if there is tag <d> first then get the data located between <d> and </d>
						if (strcmp(nodeName,"Nuclear_Features") == 0) {
							// Check first if the Class_Membership attribute exists
							if (levelThreeElement->Attribute("Class_Membership"))
								cellTypes.insert(levelThreeElement->Attribute("Class_Membership"));
							//else cerr<<"Missing Class_Membership"<<endl;
						} else {
							//cerr<<"Check XML file's format. Found element tag which is different from Nuclear_Features"<<endl;
							return false;
						}
						levelThreeElement = levelThreeElement->NextSiblingElement("Nuclear_Features");
					}
			} else {
				//cerr<<"Check XML file's format. Found element tag which is different from Nuclei"<<endl;
				return false;
			}

	  } else return false; //cerr<<"Check XML file's format. It should have Nuclei tag"<<endl;
	} else {
				//cerr<<"Check XML file's format. It should start with Farsight_Output tag"<<endl;
				return false;
			}
  }
	//levelOneElement=levelOneElement->NextSiblingElement();
  return true;
}

void CellTypes::Display() {

	int i=1;
	for (pos=cellTypes.begin();pos!=cellTypes.end();pos++) {
		if (i != cellTypes.size())
			cout<<*pos<<" ";	
		else cout<<*pos;
		i++;
	}

	//cout<<"\""<<*pos<<"\""<<" ";	
}

int main(int argc, char* argv[])
{
  	if (argc <2) {			
			//cout<<"Incorrect Usage! Enter the name of the XML file" <<endl;
			//cout<<"Usage: read_cell_types [XMLfilename]"<<endl;
			return 0;
	};

	CellTypes ct;
	
	if (! ct.ReadCellTypes(argv[1])) {
		//cerr<<"No data!"<<endl;
		return 0;
	} else {
		// Show the cell types extracted from the file
		ct.Display();
	}
	return 0;
}

