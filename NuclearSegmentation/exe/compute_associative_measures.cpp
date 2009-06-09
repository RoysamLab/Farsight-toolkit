/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/

//#include "ftkObjectAssociation.h"
#include <Nuclear_Association/ftkNuclearAssociationRules.h>
#include <iostream>
#include <stdlib.h>

using namespace std;

int main(int argc, char* argv[])
{
	string segFName;
	string xmlFname;
	
	
	//ftk::ObjectAssociation *ObjAssoc;	
	ftk::NuclearAssociationRules *ObjAssoc;
	if(argc ==2)
	{
		//instantiate the association rules object
		//ObjAssoc = new ftk::ObjectAssociation("", 0);
		ObjAssoc = new ftk::NuclearAssociationRules("",0);
		ObjAssoc->ReadRulesFromXML(argv[1]);		
	}
	else if(argc == 4)
	{		
		//instantiate the association rules object
		//ObjAssoc = new ftk::ObjectAssociation(argv[1], atoi(argv[3]));
		ObjAssoc = new ftk::NuclearAssociationRules(argv[1], atoi(argv[3]));
		
		//Now, the user will add the rules (define them) one by one
		string ruleName,targFileName;
		int outsideDistance, insideDistance, assocType, useAllObj;
		bool useAllObject; 
		for(int i=1; i<=atoi(argv[3]); i++)
		{
			cout<<"Enter Rule("<<i<<") Name:";
			cin>>ruleName;
			cout<<"Enter Target File("<<i<<") Name:";
			cin>>targFileName;
			cout<<"Enter Ouside Distance("<<i<<") value:";
			cin>>outsideDistance;
			cout<<"Use Whole cell("<<i<<") ?:(1=True,else=False)";
			cin>>useAllObj;
			if(useAllObj == 1)
			{
				useAllObject = true;
				insideDistance = 0;
			}
			else
			{
				useAllObject = false;
				cout<<"Enter Inside Distance("<<i<<") value:";
				cin>>insideDistance;
			}
			cout<<"Enter Association Type("<<i<<"):\n1=min\n2=max\n3=total\n4=average\n";
			cin>>assocType;

            //Now having all the parameters of the ith rule defined, add this rule to the list of rules
			ObjAssoc->AddAssociation(ruleName, targFileName, outsideDistance, insideDistance, useAllObject, assocType);
		}

		//Now, write the association rules into an xml file (can be used at any later time)
		ObjAssoc->WriteRulesToXML(argv[2]);
	}
	else
	{
		std::cout<<"Usage1: compute_associative_measures <XMLFileName>"<<std::endl; 
		std::cout<<"Usage2: compute_associative_measures <SegmentationResultsFileName> <AssocDefinitionXMLFileName> <NumberOfAssociationRule>"<<std::endl; 
return 0;
	}	
	
	ObjAssoc->PrintSelf();
	ObjAssoc->Compute();
	std::string inImageFName = argv[1];
	std::string outXMLFname = inImageFName.substr(0,inImageFName.find_last_of("."));
	outXMLFname.append("_AssocFeatures.XML");
	ObjAssoc->WriteAssociativeFeaturesToXML(outXMLFname.c_str());
	return 1;
}

