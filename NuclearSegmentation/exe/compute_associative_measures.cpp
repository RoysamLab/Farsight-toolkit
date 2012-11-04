/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

//#include "ftkObjectAssociation.h"
#include <Nuclear_Association/ftkNuclearAssociationRules.h>
#include <iostream>
#include <stdlib.h>

int main(int argc, char* argv[])
{
	std::string segFName;
	std::string xmlFname;


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
		std::string ruleName, targFileName;
		int outsideDistance, insideDistance, assocType, useAllObj, num_threshs, num_in_fg;
		bool useAllObject,subBkground,use_multiple_thresh;
		subBkground = false; use_multiple_thresh = false; num_threshs=1; num_in_fg=1;
		for(int i=1; i<=atoi(argv[3]); i++)
		{
			std::cout <<"Enter Rule("<<i<<") Name:";
			std::cin >> ruleName;
			std::cout << "Enter Target File("<<i<<") Name:";
			std::cin >> targFileName;
			std::cout << "Enter Ouside Distance("<<i<<") value:";
			std::cin >> outsideDistance;
			std::cout << "Use Whole cell("<<i<<") ?:(1=True,else=False)";
			std::cin >> useAllObj;
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
			cout<<"Subtract Background("<<i<<") ?:(1=True,else=False)";
			cin>>subBkground;
			if( subBkground ){
				cout<<"Use multiple level thresholding?("<<i<<") ?:(1=True,else=False)";
				cin>>use_multiple_thresh;
			}
			if( use_multiple_thresh ){
				cout<<"Number of levels to use in("<<i<<"):";
				cin>>num_threshs;
				cout<<"Number of levels to include in the foreground in("<<i<<"):";
				cin>>num_in_fg;
			}
			cout<<"Enter Association Type("<<i<<"):\n1=min\n2=max\n3=total\n4=average\n5=surroundedness\n";
			cin>>assocType;

            //Now having all the parameters of the ith rule defined, add this rule to the list of rules
			ObjAssoc->AddAssociation(ruleName, targFileName, outsideDistance, insideDistance, useAllObject, subBkground, use_multiple_thresh, num_threshs, num_in_fg, assocType, "" );
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