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

/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
#include "ftkObjectAssociation.h"
#include <iostream>
#include <fstream>
#include <iomanip>

namespace ftk 
{	

/* The constructor of the Association Rule Class */
AssociationRule::AssociationRule(std::string name)
{
	SetRuleName(name);	
	segFileName = "";
	targFileName = "";
	outsideDistance = 0;
	insideDistance = 0;
	useWholeObject = false;
	subBkground = false;
	use_multiple_thresh = false;
	grayIm_set=false;
	segIm_set=false;
	num_threshs = 1;
	num_in_fg = 1;
	assocType = ASSOC_AVERAGE;
	binary_path = "";
}

/* From here, we start defining the member functions of the ObjectAssociation class */
ObjectAssociation::ObjectAssociation(std::string AssocFName, int numOfRules)
{
	segImageName = AssocFName;
	numOfAssocRules = numOfRules;		
	assocMeasurementsList=NULL;
	numOfLabels=0;
	//added by Yousef on 10-18-2009
	//invalidObjects = NULL;
}

/* Add association rules to the list of rules */
void ObjectAssociation::AddAssociation(std::string ruleName,std::string targFileName, int outsideDistance, int insideDistance,	bool useAllObject, bool subBkground, bool use_multiple_thresh, int num_threshs, int num_in_fg, int assocType, std::string append_path)
{
	AssociationRule *assocRule = new AssociationRule(ruleName);
	assocRule->SetSegmentationFileNmae(segImageName);
	assocRule->SetTargetFileNmae(targFileName);
	assocRule->SetOutDistance(outsideDistance);
	assocRule->SetInDistance(insideDistance);
	assocRule->SetUseWholeObject(useAllObject);
	assocRule->SetUseBackgroundSubtraction(subBkground);
	if(subBkground)
		assocRule->SetUseMultiLevelThresholding(use_multiple_thresh);
	else
		assocRule->SetUseMultiLevelThresholding(false);
	if(subBkground&&use_multiple_thresh){
		assocRule->SetNumberOfThresholds(num_threshs);
		if( num_threshs >= num_in_fg )
			assocRule->SetNumberOfThresholds(num_in_fg);
		else
			assocRule->SetNumberIncludedInForeground(num_threshs);
	}
	else{
		assocRule->SetNumberOfThresholds(1);
		assocRule->SetNumberIncludedInForeground(1);
	}
	switch(assocType)
	{
	case 1:
		assocRule->SetAssocType(ASSOC_MIN);
		break;
	case 2:
		assocRule->SetAssocType(ASSOC_MAX);
		break;
	case 3:
		assocRule->SetAssocType(ASSOC_TOTAL);
		break;
	case 4:
		assocRule->SetAssocType(ASSOC_AVERAGE);
		break;
	case 5:
		assocRule->SetAssocType(ASSOC_SURROUNDEDNESS);
		break;
	case 6:
		assocRule->SetAssocType(ASSOC_DIST_OBJECT);
		break;
	default:
		assocRule->SetAssocType(ASSOC_AVERAGE);
	}
	assocRule->set_path( append_path );
	++numOfAssocRules;
	assocRulesList.push_back(*assocRule);
	
}

/* Write the defined Association Rules into an XML file */
void ObjectAssociation::WriteRulesToXML(std::string xmlFname)
{
	TiXmlDocument doc;   
 
	//Root node
	TiXmlElement * root = new TiXmlElement( "ObjectAssociationRules" );  
	doc.LinkEndChild( root );  
	root->SetAttribute("SegmentationSource", segImageName.c_str());
	root->SetAttribute("NumberOfAssociativeMeasures", numOfAssocRules);

	TiXmlComment * comment = new TiXmlComment();
	comment->SetValue(" Definition of Association Rules between different objects " );  
	root->LinkEndChild( comment );  

	//Add the Association Rules one by one
	for(int i=0; i<numOfAssocRules; i++)
	{
		TiXmlElement *element = new TiXmlElement("AssociationRule");
		element->SetAttribute("Name",assocRulesList[i].GetRuleName().c_str());
		element->SetAttribute("Target_Image",assocRulesList[i].GetTargetFileNmae().c_str());
		element->SetAttribute("Outside_Distance",assocRulesList[i].GetOutDistance());
		element->SetAttribute("Inside_Distance",assocRulesList[i].GetInDistance());
		if(assocRulesList[i].IsUseWholeObject())
			element->SetAttribute("Use_Whole_Object","True");
		else
			element->SetAttribute("Use_Whole_Object","False");
		if(assocRulesList[i].IsUseBackgroundSubtraction())
			element->SetAttribute("Use_Background_Subtraction","True");
		else
			element->SetAttribute("Use_Background_Subtraction","False");
		if(assocRulesList[i].IsUseMultiLevelThresholding())
			element->SetAttribute("Use_MultiLevel_Thresholding","True");
		else
			element->SetAttribute("Use_MultiLevel_Thresholding","False");
		element->SetAttribute("Number_Of_Thresholds",assocRulesList[i].GetNumberOfThresholds());
		element->SetAttribute("Number_Included_In_Foreground",assocRulesList[i].GetNumberIncludedInForeground());
		switch(assocRulesList[i].GetAssocType())
		{
		case ASSOC_MIN:
			element->SetAttribute("Association_Type","MIN");
			break;
		case ASSOC_MAX:
			element->SetAttribute("Association_Type","MAX");
			break;
		case ASSOC_TOTAL:
			element->SetAttribute("Association_Type","TOTAL");
			break;
		case ASSOC_AVERAGE:
			element->SetAttribute("Association_Type","AVERAGE");
			break;
		case ASSOC_SURROUNDEDNESS:
			element->SetAttribute("Association_Type","SURROUNDEDNESS");
			break;
		case ASSOC_DIST_OBJECT:
			element->SetAttribute("Association_Type","DIST_OBJECT");
			break;
		default:
			element->SetAttribute("Association_Type","AVERAGE");
		}

		root->LinkEndChild(element);
	}
	
	doc.SaveFile( xmlFname.c_str() );
}

/* Read the defined Association Rules from an XML file */
int ObjectAssociation::ReadRulesFromXML(std::string xmlFname)
{
	//open the xml file
	TiXmlDocument doc;
	if ( !doc.LoadFile( xmlFname.c_str() ) )
	{
		std::cout<<"Unable to load XML File";
		return 0;
	}

	//get the root node (exit on error)
	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "ObjectAssociationRules" ) != 0 )
	{
		std::cout<<"Incorrect XML root Element";		
		return 0;
	}
		
	//get the segmentation image name and the number of association rulese (exit on error)
	segImageName = rootElement->Attribute("SegmentationSource");
	numOfAssocRules = atoi(rootElement->Attribute("NumberOfAssociativeMeasures"));
	if(numOfAssocRules<=0)
	{
		std::cout<<"Incorrect number of association rules";		
		return 0;
	}
	

	//now get the rules one by one
	TiXmlElement* parentElement = rootElement->FirstChildElement();
	while (parentElement)
	{		
		const char * parent = parentElement->Value();

		if ( strcmp( parent, "AssociationRule" ) != 0 )
		{
			std::cout<<"The XML file format is incorrect!";
			return 0;
		}
		//get the attributes one by one
		TiXmlAttribute *atrib = parentElement->FirstAttribute();
		if(strcmp(atrib->Name(),"Name")!=0)
		{
			std::cout<<"First attribute in an Association rule must be its name";
			return 0;
		}
		int numAttribs = 9;
		//an association rule object
		AssociationRule *assocRule = new AssociationRule(atrib->ValueStr());		
		while(atrib)
		{							
			if(strcmp(atrib->Name(),"Name")==0) 
			{
				//Processed above
				//assocRulesList[i].SetRuleName(atrib->ValueStr());
				numAttribs--;
			}
			else if(strcmp(atrib->Name(),"Target_Image")==0)
			{
				assocRule->SetTargetFileNmae(atrib->ValueStr());
			}
			else if(strcmp(atrib->Name(),"Outside_Distance")==0)
			{
				assocRule->SetOutDistance(atoi(atrib->ValueStr().c_str()));
				numAttribs--;
			}
			else if(strcmp(atrib->Name(),"Inside_Distance")==0)
			{
				assocRule->SetInDistance(atoi(atrib->ValueStr().c_str()));
				numAttribs--;
			}
			else if(strcmp(atrib->Name(),"Use_Whole_Object")==0)
			{
				const char* V = atrib->Value();
				if(strcmp(V,"True")==0)
					assocRule->SetUseWholeObject(true);
				else
					assocRule->SetUseWholeObject(false);

				numAttribs--;
			}
			else if(strcmp(atrib->Name(),"Use_Background_Subtraction")==0)
			{
				const char* V = atrib->Value();
				if(strcmp(V,"True")==0)
					assocRule->SetUseBackgroundSubtraction(true);
				else
					assocRule->SetUseBackgroundSubtraction(false);

				numAttribs--;
			}
			else if(strcmp(atrib->Name(),"Use_MultiLevel_Thresholding")==0)
			{
				const char* V = atrib->Value();
				if(strcmp(V,"True")==0)
					assocRule->SetUseMultiLevelThresholding(true);
				else
					assocRule->SetUseMultiLevelThresholding(false);

				numAttribs--;
			}
			else if(strcmp(atrib->Name(),"Number_Of_Thresholds")==0)
			{
				assocRule->SetNumberOfThresholds(atoi(atrib->ValueStr().c_str()));
				numAttribs--;
			}
			else if(strcmp(atrib->Name(),"Number_Included_In_Foreground")==0)
			{
				assocRule->SetNumberIncludedInForeground(atoi(atrib->ValueStr().c_str()));
				numAttribs--;
			}
			else if(strcmp(atrib->Name(),"Association_Type")==0)
			{
				const char* V = atrib->Value();
				if(strcmp(V,"MIN")==0)
					assocRule->SetAssocType(ASSOC_MIN);
				else if(strcmp(V,"MAX")==0)
					assocRule->SetAssocType(ASSOC_MAX);
				else if(strcmp(V,"TOTAL")==0)
					assocRule->SetAssocType(ASSOC_TOTAL);
				else if(strcmp(V,"SURROUNDEDNESS")==0)
					assocRule->SetAssocType(ASSOC_SURROUNDEDNESS);
				else if(strcmp(V,"DIST_OBJECT")==0)
					assocRule->SetAssocType(ASSOC_DIST_OBJECT);
				else 
					assocRule->SetAssocType(ASSOC_AVERAGE);

				numAttribs--;
			}			
			//go to the next attribute
			atrib = atrib->Next();
		}

		//if you have a wrong number of attributes, then exit
		if(numAttribs != 0)
		{
			std::cout<<"The XML file has incorrect format";
			return 0;
		}

		//add the assocuation rule to the association rules list
		assocRulesList.push_back(*assocRule);

		//go to the next association rule element
		parentElement = parentElement->NextSiblingElement();
	} 
	
	return 1;
}

/* Write the computer Associative features of all the objects to an XML file */
void ObjectAssociation::WriteAssociativeFeaturesToXML(std::string xmlFname)
{//*************************************************************************Add code to write EC_array
	TiXmlDocument doc;   
 
	//Root node
	TiXmlElement * root = new TiXmlElement( "ObjectAssociationRules" );  
	doc.LinkEndChild( root );  
	root->SetAttribute("SegmentationSource", segImageName.c_str());
	root->SetAttribute("NumberOfAssociativeMeasures", numOfAssocRules);
	root->SetAttribute("NumberOfObjects", numOfLabels);

	TiXmlComment * comment = new TiXmlComment();
	comment->SetValue(" List Of Associative Measurements " );  
	root->LinkEndChild( comment );  

	//go over the objects one by one	
	for(unsigned int jj=1; jj<labelsList.size(); jj++)
	{
		//if the object is invalid, then just ignore it
		//if(invalidObjects[j]==1)
		//	continue;
		int j = labelsList[jj];
		TiXmlElement *element = new TiXmlElement("Object");
		element->SetAttribute("Type",objectType.c_str());		
		std::stringstream out1;
		out1 << std::setprecision(2) << std::fixed << j;//+1;
		element->SetAttribute("ID",out1.str());		
		//Add the Associative measurements one by one
		for(int i=0; i<numOfAssocRules; i++)
		{	
			TiXmlElement *element2 = new TiXmlElement("Association");
			std::stringstream out2;
			out2 << std::setprecision(2) << std::fixed << assocMeasurementsList[i][jj-1];
			//element2->SetAttribute(assocRulesList[i].GetRuleName(),out2.str());			
			element2->SetAttribute("Name",assocRulesList[i].GetRuleName());			
			element2->SetAttribute("Value",out2.str());
			element->LinkEndChild(element2);
		}
		root->LinkEndChild(element);
	}	
	
	doc.SaveFile( xmlFname.c_str() );
}

void ObjectAssociation::PrintSelf()
{
	//Print the header
	std::cout<<"\n---------------------------------------------------------------\n";	
	std::cout<<"Object Association Rules\n"; 
	std::cout<<"SegmentationSource "<<segImageName.c_str()<<std::endl;
	std::cout<<"NumberOfAssociativeMeasures "<<numOfAssocRules<<std::endl;
	std::cout<<".................................................................\n";	
	//Print the Association Rules one by one
	for(int i=0; i<numOfAssocRules; i++)
	{		
		std::cout<<"Name("<<i<<"): "<<assocRulesList[i].GetRuleName().c_str()<<std::endl;
		std::cout<<"Target_Image("<<i<<"): "<<assocRulesList[i].GetTargetFileNmae().c_str()<<std::endl;
		std::cout<<"Outside_Distance("<<i<<"): "<<assocRulesList[i].GetOutDistance()<<std::endl;
		std::cout<<"Inside_Distance("<<i<<"): "<<assocRulesList[i].GetInDistance()<<std::endl;		
		if(assocRulesList[i].IsUseWholeObject())
			std::cout<<"Use_Whole_Object("<<i<<"): True"<<std::endl;			
		else
			std::cout<<"Use_Whole_Object("<<i<<"): False"<<std::endl;
		if(assocRulesList[i].IsUseBackgroundSubtraction())
			std::cout<<"Use_Background_Subtraction("<<i<<"): True"<<std::endl;			
		else
			std::cout<<"Use_Background_Subtraction("<<i<<"): False"<<std::endl;
		if(assocRulesList[i].IsUseMultiLevelThresholding())
			std::cout<<"Use_MultiLevel_Thresholding("<<i<<"): True"<<std::endl;			
		else
			std::cout<<"Use_MultiLevel_Thresholding("<<i<<"): False"<<std::endl;
		std::cout<<"Number_Of_Thresholds("<<i<<"): "<<assocRulesList[i].GetNumberOfThresholds()<<std::endl;		
		std::cout<<"Number_Included_In_Foreground("<<i<<"): "<<assocRulesList[i].GetNumberIncludedInForeground()<<std::endl;		
		switch(assocRulesList[i].GetAssocType())
		{
		case ASSOC_MIN:
			std::cout<<"Association_Type("<<i<<"): MIN"<<std::endl;			
			break;
		case ASSOC_MAX:
			std::cout<<"Association_Type("<<i<<"): MAX"<<std::endl;
			break;
		case ASSOC_TOTAL:
			std::cout<<"Association_Type("<<i<<"): TOTAL"<<std::endl;
			break;
		case ASSOC_AVERAGE:
			std::cout<<"Association_Type("<<i<<"): AVERAGE"<<std::endl;
			break;
		case ASSOC_SURROUNDEDNESS:
			std::cout<<"Association_Type("<<i<<"): SURROUNDEDNESS"<<std::endl;
			break;
		case ASSOC_DIST_OBJECT:
			std::cout<<"Association_Type("<<i<<"): DIST_OBJECT"<<std::endl;
			break;
		default:
			std::cout<<"Association_Type("<<i<<"): AVERAGE"<<std::endl;
		}		
	}
	std::cout<<"---------------------------------------------------------------\n";
}

} //end namespace ftk

