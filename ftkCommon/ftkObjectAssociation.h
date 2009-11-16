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
  Module:    $RCSfile: ftkObjectAssociation.h,v $
  Language:  C++
  Date:      $Date: 2008/11/26 12:30:00 $
  Version:   $Revision: 1 $
 
=========================================================================*/
#ifndef __ftkObjectAssociation_h
#define __ftkObjectAssociation_h

#include <tinyxml/tinyxml.h> //needed to read from and write to XML files
#include <ftkObject.h> //don't remove this line
//need to include some other headers here
//I think ftkLableImageFeatures.h should be included


namespace ftk
{ 
/** \class AssociationRule
 *  \brief To define an association rule between two objects
 * \author Yousef Al-Kofahi, Rensselear Polytechnic Institute (RPI), Troy, NY USA
 */
enum AssociationType {
      ASSOC_MIN = 1,
      ASSOC_MAX = 2,
      ASSOC_TOTAL = 3,
      ASSOC_AVERAGE = 4
};

class AssociationRule
{
public:		
	/* Constructor */
	AssociationRule(std::string ruleName);	
	/* Set and Get functions */
	void SetRuleName(std::string rName) {ruleName = rName;};
	void SetSegmentationFileNmae(std::string name){ segFileName = name; };
	void SetTargetFileNmae(std::string name){ targFileName = name; };
	void SetOutDistance(int dist){outsideDistance = dist;};
	void SetInDistance(int dist){insideDistance = dist;};
	void SetUseWholeObject(bool useAll){useWholeObject = useAll;};
	void SetAssocType(AssociationType tp ){assocType = tp;};
	std::string GetRuleName() { return ruleName;};
	std::string GetSegmentationFileNmae(){ return segFileName;};
	std::string GetTargetFileNmae(){ return targFileName; };
	int GetOutDistance(){ return outsideDistance; };
	int GetInDistance(){ return insideDistance; };
	bool IsUseWholeObject() {return useWholeObject; };
	AssociationType GetAssocType() {return assocType; };
private:
	std::string ruleName;
	std::string segFileName;
	std::string targFileName;
	int outsideDistance;
	int insideDistance;
	bool useWholeObject;
	AssociationType assocType;
}; // end AssociationRule


/** \class ObjectAssociation
 *  \brief To define a set of association rules between different objects and to compute the associative measures. 
 *  Inherit from this class to extend the functionality for a specific Segmentation type 
 * \author Yousef Al-Kofahi, Rensselear Polytechnic Institute (RPI), Troy, NY USA
 */
class ObjectAssociation
{
public:
	/* Contsructor */
	ObjectAssociation(std::string segImageName, int numOfAssocRules);
	/* Used to add a new association rule to the rules list */
	void AddAssociation(std::string ruleName,std::string targFileName, int outsideDistance, int insideDistance,	bool useAllObject, int assocType);
	/* I/O	*/
	void WriteRulesToXML(std::string xmlFname);
	int ReadRulesFromXML(std::string xmlFname);
	void WriteAssociativeFeaturesToXML(std::string xmlFname);

	/* used for validation. In other words, it can be used to make sure that the XML reader works fine */
	void PrintSelf(); 

	/* Get the values of private member variables */
	std::string GetSegImgName() {return segImageName;};
	int GetNumofAssocRules() {return numOfAssocRules;};

	/* Get the features list*/
	float** GetAssocFeaturesList() {return assocMeasurementsList;};
	
protected:
	/* This is the list of association rules */
	std::vector<AssociationRule> assocRulesList;
	/* This is the list of associative measurements */
	float** assocMeasurementsList;
	/* This is the total number of objects/labels */
	int numOfLabels;
	/* This is the object type (ex. nucleus, spine, etc..) */
	std::string objectType;

	//added by Yousef on 10-18-2009
	//unsigned short* invalidObjects;

	//added by Yousef on 10/20/2009
	std::vector< unsigned short > labelsList;

private:
	std::string segImageName;
	int numOfAssocRules;	
			
}; // end ObjectAssociation

}  // end namespace ftk

#endif	// end __ftkObjectAssociation_h
