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
#ifndef __ftkObject_h
#define __ftkObject_h

#include <vector>
#include <string>

#include "ftkIntrinsicFeatures.h"

namespace ftk
{
/** \class Object
 *  \brief Information stored about individual object
 * 
 * (for valid/duplicated/class use -1 for not used) 
 *  0=false, 1=true 
 * (all classes are positive)
 */
class Object
{
public:
	enum ValidityTypes { VALID, EXCLUDED, DELETED, MERGED, SPLIT };
	typedef struct { int x; int y; int z; int t; } Point;
	typedef struct { Point min; Point max; } Box;
	typedef struct { std::string date; std::string description; } EditRecord;

	Object(std::string type);			//USE type = "null" for an empty object.

	void SetId(int id) { myId = id; };
	void SetValidity( char val) { valid = val; };
	void SetDuplicated( char dup) { duplicated = dup; };
	void SetClass(char c) { myClass = c; };
	void SetCentroid( Point p ) { myCentroid = p; };
	void SetBoundingBox( Box b ) { myBoundingBox = b; };
	void SetFeatures( std::vector<float> f ){ myFeatures = f; };
	void AddEditRecord( EditRecord record ) { myHistory.push_back(record); };

	std::string GetType() { return myType; };
    int GetId () { return myId; };
	char GetValidity() { return valid; };
	char GetDuplicated() { return duplicated; };
	char GetClass() { return myClass; };
	Point GetCentroid() { return myCentroid; };
	Box GetBoundingBox() { return myBoundingBox; };
	std::vector<float> GetFeatures() { return myFeatures; };
	std::vector<EditRecord> getHistory() { return myHistory; };

private:
	std::string myType;				//The type of object this is 
	int myId;						//Must have an id
	char valid;						//May be marked as invalid 
	char duplicated;				//Useful for nuclei
	char myClass;					//May have class identifier
	Point myCentroid;				//Should have one centroid (splines have many)
	Box myBoundingBox;				//should have one bounding box
	std::vector<float> myFeatures;	//May have a variable list of features
	std::vector<EditRecord> myHistory;   //May have a list of edits/other modifications
	
}; // end Object

}  // end namespace ftk

#endif	// end __ftkObject_h
