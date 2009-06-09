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

using namespace std;

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
	typedef struct { int x; int y; int z; int t; } Point;
	typedef struct { Point min; Point max; } Box;
	typedef struct {string date; string description; } EditRecord;

	Object(string type);			//USE type = "null" for an empty object.

	void SetId(int id) { myId = id; };
	void SetValidity( char val) { valid = val; };
	void SetDuplicated( char dup) { duplicated = dup; };
	void SetClass(char c) { myClass = c; };
	void AddCenter( Point p ) { myCenters.push_back(p); };
	void ClearCenters(){ myCenters.clear(); };
	void AddBound( Box b ) { myBounds.push_back(b); };
	void ClearBounds(){ myBounds.clear(); };
	void SetFeatures( vector<float> f ){ myFeatures = f; };
	void AddEditRecord( EditRecord record ) { myHistory.push_back(record); };

	string GetType() { return myType; };
    int GetId () { return myId; };
	char GetValidity() { return valid; };
	char GetDuplicated() { return duplicated; };
	char GetClass() { return myClass; };
	vector<Point> GetCenters() { return myCenters; };
	vector<Box> GetBounds() { return myBounds; };
	vector<float> GetFeatures() { return myFeatures; };
	vector<EditRecord> getHistory() { return myHistory; };

private:
	string myType;					//The type of object this is 
	int myId;						//Must have an id
	char valid;						//May be marked as invalid 
	char duplicated;				//Useful for nuclei
	char myClass;					//May have class identifier
	vector<Point> myCenters;		//Should have atleast one center (splines have many)
	vector<Box> myBounds;			//Could have several boxes that determine bounds
	vector<float> myFeatures;		//May have a variable list of features
	vector<EditRecord> myHistory;   //May have a list of edits/other modifications
	
}; // end Object

}  // end namespace ftk

#endif	// end __ftkObject_h
