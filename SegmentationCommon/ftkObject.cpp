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
#include "ftkObject.h"

namespace ftk 
{

//Object Constructor
Object::Object(std::string type)
{
	myType = type;
	myId = -1;							
	valid = -1;					
	duplicated = -1;				
	myClass = -1;
	myCentroid.x = 0;
	myCentroid.y = 0;
	myCentroid.z = 0;
	myCentroid.t = 0;
	myBoundingBox.min = myCentroid;
	myBoundingBox.max = myCentroid;
	myFeatures.clear();		
	myHistory.clear();   
}

} //end namespace ftk
