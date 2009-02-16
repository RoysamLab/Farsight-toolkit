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
Object::Object(string type)
{
	myType = type;
	myId = -1;							
	valid = -1;					
	duplicated = -1;				
	myClass = -1;					
	myCenters.clear();		
	myBounds.clear();			
	myFeatures.clear();		
	myHistory.clear();   
}

} //end namespace ftk
