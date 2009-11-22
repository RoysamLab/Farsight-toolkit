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
#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif
#include "ftkUtils.h"

namespace ftk
{

bool FileExists(std::string filename)
{
	FILE * pFile = fopen (filename.c_str(),"r");
	if (pFile==NULL)
		return false;
	fclose (pFile);
	return true;
}

std::string NumToString(double d)
{
	std::stringstream out;
	out << std::setprecision(2) << std::fixed << d;	//Default is to use 2 decimal places
	return out.str();
}

std::string NumToString(int i)
{
	std::stringstream out;
	out << i ;	 
	return out.str();
}

std::string NumToString(double d, int p)
{
	std::stringstream out;
	out << std::setprecision(p) << std::fixed << d;	
	return out.str();
}

std::string TimeStamp()
{
	time_t rawtime;
	struct tm *timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	std::string dt = asctime(timeinfo);
	size_t end = dt.find('\n');
	dt.erase(end);
	return dt;
}

}  // end namespace ftk
