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

/** @file pix_t.h
*   @brief class for pix
*   This is the class that represents the pix in the image
*
*   @author Maciej Wotjon
*/

#ifndef __PIX_T_H_
#define __PIX_T_H_

/** @brief pix class
*/
class pix_t
{
public:
	//coordinates
	int x_, y_, z_;

	pix_t( int x= -1, int y= -1, int z= -1)
		:x_(x)
		, y_(y)
		, z_(z)
	{
	};

	pix_t(const pix_t& p )
		:x_(p.x_)
		, y_(p.y_)
		, z_(p.z_)
	{
	};

	~pix_t(void)
	{
	};
};

#endif
