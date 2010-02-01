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
/**************************************************************************  
 *  Volume Density Iso- Anisotropic Diffusion
 *   Author: Xiaosong Yuan, RPI
 *  Created on Nov. 16, 2005         
 *************************************************************************/
#ifndef __mdlAnisoDiffuse_h
#define __mdlAnisoDiffuse_h

#include "mdlTypes.h"

#include <string>
#include <vector>
#include <iostream>

#include "itkVector.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"


namespace mdl
{

class AnisoDiffuse
{
public:
	AnisoDiffuse();
	void SetDebug(bool inp = true){ debug = inp; };
	void SetOverwrite(bool inp = true){ overwriteInput = inp; };
	void SetUseIso(bool inp = true){ useIso = inp; };
	void SetInput(ImageType::Pointer inImage);
	bool Update();
	ImageType::Pointer GetOutput();

private:
	struct mdlVector{float xd; float yd; float zd;};

	//Images
	ImageType::Pointer m_inputImage;
	ImageType::Pointer m_outputImage;

	bool debug;				//If debug is true, process in steps and print stuff
	bool overwriteInput;	//Replace the original image with processed one
	bool useIso;			//true Isotropic, false = Anisotropic

	int timesDiffuse;

	bool DoMDLAniso();
};

}  // end namespace mdl

#endif
