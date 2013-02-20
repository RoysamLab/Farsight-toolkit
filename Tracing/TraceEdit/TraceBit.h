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

#ifndef __TRACEBIT_H
#define __TRACEBIT_H

#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkAbstractArray.h"
#include "vtkVariantArray.h"

#include <iostream>
/**
 * A TraceBit has the x,y,z and id of a point on a trace
 **/
class TraceBit
{
  public:
    TraceBit();
	TraceBit(const TraceBit& T);
    ~TraceBit();
	double GetCoordinateByRef(int asint);
	void setCoordinateByRef(int asint,double value);
    double x,y,z,r, I;		//
	float dx,dy,dz;
	//unsigned char I;
    unsigned int id;
    unsigned int marker;
    void Print(std::ostream &c);

	vtkSmartPointer<vtkVariantArray> DataRow();

	bool modified;

private:
	vtkSmartPointer<vtkVariantArray> CellData;
};
#endif

