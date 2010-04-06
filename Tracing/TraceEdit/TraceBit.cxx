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

#include "TraceBit.h"

TraceBit::TraceBit()
{
}
TraceBit::TraceBit(const TraceBit& T)
{
	this->x = T.x;
	this->y = T.y;
	this->z = T.z;
	this->r = T.r;
	this->I = T.I; 		
	this->id = T.id;
    this->marker = T.marker;
	this->dx = T.dx;
	this->dy = T.dy;
	this->dz = T.dz;
}

TraceBit::~TraceBit()
{
}

//0 = x, 1 = y , 2 = z
double TraceBit::GetCoordinateByRef(int asint){
	switch(asint){
		case 0: return this->x; break;
		case 1: return this->y; break;
		case 2: return this->z; break;
	}
	return 0.0;
}
void TraceBit::setCoordinateByRef(int asint,double value){
	switch(asint){
		case 0: this->x = value; break;
		case 1:  this->y = value; break;
		case 2:  this->z = value; break;
	}
}
void TraceBit::Print(std::ostream &c)
{
  c<<"\t\tTraceBit:"<<std::endl;
  c<<"\t\tx:"<<x<<" y:"<<y<<" z:"<<z<<" r:"<<r<<" id:"<<id<<" marker:"<<marker;
}
