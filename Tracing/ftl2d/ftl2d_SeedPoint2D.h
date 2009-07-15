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

#ifndef SEEDPOINT2D_H
#define SEEDPOINT2D_H

class SeedPoint2D	{
	public :
	unsigned int getx() {return(x);}
	unsigned int gety() {return(y);}
	float getScale() {return(Scale);}
	float getIntensity() {return(Intensity);}
	
	void setx(unsigned int i) {x = i;}
	void sety(unsigned int i) {y = i;}
	void setScale(float i) {Scale = i;}
	void setIntensity(float i) {Intensity = i;}
	void PrintSelf();
	
	private:
	unsigned int x;
    unsigned int y;
    float Intensity;
    float Scale;
 };
 
 #endif
