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


#include "face.h"

inline float DistanceBetweenTwoPoints(const Vec3f &p1, const Vec3f &p2) {
  Vec3f v(p1,p2);
  return v.Length();
}

float AreaOfTriangle(float a, float b, float c) {
  // Area of Triangle =  (using Heron's Formula)
  //  sqrt[s*(s-a)*(s-b)*(s-c)]
  //    where s = (a+b+c)/2
  float s = (a+b+c) / (float)2;
  return sqrt(s*(s-a)*(s-b)*(s-c));
}

float Face::getArea() const {
  Vec3f a = (*this)[0]->get();
  Vec3f b = (*this)[1]->get();
  Vec3f c = (*this)[2]->get();
  //Vec3f d = (*this)[3]->get();
  return AreaOfTriangle(DistanceBetweenTwoPoints(a,b),DistanceBetweenTwoPoints(a,c),DistanceBetweenTwoPoints(b,c));
}

