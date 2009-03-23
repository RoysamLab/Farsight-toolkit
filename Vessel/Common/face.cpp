
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

