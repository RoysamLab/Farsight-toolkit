#ifndef _BOUNDING_BOX_H_
#define _BOUNDING_BOX_H_

#include <assert.h>
#include "vectors.h"
#include "utils.h"


// ====================================================================
// ====================================================================

class BoundingBox {

public:

  // CONSTRUCTOR & DESTRUCTOR
  BoundingBox(Vec3f _min, Vec3f _max) {
    Set(_min,_max); }
  ~BoundingBox() {}

  // ACCESSORS
  void Get(Vec3f &_min, Vec3f &_max) const {
    _min = min;
    _max = max; }
  Vec3f getMin() const { return min; }
  Vec3f getMax() const { return max; }

  void getCenter(Vec3f &c) const {
    c = max; 
    c -= min;
    c *= 0.5f;
    c += min;
  }

  double maxDim() const {
    double x = max.x() - min.x();
    double y = max.y() - min.y();
    double z = max.z() - min.z();
    return max3(x,y,z);
  }

  // MODIFIERS
  void Set(BoundingBox *bb) {
    assert (bb != NULL);
    min = bb->min;
    max = bb->max; }
  void Set(Vec3f _min, Vec3f _max) {
    assert (min.x() <= max.x() &&
	    min.y() <= max.y() &&
	    min.z() <= max.z());
    min = _min;
    max = _max; }
  void Extend(const Vec3f v) {
    min = Vec3f(min2(min.x(),v.x()),
		min2(min.y(),v.y()),
		min2(min.z(),v.z()));
    max = Vec3f(max2(max.x(),v.x()),
		max2(max.y(),v.y()),
		max2(max.z(),v.z())); }
  void Extend(BoundingBox *bb) {
    assert (bb != NULL);
    Extend(bb->min);
    Extend(bb->max); }

  // DEBUGGING 
  void Print(const char *s="") const {
  /*  printf ("BOUNDING BOX %s: %f %f %f  -> %f %f %f\n", s,
            min.x(),min.y(),min.z(),
            max.x(),max.y(),max.z());*/ }
  void paint() const;

private:
  BoundingBox() { assert(0); } // don't use this constructor

  // REPRESENTATION
  Vec3f min;
  Vec3f max;
};

// ====================================================================
// ====================================================================

#endif
