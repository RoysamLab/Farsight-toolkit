#include <vbl/vbl_array_2d.txx>
#include "itkAffineTransform.h"

typedef itk::AffineTransform< double, 3>   TransformType;

VBL_ARRAY_2D_INSTANTIATE(TransformType::Pointer);

