#ifndef COMMONTYPEDEFS_H
#define COMMONTYPEDEFS_H

typedef unsigned short                                    ImagePixelType; 
typedef itk::Image< ImagePixelType, spr_SPIMGDIMENSION>  SpineImageType;
typedef itk::Image <ImagePixelType, 2>           SpineImage2DType;

typedef   unsigned char   BinaryPixelType; 
typedef itk::Image< BinaryPixelType, spr_SPIMGDIMENSION>   BinaryImageType;
typedef itk::Image< BinaryPixelType, 2>   BinaryImage2DType;

typedef itk::PointSet <ImagePixelType, spr_SPIMGDIMENSION>	PointSetType;
typedef PointSetType::PointsContainer				PointSetContainerType;
typedef PointSetType::PointType					PointType;
typedef PointSetContainerType::Iterator					PointSetContIterType;

class TraceSegNode;
typedef std::vector<TraceSegNode*> TraceSegNodeVecType;
typedef std::vector<unsigned short> TraceIDVecType;
typedef std::map<unsigned short, PointSetType::Pointer> PointSetPointerMap;

typedef struct _tagMatrix
	{
	    double r[3][3];
}TRMatrix;


// FTK FIXME: image reading needs to be modified when part of FTK
//spine image type stuff
typedef SpineImageType::IndexType                SpineImageIndexType;
typedef std::list<SpineImageIndexType>           IndexListType;
typedef std::vector<SpineImageIndexType>		 IndexVecType;

#endif
