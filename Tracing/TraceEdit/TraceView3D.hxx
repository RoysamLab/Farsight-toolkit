#ifndef _TraceView3DTemplates_h_
#define _TraceView3DTemplates_h_

template<typename TINPUT >
typename TINPUT::Pointer cropImages( typename TINPUT::Pointer inputImage, int x, int y, int z)
{
	typename TINPUT::IndexType start;
	start[0] = ((x - _xTile/2)>0) ? (x - _xTile/2):0;
	start[1] = ((y - _yTile/2)>0) ? (y - _yTile/2):0;
	start[2] = ((z - _zTile/2)>0) ? (z - _zTile/2):0;
	typename TINPUT::SizeType size;
	size[0] = ((x+_xTile/2)>_xSize) ? (_xSize-x) : _xTile/2;
	size[1] = ((y+_yTile/2)>_ySize) ? (_ySize-y) : _yTile/2; 
	size[2] = ((z+_zTile/2)>_zSize) ? (_zSize-z) : _zTile/2; 
	typename TINPUT::SizeType size2;
	size2[0] = ((x-_xTile/2)>0) ? _xTile/2 : x; 
	size2[1] = ((y-_yTile/2)>0) ? _yTile/2 : y;
	size2[2] = ((z-_zTile/2)>0) ? _zTile/2 : z;
	size[0] = size2[0] + size[0];
	size[1] = size2[1] + size[1];
	size[2] = size2[2] + size[2];

	typename TINPUT::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);
	
	typedef itk::RegionOfInterestImageFilter< TINPUT, TINPUT > ROIFilterType_5;
	typename ROIFilterType_5::Pointer ROIfilter3 = ROIFilterType_5::New();
	ROIfilter3->SetRegionOfInterest(desiredRegion);
	ROIfilter3->SetInput(inputImage);
// #pragma omp critical
	ROIfilter3->Update();
	
	typedef itk::ImageDuplicator< TINPUT >  DuplicatorType_5;
	typename DuplicatorType_5::Pointer LabelDuplicator = DuplicatorType_5::New();
	LabelDuplicator->SetInputImage(ROIfilter3->GetOutput());
// #pragma omp critical
	LabelDuplicator->Update();
	
	return LabelDuplicator->GetOutput();

}
#endif