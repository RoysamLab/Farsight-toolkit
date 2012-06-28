
template<typename TINPUT >
typename TINPUT::Pointer ftkMainDarpaTrace::cropImages( typename TINPUT::Pointer inputImage, int x, int y, int z)
{
	typename TINPUT::IndexType start;
	start[0] = ((x - _xTile/2)>0) ? (x - _xTile/2):0;
	start[1] = ((y - _yTile/2)>0) ? (y - _yTile/2):0;
	start[2] = ((z - _zTile/2)>0) ? (z - _zTile/2):0;
// 				if(_zSize <= _zTile)
// 					start[2] = 0;
// 				else
// 					start[2] = ((z - _zTile/2)>0) ? (z - _zTile/2):0;

// 	rawImageType_8bit::IndexType start2;
// 	start2[0] = ((x - _xTile/2)>0) ? (x - _xTile/2):0;
// 	start2[1] = ((y - _yTile/2)>0) ? (y - _yTile/2):0;
// 	start2[2] = ((z - _zTile/2)>0) ? (z - _zTile/2):0;
// // 				if(_zSize <= _zTile)
// // 					start2[2] = 0;
// // 				else
// // 					start2[2] = ((z - _zTile/2)>0) ? (z - _zTile/2):0;

// 	std::cout << std::endl << "x: " << x << ", y: " << y << ", z: " << z;
	
	typename TINPUT::SizeType size;
	size[0] = ((x+_xTile/2)>_xSize) ? (_xSize-x) : _xTile/2;
	size[1] = ((y+_yTile/2)>_ySize) ? (_ySize-y) : _yTile/2; 
	size[2] = ((z+_zTile/2)>_zSize) ? (_zSize-z) : _zTile/2; 
	
// 	std::cout << std::endl << "SIE O: " << size;
	
	typename TINPUT::SizeType size2;
	size2[0] = ((x-_xTile/2)>0) ? _xTile/2 : x; 
	size2[1] = ((y-_yTile/2)>0) ? _yTile/2 : y;
	size2[2] = ((z-_zTile/2)>0) ? _zTile/2 : z;
	
// 	std::cout << std::endl << "SIE ----O2: " << size2 <<" " <<((z-_zTile/2)>0) ? _zTile/2 : z;
	
	size[0] = size2[0] + size[0];
	size[1] = size2[1] + size[1];
	size[2] = size2[2] + size[2];
	

// 	std::cout << std::endl << "SIE ----O: " << size <<" " <<((z-_zTile/2)>0) ? _zTile/2 : z;
	
// 	if( _xTile/2 > x )
// 		size[0] = size[0] + x;
// 	else
// 		size[0] = size[0] + _xTile/2;
// 	if( _yTile/2 > y )
// 		size[1] = size[1] + y;
// 	else
// 		size[1] = size[1] + _yTile/2;
// 	if( _zTile/2 > z )
// 		size[2] = size[2] + z;
// 	else
// 		size[2] = size[2] + _zTile/2;
// 	
// 	std::cout << std::endl << "SIE -----O: " << size;
	
// 				if(_zSize <= _zTile)
// 					size[2] = _zSize;
// 				else
// 					size[2] = ((z+_zTile/2)<_zSize) ? _zTile : (_zTile/2+_zSize-z-1);

// 	rawImageType_8bit::SizeType size2;
// 	size2[0] = ((x+_xTile/2)<_xSize) ? _xTile : (_xTile/2+_xSize-x-1);
// 	size2[1] = ((y+_yTile/2)<_ySize) ? _yTile : (_yTile/2+_ySize-y-1);
// 	size2[2] = ((z+_zTile/2)<_zSize) ? _zTile : (_zTile/2+_zSize-z-1);
// // 				if(_zSize <= _zTile)
// // 					size2[2] = _zSize;
// // 				else
// // 					size2[2] = ((z+_zTile/2)<_zSize) ? _zTile : (_zTile/2+_zSize-z-1);

	typename TINPUT::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);

	
// 	std::cout << std::endl << "DESIRE REGION: " << desiredRegion << std::flush;
// 	std::cout << std::endl << "INPUT REGION: " << inputImage->GetLargestPossibleRegion() << std::flush;
// 	rawImageType_8bit::RegionType desiredRegion2;
// 	desiredRegion2.SetSize(size2);
// 	desiredRegion2.SetIndex(start2);
// 	
// 	rawImageType_16bit::RegionType desiredRegion2_16bits;
// 	desiredRegion2_16bits.SetSize(size2);
// 	desiredRegion2_16bits.SetIndex(start2);
	
// 	rawImageType_flo::RegionType desiredRegion2_float;
// 	desiredRegion2_float.SetSize(size2);
// 	desiredRegion2_float.SetIndex(start2);
	
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