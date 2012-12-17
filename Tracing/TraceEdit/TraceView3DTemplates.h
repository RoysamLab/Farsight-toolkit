
template<typename TINPUT >
typename TINPUT::Pointer cropImages( typename TINPUT::Pointer inputImage, int x, int y, int z, int xTile, int yTile, int zTile,int xSize,int ySize, int zSize)
{
	typename TINPUT::IndexType start;
	start[0] = ((x - xTile/2)>0) ? (x - xTile/2):0;
	start[1] = ((y - yTile/2)>0) ? (y - yTile/2):0;
	start[2] = ((z - zTile/2)>0) ? (z - zTile/2):0;
	typename TINPUT::SizeType size;
	size[0] = ((x+xTile/2)>xSize) ? (xSize-x) : xTile/2;
	size[1] = ((y+yTile/2)>ySize) ? (ySize-y) : yTile/2; 
	size[2] = ((z+zTile/2)>zSize) ? (zSize-z) : zTile/2; 
	typename TINPUT::SizeType size2;
	size2[0] = ((x-xTile/2)>0) ? xTile/2 : x; 
	size2[1] = ((y-yTile/2)>0) ? yTile/2 : y;
	size2[2] = ((z-zTile/2)>0) ? zTile/2 : z;
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

template <typename T>
int writeImage(typename T::Pointer im, const char* filename)
{
// 	printf("Writing %s ... \n",filename);
	std::cout << std::endl << "Writing ... " << filename;
	typedef typename itk::ImageFileWriter<T> WriterType;

	typename WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(filename);
	writer->SetInput(im);
	try
	{
		writer->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}
	itk::Size<3> inputImageSize = im->GetLargestPossibleRegion().GetSize();
	std::cout<<" done: Image size: "<<inputImageSize[0]<<"x"<<inputImageSize[1]<<"x"<<inputImageSize[2];
	return EXIT_SUCCESS;
}

