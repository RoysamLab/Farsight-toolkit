namespace ftk
{
template <typename InputPixelType, typename LabelPixelType>  void
	ProjectProcessor::ComputeMontageIntrinsicFeatures( int nucChannel )
{
	typedef typename itk::Image< InputPixelType, 3 > InputImageType;
	//The computation engine assumes 8-bit images and till it is fixed
	typedef typename itk::Image< unsigned char,  3 > InputImageType1;
	typedef typename itk::Image< LabelPixelType, 3 > LabelImageType1;
	typedef typename itk::RegionOfInterestImageFilter
				< InputImageType, InputImageType1 > ROIFilterType;
	typedef typename itk::RegionOfInterestImageFilter
				< LabelImageType1, LabelImageType1 > LabelROIFilterType;
	typedef typename itk::ImageRegionConstIterator< InputImageType1 > ConstIteratorType;
	typedef typename itk::ImageRegionConstIterator< LabelImageType1 > UIntConstIteratorType;
	typedef itk::ImageRegionConstIterator< LabelImageType > UShrtConstIteratorType;
	typedef typename itk::ImageRegionIteratorWithIndex< InputImageType1 > IteratorType;
	typedef typename itk::ImageRegionIteratorWithIndex< LabelImageType1 > UIntIteratorType;
	typedef itk::ImageRegionIteratorWithIndex< LabelImageType > UShrtIteratorType;
	typedef typename itk::LabelStatisticsImageFilter< InputImageType, LabelImageType1 > 
								LabelStatisticsImageFilterType;
	typedef typename LabelStatisticsImageFilterType::ValidLabelValuesContainerType
									ValidLabelValuesType;

	typename LabelImageType1::Pointer LabelImage = outputImage->
		GetItkPtr< typename LabelImageType1::PixelType >
			( 0, 0,	    ftk::Image::DEFAULT );
	typename InputImageType::Pointer  InputImage = inputImage->
		GetItkPtr< typename InputImageType::PixelType >
			( 0, nucChannel, ftk::Image::DEFAULT );

	std::vector< std::string > FeatureNames;
	typename LabelStatisticsImageFilterType::Pointer LabelStatisticsImageFilter =
							LabelStatisticsImageFilterType::New();
	LabelStatisticsImageFilter->SetLabelInput( LabelImage );
	LabelStatisticsImageFilter->SetInput( InputImage );
	LabelStatisticsImageFilter->UseHistogramsOff();
	try
	{
		LabelStatisticsImageFilter->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Error in montage compute stats:"<< excp << std::endl;
	}

	std::vector< std::string > SNames;
	//Compute features for one image and get the names of the features that are turned on
	{
		LabelPixelType Index = 1;
		while( LabelStatisticsImageFilter->GetCount(Index)<10 ) ++Index;
		typename LabelStatisticsImageFilterType::BoundingBoxType BoundBox;
		BoundBox = LabelStatisticsImageFilter->GetBoundingBox( Index );
		typename LabelImageType1::IndexType CCStart;
		CCStart[0] = BoundBox[0];
		CCStart[1] = BoundBox[2];
		CCStart[2] = BoundBox[4];
		typename LabelImageType1::SizeType CCSize;
		CCSize[0] = BoundBox[1]-CCStart[0]+1;
		CCSize[1] = BoundBox[3]-CCStart[1]+1;
		CCSize[2] = BoundBox[5]-CCStart[2]+1;
		typename LabelImageType1::RegionType CCRegion;
		CCRegion.SetSize ( CCSize  );
		CCRegion.SetIndex( CCStart );
		LabelImageType::Pointer CroppedImage = LabelImageType::New();
		CroppedImage->SetRegions( CCRegion );
		CroppedImage->Allocate();
		CroppedImage->FillBuffer(0);
		CroppedImage->Update();

		typename LabelROIFilterType::Pointer LabelCropFilter = LabelROIFilterType::New();
		LabelCropFilter->SetInput( LabelImage );
		LabelCropFilter->SetRegionOfInterest( CCRegion );
		try{
			LabelCropFilter->Update();
		}
		catch( itk::ExceptionObject & excp )
		{
			std::cerr << "Error in crop filter while computing features:"
				  << excp << std::endl;
		}

		UIntConstIteratorType CroppedIter ( LabelCropFilter->GetOutput(),
			LabelCropFilter->GetOutput()->GetLargestPossibleRegion() );
		UShrtIteratorType CopyCropIter( CroppedImage,
			CroppedImage->GetLargestPossibleRegion() );
			CroppedIter.GoToBegin();
		for( CopyCropIter.GoToBegin(); !CopyCropIter.IsAtEnd();
		     ++CopyCropIter,++CroppedIter )
			if( CroppedIter.Get()==Index )
				CopyCropIter.Set( 1 );

		ftk::IntrinsicFeatureCalculator *iCalc = new ftk::IntrinsicFeatureCalculator();
		ftk::Image::Pointer IImage = ftk::Image::New();
		ftk::Image::Pointer OImage = ftk::Image::New();
		std::vector<unsigned char> color;
		color.assign(3,255);
		IImage->AppendChannelFromData3D( (void*)CroppedImage->GetBufferPointer(),
				 itk::ImageIOBase::UCHAR, sizeof(unsigned char),
				 CCSize[0], CCSize[1], CCSize[2], "nuc", color, true); 
		OImage->AppendChannelFromData3D( (void*)CroppedImage->GetBufferPointer(),
				itk::ImageIOBase::USHORT, sizeof(unsigned short),
				CCSize[0], CCSize[1], CCSize[2], "nuc", color, true); 
		iCalc->SetInputImages(IImage,OImage,0,0);
		if(definition->intrinsicFeatures.size() > 0)
			iCalc->SetFeaturesOn( GetOnIntrinsicFeatures() );
		vtkSmartPointer< vtkTable > CCtable = iCalc->Compute();	//Create a new temp table
		for( vtkIdType i=0; i<CCtable->GetNumberOfColumns(); ++i )
			SNames.push_back( CCtable->GetColumnName(i) );
		delete iCalc;
	}

	typename std::vector<LabelPixelType> labelsList;
	for( typename ValidLabelValuesType::const_iterator
		vIt=LabelStatisticsImageFilter->GetValidLabelValues().begin();
		vIt != LabelStatisticsImageFilter->GetValidLabelValues().end();	++vIt )
	{
		if ( LabelStatisticsImageFilter->HasLabel(*vIt) )
		{
			if( *vIt ) labelsList.push_back(*vIt); //Skips 0 if present
		}
		else
			std::cout<<"The id "<<*vIt<<" is present but was discarded\n";
	}
	std::vector< std::vector<double> > FVector;
	FVector.resize( labelsList.size() );
	for( LabelPixelType LabelIndex=0; LabelIndex<labelsList.size(); ++LabelIndex )
			FVector.at(LabelIndex).resize(SNames.size());
	std::sort( labelsList.begin(), labelsList.end() );
	std::cout<<"Found "<<labelsList.size()<<" labels for computing features\n"<<std::flush;

#ifdef _OPENMP
	//Use 95% of the cores by default n save a little for the OS
	if(!numThreadsSet) n_thr = 0.95*omp_get_max_threads();
	itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);
	std::cout<<"Using "<<n_thr<<" threads\n"<<std::flush;
#if _OPENMP > 200805L
	omp_set_max_active_levels(1);
#endif
	#pragma omp parallel for num_threads(n_thr)
#if _OPENMP > 200805L
	for( LabelPixelType i=0; i<labelsList.size(); ++i )
#else
	for( itk::IndexValueType i=0; i<labelsList.size(); ++i )
#endif
#else
	for( LabelPixelType i=0; i<labelsList.size(); ++i )
#endif
	{
		LabelPixelType LabelIndex = labelsList.at(i);
		typename LabelStatisticsImageFilterType::BoundingBoxType BoundBox;
		BoundBox = LabelStatisticsImageFilter->GetBoundingBox( LabelIndex );
		typename InputImageType1::IndexType Start;
		typename LabelImageType::IndexType CCStart;
		CCStart[0] = Start[0] = BoundBox[0];
		CCStart[1] = Start[1] = BoundBox[2];
		CCStart[2] = Start[2] = BoundBox[4];
		InputImageType1::SizeType Size;
		LabelImageType::SizeType CCSize;
		CCSize[0] = Size[0] = BoundBox[1]-BoundBox[0]+1;
		CCSize[1] = Size[1] = BoundBox[3]-BoundBox[2]+1;
		CCSize[2] = Size[2] = BoundBox[5]-BoundBox[4]+1;
		InputImageType1::RegionType CroppedRegion;
		CroppedRegion.SetSize ( Size  );
		CroppedRegion.SetIndex( Start );
		LabelImageType::RegionType CCRegion;
		CCRegion.SetSize ( CCSize  );
		CCRegion.SetIndex( CCStart );

		//Need these images copied and the fudge for the centroids
		//because the crop filter does not preserve the origin
		LabelImageType::Pointer CroppedImage = LabelImageType::New();
		CroppedImage->SetRegions( CCRegion );
		CroppedImage->Allocate();
		CroppedImage->FillBuffer(0);
		CroppedImage->Update();
		InputImageType1::Pointer CIImage = InputImageType1::New();
		CIImage->SetRegions( CroppedRegion );
		CIImage->Allocate();
		CIImage->FillBuffer(0);
		CIImage->Update();

		typename ROIFilterType::Pointer CropImageFilter = ROIFilterType::New();
		typename LabelROIFilterType::Pointer LabelCropFilter = LabelROIFilterType::New();
		LabelCropFilter->SetInput( LabelImage );
		CropImageFilter->SetInput( InputImage );
		LabelCropFilter->SetRegionOfInterest( CCRegion );
		CropImageFilter->SetRegionOfInterest( CroppedRegion );

	#pragma omp critical //May be able to elimiate this
	{
		try
		{
			CropImageFilter->Update();
			LabelCropFilter->Update();
		}
		catch( itk::ExceptionObject & excp )
		{
			std::cerr << "Error in crop filter " << excp << std::endl;
		}
	}

		UIntConstIteratorType ConstLabelIter ( LabelCropFilter->GetOutput(),
				LabelCropFilter->GetOutput()->GetLargestPossibleRegion() );
		ConstIteratorType ConstIntensityIter ( CropImageFilter->GetOutput(),
				CropImageFilter->GetOutput()->GetLargestPossibleRegion() );
		UShrtIteratorType CopyLabelIter( CroppedImage,
				CroppedImage->GetLargestPossibleRegion() );
		IteratorType CopyIntensityIter( CIImage, CIImage->GetLargestPossibleRegion() );
		ConstLabelIter.GoToBegin(); ConstIntensityIter.GoToBegin();
		CopyIntensityIter.GoToBegin();
		for( CopyLabelIter.GoToBegin(); !CopyLabelIter.IsAtEnd();
		     ++CopyLabelIter,++ConstLabelIter,++ConstIntensityIter,++CopyIntensityIter )
		{
			if( ConstLabelIter.Get()==LabelIndex )
			{
				CopyLabelIter.Set( 1 );
			}
			else
				CopyLabelIter.Set( 0 );
			CopyIntensityIter.Set( ConstIntensityIter.Get() );
		}

		ftk::Image::Pointer IImage = ftk::Image::New();
		ftk::Image::Pointer OImage = ftk::Image::New();
		std::vector<unsigned char> color;
		color.assign(3,255);
		IImage->AppendChannelFromData3D( (void*)CIImage->GetBufferPointer(),
			itk::ImageIOBase::UCHAR, sizeof(unsigned char),
			CCSize[0], CCSize[1], CCSize[2], "nuc", color, true); 
		OImage->AppendChannelFromData3D( (void*)CroppedImage->GetBufferPointer(),
			itk::ImageIOBase::USHORT, sizeof(unsigned short),
			CCSize[0], CCSize[1], CCSize[2], "nuc", color, true); 
		ftk::IntrinsicFeatureCalculator *iCalc = new ftk::IntrinsicFeatureCalculator();
		iCalc->SetInputImages(IImage,OImage,0,0);
		if(definition->intrinsicFeatures.size() > 0)
			iCalc->SetFeaturesOn( GetOnIntrinsicFeatures() );
		vtkSmartPointer< vtkTable > CCtable = iCalc->Compute();
		for( vtkIdType j=1; j<CCtable->GetNumberOfColumns(); ++j )
			FVector.at(i).at(j) = CCtable->GetValue(0,j).ToDouble();

		//Change the id, centroid_x, centroid_y and centroid_z
		FVector.at(i).at(0) = LabelIndex;
		FVector.at(i).at(1) += CCStart[0];
		FVector.at(i).at(2) += CCStart[1];
		FVector.at(i).at(3) += CCStart[2];
		delete iCalc;
	}
	std::cout<<"Features done. Populating table...\n"<<std::flush;
	table = vtkSmartPointer<vtkTable>::New();
	for( unsigned i=0; i<SNames.size(); ++i )
	{
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( SNames.at(i).c_str() );
		table->AddColumn( column );
	}
	for( LabelPixelType LabelIndex=0; LabelIndex<labelsList.size(); ++LabelIndex )
	{
		vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
		for( unsigned i=0; i<SNames.size(); ++i )
			row->InsertNextValue( vtkVariant( FVector.at(LabelIndex).at(i) ) );
		table->InsertNextRow(row);
	}

#ifdef _OPENMP
	itk::MultiThreader::SetGlobalMaximumNumberOfThreads(n_thr);
#endif

	std::cout<<"Features for the montage done with "<<table->GetNumberOfRows()<<" labels\n";

	return;
}
}
