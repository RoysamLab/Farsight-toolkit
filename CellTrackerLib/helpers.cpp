#include "helpers.h"

#if defined(_MSC_VER)
#pragma warning(disable: 4018)
#pragma warning(disable: 4996)
#pragma warning(disable: 4101)
#endif
extern int rank, npes;

#define DEBUG
#ifdef DEBUG
#define SHORT(x) (strrchr(x,'\\') ? strrchr(x,'\\')+1: x)
#define _TRACE {printf("In-%s:%d:%s:\n",SHORT(__FILE__),__LINE__,__FUNCTION__);}
#define _ETRACE {printf("Entering-%s:%d:%s:\n",SHORT(__FILE__),__LINE__,__FUNCTION__);}
#define _LTRACE {printf("Leaving-%s:%d:%s:\n",SHORT(__FILE__),__LINE__,__FUNCTION__);}
#else
#define _TRACE
#define _ETRACE
#define _LTRACE
#endif


namespace helpers{
  //---------------------------------------------------------------------------------------------------------------------
  Color2DImageType::Pointer getColorFromGrayScale(Input2DImageType::Pointer im)
  {
    Color2DImageType::Pointer out = Color2DImageType::New();
    out->SetRegions(im->GetLargestPossibleRegion());
    out->Allocate();
    typedef itk::ImageRegionIterator<Color2DImageType> Color2DIteratorType;
    Color2DIteratorType coloriter(out,out->GetLargestPossibleRegion());
    twoDIteratorType grayiter(im,im->GetLargestPossibleRegion());

    VectorPixelType pixel;
    for(grayiter.GoToBegin(),coloriter.GoToBegin();!coloriter.IsAtEnd();++coloriter,++grayiter)
    {
      unsigned char val = grayiter.Get();
      pixel[0] = val; pixel[1] = val; pixel[2] = val;
      coloriter.Set(pixel);
    }
    return out;
  }
  //---------------------------------------------------------------------------------------------------------------------
  Input2DImageType::Pointer getProjection(InputImageType::Pointer im)
  {
    Input2DImageType::Pointer output = Input2DImageType::New();
    Input2DImageType::RegionType region;
    Input2DImageType::SizeType size;
    Input2DImageType::IndexType index;
    index[0]=0;
    index[1]=0;
    size[0] = im->GetLargestPossibleRegion().GetSize()[0];
    size[1] = im->GetLargestPossibleRegion().GetSize()[1];
    region.SetSize(size);
    region.SetIndex(index);
    output->SetRegions(region);
    output->Allocate();

    SliceIteratorType inputIt(im,im->GetLargestPossibleRegion());
    LinearIteratorType outputIt(output,output->GetLargestPossibleRegion());

    inputIt.SetFirstDirection(0);
    inputIt.SetSecondDirection(1);
    outputIt.SetDirection(0);


    outputIt.GoToBegin();
    while ( ! outputIt.IsAtEnd() )
    {
      while ( ! outputIt.IsAtEndOfLine() )
      {
        outputIt.Set( itk::NumericTraits<unsigned short>::NonpositiveMin() );
        ++outputIt;
      }
      outputIt.NextLine();
    }

    inputIt.GoToBegin();
    outputIt.GoToBegin();

    while( !inputIt.IsAtEnd() )
    {
      while ( !inputIt.IsAtEndOfSlice() )
      {
        while ( !inputIt.IsAtEndOfLine() )
        {
          outputIt.Set( MAX( outputIt.Get(), inputIt.Get() ));
          ++inputIt;
          ++outputIt;
        }
        outputIt.NextLine();
        inputIt.NextLine();
      }
      outputIt.GoToBegin();
      inputIt.NextSlice();

    }
    return output;
  }

  //---------------------------------------------------------------------------------------------------------------------
  Color2DImageType::Pointer getColorProjection(ColorImageType::Pointer im)
  {
    Color2DImageType::Pointer output = Color2DImageType::New();

    Color2DImageType::RegionType region;
    Color2DImageType::SizeType size;
    Color2DImageType::IndexType index;
    index[0]=0;
    index[1]=0;
    size[0] = im->GetLargestPossibleRegion().GetSize()[0];
    size[1] = im->GetLargestPossibleRegion().GetSize()[1];
    region.SetSize(size);
    region.SetIndex(index);
    output->SetRegions(region);
    output->Allocate();

    SliceColorIteratorType inputIt(im,im->GetLargestPossibleRegion());
    LinearColorIteratorType outputIt(output,output->GetLargestPossibleRegion());

    inputIt.SetFirstDirection(0);
    inputIt.SetSecondDirection(1);
    outputIt.SetDirection(0);

    VectorPixelType nonpositivemin_vector;
    nonpositivemin_vector[0]=0;
    nonpositivemin_vector[1]=0;
    nonpositivemin_vector[2]=0;

    outputIt.GoToBegin();
    while ( ! outputIt.IsAtEnd() )
    {
      while ( ! outputIt.IsAtEndOfLine() )
      {
        outputIt.Set( nonpositivemin_vector );
        ++outputIt;
      }
      outputIt.NextLine();
    }

    inputIt.GoToBegin();
    outputIt.GoToBegin();


    while( !inputIt.IsAtEnd() )
    {
      while ( !inputIt.IsAtEndOfSlice() )
      {
        while ( !inputIt.IsAtEndOfLine() )
        {
          if(outputIt.Get().GetNorm()==0)
          {
            if(inputIt.Get().GetNorm()!=0)
            {
              outputIt.Set(inputIt.Get());
            }
          }
          ++inputIt;
          ++outputIt;
        }
        outputIt.NextLine();
        inputIt.NextLine();
      }
      outputIt.GoToBegin();
      inputIt.NextSlice();

    }
    return output;
  }

  //---------------------------------------------------------------------------------------------------------------------
  Input2DImageType::Pointer getCollage(InputImageType::Pointer im[4])
  {
    InputImageType::SizeType size = im[0]->GetLargestPossibleRegion().GetSize();
    Input2DImageType::Pointer imcollage = Input2DImageType::New();
    Input2DImageType::SizeType size2d;
    Input2DImageType::RegionType region;
    Input2DImageType::IndexType index;
    size2d[0] = size[0]*2;
    size2d[1] = size[1]*2;
    index[0]=0;
    index[1]=0;

    region.SetIndex(index);
    region.SetSize(size2d);
    imcollage->SetRegions(region);
    imcollage->Allocate(); 

    int startpoints[][2]={ {0,0},{size[0],0},{0,size[1]},{size[0],size[1]}};
    size2d[0]=size[0];
    size2d[1]=size[1];
    region.SetSize(size2d);
    for(int co = 0; co<4; co++)
    {
      index[0]=startpoints[co][0];
      index[1]=startpoints[co][1];
      region.SetIndex(index);
      Input2DImageType::Pointer proj = getProjection(im[co]);
      Const2DIteratorType inputIt(proj,proj->GetLargestPossibleRegion());
      twoDIteratorType outputIt(imcollage,region);
      for(inputIt.GoToBegin(),outputIt.GoToBegin();!inputIt.IsAtEnd();++inputIt,++outputIt)
      {
        outputIt.Set(inputIt.Get());
      }
    }

    return imcollage;

  }

  //---------------------------------------------------------------------------------------------------------------------
  unsigned char getMedianValue(NeighborhoodIteratorType it,int size)
  {
    static std::vector<unsigned char> data(size);
    int pc = 0;
    for(unsigned int i = 0; i < it.Size(); i++,pc++)
      data[pc]=it.GetPixel(i);
    std::sort(data.begin(),data.end());
    return data[data.size()/2];
  }

  //---------------------------------------------------------------------------------------------------------------------
  void unmix_median(InputImageType::Pointer im[4],InputImageType::Pointer om[4],InputImageType::Pointer assignment[4])
  {
    printf("Performing unmixing ...\n");
    InputImageType::SizeType size = im[0]->GetLargestPossibleRegion().GetSize();

    MedianFilterType::Pointer filt[4];
    IteratorType iterator[4];
    IteratorType assigniter[4];

    InputImageType::SizeType radius;
    radius[0]=1;
    radius[1]=1;
    radius[2]=1;

    InputImageType::SizeType imagesize = im[0]->GetLargestPossibleRegion().GetSize();
    InputImageType::IndexType imageindex;
    imageindex.Fill(0);
    InputImageType::RegionType region;
    region.SetSize(imagesize);
    region.SetIndex(imageindex);


    for(int counter=0; counter<4; counter++)
    {
      printf("\tPerforming median filtering on channel %d ...",counter+1);
      filt[counter]=MedianFilterType::New();
      filt[counter]->SetRadius(radius);
      filt[counter]->SetInput(im[counter]);
      filt[counter]->Update();
      om[counter]=filt[counter]->GetOutput();
      assignment[counter]=InputImageType::New();
      assignment[counter]->SetRegions(region);
      assignment[counter]->Allocate();
      assigniter[counter]=IteratorType(assignment[counter],assignment[counter]->GetLargestPossibleRegion());
      assigniter[counter].GoToBegin();
      iterator[counter]=IteratorType(om[counter],om[counter]->GetLargestPossibleRegion());
      iterator[counter].GoToBegin();
      printf(" Done %d.\n",counter+1);
    }

    //int total_voxels = size[0]*size[1]*size[2];
    int num_processed = 0;
    printf("\tComputing maximum among channels ... ");
    for(;!iterator[0].IsAtEnd();)
    {
      num_processed++;
      /*
         if(num_processed % 100 ==0)
         printf("Processed %0.2lf%% voxels\r",100.0/total_voxels*num_processed);
         */
      //		if(100.0/total_voxels*num_processed> 100)
      //			break;
      double max = -1;
      int maxpos = -1;
      for(int co = 0; co < 4; co++)
      {
        unsigned char temp = iterator[co].Value();
        if(max < temp)
        {
          max = temp;
          maxpos = co;
        }
      }
      for(int co = 0; co < 4; co++)
      {
        if(maxpos != co)
        {
          iterator[co].Set(0);
          assigniter[co].Set(0);
        }
        else
        {
          assigniter[co].Set(255);
        }
      }
      for(int co = 0; co<4; co++)
      {
        ++iterator[co];
        ++assigniter[co];
      }
    }
    printf(" Done.\n");
  }

  //---------------------------------------------------------------------------------------------------------------------
  void unmix_neighborhood(InputImageType::Pointer im[4],InputImageType::Pointer om[4])
  {

    printf("Beginning iterators\n");
    IteratorType iterator[]={IteratorType(om[0],om[0]->GetLargestPossibleRegion()),
      IteratorType(om[1],om[1]->GetLargestPossibleRegion()),
      IteratorType(om[2],om[2]->GetLargestPossibleRegion()),
      IteratorType(om[3],om[3]->GetLargestPossibleRegion()),
    };

    iterator[0].GoToBegin();
    iterator[1].GoToBegin();
    iterator[2].GoToBegin();
    iterator[3].GoToBegin();



    printf("Beginning Neighborhood iterators\n");
    NeighborhoodIteratorType::RadiusType radius;
    radius.Fill(1);
    NeighborhoodIteratorType nit[]={NeighborhoodIteratorType(radius,im[0],im[0]->GetLargestPossibleRegion()),
      NeighborhoodIteratorType(radius,im[1],im[1]->GetLargestPossibleRegion()),
      NeighborhoodIteratorType(radius,im[2],im[2]->GetLargestPossibleRegion()),
      NeighborhoodIteratorType(radius,im[3],im[3]->GetLargestPossibleRegion())
    };

    nit[0].GoToBegin();
    nit[1].GoToBegin();
    nit[2].GoToBegin();
    nit[3].GoToBegin();
    InputImageType::SizeType size = im[0]->GetLargestPossibleRegion().GetSize();
    int total_voxels = size[0]*size[1]*size[2];
    int num_processed = 0;

    for(;!iterator[0].IsAtEnd();++iterator[0],++iterator[1],++iterator[2],++iterator[3],++nit[0],++nit[1],++nit[2],++nit[3])
    {
      num_processed++;
      if(num_processed % 100 ==0)
        printf("Processed %0.2lf%% voxels\r",100.0/total_voxels*num_processed);
      //		if(100.0/total_voxels*num_processed> 100)
      //			break;
      double max = -1;
      int maxpos = -1;
      for(int co = 0; co < 4; co++)
      {
        unsigned char temp = getMedianValue(nit[co],27);
        if(max < temp)
        {
          max = temp;
          maxpos = co;
        }
      }
      for(int co = 0; co < 4; co++)
      {
        if(maxpos != co)
        {
          iterator[co].Set(0);
        }
        else
        {
          iterator[co].Set((unsigned char)max);
        }
      }
    }
  }

  //---------------------------------------------------------------------------------------------------------------------
  OutputImageType::Pointer getThresholded(InputImageType::Pointer im,int n)
  {
    printf("Performing binary thresholding ...");
    ThresholdFilterType::Pointer tfilt = ThresholdFilterType::New();
    tfilt->SetInput(im);
    tfilt->SetInsideValue(255);
    tfilt->SetOutsideValue(0);
    tfilt->SetLowerThreshold(n);
    tfilt->SetUpperThreshold(255);
    tfilt->Update();
    printf(" Done.\n");
    return tfilt->GetOutput();
  }

  //---------------------------------------------------------------------------------------------------------------------
  OutputImageType::Pointer getOtsuThresholded(InputImageType::Pointer im)
  {
    printf("Performing binary thresholding ...");
    OtsuThresholdFilterType::Pointer tfilt = OtsuThresholdFilterType::New();
    tfilt->SetInput(im);
    tfilt->SetInsideValue(0);
    tfilt->SetOutsideValue(255);
    tfilt->SetNumberOfHistogramBins(256);
    tfilt->Update();
    printf(" Done.\n");
    return tfilt->GetOutput();
  }

  //---------------------------------------------------------------------------------------------------------------------
  OutputImageType::Pointer getBinaryMedianFiltered(InputImageType::Pointer im, InputImageType::SizeType radius)
  {
    printf("Performing binary median filtering  ...");
    BinaryMedianFilterType::Pointer tfilt = BinaryMedianFilterType::New();
    tfilt->SetInput(im);
    tfilt->SetForegroundValue(255);
    tfilt->SetBackgroundValue(0);
    tfilt->SetRadius(radius);
    tfilt->Update();
    printf(" Done.\n");
    return tfilt->GetOutput();
  }

  //---------------------------------------------------------------------------------------------------------------------
  OutputImageType::Pointer getScaledFromBool(BoolImageType::Pointer im)
  {
    typedef itk::ShiftScaleImageFilter<BoolImageType,OutputImageType> ScaleFilterType;

    ScaleFilterType::Pointer scalef = ScaleFilterType::New();

    scalef->SetInput(im);
    scalef->SetScale(255);
    scalef->Update();
    return scalef->GetOutput();
  }

  //---------------------------------------------------------------------------------------------------------------------
  InputImageType::Pointer getLargeComponents(InputImageType::Pointer im, int n)
  {
    printf("Removing small connected components ...");
    //	typedef itk::Image<short int,3> LabelImageType;
    typedef itk::Image<LabelPixelType,3> LabelImageType;
    typedef itk::ConnectedComponentImageFilter<InputImageType,LabelImageType> ConnectedFilterType;
    typedef itk::RelabelComponentImageFilter<LabelImageType,LabelImageType> RelabelFilterType;

    ConnectedFilterType::Pointer cfilter = ConnectedFilterType::New();
    cfilter->SetFullyConnected(1);
    cfilter->SetInput(im);
    cfilter->Update();

    RelabelFilterType::Pointer rfilter = RelabelFilterType::New();
    rfilter->SetInput(cfilter->GetOutput());
    rfilter->InPlaceOn();

    rfilter->Update();
    RelabelFilterType::ObjectSizeInPixelsContainerType sizes = rfilter->GetSizeOfObjectsInPixels();
    //std::vector<long unsigned int> sizes = rfilter->GetSizeOfObjectsInPixels();
    int threshsize = -1;
    for(unsigned int counter=0; counter<sizes.size(); counter++)
    {
      if(sizes[counter] < (unsigned int)n)
      {
        threshsize = counter;
        break;
      }
    }

    typedef itk::BinaryThresholdImageFilter<LabelImageType,InputImageType> BinFilterType;
    BinFilterType::Pointer bfilter = BinFilterType::New();

    bfilter->SetInput(rfilter->GetOutput());
    bfilter->SetInsideValue(255);
    bfilter->SetOutsideValue(0);
    bfilter->SetLowerThreshold(1);
    if(threshsize>1)
      bfilter->SetUpperThreshold(threshsize-1);
    else
      bfilter->SetUpperThreshold(1);
    bfilter->Update();
    printf(" Done.\n");

    return bfilter->GetOutput();

  }

  //---------------------------------------------------------------------------------------------------------------------
  LabelImageType::Pointer getLargeLabels(LabelImageType::Pointer im, int n)
  {
    printf("getLargeLabels called with input n = %d\n",n);
    printf("Removing small connected components ...\n");
    //	typedef itk::Image<short int,3> LabelImageType;
    typedef itk::Image<LabelPixelType,3> LabelImageType;
    typedef itk::RelabelComponentImageFilter<LabelImageType,LabelImageType> RelabelFilterType;
    typedef itk::ScalarConnectedComponentImageFilter<LabelImageType,LabelImageType> ConnectedFilterType;

    ConnectedFilterType::Pointer cfilter = ConnectedFilterType::New();
    cfilter->SetInput(im);
    cfilter->SetFullyConnected(1);
    cfilter->SetDistanceThreshold(0);
    cfilter->Update();

    LabelIteratorType it(cfilter->GetOutput(),cfilter->GetOutput()->GetLargestPossibleRegion());
    LabelIteratorType it1(im,im->GetLargestPossibleRegion());
    for(it.GoToBegin(),it1.GoToBegin();!it.IsAtEnd(); ++it,++it1)
    {
      if(it.Get()==1)
        it.Set(0);
      if(it1.Get()==0)
        it.Set(0);
    }
    //	cfilter->SetBackgroundValue(0);


    RelabelFilterType::Pointer rfilter = RelabelFilterType::New();
    rfilter->SetInput(cfilter->GetOutput());
    rfilter->InPlaceOn();

    rfilter->Update();
    //std::vector<long unsigned int> sizes = rfilter->GetSizeOfObjectsInPixels();
    RelabelFilterType::ObjectSizeInPixelsContainerType sizes = rfilter->GetSizeOfObjectsInPixels();
    int threshsize = sizes.size()+1;// just to be safe
    for(unsigned int counter=0; counter<sizes.size(); counter++)
    {
      if(sizes[counter] < (unsigned int)n)
      {
        threshsize = counter;
        break;
      }
    }

    LabelIteratorType iter = LabelIteratorType(rfilter->GetOutput(),rfilter->GetOutput()->GetLargestPossibleRegion());
    for(iter.GoToBegin();!iter.IsAtEnd();++iter)
    {
      if(iter.Get()>threshsize)
        iter.Set(0);
    }

    //printf("%d/%d: Done removing small connected components.\n",rank,npes);
    return rfilter->GetOutput();


  }

  //---------------------------------------------------------------------------------------------------------------------
  double features_box_overlap(FeatureType &f1, FeatureType &f2)
  {
    float sx,sy,sz;
    float ex,ey,ez;
    sx = MAX(f1.BoundingBox[0],f2.BoundingBox[0]);
    sy = MAX(f1.BoundingBox[2],f2.BoundingBox[2]);
    sz = MAX(f1.BoundingBox[4],f2.BoundingBox[4]);
    ex = MIN(f1.BoundingBox[1],f2.BoundingBox[1]);
    ey = MIN(f1.BoundingBox[3],f2.BoundingBox[3]);
    ez = MIN(f1.BoundingBox[5],f2.BoundingBox[5]);

    double overlap=0;
    if((sx<ex) && (sy<ey) && (sz<ez))
    {
      overlap = (ex-sx)*(ey-sy)*(ez-sz);
    }

    return overlap;
  }

  //---------------------------------------------------------------------------------------------------------------------
  double features_diff(FeatureType &f1, FeatureType &f2,bool overlap)
  {
    //parameters
    double distance_threshold = 12;
    double distance_variance = 100;
    double intensity_variance = 100;
    double volume_variance = 100000;

    typedef ftk::IntrinsicFeatures FT;
    float spacing[3] = {0.357,0.357,2.0};
    float x1=f1.Centroid[0]*spacing[0];
    float y1=f1.Centroid[1]*spacing[1];
    float z1=f1.Centroid[2]*spacing[2];
    float x2=f2.Centroid[0]*spacing[0]; 
    float y2=f2.Centroid[1]*spacing[1]; 
    float z2=f2.Centroid[2]*spacing[2];
    double dist = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
    if(dist>distance_threshold)
    {
      //printf("I got distance of %lf and I rejected it %f %f %f %f %f %f\n",dist,f1.centroid[0],f1.centroid[1],f1.centroid[2],f2.centroid[0],f2.centroid[1],f2.centroid[2]);
      return 1e10;
    }

    double overlap_measure = 0;
    if(overlap)
    {
      float sx,sy,sz;
      float ex,ey,ez;
      sx = MAX(f1.BoundingBox[0],f2.BoundingBox[0]);
      sy = MAX(f1.BoundingBox[2],f2.BoundingBox[2]);
      sz = MAX(f1.BoundingBox[4],f2.BoundingBox[4]);
      ex = MIN(f1.BoundingBox[1],f2.BoundingBox[1]);
      ey = MIN(f1.BoundingBox[3],f2.BoundingBox[3]);
      ez = MIN(f1.BoundingBox[5],f2.BoundingBox[5]);
      if((sx<ex) && (sy<ey) && (sz<ez))
        overlap_measure = (ex-sx)*(ey-sy)*(ez-sz)/(0.5*(f1.ScalarFeatures[FT::BBOX_VOLUME]+f2.ScalarFeatures[FT::BBOX_VOLUME]));
    }
    if(overlap_measure>1)
    {
      printf("I got >1 overlap_measure\n");
    }
    double cost = (1-overlap_measure)*(1-exp(-(dist*dist/2/distance_variance+(f1.ScalarFeatures[FT::MEAN]-f2.ScalarFeatures[FT::MEAN])*(f1.ScalarFeatures[FT::MEAN]-f2.ScalarFeatures[FT::MEAN])/2/intensity_variance+(f1.ScalarFeatures[FT::VOLUME]-f2.ScalarFeatures[FT::VOLUME])*(f1.ScalarFeatures[FT::VOLUME]-f2.ScalarFeatures[FT::VOLUME])/2/volume_variance)));
    return cost;
  }

  //---------------------------------------------------------------------------------------------------------------------
  void getFeatureVectorsFarsight(LabelImageType::Pointer im, InputImageType::Pointer in_image, std::vector<ftk::IntrinsicFeatures> & feature_vector, int time, int tag)
  {


    // std::stringstream testName2;
    // testName2<<time;
    // std::string testImage2 = "/data/nicolas/test/Image2_" + testName2.str() + ".tif";
    // std::string testLabel2 = "/data/nicolas/test/Label2_" + testName2.str() + ".tif";

    // typedef itk::ImageFileWriter<LabelImageType> WriterTypeLabel2;
    //             WriterTypeLabel2::Pointer writerLabel2 = WriterTypeLabel2::New();
    //                 writerLabel2->SetInput(im);
    //                 writerLabel2->SetFileName(testLabel2.c_str());
    //                  writerLabel2->Update();
    // 
    // typedef itk::ImageFileWriter<InputImageType> WriterTypeImage2;
    //                   WriterTypeImage2::Pointer writerImage2 = WriterTypeImage2::New();
    //                    writerImage2->SetInput(in_image);
    //                     writerImage2->SetFileName(testImage2.c_str());
    //                       writerImage2->Update();


    //printf("Started feature calculation\n");
    if(im->GetLargestPossibleRegion().GetSize()[2]==1)
    {
      //convert it to a 2D Image
      LabelImageType::SizeType linsize = im->GetLargestPossibleRegion().GetSize();
      Label2DImageType::Pointer l2d = Label2DImageType::New();
      Label2DImageType::SizeType l2dsize; l2dsize[0] = linsize[0]; l2dsize[1] = linsize[1];
      Label2DImageType::IndexType l2dindex; l2dindex.Fill(0);
      Label2DImageType::RegionType l2dregion; l2dregion.SetSize(l2dsize); l2dregion.SetIndex(l2dindex);
      l2d->SetRegions(l2dregion);
      l2d->Allocate();
      l2d->Update();

      LabelImageType::PixelType *iterLabelCopy = l2d->GetBufferPointer();
      LabelImageType::PixelType *iterLabel = im->GetBufferPointer();

#pragma omp parallel for
      for( int ii=0; ii<l2dsize[0]*l2dsize[1]; ++ii )
      {
        iterLabelCopy[ii] = iterLabel[ii];
      }

      //memcpy(im->GetBufferPointer(),l2d->GetBufferPointer(),sizeof(LabelImageType::PixelType)*l2dsize[0]*l2dsize[1]);

      Input2DImageType::Pointer i2d = Input2DImageType::New();
      Input2DImageType::RegionType l2dregionInput; l2dregionInput.SetSize(l2dsize); l2dregionInput.SetIndex(l2dindex);
      i2d->SetRegions(l2dregionInput);
      i2d->Allocate();
      i2d->Update();

      InputImageType::PixelType *iterInputCopy = i2d->GetBufferPointer();
      InputImageType::PixelType *iterInput = in_image->GetBufferPointer();

#pragma omp parallel for
      for( int ii=0; ii<l2dsize[0]*l2dsize[1]; ++ii )
      {
        iterInputCopy[ii] = iterInput[ii];
      }
      //memcpy(in_image->GetBufferPointer(),i2d->GetBufferPointer(),sizeof(Input2DImageType::PixelType)*l2dsize[0]*l2dsize[1]);


      // 		std::stringstream testName;
      // 		testName<<time;
      // 		std::string testImage = "/data/nicolas/test/Image_" + testName.str() + ".tif";
      // 		std::string testLabel = "/data/nicolas/test/Label_" + testName.str() + ".tif";


      // 		typedef itk::ImageFileWriter<Label2DImageType> WriterTypeLabel;
      // 		WriterTypeLabel::Pointer writerLabel = WriterTypeLabel::New();
      // 		writerLabel->SetInput(l2d);
      // 		writerLabel->SetFileName(testLabel.c_str());
      // 		writerLabel->Update();
      // 
      //                 typedef itk::ImageFileWriter<Input2DImageType> WriterTypeImage;
      //                 WriterTypeImage::Pointer writerImage = WriterTypeImage::New();
      //                 writerImage->SetInput(i2d);
      //                 writerImage->SetFileName(testImage.c_str());
      // 
      // 		writerImage->Update();


      typedef ftk::LabelImageToFeatures<Input2DImageType::PixelType, Label2DImageType::PixelType, 2> FeatureCalculator2DType;
      FeatureCalculator2DType::Pointer fc2d = FeatureCalculator2DType::New();
      fc2d->SetImageInputs(i2d,l2d);

      fc2d->SetLevel(1);
      //fc2d->ComputeTexturesOn();
      //	//fc2d->ComputeHistogramOn();		Amin comment
      fc2d->Update();
      std::vector<Label2DImageType::PixelType> labels = fc2d->GetLabels();
      //	printf("label size:%d\n",labels.size());
      for(unsigned int counter=0; counter<labels.size();counter++)
      {
        //		printf("label:%d\n",labels[counter]);
        if(labels[counter]==0)
          continue;
        feature_vector.push_back(*(fc2d->GetFeatures(labels[counter])));
        feature_vector.back().num=labels[counter];
        feature_vector.back().tag = tag;
        feature_vector.back().time = time;
        //		printf("added:%d\n",feature_vector.back().num);
      }

    }
    else
    {
      typedef ftk::LabelImageToFeatures<InputImageType::PixelType,LabelImageType::PixelType,3> FeatureCalculatorType;
      FeatureCalculatorType::Pointer fc = FeatureCalculatorType::New();
      fc->SetImageInputs(in_image,im);
      fc->SetLevel(3);
      //	fc->ComputeTexturesOn();
      //	fc->ComputeHistogramOn();
      fc->Update();
      std::vector<LabelImageType::PixelType> labels = fc->GetLabels();
      for(unsigned int counter=0; counter< labels.size(); counter++)
      {
        if(labels[counter]==0)
          continue;
        //printf("label:%d\n",labels[counter]);
        feature_vector.push_back(*(fc->GetFeatures(labels[counter])));
        feature_vector.back().num = labels[counter];
        feature_vector.back().tag = tag;
        feature_vector.back().time = time;
      }
    }
    //printf("Ended feature calculation\n");
  }

  //---------------------------------------------------------------------------------------------------------------------
  InputImageType::Pointer getDilated(InputImageType::Pointer im, int n)
  {
    printf("Performing morphological dilation ...");
    typedef itk::BinaryBallStructuringElement<InputPixelType,3> StructuringElementType;
    typedef itk::Neighborhood<InputPixelType,3> NeighborhoodElementType;
    typedef itk::BinaryDilateImageFilter<InputImageType,InputImageType,NeighborhoodElementType> DilateFilterType;

    StructuringElementType selement;
    NeighborhoodElementType::SizeType size;
    size[0]=n;
    size[1]=n;
    size[2]=n/3;//FIXME
    selement.SetRadius(size);
    selement.CreateStructuringElement();
    DilateFilterType::Pointer dfilter = DilateFilterType::New();
    dfilter->SetKernel(selement);
    dfilter->SetInput(im);
    dfilter->Update();
    printf(" Done.\n");
    return dfilter->GetOutput();
  }

  //---------------------------------------------------------------------------------------------------------------------
  InputImageType::Pointer getEroded(InputImageType::Pointer im, int n)
  {
    printf("Performing morphological Erosion ...");

    typedef itk::BinaryBallStructuringElement<InputPixelType,3> StructuringElementType;
    typedef itk::Neighborhood<InputPixelType,3> NeighborhoodElementType;
    typedef itk::BinaryErodeImageFilter<InputImageType,InputImageType,NeighborhoodElementType> ErodeFilterType;

    StructuringElementType selement;
    NeighborhoodElementType::SizeType size;
    size[0]=n;
    size[1]=n;
    size[2]=n/5;
    selement.SetRadius(size);
    selement.CreateStructuringElement();
    ErodeFilterType::Pointer efilter = ErodeFilterType::New();
    efilter->SetInput(im);
    efilter->SetKernel(selement);
    efilter->Update();

    printf(" Done.\n");
    return efilter->GetOutput();
  }

  //---------------------------------------------------------------------------------------------------------------------
  DistanceMapFilterType::Pointer getDistanceMap(InputImageType::Pointer im)
  {
    printf("Computing Danielsson distance map ... ");
    DistanceMapFilterType::Pointer distfilter = DistanceMapFilterType::New();
    distfilter->SetInput(im);
    distfilter->InputIsBinaryOn();
    distfilter->Update();
    printf("Done.\n");
    return distfilter;
  }
  /* obselete */
  //InputImageType::Pointer getOldSegmented(InputImageType::Pointer im_input,int threshold, int min_component_size, int morph_opening_depth)
  //{
  //	InputImageType::Pointer threshVessel = getThresholded(im_input,threshold);
  //	threshVessel = getDilated(threshVessel,morph_opening_depth);
  //	threshVessel = getEroded(threshVessel,morph_opening_depth);
  //	threshVessel = getLargeComponents(threshVessel,min_component_size);
  //	return threshVessel;
  //}

  //---------------------------------------------------------------------------------------------------------------------
  LabelImageType::Pointer getLabelled(InputImageType::Pointer im_input,int threshold, int min_component_size, int morph_opening_depth)
  {
    InputImageType::Pointer threshVessel = getThresholded(im_input,threshold);
    threshVessel = getDilated(threshVessel,morph_opening_depth);
    threshVessel = getEroded(threshVessel,morph_opening_depth);
    threshVessel = getLargeComponents(threshVessel,min_component_size);



    ConnectedFilterType::Pointer cfilter = ConnectedFilterType::New();
    cfilter->SetFullyConnected(1);
    cfilter->SetInput(threshVessel);
    cfilter->Update();

    RelabelFilterType::Pointer rfilter = RelabelFilterType::New();
    rfilter->SetInput(cfilter->GetOutput());
    rfilter->InPlaceOn();

    rfilter->Update();
    return rfilter->GetOutput();
  }

  //---------------------------------------------------------------------------------------------------------------------
  InputImageType::Pointer getLabelToBinary(LabelImageType::Pointer l)
  {
    InputImageType::Pointer im = InputImageType::New();
    im->SetRegions(l->GetLargestPossibleRegion());
    im->Allocate();

    IteratorType iter(im,im->GetLargestPossibleRegion());
    ConstLabelIteratorType liter(l,l->GetLargestPossibleRegion());
    for(iter.GoToBegin(),liter.GoToBegin();!liter.IsAtEnd();++liter,++iter)
    {
      if(liter.Get()>0)
        iter.Set(255);
    }
    return im;

  }

  //---------------------------------------------------------------------------------------------------------------------
  ColorImageType::Pointer getColorCompositeImage(InputImageType::Pointer im[4],VectorPixelType colors[4])
  {

    ColorImageType::Pointer output = ColorImageType::New();
    ColorImageType::SizeType size = im[0]->GetLargestPossibleRegion().GetSize();
    ColorImageType::IndexType index;
    index.Fill(0);
    ColorImageType::RegionType region;
    region.SetSize(size);
    region.SetIndex(index);
    output->SetRegions(region);
    output->Allocate();

    ColorIteratorType oiter(output,output->GetLargestPossibleRegion());
    ConstIteratorType iiter[4];
    for(int counter=0; counter <4; counter++)
    {
      iiter[counter] = ConstIteratorType(im[counter],im[counter]->GetLargestPossibleRegion());
      iiter[counter].GoToBegin();
    }

    for(oiter.GoToBegin();!oiter.IsAtEnd(); ++oiter,++iiter[0],++iiter[1],++iiter[2],++iiter[3])
    {
      for(int counter=0; counter<4; counter++)
      {
        if(iiter[counter].Get()!=0)
        {
          oiter.Set(colors[counter]);
          break;
        }
      }
    }
    return output;

  }

  //---------------------------------------------------------------------------------------------------------------------
  ColorImageType::Pointer getColorComposite(InputImageType::Pointer im[],int n, VectorPixelType colors[])
  {

    ColorImageType::Pointer output = ColorImageType::New();
    ColorImageType::SizeType size = im[0]->GetLargestPossibleRegion().GetSize();
    ColorImageType::IndexType index;
    index.Fill(0);
    ColorImageType::RegionType region;
    region.SetSize(size);
    region.SetIndex(index);
    output->SetRegions(region);
    output->Allocate();

    ColorIteratorType oiter(output,output->GetLargestPossibleRegion());
    ConstIteratorType iiter[4];
    for(int counter=0; counter <n; counter++)
    {
      iiter[counter] = ConstIteratorType(im[counter],im[counter]->GetLargestPossibleRegion());
      iiter[counter].GoToBegin();
    }

    printf("Starting Iteration\n");
    for(oiter.GoToBegin();!oiter.IsAtEnd(); ++oiter)
    {
      for(int counter=0; counter<n; counter++)
      {
        if(iiter[counter].Get()!=0)
        {
          oiter.Set(colors[counter]);
          break;
        }
      }
      for(int c =0; c<n; c++)
        ++iiter[c];
    }
    return output;

  }

  //---------------------------------------------------------------------------------------------------------------------
  InputImageType::Pointer getImageFromNPTS(char *filename_npts,int imagesize[])
  {
    InputImageType::Pointer output = InputImageType::New();
    InputImageType::RegionType region;
    InputImageType::SizeType size;
    InputImageType::IndexType index;
    index.Fill(0);
    size[0]=imagesize[0];
    size[1]=imagesize[1];
    size[2]=imagesize[2];
    region.SetSize(size);
    region.SetIndex(index);
    output->SetRegions(region);
    output->Allocate();

    FILE *fp = fopen(filename_npts,"r");
    int x,y,z;
    while(fscanf(fp,"%d %d %d %*f %*f",&z,&y,&x)>0)
    {
      index[0]=x-1;
      index[1]=y-1;
      index[2]=z-1;
      output->SetPixel(index,255);
    }
    fclose(fp);
    return output;


  }

  //---------------------------------------------------------------------------------------------------------------------
  //
  void getClassified(DistanceImageType::Pointer dist, InputImageType::Pointer micro, InputImageType::Pointer &p1, InputImageType::Pointer &p2)
  {

    //	typedef itk::Image<short int,3> LabelImageType;
    typedef itk::Image<LabelPixelType,3> LabelImageType;
    typedef itk::ConnectedComponentImageFilter<InputImageType,LabelImageType> ConnectedFilterType;
    typedef itk::RelabelComponentImageFilter<LabelImageType,LabelImageType> RelabelFilterType;

    ConnectedFilterType::Pointer cfilter = ConnectedFilterType::New();
    cfilter->SetFullyConnected(1);
    cfilter->SetInput(micro);
    cfilter->Update();

    RelabelFilterType::Pointer rfilter = RelabelFilterType::New();
    rfilter->SetInput(cfilter->GetOutput());
    rfilter->InPlaceOn();

    rfilter->Update();

    int num_objects = rfilter->GetNumberOfObjects();
    typedef itk::ImageRegionIterator<LabelImageType> LabelIteratorType;
    LabelIteratorType liter(rfilter->GetOutput(),rfilter->GetOutput()->GetLargestPossibleRegion());

    typedef itk::ImageRegionIterator<DistanceImageType> DistanceIteratorType;

    DistanceIteratorType diter(dist,dist->GetLargestPossibleRegion());


    p1=InputImageType::New();
    p2=InputImageType::New();
    InputImageType::SizeType size=micro->GetLargestPossibleRegion().GetSize();
    InputImageType::IndexType index1;
    index1.Fill(0);
    InputImageType::RegionType region;
    region.SetSize(size);
    region.SetIndex(index1);
    p1->SetRegions(region);
    p2->SetRegions(region);
    p1->Allocate();
    p2->Allocate();

    IteratorType iter1(p1,p1->GetLargestPossibleRegion()),iter2(p2,p2->GetLargestPossibleRegion());

    //	std::vector<short int> distance;
    std::vector<LabelPixelType> distance;
    for(int counter=0; counter<num_objects; counter++)
      distance.push_back(1000);
    printf("About to compute distances\n");
    for(liter.GoToBegin(),diter.GoToBegin(),iter1.GoToBegin(),iter2.GoToBegin();!liter.IsAtEnd(); ++liter,++diter)
    {
      if(liter.Get()>0)
      {
        if(distance[liter.Get()-1]>diter.Get())
        {
          if(diter.Get()<0)
            printf("What the hell!\n");
          distance[liter.Get()-1]=diter.Get();
        }
      }
    }
    for(int counter=0; counter<num_objects;counter++)
    {
      printf("%d \t",(int)distance[counter]);
    }
    for(liter.GoToBegin();!liter.IsAtEnd();++iter1,++iter2,++liter)
    {
      if(liter.Get()!=0)
      {
        if(distance[liter.Get()-1]<4)
        {
          iter1.Set(255);
          iter2.Set(0);
        }
        else
        {
          iter2.Set(255);
          iter1.Set(0);
        }
      }
    }
  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  void writeNPTS(InputImageType::Pointer im, char*filename)
  {
    //	printf("About to write npts file %s\n",filename);
    FILE*fp = fopen(filename,"w");
    InputImageType::SizeType size = im->GetLargestPossibleRegion().GetSize();
    InputImageType::IndexType index;
    for(unsigned int cx=0; cx<size[0];cx++)
    {
      index[0]=cx;
      for(unsigned int cy=0; cy<size[1];cy++)
      {
        index[1]=cy;
        for(unsigned int cz=0; cz<size[2];cz++)
        {
          index[2]=cz;
          if(im->GetPixel(index)>0)
          {
            fprintf(fp,"%d %d %d 1 1\n",cz+1,cy+1,cx+1);
          }
        }
      }
    }
    fclose(fp);
  }

  /*void MSA_classifymicroglia(char*filename_vesselnpts,char*filename_microglianpts)
    {

    int sizes[3]={1024,1024,77};
    InputImageType::Pointer microglia = getImageFromNPTS(filename_microglianpts,sizes);
    InputImageType::Pointer threshVessel = getImageFromNPTS(filename_vesselnpts,sizes);

    DistanceMapFilterType::Pointer distfilt ;//= getDistanceMap(threshVessel);

    DistanceImageType::Pointer distImage=readImage<DistanceImageType>("G:\\Arun\\MSA paper\\distance_map.tif");

    InputImageType::Pointer mic1,mic2;
    getClassified(distImage,microglia,mic1,mic2);
    writeNPTS(mic1,"G:\\Arun\\MSA paper\\mic1.npts");
    writeNPTS(mic2,"G:\\Arun\\MSA paper\\mic2.npts");

  //	writeImage<DistanceImageType>(distfilt->GetOutput(),"G:\\Arun\\MSA paper\\distance_map.tif");
  }
  */

  //---------------------------------------------------------------------------------------------------------------------
  // 
  InputImageType::Pointer getPreprocessed(InputImageType::Pointer im)
  {
    MedianFilterType::Pointer filter = MedianFilterType::New();
    filter->SetInput(im);
    InputImageType::SizeType radius;
    radius[0]=1;
    radius[1]=1;
    radius[2]=1;
    filter->SetRadius(radius);
    filter->Update();
    return filter->GetOutput();
  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  ColorImageType::Pointer getColorCompositeImageFromLabelled(LabelImageType::Pointer im,VectorPixelType color)
  {
    ColorImageType::Pointer output = ColorImageType::New();
    ColorImageType::SizeType size = im->GetLargestPossibleRegion().GetSize();
    ColorImageType::IndexType index;
    index.Fill(0);
    ColorImageType::RegionType region;
    region.SetSize(size);
    region.SetIndex(index);
    output->SetRegions(region);
    output->Allocate();
    VectorPixelType black;
    black[0]=0;black[1]=0;black[2]=0;

    ColorIteratorType oiter(output,output->GetLargestPossibleRegion());
    ConstLabelIteratorType iiter;
    iiter = ConstLabelIteratorType(im,im->GetLargestPossibleRegion());
    iiter.GoToBegin();

    for(oiter.GoToBegin();!oiter.IsAtEnd(); ++oiter,++iiter)
    {
      if(iiter.Get())
        oiter.Set(color);
      else
        oiter.Set(black);
    }
    return output;
  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  Color2DImageType::Pointer getColor2DImage(LabelImageType::Pointer labelled, int channel)
  {
    unsigned char colorarray[][3]={
      {255,0,0},
      {0,154,25},
      {207,141,0},
      {255,0,0} };
    VectorPixelType colorcodes[4];

    for(int counter=0; counter<4; counter++)
      colorcodes[counter]=colorarray[counter];
    ColorImageType::Pointer output_segmentation = getColorCompositeImageFromLabelled(labelled,colorcodes[channel-1]);
    return getColorProjection(output_segmentation);	
  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  InputImageType::Pointer getEmpty(int s1,int s2, int s3)
  {
    InputImageType::Pointer p = InputImageType::New();
    InputImageType::SizeType size;
    InputImageType::IndexType index;
    InputImageType::RegionType region;
    size[0] = s1; size[1] = s2; size[2] = s3;
    index.Fill(0);
    region.SetSize(size);
    region.SetIndex(index);
    p->SetRegions(region);
    p->Allocate();
    return p;
  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  Input2DImageType::Pointer get2DEmpty(int s1, int s2)
  {
    Input2DImageType::Pointer p = Input2DImageType::New();
    Input2DImageType::SizeType size;
    Input2DImageType::IndexType index;
    Input2DImageType::RegionType region;
    size[0] = s1;size[1] = s2;
    index.Fill(0);
    region.SetSize(size);
    region.SetIndex(index);
    p->SetRegions(region);
    p->Allocate();
    return p;
  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  Input2DImageType::Pointer get2DBoundary(LabelImageType::Pointer label)
  {
    LabelIteratorType liter = LabelIteratorType(label,label->GetLargestPossibleRegion());
    liter.GoToBegin();

    //find the maximum number of cells
    unsigned short max1 = 0;
    for(liter.GoToBegin();!liter.IsAtEnd();++liter)
      max1 = MAX(max1,liter.Get());

    //find all the cubes in which cells lie
    std::vector<cubecoord> carray(max1+1);
    for(int counter=0; counter<=max1; counter++)
    {
      carray[counter].sx=60000;carray[counter].sy=60000;carray[counter].sz=60000;
      carray[counter].ex=0;carray[counter].ey=0;carray[counter].ez=0;
    }

    typedef itk::ImageRegionConstIteratorWithIndex<LabelImageType> ConstLabelIteratorWithIndex;
    ConstLabelIteratorWithIndex cliter = ConstLabelIteratorWithIndex(label,label->GetLargestPossibleRegion());
    InputImageType::IndexType index;
    for(cliter.GoToBegin();!cliter.IsAtEnd();++cliter)
    {
      int cur = cliter.Get();
      if(cur!=0)
      {
        index = cliter.GetIndex();
        carray[cur].sx= MIN(index[0],carray[cur].sx);
        carray[cur].sy= MIN(index[1],carray[cur].sy);
        carray[cur].sz= MIN(index[2],carray[cur].sz);
        carray[cur].ex= MAX(index[0],carray[cur].ex);
        carray[cur].ey= MAX(index[1],carray[cur].ey);
        carray[cur].ez= MAX(index[2],carray[cur].ez);
      }
    }

    //find the largest image size we need
    unsigned short wx=0,wy=0,wz=0;
    for(int counter=1; counter<=max1; counter++)
    {
      wx = MAX(carray[counter].ex-carray[counter].sx+1,wx);
      wy = MAX(carray[counter].ey-carray[counter].sy+1,wy);
      wz = MAX(carray[counter].ez-carray[counter].sz+1,wz);
    }
    // accommodate padding
    wx = wx+2;wy = wy +2; wz = wz+2;
    // create a tiny image of maximum size

    LabelImageType::SizeType globalsize = label->GetLargestPossibleRegion().GetSize();
    Input2DImageType::Pointer t2d = get2DEmpty(wx,wy);
    Input2DImageType::Pointer bound_im = get2DEmpty(globalsize[0],globalsize[1]);

    bound_im->FillBuffer(0);
    typedef itk::ImageRegionIteratorWithIndex<Input2DImageType> Iterator2DWithIndexType;

    for(int counter=1; counter<=max1; counter++)
    {

      Input2DImageType::SizeType size;
      Input2DImageType::IndexType index2d;
      Input2DImageType::RegionType region;
      index2d.Fill(1);

      region.SetIndex(index2d);

      LabelImageType::SizeType lsize;
      LabelImageType::IndexType lindex;
      LabelImageType::RegionType lregion;

      t2d->FillBuffer(0);
      lsize[0] = carray[counter].ex-carray[counter].sx+1;
      lsize[1] = carray[counter].ey-carray[counter].sy+1;
      lsize[2] = carray[counter].ez-carray[counter].sz+1;

      lindex[0] = carray[counter].sx;
      lindex[1] = carray[counter].sy;
      lindex[2] = carray[counter].sz;

      lregion.SetIndex(lindex);
      lregion.SetSize(lsize);

      ConstLabelIteratorWithIndex localiter = ConstLabelIteratorWithIndex(label,lregion);

      size[0] = lsize[0];size[1] = lsize[1];
      region.SetSize(size);
      twoDIteratorType iter = twoDIteratorType(t2d,region);


      t2d->FillBuffer(0);
      for(localiter.GoToBegin(),iter.GoToBegin();!localiter.IsAtEnd();++localiter,++iter)
      {
        if(iter.IsAtEnd())
          iter.GoToBegin();
        if(iter.Get()==0)
        {
          iter.Set(255*(localiter.Get()!=0));
        }
      }

      Iterator2DWithIndexType iter2d = Iterator2DWithIndexType(t2d,region);
      Input2DImageType::IndexType idx;
      idx[0] = lindex[0];idx[1] = lindex[1];
      region.SetIndex(idx);
      twoDIteratorType global2diter = twoDIteratorType(bound_im,region);
      global2diter.GoToBegin();
      for(iter2d.GoToBegin(); !iter2d.IsAtEnd();++iter2d,++global2diter)
      {
        unsigned char cur = iter2d.Get();
        if(cur==0)
          continue;
        index2d = iter2d.GetIndex();
        //check the four neighbors to see if this pixel is in boundary
        index2d[0]--;

        if(t2d->GetPixel(index2d)==0)
        {
          global2diter.Set(255);
          continue;
        }
        index2d[0]++;
        index2d[1]--;
        if(t2d->GetPixel(index2d)==0)
        {
          global2diter.Set(255);
          continue;
        }
        index2d[1]++;
        index2d[1]++;
        if(t2d->GetPixel(index2d)==0)
        {
          global2diter.Set(255);
          continue;
        }
        index2d[1]--;
        index2d[0]++;
        if(t2d->GetPixel(index2d)==0)
        {
          global2diter.Set(255);
          continue;
        }
        index2d[0]--;
      }

    }
    return bound_im;

  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  Color2DImageType::Pointer getColorBoundaryImage(LabelImageType::Pointer labelled, InputImageType::Pointer im, int channel)
  {
    double multiplier;
    if(channel==1)
      multiplier = 6;
    else
      multiplier = 1.5;
    unsigned char colorarray[][3]={
      {255,0,0},
      {0,154,25},
      {207,141,0},
      {255,0,0} };
    VectorPixelType colorcodes[4];
    VectorPixelType black;black[0]=0;black[1]=0;black[2]=0;

    for(int counter=0; counter<4; counter++)
      colorcodes[counter]=colorarray[counter];

    Input2DImageType::Pointer boundary = get2DBoundary(labelled);
    Color2DImageType::Pointer cimage = Color2DImageType::New();
    cimage->SetRegions(boundary->GetLargestPossibleRegion());
    cimage->Allocate();

    Input2DImageType::Pointer im2d = getProjection(im);
    VectorPixelType pix;
    typedef itk::ImageRegionIterator<Color2DImageType> Color2DIteratorType;
    Color2DIteratorType iter = Color2DIteratorType(cimage,cimage->GetLargestPossibleRegion());
    twoDIteratorType biter = twoDIteratorType(boundary,boundary->GetLargestPossibleRegion());
    twoDIteratorType inputiter = twoDIteratorType(im2d,im2d->GetLargestPossibleRegion());
    //generate a mask image 
    InputImageType::Pointer mask=InputImageType::New();
    mask->SetRegions(labelled->GetLargestPossibleRegion());
    mask->Allocate();
    LabelIteratorType liter = LabelIteratorType(labelled,labelled->GetLargestPossibleRegion());
    IteratorType maskiter = IteratorType(mask,mask->GetLargestPossibleRegion());
    for(liter.GoToBegin(),maskiter.GoToBegin();!maskiter.IsAtEnd();++maskiter,++liter)
    {
      maskiter.Set(255*(liter.Get()!=0));
    }
    Input2DImageType::Pointer mask2d = getProjection(mask);
    twoDIteratorType mask2diter = twoDIteratorType(mask2d,mask2d->GetLargestPossibleRegion());
    //initialize all iterators
    iter.GoToBegin();biter.GoToBegin();inputiter.GoToBegin();mask2diter.GoToBegin();
    for(;!iter.IsAtEnd();++iter,++biter,++inputiter,++liter,++mask2diter)
    {
      if(biter.Get()!=0)
      {
        pix = colorcodes[channel-1];
      }
      else if(mask2diter.Get()!=0)
      {
        pix[0] = (unsigned char)MIN(inputiter.Get()*multiplier,255);
        pix[1] = pix[0];
        pix[2] = pix[0];
      }
      else
      {
        pix = black;
      }
      iter.Set(pix);
    }

    return cimage;
  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  LabelImageType::Pointer getLabelsMapped(LabelImageType::Pointer label, std::vector<FeatureType> &fvec, unsigned int * indices)
  {
    //ideally .num should have the actual cell id which we want to map to track id, but its not so :-/
    printf("Entering getLabelsMapped\n");

    LabelIteratorType liter = LabelIteratorType(label,label->GetLargestPossibleRegion());
    LabelImageType::Pointer lout = LabelImageType::New();
    lout->SetRegions(label->GetLargestPossibleRegion());
    lout->Allocate();
    LabelIteratorType oiter = LabelIteratorType(lout,lout->GetLargestPossibleRegion());
    liter.GoToBegin();
    oiter.GoToBegin();
    for(;!liter.IsAtEnd(); ++liter,++oiter)
    {
      if(liter.Get()>0)
        oiter.Set(indices[liter.Get()-1]+1);
      else
        oiter.Set(0);
    }
    for(unsigned int counter=0; counter<fvec.size(); counter++)
    {
      fvec[counter].num = indices[fvec[counter].num-1]+1;
    }
    printf("Exiting getLabelsMapped\n");
    return lout;
  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  void CLAMP(InputImageType::IndexType &i,InputImageType::SizeType size)
  {
    if((unsigned int)i[0]>size[0]-1)
    {
      i[0]=size[0]-1;
    }
    if((unsigned int)i[1]>size[1]-1)
    {
      i[1]=size[1]-1;
    }
    if((unsigned int)i[2]>size[2]-1)
    {
      i[2]=size[2]-1;
    }
  }

  std::vector<float> traverseCenterline(itk::ImageRegionIteratorWithIndex<InputImageType> iter,InputImageType::Pointer im,char neighbors[26][3],int n)
  {
    iter.Set(2);
    InputImageType::IndexType tindex;
    std::queue<PQdata> q;
    PQdata temp;
    temp.index = iter.GetIndex();
    temp.depth = n;
    q.push(temp);
    std::vector<float> output;
    std::vector<InputImageType::IndexType> coordinates;
    PQdata newelement;
    InputImageType::SizeType size = im->GetLargestPossibleRegion().GetSize();
    while(!q.empty())
    {
      PQdata top = q.front();
      q.pop();
      coordinates.push_back(top.index);
      if(top.depth==0)
      {
        output.push_back(top.index[0]);
        output.push_back(top.index[1]);
        output.push_back(top.index[2]);
        continue;
      }
      bool done_once = false;
      int pushed_count=0;
      for(int counter=0; counter<26;counter++)
      {
        tindex[0]=top.index[0]+neighbors[counter][0];
        tindex[1]=top.index[1]+neighbors[counter][1];
        tindex[2]=top.index[2]+neighbors[counter][2];
        if(tindex[0]>=0 && tindex[1] >= 0 && tindex[2] >=0 &&
            (unsigned int)tindex[0] <size[0] && (unsigned int)tindex[1] < size[1]
            && (unsigned int)tindex[2] < size[2])
        {
          if(im->GetPixel(tindex)==1)
          {
            done_once = true;
            pushed_count++;
            newelement.index = tindex;
            newelement.depth = top.depth-1;
            im->SetPixel(tindex,2);
            q.push(newelement);
          }
        }
      }
      //printf("I pushed %d values top.depth = %d qsize = %d\n",pushed_count,top.depth,q.size());
      if(!done_once)
      {
        output.push_back(top.index[0]);
        output.push_back(top.index[1]);
        output.push_back(top.index[2]);
      }
    }
    if(output.size()==3)
    {
      output.push_back(temp.index[0]);
      output.push_back(temp.index[1]);
      output.push_back(temp.index[2]);
    }
    //printf("coordinates.size() = %d\n",coordinates.size());
    for(unsigned int counter=0; counter<coordinates.size(); counter++)
    {
      im->SetPixel(coordinates[counter],1);
      //	printf("Coordinates: %d %d %d\n",coordinates[counter][0],coordinates[counter][1],coordinates[counter][2]);
    }
    /*
       for(int counter=0; counter<output.size()/3; counter++)
       {
       printf("Output: %0.3f %0.3f %0.3f\n",output[3*counter],output[3*counter+1],output[3*counter+2]);
       }*/
    //scanf("%*d");
    assert(output.size()>=6);
    return output;
  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  FloatImageType::IndexType searchNearestVesselDirection(FloatImageType::Pointer dir_image[3],FloatImageType::IndexType index,InputImageType::Pointer vesselim)
  {
    const int sdepth = 5;
    //int cur_depth = 0;
    InputImageType::IndexType tindex;
    tindex[0] = -1;
    tindex[1] = -1;
    tindex[2] = -1;
    int xind, yind,zind;
    for(int cur_depth = 0; cur_depth<sdepth ; cur_depth++)
    {
      zind = index[2]+cur_depth;
      tindex[2] = zind;
      for(xind = index[0]-cur_depth; xind <= index[0]+cur_depth; xind++)
      {
        tindex[0] = xind;
        for(yind = index[1]-cur_depth; yind <= index[1]+cur_depth; yind++)
        {
          tindex[1] = yind;
          if(vesselim->GetPixel(tindex)!=0)
          {
            return tindex;
          }
        }
      }
      zind = index[2]-cur_depth;
      tindex[2] = zind;
      for(xind = index[0]-cur_depth; xind <= index[0]+cur_depth; xind++)
      {
        tindex[0] = xind;
        for(yind = index[1]-cur_depth; yind <= index[1]+cur_depth; yind++)
        {
          tindex[1] = yind;
          if(vesselim->GetPixel(tindex)!=0)
          {
            return tindex;
          }
        }
      }
      xind = index[0]+cur_depth;
      tindex[0] = xind;
      for(zind = index[2]-cur_depth; zind <= index[2]+cur_depth; zind++)
      {
        tindex[2] = zind;
        for(yind = index[1]-cur_depth; yind <= index[1]+cur_depth; yind++)
        {
          tindex[1] = yind;
          if(vesselim->GetPixel(tindex)!=0)
          {
            return tindex;
          }
        }
      }
      xind = index[0]-cur_depth;
      tindex[0] = xind;
      for(zind = index[2]-cur_depth; zind <= index[2]+cur_depth; zind++)
      {
        tindex[2] = zind;
        for(yind = index[1]-cur_depth; yind <= index[1]+cur_depth; yind++)
        {
          tindex[1] = yind;
          if(vesselim->GetPixel(tindex)!=0)
          {
            return tindex;
          }
        }
      }
      yind = index[1]+cur_depth;
      tindex[1] = yind;
      for(xind = index[0]-cur_depth; xind <= index[0]+cur_depth; xind++)
      {
        tindex[0] = xind;
        for(zind = index[2]-cur_depth; yind <= index[2]+cur_depth; zind++)
        {
          tindex[2] = zind;
          if(vesselim->GetPixel(tindex)!=0)
          {
            return tindex;
          }
        }
      }
      yind = index[1]-cur_depth;
      tindex[1] = yind;
      for(xind = index[0]-cur_depth; xind <= index[0]+cur_depth; xind++)
      {
        tindex[0] = xind;
        for(zind = index[2]-cur_depth; yind <= index[2]+cur_depth; zind++)
        {
          tindex[2] = zind;
          if(vesselim->GetPixel(tindex)!=0)
          {
            return tindex;
          }
        }
      }
    }
    return tindex;
  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  void AnalyzeTimeFeatures(std::vector<ftk::TrackFeatures> &tfs, float spacing[3])
  {
    for(unsigned int tcounter=0; tcounter < tfs.size(); tcounter++) // looping over labels (tracks)
    {
      printf("Beginning to read track no: %d/%d\n",tcounter+1, (int)tfs.size());
      //std::vector<TrackPoint> track = total_tracks[tcounter];
      if(tfs[tcounter].intrinsic_features.size()<2)
      {
        printf("Ignored a tiny track of size %d track.size()\n",(int)tfs[tcounter].intrinsic_features.size());
        continue;
      }


      printf("I'm working on a track of size %d\n",(int)tfs[tcounter].intrinsic_features.size());
      ftk::TrackFeatures t = tfs[tcounter];
      typedef ftk::TrackPointFeatures TPF;
      typedef ftk::TrackFeatures TF;
      typedef ftk::IntrinsicFeatures FeatureType;

      for(unsigned int counter=0; counter< t.intrinsic_features.size(); counter++) // looping over time for each labeled cell
      {
        ftk::TrackPointFeatures tpf;
        Vec3f dir;
        Vec3f dirnext;
        if(counter>0)
        {

          dir.x=t.intrinsic_features[counter].Centroid[0]*spacing[0]-t.intrinsic_features[counter-1].Centroid[0]*spacing[0];
          dir.y=t.intrinsic_features[counter].Centroid[1]*spacing[1]-t.intrinsic_features[counter-1].Centroid[1]*spacing[1];
          dir.z=t.intrinsic_features[counter].Centroid[2]*spacing[2]-t.intrinsic_features[counter-1].Centroid[2]*spacing[2];
          double dirdist = sqrt(dir.x*dir.x+dir.y*dir.y+dir.z*dir.z);
          tpf.scalars[TPF::DISPLACEMENT_VEC_X]=dir.x;
          tpf.scalars[TPF::DISPLACEMENT_VEC_Y]=dir.y;
          tpf.scalars[TPF::DISPLACEMENT_VEC_Z]=dir.z;
          dir.Normalize();
          tpf.scalars[TPF::INST_SPEED] = dirdist/(t.intrinsic_features[counter].time-t.intrinsic_features[counter-1].time);
          tpf.scalars[TPF::DISTANCE] = dirdist;
        }
        else
        {
          for (int i = 0; i<tpf.M; ++i)
            tpf.scalars[i]=0;
        }
        if(counter +1 > t.tfeatures.size())
          t.tfeatures.push_back(tpf);
        else
          t.tfeatures[counter] = tpf;

      }

      float avg_speed = 0;
      float pathlength = 0;
      float max_speed = 0;
      float min_speed = std::numeric_limits<float>::max();
      for(unsigned int counter =0; counter < t.tfeatures.size(); counter++)
      {
        pathlength += t.tfeatures[counter].scalars[TPF::DISTANCE];
        // Compute stats of track point features:
        max_speed = MAX(t.tfeatures[counter].scalars[TPF::INST_SPEED],max_speed);
        avg_speed += t.tfeatures[counter].scalars[TPF::INST_SPEED];
        if (counter != 0)
          min_speed = MIN(t.tfeatures[counter].scalars[TPF::INST_SPEED],min_speed);
      }
      if (min_speed == std::numeric_limits<float>::max())
        min_speed =0;
      int tnum = t.tfeatures.size();
      t.scalars[TF::AVG_SPEED] = avg_speed/(float)(t.tfeatures.size()-1);
      t.scalars[TF::MAX_SPEED] = max_speed;
      t.scalars[TF::MIN_SPEED] = min_speed;

      t.scalars[TF::PATHLENGTH] = pathlength;
      t.scalars[TF::DISPLACEMENT_VEC_X] = t.intrinsic_features[tnum-1].Centroid[0]*spacing[0]-t.intrinsic_features[0].Centroid[0]*spacing[0];
      t.scalars[TF::DISPLACEMENT_VEC_Y] = t.intrinsic_features[tnum-1].Centroid[1]*spacing[1]-t.intrinsic_features[0].Centroid[1]*spacing[1];
      t.scalars[TF::DISPLACEMENT_VEC_Z] = t.intrinsic_features[tnum-1].Centroid[2]*spacing[2]-t.intrinsic_features[0].Centroid[2]*spacing[2];
      float total_distance = sqrt(t.scalars[TF::DISPLACEMENT_VEC_X]*t.scalars[TF::DISPLACEMENT_VEC_X]+t.scalars[TF::DISPLACEMENT_VEC_Y]*t.scalars[TF::DISPLACEMENT_VEC_Y]+t.scalars[TF::DISPLACEMENT_VEC_Z]*t.scalars[TF::DISPLACEMENT_VEC_Z]);
      t.scalars[TF::CONFINEMENT_RATIO] = pathlength/total_distance;
      //t.scalars[TF::PATHLENGTH] /= (float)(t.tfeatures.size()-1);
      t.scalars[TF::TOTAL_DISTANCE] = total_distance;
      //t.scalars[TF::TOTAL_DISTANCE] = total_distance/(float)(t.tfeatures.size()-1);
      tfs[tcounter] = t;
    }
  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  void PrintTrackFeatures(std::vector<ftk::TrackFeatures> &tfs,std::string path)
  {

    std::string filename = path+"\\TrackFeatures.txt";
    FILE *fp = fopen(filename.c_str(),"w");
    fprintf(fp,"ID\t avg_speed\t max_speed\t min_speed\t displacement_vec_x\t displacement_vec_y\t displacement_vec_z\t pathlength\t total_distance\t confinement_ratio\n");

    for(unsigned int tcounter=0; tcounter < tfs.size(); tcounter++) // looping over labels (tracks)
    {
      printf("Beginning to read track no: %d/%d\n",tcounter+1, (int)tfs.size());
      //std::vector<TrackPoint> track = total_tracks[tcounter];
      if(tfs[tcounter].intrinsic_features.size()<2)
      {
        printf("Ignored a tiny track of size %d track.size()\n",(int)tfs[tcounter].intrinsic_features.size());
        continue;
      }
      printf("I'm working on a track of size %d\n",(int)tfs[tcounter].intrinsic_features.size());
      ftk::TrackFeatures t = tfs[tcounter];
      typedef ftk::TrackPointFeatures TPF;
      typedef ftk::TrackFeatures TF;
      typedef ftk::IntrinsicFeatures FeatureType;

      fprintf(fp,"%d\t %0.3f\t %0.3f\t %0.3f\t %0.3f\t %0.3f\t %0.3f\t %0.3f\t %0.3f\t %0.3f\n",\
          t.intrinsic_features[0].num,\
          t.scalars[TF::AVG_SPEED],t.scalars[TF::MAX_SPEED],t.scalars[TF::MIN_SPEED], t.scalars[TF::DISPLACEMENT_VEC_X],t.scalars[TF::DISPLACEMENT_VEC_Y],t.scalars[TF::DISPLACEMENT_VEC_Z],t.scalars[TF::PATHLENGTH],\
          t.scalars[TF::TOTAL_DISTANCE],t.scalars[TF::CONFINEMENT_RATIO]);

    }
    fclose(fp);

  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  void AnalyzeVesselCenterlines(InputImageType::Pointer cline, std::vector<ftk::TrackFeatures> &tfs,float spacing[3])
  {
    InputImageType::SizeType clinesize = cline->GetLargestPossibleRegion().GetSize();
    DistanceMapFilterType::Pointer distfilt = getDistanceMap(cline);
    distfilt->Update();
    DistanceImageType::Pointer dmap = distfilt->GetOutput();
    OffsetImageType::Pointer vectmap = distfilt->GetVectorDistanceMap();

    FloatImageType::Pointer dir_image[3];

    dir_image[0]=FloatImageType::New();
    dir_image[0]->SetRegions(cline->GetLargestPossibleRegion());
    dir_image[0]->Allocate();
    dir_image[1]=FloatImageType::New();
    dir_image[1]->SetRegions(cline->GetLargestPossibleRegion());
    dir_image[1]->Allocate();
    dir_image[2]=FloatImageType::New();
    dir_image[2]->SetRegions(cline->GetLargestPossibleRegion());
    dir_image[2]->Allocate();

    typedef itk::ImageRegionIterator<FloatImageType> FloatIter;
    FloatIter fiter1 = FloatIter(dir_image[0],dir_image[0]->GetLargestPossibleRegion());
    FloatIter fiter2 = FloatIter(dir_image[1],dir_image[1]->GetLargestPossibleRegion());
    FloatIter fiter3 = FloatIter(dir_image[2],dir_image[2]->GetLargestPossibleRegion());

    typedef itk::ImageRegionIteratorWithIndex<InputImageType> IterWithIndexType;
    IterWithIndexType iter = IterWithIndexType(cline,cline->GetLargestPossibleRegion());
    //float points[6];

    char neighbors[26][3];
    int pc = 0;
    for(int x=-1;x<=1;x++)
    {
      for(int y=-1;y<=1;y++)
      {
        for(int z=-1;z<=1;z++)
        {
          if((x!=0)||(y!=0)||(z!=0))
          {
            neighbors[pc][0]=x;
            neighbors[pc][1]=y;
            neighbors[pc][2]=z;
            pc++;
          }
        }
      }
    }
    printf("Done loading neighbors\n");

    InputImageType::Pointer debug_image = InputImageType::New();
    debug_image->SetRegions(cline->GetLargestPossibleRegion());
    debug_image->Allocate();

    fiter1.GoToBegin();fiter2.GoToBegin();fiter3.GoToBegin();
    float dx,dy,dz,dsum;
    InputImageType::IndexType idx;//,idx2;
    for(iter.GoToBegin();!iter.IsAtEnd();++iter,++fiter1,++fiter2,++fiter3)
    {

      if(iter.Get()!=0)
      {
        idx = iter.GetIndex();
        std::vector<float> fpoints = traverseCenterline(iter,cline,neighbors,5);
        //printf("I finished once fpoints.size() = %d\n",fpoints.size());
        if(fpoints.size()>=6)
        {

          dx = fpoints[3]-fpoints[0];
          dy = fpoints[4]-fpoints[1];
          dz = fpoints[5]-fpoints[2];
          dsum = sqrt(dx*dx+dy*dy+dz*dz);
          if(dsum>0.0001)
          {
            dx/=dsum;dy/=dsum;dz/=dsum;
          }
          else
          {
            printf("Error in dx,dy,dz\n");
          }
          fiter1.Set(dx);fiter2.Set(dy);fiter3.Set(dz);
        }
        else
        {

          dx = idx[0]-fpoints[0];
          dy = idx[1]-fpoints[1];
          dz = idx[2]-fpoints[2];
          dsum = sqrt(dx*dx+dy*dy+dz*dz);
          if(dsum>0.0001)
          {
            dx/=dsum;dy/=dsum;dz/=dsum;
          }
          fiter1.Set(dx);fiter2.Set(dy);fiter3.Set(dz);
          printf("I got fpoints.size() = %d\n",(int)fpoints.size());
          //debug_image->SetPixel(iter.GetIndex(),255);;
        }
        //	idx2[0] = idx[0]+dx*5;idx2[1] = idx[1]+dy*5; idx2[2]= idx[2]+dz*5;

        //drawLine(debug_image,idx,idx2);
        //printf("finished setting fiters %f\n", dsum);
      }
      else
      {

        fiter1.Set(1);fiter2.Set(0);fiter3.Set(0);
      }
    }

    printf("No segfault\n");



    //OffsetImageType::IndexType offset;
    InputImageType::IndexType cellindex;
    InputImageType::IndexType vesselindex;
    //InputImageType::IndexType celldirindex;
    for(unsigned int tcounter=0; tcounter < tfs.size(); tcounter++)
    {
      printf("Beginning to read track no: %d/%d\n",tcounter+1,(int)tfs.size());
      //std::vector<TrackPoint> track = total_tracks[tcounter];
      if(tfs[tcounter].intrinsic_features.size()<3)
      {
        printf("Ignored a tiny track of size %d track.size()\n",(int)tfs[tcounter].intrinsic_features.size());
        continue;
      }
      printf("I'm working on a track of size %d\n",(int)tfs[tcounter].intrinsic_features.size());
      ftk::TrackFeatures t = tfs[tcounter];
      //std::vector<float> spacing(3);
      //spacing[0] = spacing[1] = 0.357;
      //spacing[2] = 2.0;
      typedef ftk::TrackPointFeatures TPF;
      typedef ftk::TrackFeatures TF;
      for(unsigned int counter=0; counter< t.intrinsic_features.size(); counter++)
      {
        ftk::TrackPointFeatures tpf = t.tfeatures[counter];
        Vec3f dir;

        if(counter>0)
        {

          dir.x = tpf.scalars[TPF::DISPLACEMENT_VEC_X];
          dir.y = tpf.scalars[TPF::DISPLACEMENT_VEC_Y];
          dir.z = tpf.scalars[TPF::DISPLACEMENT_VEC_Z];
          dir.Normalize();
        }

        cellindex[0] = static_cast<int>(t.intrinsic_features[counter].Centroid[0]+0.5);
        cellindex[1] = static_cast<int>(t.intrinsic_features[counter].Centroid[1]+0.5);
        cellindex[2] = static_cast<int>(t.intrinsic_features[counter].Centroid[2]+0.5);

        CLAMP(cellindex,cline->GetLargestPossibleRegion().GetSize());
        //			printf("Next Cell Index = %d %d %d\n",cellindex[0],cellindex[1],cellindex[2]);
        OffsetImageType::PixelType offpix = vectmap->GetPixel(cellindex);

        tpf.scalars[TPF::DISTANCE_TO_1] = sqrt(offpix[0]*offpix[0]*spacing[0]*spacing[0]+offpix[1]*offpix[1]*spacing[0]*spacing[0]+offpix[2]*offpix[2]*spacing[1]*spacing[1]);

        if(counter>0)
        {
          vesselindex[0] = offpix[0]+cellindex[0];
          vesselindex[1] = offpix[1]+cellindex[1];
          vesselindex[2] = offpix[2]+cellindex[2];


          Vec3f dir2;
          if(cline->GetPixel(vesselindex)==0)
          {
            printf("I'm calling searchNearestVesselDirection\n");
            vesselindex = searchNearestVesselDirection(dir_image,vesselindex,cline);
          }
          dir2.x = spacing[0]*dir_image[0]->GetPixel(vesselindex);
          dir2.y = spacing[1]*dir_image[1]->GetPixel(vesselindex);
          dir2.z = spacing[2]*dir_image[2]->GetPixel(vesselindex);
          if((dir2.x*dir2.x+dir2.y*dir2.y+dir2.z*dir2.z)<1e-6)
          {
            printf("Error! vessel direction is zero\n");
            tpf.scalars[TPF::ANGLE_REL_TO_1] = 0.0;
            continue;
          }
          dir2.Normalize();

          tpf.scalars[TPF::ANGLE_REL_TO_1] = 90.0/acosf(0)*acosf(fabs(dir.x*dir2.x+dir.y*dir2.y+dir.z*dir2.z));
        }
        else
        {
          tpf.scalars[TPF::ANGLE_REL_TO_1] = 0.0;
        }
        printf("Just before\n");
        if(counter > t.tfeatures.size()-1)
        {
          t.tfeatures.push_back(tpf);
        }
        else
        {
          t.tfeatures[counter] = tpf;
        }
        printf("Just After\n");
      }

      float avg_distance_to_1 = 0;
      float avg_angle_rel_to_1 = 0;
      float contact_to_2 = 0;
      //float change_distance_to_1=0.0;
      for(unsigned int counter =0; counter < t.tfeatures.size(); counter++)
      {
        if(counter>0)
        {
          t.tfeatures[counter].scalars[TPF::CHANGE_DISTANCE_TO_1] = t.tfeatures[counter].scalars[TPF::DISTANCE_TO_1] - t.tfeatures[counter-1].scalars[TPF::DISTANCE_TO_1];
        }
        else
        {
          t.tfeatures[counter].scalars[TPF::CHANGE_DISTANCE_TO_1] = t.tfeatures[counter].scalars[TPF::DISTANCE_TO_1];
        }
        avg_distance_to_1 += t.tfeatures[counter].scalars[TPF::DISTANCE_TO_1];
        avg_angle_rel_to_1 += t.tfeatures[counter].scalars[TPF::ANGLE_REL_TO_1];
        contact_to_2 += t.tfeatures[counter].scalars[TPF::HAS_CONTACT_TO_2];
      }
      int tnum = t.tfeatures.size();

      t.scalars[TF::AVG_DIST_TO_1] = avg_distance_to_1/tnum;
      t.scalars[TF::AVG_ANGLE_REL_TO_1] = avg_angle_rel_to_1/(tnum-1);
      t.scalars[TF::CONTACT_TO_2] = contact_to_2/tnum;
      t.scalars[TF::CHANGE_DISTANCE_TO_1] = t.tfeatures[tnum-1].scalars[TPF::DISTANCE_TO_1] - t.tfeatures[0].scalars[TPF::DISTANCE_TO_1];
      tfs[tcounter] = t;
    }

  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  InputImageType::Pointer getMaxImage(InputImageType::Pointer im1, InputImageType::Pointer im2)
  {
    InputImageType::Pointer im3 = InputImageType::New();
    im3->SetRegions(im1->GetLargestPossibleRegion());
    im3->Allocate();

    IteratorType iter3(im3,im3->GetLargestPossibleRegion());
    ConstIteratorType iter1(im1,im1->GetLargestPossibleRegion());
    ConstIteratorType iter2(im2,im2->GetLargestPossibleRegion());

    for(iter3.GoToBegin(),iter1.GoToBegin(),iter2.GoToBegin(); !iter3.IsAtEnd(); ++iter1,++iter2,++iter3)
    {
      iter3.Set(MAX(iter1.Get(),iter2.Get()));
    }
    return im3;
  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  void AnalyzeDCContact(LabelImageType::Pointer segmented[][4], std::vector<ftk::TrackFeatures> &tfs, int c, int num_t, float spacing[3])
  {
    //int min_voxels = 1000;
    //float perc = 5.0f;
    float inc_factor = 1.25;
    int wsize = 20;
    float radius = -1;
    float factor = 1;
    float ratio_threshold = 1.7;
    LabelImageType::SizeType bound = segmented[0][0]->GetLargestPossibleRegion().GetSize();

    typedef ftk::TrackPointFeatures TPF;
    typedef ftk::TrackFeatures TF;
    std::vector<int> distances(10000);
    //std::vector<float> spacing(3);
    //spacing[0] = spacing[1] = 0.357;
    //spacing[2] = 2.0;
    for(unsigned int tc=0; tc < tfs.size(); tc++)
    {
      //	printf("running for track :%d\n",tc);
      // for this track = counter find dc contact feature
      while(tfs[tc].tfeatures.size() < tfs[tc].intrinsic_features.size())
      {
        TPF temptpf;
        tfs[tc].tfeatures.push_back(temptpf);
      }
      for(unsigned int tp = 0; tp < tfs[tc].intrinsic_features.size(); tp++)
      {
        //	printf("running for trackpoint : %d/%d\n",tc,tp);
        int t = tfs[tc].intrinsic_features[tp].time;

        int x,y,z;
        x = static_cast<int>(tfs[tc].intrinsic_features[tp].Centroid[0]);//*spacing[0]);
        y = static_cast<int>(tfs[tc].intrinsic_features[tp].Centroid[1]);//*spacing[1]);
        z = static_cast<int>(tfs[tc].intrinsic_features[tp].Centroid[2]);//*spacing[2]);

        while(true)
        {
          //	printf("Staring while\n");
          LabelImageType::SizeType lsize;
          LabelImageType::IndexType lindex;
          LabelImageType::IndexType lend;
          LabelImageType::RegionType lregion;

          lindex[0] = MAX(MIN(x-wsize*factor, bound[0]-1),0);
          lindex[1] = MAX(MIN(y-wsize*factor, bound[1]-1),0);
          lindex[2] = MAX(MIN(z-wsize*factor, bound[2]-1),0);

          lend[0] = MAX(MIN(x+wsize*factor, bound[0]-1),0);
          lend[1] = MAX(MIN(y+wsize*factor, bound[1]-1),0);
          lend[2] = MAX(MIN(z+wsize*factor, bound[2]-1),0);


          lsize[0] = lend[0] - lindex[0] + 1;				// 40x40 window centered at the centroid of the cell
          lsize[1] = lend[1] - lindex[1] + 1;
          lsize[2] = lend[2] - lindex[2] + 1;

          lregion.SetSize(lsize);
          lregion.SetIndex(lindex);

          //lregion.Print(std::cout);
          typedef itk::ImageRegionIteratorWithIndex<LabelImageType> LIWI; //Label Image Window Iterator 
          LIWI liter(segmented[t][c-1], lregion);
          int count = 0;
          for(liter.GoToBegin(); !liter.IsAtEnd();++liter)
          {
            if(liter.Get()>0)
              count++;									// sum up the dendritic cell pixels					
          }
          if(count>1000)
          {
            distances.clear();
            for(liter.GoToBegin(); !liter.IsAtEnd();++liter)
            {
              if(liter.Get()>0)
              {
                LabelImageType::IndexType tindex = liter.GetIndex();
                distances.push_back(sqrt(float(x-tindex[0])*(x-tindex[0])*spacing[0]*spacing[0]+(y - tindex[1])*(y-tindex[1])*spacing[1]*spacing[1]+(z-tindex[2])*(z-tindex[2])*spacing[2]*spacing[2]));
              }
            }	// finshied computing the distances of the DCs to the centers of the TCs
            sort(distances.begin(),distances.end());
            float sum = 0;
            int num_num = distances.size()*0.05;
            for(int counter = 0; counter< num_num; counter++)	// take only the 5% closest DC pixels and averge them up.
            {
              sum = sum + distances[counter];
            }
            sum /= num_num;
            radius = pow(static_cast<float>(tfs[tc].intrinsic_features[tp].ScalarFeatures[FeatureType::VOLUME]*spacing[0]*spacing[1]*spacing[2]*3.0/8.0/acos(double(0))),1/3.0f); // approximate the radius of the cell
            printf("Sum = %0.2f radius = %0.2f ratio = %0.3f\n",sum,radius,sum/radius);
            if(sum/radius<= ratio_threshold)
            {
              // it is making contact with a DC
              printf("I set it to true\n");
              tfs[tc].tfeatures[tp].scalars[TPF::HAS_CONTACT_TO_2]=1.0;
            }
            else
            {
              //	printf("I set it to false\n");
              tfs[tc].tfeatures[tp].scalars[TPF::HAS_CONTACT_TO_2]=0.0;
            }
            break;
          }
          else
          {
            if(factor >100)
            {
              //	printf("Error: I still didnt find enough DC points\n");
              tfs[tc].tfeatures[tp].scalars[TPF::HAS_CONTACT_TO_2]=0.0;
              break;
            }
            factor = factor * inc_factor;
            continue;
          }
        }
      }
      tfs[tc].scalars[TF::CONTACT_TO_2] = 0;
      //printf("starting to add dc contact numbers\n");
      for(unsigned int tp = 0; tp < tfs[tc].tfeatures.size(); tp++)
      {
        tfs[tc].scalars[TF::CONTACT_TO_2] += tfs[tc].tfeatures[tp].scalars[TPF::HAS_CONTACT_TO_2];
      }


      //	printf("Finished adding dc contact numbers\n");
      tfs[tc].scalars[TF::CONTACT_TO_2] /=tfs[tc].intrinsic_features.size();
    }


    printf("tfs.size()\n");


  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  InputImageType::Pointer extract_raw_image(float bbox[6],InputImageType::Pointer r)
  {
    int bb[6];
    for(int co = 0; co < 6; co++)
    {
      bb[co] = int(bbox[co]+0.5);
    }
    InputImageType::Pointer lp = InputImageType::New();
    InputImageType::IndexType lindex;
    InputImageType::SizeType lsize;
    InputImageType::RegionType lregion;
    lindex.Fill(0);
    lsize[0] = bb[1]-bb[0]+1;
    lsize[1] = bb[3]-bb[2]+1;
    lsize[2] = bb[5]-bb[4]+1;
    lregion.SetIndex(lindex);
    lregion.SetSize(lsize);
    lp->SetRegions(lregion);
    lp->Allocate();
    lp->FillBuffer(0);
    IteratorType lpiter(lp,lp->GetLargestPossibleRegion());
    lindex[0] = bb[0];
    lindex[1] = bb[2];
    lindex[2] = bb[4];
    lregion.SetIndex(lindex);
    IteratorType liter(r,lregion);
    for(;!lpiter.IsAtEnd(); ++lpiter,++liter)
    {
      lpiter.Set(liter.Get());
    }
    return lp;

  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  LabelImageType::Pointer extract_label_image(int label, float bbox[6],LabelImageType::Pointer l)
  {

    // 	LabelImageType::RegionType test = l->GetLargestPossibleRegion();
    // 	std::cout << std::endl << test;

    int bb[6];
    for(int co = 0; co < 6; co++)
    {
      bb[co] = int(bbox[co]+0.5);
      // 		std::cout << std::endl << bbox[co];
    }
    LabelImageType::Pointer lp = LabelImageType::New();
    LabelImageType::IndexType lindex;
    LabelImageType::SizeType lsize;
    LabelImageType::RegionType lregion;
    lindex.Fill(0);
    lsize[0] = bb[1]-bb[0]+1;
    lsize[1] = bb[3]-bb[2]+1;
    lsize[2] = bb[5]-bb[4]+1;
    lregion.SetIndex(lindex);
    lregion.SetSize(lsize);
    lp->SetRegions(lregion);
    lp->Allocate();
    lp->FillBuffer(0);
    LabelIteratorType lpiter(lp,lp->GetLargestPossibleRegion());
    lindex[0] = bb[0];
    lindex[1] = bb[2];
    lindex[2] = bb[4];
    lregion.SetIndex(lindex);

    // 	std::cout << std::endl << lregion;
    // 	std::cout << std::flush;

    LabelIteratorType liter(l,lregion);
    for(;!lpiter.IsAtEnd(); ++lpiter,++liter)
    {
      if(liter.Get()==label)
        lpiter.Set(liter.Get());
    }
    return lp;
  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  void annotateImage(Color2DImageType::Pointer number,Color2DImageType::Pointer orig, int n, int x, int y)
  {
    //	printf("annotateImage called with n = %d x = %d y = %d ... ", n,x,y);
    typedef itk::ImageRegionConstIterator<Color2DImageType> ConstColor2DIteratorType;
    typedef itk::ImageRegionIterator<Color2DImageType> Color2DIteratorType;

    Color2DImageType::RegionType region;
    Color2DImageType::IndexType index;
    Color2DImageType::SizeType size;

    //printf("n = %d xsize : %d ysize: %d x: %d y = %d\n",n,number->GetLargestPossibleRegion().GetSize()[0],number->GetLargestPossibleRegion().GetSize()[1],x,y);
    double p;
    if(n<10)
      p=3;
    else if (n<100)
      p=3.0/2;
    else 
      p=1;
    size[0]=25/p;
    size[1]=11;
    index[0] = (int(n/20)*34);
    index[1] = number->GetLargestPossibleRegion().GetSize()[1]-((n%20)*31.5)-size[1]+1;
    index[0]= MIN(MAX(index[0],0),number->GetLargestPossibleRegion().GetSize()[0]-size[0]-1);
    index[1]= MIN(MAX(index[1],0),number->GetLargestPossibleRegion().GetSize()[1]-size[1]-1);

    VectorPixelType white; white[0]=255; white[1]=255; white[2]=255;
    VectorPixelType blue; blue[0]=0; blue[1]=0; blue[2]=255;

    region.SetSize(size);
    region.SetIndex(index);

    //printf("index : %d %d\n",index[0],index[1]);
    ConstColor2DIteratorType numberiter(number,region);

    if(x<size[0]/2)
      index[0]=0;
    else
      index[0] = MIN(MAX(x-size[0]/2,0),orig->GetLargestPossibleRegion().GetSize()[0]-size[0]);
    if(y<size[1]/2)
      index[1]=0;
    else
      index[1] = MIN(MAX(y-size[1]/2,0),orig->GetLargestPossibleRegion().GetSize()[1]-size[1]);
    printf("index : %d %d\n",index[0],index[1]);
    region.SetIndex(index);

    Color2DIteratorType origiter(orig,region);


    for(origiter.GoToBegin(),numberiter.GoToBegin();!origiter.IsAtEnd() && !numberiter.IsAtEnd(); ++numberiter,++origiter)
    {
      if(numberiter.Get()!=white)
        origiter.Set(blue);
    }
    //printf("Done\n");
  }

  //---------------------------------------------------------------------------------------------------------------------
  //
  ColorImageType::Pointer getColorImageFromColor2DImages(std::vector<Color2DImageType::Pointer> input)
  {
    ColorImageType::Pointer col = ColorImageType::New();
    ColorImageType::SizeType csize;
    ColorImageType::IndexType cindex;
    ColorImageType::RegionType cregion;
    Color2DImageType::SizeType c2dsize = input[0]->GetLargestPossibleRegion().GetSize();
    csize[0] = c2dsize[0];
    csize[1] = c2dsize[1];
    csize[2] = input.size();
    cindex.Fill(0);
    cregion.SetSize(csize);
    cregion.SetIndex(cindex);
    col->SetRegions(cregion);
    col->Allocate();

    ColorIteratorType citer(col,col->GetLargestPossibleRegion());
    citer.GoToBegin();
    typedef itk::ImageRegionIterator<Color2DImageType> Color2DIteratorType;
    for(int counter=0; counter< input.size(); counter++)
    {
      Color2DIteratorType c2diter(input[counter],input[counter]->GetLargestPossibleRegion());
      c2diter.GoToBegin();
      for(;!c2diter.IsAtEnd(); ++c2diter, ++citer)
      {
        citer.Set(c2diter.Get());
      }
    }
    return col;
  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  void drawLine(ColorImageType::Pointer input, VectorPixelType color1, VectorPixelType color2, int x1, int y1, int z1, int x2, int y2, int z2)
  {
    //z1 has to be = z2
    int upscale = 1;
    ColorImageType::IndexType index1, index2;
    ColorImageType::SizeType size = input->GetLargestPossibleRegion().GetSize();
    index1[0] = MAX(MIN(upscale*x1,size[0]-1),0);
    index1[1] = MAX(MIN(upscale*y1,size[1]-1),0);
    index1[2] = MAX(MIN(z1,size[2]-1),0);
    index2[0] = MAX(MIN(upscale*x2,size[0]-1),0);
    index2[1] = MAX(MIN(upscale*y2,size[1]-1),0);
    index2[2] = MAX(MIN(z2,size[2]-1),0);

    //printf("drawing line...");
    typedef itk::LineIterator<ColorImageType> LineIteratorType;
    LineIteratorType li(input,index1,index2);

    li.GoToBegin();
    int pc=0,pc1 = 0;
    for(;!li.IsAtEnd(); ++li)
    {
      pc++;
    }

    for(li.GoToBegin();!li.IsAtEnd();++li)
    {
      float weight = pc1*1.0/pc;
      VectorPixelType color;
      color[0] = weight*float(color1[0]) + (1-weight)*float(color2[0]);
      color[1] = weight*float(color1[1]) + (1-weight)*float(color2[1]);
      color[2] = weight*float(color2[2]) + (1-weight)*float(color2[2]);
      //printf("color = [%d %d %d] \n",color[0],color[1],color[2]);
      VectorPixelType getcolor = li.Get();
      color[0] = MAX(getcolor[0],color[0]);
      color[1] = MAX(getcolor[1],color[1]);
      color[2] = MAX(getcolor[2],color[2]);
      li.Set(color);
      pc1++;
    }
    //scanf("%*d");
    //printf("\n");
  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  std::vector<FeatureType> get_all_connected_components(LabelImageType::Pointer labelim,FeatureType f)
  {
    int id = f.num;
    std::vector<FeatureType> comps;
    LabelImageType::RegionType lregion;
    LabelImageType::IndexType index;
    index[0] = f.BoundingBox[0];
    index[1] = f.BoundingBox[2];
    index[2] = f.BoundingBox[4];
    LabelImageType::SizeType size;
    size[0] = f.BoundingBox[1] - f.BoundingBox[0] + 1;
    size[1] = f.BoundingBox[3] - f.BoundingBox[2] + 1;
    size[2] = f.BoundingBox[5] - f.BoundingBox[4] + 1;
    lregion.SetSize(size);

    InputImageType::Pointer im = InputImageType::New();
    LabelImageType::IndexType empty;
    empty.Fill(0);
    lregion.SetIndex(empty);
    im->SetRegions(lregion);
    im->Allocate();
    im->FillBuffer(0);

    lregion.SetIndex(index);
    LabelIteratorType liter(labelim,lregion);
    IteratorType iter(im,im->GetLargestPossibleRegion());
    for(iter.GoToBegin(),liter.GoToBegin();!liter.IsAtEnd();++iter,++liter)
    {
      if(liter.Get()==id)
        iter.Set(id);
    }
    ConnectedFilterType::Pointer cfilter = ConnectedFilterType::New();
    cfilter->SetFullyConnected(1);
    cfilter->SetInput(im);
    cfilter->Update();

    getFeatureVectorsFarsight(cfilter->GetOutput(),im,comps,1,1);
    for(int counter=0; counter < comps.size(); counter++)
    {
      comps[counter].Centroid[0] += index[0];
      comps[counter].Centroid[1] += index[1];
      comps[counter].Centroid[2] += index[2];
    }
    return comps;

  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  void MergeCells(std::vector<LabelImageType::Pointer> lin, std::vector<InputImageType::Pointer> imin, std::vector<FeatureType> fin, FeatureVariances fvar, LabelImageType::Pointer &lout, InputImageType::Pointer &rout, FeatureType &fout)
  {
    printf("In MergeCell\n");
    LabelImageType::Pointer p1,p2;
    InputImageType::Pointer r1,r2;
    p1 = lin[0];
    p2 = lin[1];
    r1 = imin[0];
    r2 = imin[1];


    LabelImageType::SizeType ls;

    int lbounds[6];

    lbounds[0] = MIN(fin[0].BoundingBox[0],fin[1].BoundingBox[0]);
    lbounds[2] = MIN(fin[0].BoundingBox[2],fin[1].BoundingBox[2]);
    lbounds[4] = MIN(fin[0].BoundingBox[4],fin[1].BoundingBox[4]);
    lbounds[1] = MAX(fin[0].BoundingBox[1],fin[1].BoundingBox[1]);
    lbounds[3] = MAX(fin[0].BoundingBox[3],fin[1].BoundingBox[3]);
    lbounds[5] = MAX(fin[0].BoundingBox[5],fin[1].BoundingBox[5]);

    ls[0] = lbounds[1]-lbounds[0]+1;
    ls[1] = lbounds[3]-lbounds[2]+1;
    ls[2] = lbounds[5]-lbounds[4]+1;

    LabelImageType::Pointer p = LabelImageType::New();
    InputImageType::Pointer r = InputImageType::New();
    LabelImageType::IndexType lindex;
    lindex.Fill(0);
    LabelImageType::RegionType lregion;
    lregion.SetIndex(lindex);
    lregion.SetSize(ls);
    p->SetRegions(lregion);
    p->Allocate();
    p->FillBuffer(0);
    r->SetRegions(lregion);
    r->Allocate();
    r->FillBuffer(0);
    LabelIteratorType liter1(p1,p1->GetLargestPossibleRegion());
    IteratorType riter1(r1,r1->GetLargestPossibleRegion());


    lindex[0] = fin[0].BoundingBox[0]-lbounds[0];
    lindex[1] = fin[0].BoundingBox[2]-lbounds[2];
    lindex[2] = fin[0].BoundingBox[4]-lbounds[4];

    lregion.SetSize(p1->GetLargestPossibleRegion().GetSize());
    lregion.SetIndex(lindex);

    LabelIteratorType liter(p,lregion);
    IteratorType riter(r,lregion);
    for(liter1.GoToBegin(),riter1.GoToBegin(),liter.GoToBegin(),riter.GoToBegin();!liter1.IsAtEnd(); ++liter1,++riter1,++liter,++riter)
    {
      if(liter1.Get()==fin[0].num)
        liter.Set(255);
      riter.Set(riter1.Get());
    }

    LabelIteratorType liter2(p2,p2->GetLargestPossibleRegion());
    IteratorType riter2(r2,r2->GetLargestPossibleRegion());

    lindex[0] = fin[1].BoundingBox[0]-lbounds[0];
    lindex[1] = fin[1].BoundingBox[2]-lbounds[2];
    lindex[2] = fin[1].BoundingBox[4]-lbounds[4];
    lregion.SetIndex(lindex);
    lregion.SetSize(p2->GetLargestPossibleRegion().GetSize());

    liter = LabelIteratorType(p,lregion);
    riter = IteratorType(r,lregion);

    for(liter2.GoToBegin(),riter2.GoToBegin(),liter.GoToBegin(),riter.GoToBegin();!liter2.IsAtEnd(); ++liter2,++liter,++riter,++riter2)
    {
      if(liter2.Get()==fin[1].num)
        liter.Set(255);
      riter.Set(riter2.Get());
    }


    std::vector<FeatureType> f1;
    getFeatureVectorsFarsight(p,r,f1,fin[0].time,fin[0].tag);

    FeatureType f = f1[0];
    f.Centroid[0]+=lbounds[0];
    f.Centroid[1]+=lbounds[2];
    f.Centroid[2]+=lbounds[4];
    f.WeightedCentroid[0]+=lbounds[0];
    f.WeightedCentroid[2]+=lbounds[2];
    f.WeightedCentroid[4]+=lbounds[4];
    f.BoundingBox[0]+=lbounds[0];
    f.BoundingBox[1]+=lbounds[0];
    f.BoundingBox[2]+=lbounds[2];
    f.BoundingBox[3]+=lbounds[2];
    f.BoundingBox[4]+=lbounds[4];
    f.BoundingBox[5]+=lbounds[5];
    lout = p;
    rout = r;
    fout = f;
    //return f;
    printf("End mergecell\n");
  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  void SplitCell(LabelImageType::Pointer lin, InputImageType::Pointer imin,FeatureType fin, FeatureVariances fvar,std::vector<LabelImageType::Pointer> &lout,std::vector<InputImageType::Pointer> &rout,std::vector<FeatureType> &fvecout)
  {

    // Amin debugging k-means:
    typedef  itk::ImageFileWriter<LabelImageType> WriterType;



    printf("In SplitCell:\n");
    float c1[3],c2[3];
    c1[0] = fin.Centroid[0]-3*(1.0*rand()/RAND_MAX)-fin.BoundingBox[0];
    c1[1] = fin.Centroid[1]-3*(1.0*rand()/RAND_MAX)-fin.BoundingBox[2];
    c1[2] = fin.Centroid[2]-0*(1.0*rand()/RAND_MAX)-fin.BoundingBox[4];
    c2[0] = fin.Centroid[0]+3*(1.0*rand()/RAND_MAX)-fin.BoundingBox[0];
    c2[1] = fin.Centroid[1]+3*(1.0*rand()/RAND_MAX)-fin.BoundingBox[2];
    c2[2] = fin.Centroid[2]+0*(1.0*rand()/RAND_MAX)-fin.BoundingBox[4];

    bool converged = false;
    LabelImageType::Pointer lcopy = LabelImageType::New();
    lcopy->SetRegions(lin->GetLargestPossibleRegion());
    lcopy->Allocate();
    lcopy->FillBuffer(0);
    while(!converged)
    {
      //	printf("In loop\t");
      int num1 =0, num2 = 0;
      LabelImageType::IndexType index1, index2;
      index1.Fill(0);index2.Fill(0);
      typedef itk::ImageRegionIteratorWithIndex<LabelImageType> LabelIteratorWithIndex;
      LabelIteratorWithIndex liter(lin,lin->GetLargestPossibleRegion());
      LabelIteratorType lcopyiter(lcopy,lcopy->GetLargestPossibleRegion());
      for(lcopyiter.GoToBegin(),liter.GoToBegin();!liter.IsAtEnd();++liter,++lcopyiter)
      {
        if(liter.Get()!=0)
        {
          LabelImageType::IndexType index = liter.GetIndex();
          float dist1, dist2;
          dist1 = sqrt((fvar.spacing[0]*(index[0]-c1[0]))*(fvar.spacing[0]*(index[0]-c1[0]))+(fvar.spacing[1]*(index[1]-c1[1]))*(fvar.spacing[1]*(index[1]-c1[1]))+(fvar.spacing[2]*(index[2]-c1[2]))*(fvar.spacing[2]*(index[2]-c1[2])));
          dist2 = sqrt((fvar.spacing[0]*(index[0]-c2[0]))*(fvar.spacing[0]*(index[0]-c2[0]))+(fvar.spacing[1]*(index[1]-c2[1]))*(fvar.spacing[1]*(index[1]-c2[1]))+(fvar.spacing[2]*(index[2]-c2[2]))*(fvar.spacing[2]*(index[2]-c2[2])));
          if(dist1< dist2)
          {
            lcopyiter.Set(1);
            index1[0] = index1[0] + index[0];
            index1[1] = index1[1] + index[1];
            index1[2] = index1[2] + index[2];
            num1++;
          }
          else
          {
            lcopyiter.Set(2);
            index2[0] = index2[0] + index[0];
            index2[1] = index2[1] + index[1];
            index2[2] = index2[2] + index[2];
            num2++;
          }
        }
      }
      LabelImageType::SizeType lsize1 = lin->GetLargestPossibleRegion().GetSize();
      if(num1+num2 == lsize1[0]*lsize1[1]*lsize1[2])
      {
        printf("num1 = %d num2 = %d volume = %d\n",num1,num2,lsize1[0]*lsize1[1]*lsize1[2]);
        // 			scanf("%*d");
      }
      if(num1==0 || num2 == 0)
      {			
        c1[0] = fin.Centroid[0]-3*(1.0*rand()/RAND_MAX)-fin.BoundingBox[0];
        c1[1] = fin.Centroid[1]-3*(1.0*rand()/RAND_MAX)-fin.BoundingBox[2];
        c1[2] = fin.Centroid[2]-0*(1.0*rand()/RAND_MAX)-fin.BoundingBox[4];
        c2[0] = fin.Centroid[0]+3*(1.0*rand()/RAND_MAX)-fin.BoundingBox[0];
        c2[1] = fin.Centroid[1]+3*(1.0*rand()/RAND_MAX)-fin.BoundingBox[2];
        c2[2] = fin.Centroid[2]+0*(1.0*rand()/RAND_MAX)-fin.BoundingBox[4];
        continue;
      }

      float change = sqrt((c1[0] - index1[0]*1.0/num1)*(c1[0] - index1[0]*1.0/num1)+(c1[1] - index1[1]*1.0/num1)*(c1[1] - index1[1]*1.0/num1)+(c1[2] - index1[2]*1.0/num1)*(c1[2] - index1[2]*1.0/num1));
      float change1 = sqrt((c2[0] - index2[0]*1.0/num2)*(c2[0] - index2[0]*1.0/num2)+(c2[1] - index2[1]*1.0/num2)*(c2[1] - index2[1]*1.0/num2)+(c2[2] - index2[2]*1.0/num2)*(c2[2] - index2[2]*1.0/num2));

      printf("change = %f\n", change);
      printf("change1 = %f\n", change1);
      printf("time = %d\n", fin.time);
      printf("id = %d\n", fin.num);


      if(change < 0.02)

        //if(change < 0.02 && change1 < 0.02)
      {
        //std::stringstream ss1,ss2;//create a stringstream
        //ss1 << fin.time;
        //ss2 << fin.num;
        //std::string outfilenametmp = "C:\\Users\\amerouan\\Desktop\\FeaturesTests\\small_im_"+ss1.str()+"_"+ss2.str()+".tiff";
        //WriterType::Pointer writer = WriterType::New();
        //writer->SetFileName(outfilenametmp.c_str());
        //writer->SetInput(lcopy);
        //writer->Update();

        converged = true;			

        std::vector<FeatureType> ftemp;
        _TRACE;
        getFeatureVectorsFarsight(lcopy,imin,ftemp,fin.time,0);
        _TRACE;
        printf("ftemp.size() = %d\n",ftemp.size());
        LabelImageType::Pointer l1,l2;
        l1 = LabelImageType::New();
        l2 = LabelImageType::New();
        LabelImageType::RegionType lregion;
        LabelImageType::SizeType lsize;
        LabelImageType::IndexType lindex;
        lindex.Fill(0);
        lsize[0] = ftemp[0].BoundingBox[1]-ftemp[0].BoundingBox[0]+1;
        lsize[1] = ftemp[0].BoundingBox[3]-ftemp[0].BoundingBox[2]+1;
        lsize[2] = ftemp[0].BoundingBox[5]-ftemp[0].BoundingBox[4]+1;
        lregion.SetSize(lsize);
        lregion.SetIndex(lindex);
        l1->SetRegions(lregion);

        l1->Allocate();
        l1->FillBuffer(0);
        LabelIteratorType loutiter(l1,l1->GetLargestPossibleRegion());

        lindex[0] = ftemp[0].BoundingBox[0];
        lindex[1] = ftemp[0].BoundingBox[2];
        lindex[2] = ftemp[0].BoundingBox[4];

        lregion.SetIndex(lindex);

        //lcopy->Print(std::cout);
        //lregion.Print(std::cout);
        LabelIteratorType liniter(lcopy,lregion);

        for(loutiter.GoToBegin(),liniter.GoToBegin();!liniter.IsAtEnd();++liniter,++loutiter)
        {
          if(liniter.Get()==1)
            loutiter.Set(255);
        }


        lindex.Fill(0);
        lsize[0] = ftemp[1].BoundingBox[1]-ftemp[1].BoundingBox[0]+1;
        lsize[1] = ftemp[1].BoundingBox[3]-ftemp[1].BoundingBox[2]+1;
        lsize[2] = ftemp[1].BoundingBox[5]-ftemp[1].BoundingBox[4]+1;
        lregion.SetSize(lsize);
        lregion.SetIndex(lindex);
        l2->SetRegions(lregion);

        l2->Allocate();
        l2->FillBuffer(0);

        loutiter = LabelIteratorType(l2,l2->GetLargestPossibleRegion());

        lindex[0] = ftemp[1].BoundingBox[0];
        lindex[1] = ftemp[1].BoundingBox[2];
        lindex[2] = ftemp[1].BoundingBox[4];
        lregion.SetIndex(lindex);

        //lcopy->Print(std::cout);
        //lregion.Print(std::cout);
        liniter = LabelIteratorType(lcopy,lregion);
        _TRACE;
        for(loutiter.GoToBegin(),liniter.GoToBegin();!liniter.IsAtEnd();++liniter,++loutiter)
        {
          if(liniter.Get()==2)
            loutiter.Set(255);
        }


        //add lower bounds to centroid and bounding box and extract two images based on bounding box TODO

        for(int counter = 0; counter < 2; counter++)
        {
          ftemp[counter].BoundingBox[0]+=fin.BoundingBox[0];
          ftemp[counter].BoundingBox[1]+=fin.BoundingBox[0];
          ftemp[counter].BoundingBox[2]+=fin.BoundingBox[2];
          ftemp[counter].BoundingBox[3]+=fin.BoundingBox[2];
          ftemp[counter].BoundingBox[4]+=fin.BoundingBox[4];
          ftemp[counter].BoundingBox[5]+=fin.BoundingBox[4];
          ftemp[counter].Centroid[0] +=fin.BoundingBox[0];
          ftemp[counter].Centroid[1] +=fin.BoundingBox[2];
          ftemp[counter].Centroid[2] +=fin.BoundingBox[4];
        }

        lout.push_back(l1);
        lout.push_back(l2);
        InputImageType::Pointer r1,r2;
        r1 = InputImageType::New();
        r2 = InputImageType::New();
        r1->SetRegions(l1->GetLargestPossibleRegion());
        r2->SetRegions(l2->GetLargestPossibleRegion());
        r1->Allocate();
        r2->Allocate();
        rout.push_back(r1);
        rout.push_back(r2);
        fvecout.push_back(ftemp[0]);
        fvecout.push_back(ftemp[1]);
        _TRACE;
        return;
      }
      c1[0] = index1[0]*1.0/num1;
      c1[1] = index1[1]*1.0/num1;
      c1[2] = index1[2]*1.0/num1;

      c2[0] = index2[0]*1.0/num2;
      c2[1] = index2[1]*1.0/num2;
      c2[2] = index2[2]*1.0/num2;
    }


    _TRACE;


  }

  //---------------------------------------------------------------------------------------------------------------------
  // 
  LabelImageType::Pointer fillHoles(LabelImageType::Pointer im, int n)
  {
    //InputImageType::Pointer bin = InputImageType::New();
    //bin->SetRegions(im->GetLargestPossibleRegion());
    //bin->Allocate();
    //bin->FillBuffer(0);

    //IteratorType iter(bin,bin->GetLargestPossibleRegion());
    LabelIteratorType liter(im,im->GetLargestPossibleRegion());

    //for(iter.GoToBegin(),liter.GoToBegin();!iter.IsAtEnd();++iter,++liter)
    //{
    //	if(liter.Get()==0)
    //	{
    //		iter.Set(255);
    //	}
    //}

    //bin = getLargeComponents(bin,n);

    typedef itk::BinaryBallStructuringElement<InputPixelType,3> StructuringElementType;
    typedef itk::Neighborhood<InputPixelType,3> NeighborhoodElementType;
    typedef itk::GrayscaleDilateImageFilter<LabelImageType,LabelImageType,NeighborhoodElementType> DilateFilterType;
    typedef itk::GrayscaleErodeImageFilter<LabelImageType,LabelImageType,NeighborhoodElementType> ErodeFilterType;

    StructuringElementType selement1,selement2;
    NeighborhoodElementType::SizeType size;
    size[0]=n;
    size[1]=n;
    size[2]=1;//FIXME
    selement1.SetRadius(size);
    selement1.CreateStructuringElement();
    selement2.SetRadius(size);
    selement2.CreateStructuringElement();
    DilateFilterType::Pointer dfilter = DilateFilterType::New();
    dfilter->SetKernel(selement1);
    dfilter->SetInput(im);
    dfilter->Update();
    ErodeFilterType::Pointer efilter = ErodeFilterType::New();
    efilter->SetKernel(selement2);
    efilter->SetInput(dfilter->GetOutput());
    efilter->Update();

    LabelImageType::Pointer dilated = efilter->GetOutput();

    LabelImageType::Pointer out = LabelImageType::New();
    out->SetRegions(im->GetLargestPossibleRegion());
    out->Allocate();
    out->FillBuffer(0);

    LabelIteratorType liter1(out,out->GetLargestPossibleRegion());
    LabelIteratorType liter2(dilated,dilated->GetLargestPossibleRegion());
    liter.GoToBegin();
    liter1.GoToBegin();
    liter2.GoToBegin();

    //IteratorType iter1(bin,bin->GetLargestPossibleRegion());
    //iter1.GoToBegin();
    for(;!liter.IsAtEnd(); ++liter,++liter1,++liter2)
    {
      if(liter.Get()!=0)
        liter1.Set(liter.Get());
      else
        liter1.Set(liter2.Get());
    }

    return out;

  }



  //LabelImageType::Pointer removeHoles(LabelImageType::Pointer im,int size)
  //{
  //	LabelImageType::Pointer out = LabelImageType::New();
  //	out->SetRegions(im,im->GetLargestPossibleRegion());
  //	out->Allocate();
  //
  //	LabelIteratorType iter1(im,im->GetLargestPossibleRegion());
  //	LabelIteratorType iter2(out,out->GetLargestPossibleRegion());
  //
  //	for(iter1.GoToBegin(),iter2.GoToBegin();!iter1.IsAtEnd(); ++iter1,++iter2)
  //	{
  //		if(iter1.Get()==0)
  //			iter2.Set(255);
  //	}
  //
  //	out = getLargeLabels(out,size);
  //	iter2 = LabelIteratorType(out,out->GetLargestPossibleRegion());
  //	for(iter1.GoToBegin(),iter2.GoToBegin();!iter1.IsAtEnd();++iter1,++iter2)
  //	{
  //		
  //		iter2.Set(
  //	}
  //}

  //---------------------------------------------------------------------------------------------------------------------
  // 
  int relabelWells(std::vector<LabelImageType::Pointer> & tracked_images, int maxPreviousLabel )
  {
    int maxValue, minValue;
    minValue = 2000000;
    maxValue = 0;
    bool flagMinValue = 0; // In case the imag is blank

    for( int ii=0;ii<tracked_images.size();++ii )
    {
      LabelIteratorType it(tracked_images.at(ii),tracked_images.at(ii)->GetLargestPossibleRegion());
      for(it.GoToBegin();!it.IsAtEnd(); ++it)
      {
        if(it.Get() > 0 )
        {
          minValue = MIN(it.Get(),minValue);
          flagMinValue = 1;
        }
        maxValue = MAX(it.Get(),maxValue);
      }
    }
    if(!flagMinValue)
      return maxPreviousLabel;
    if( maxValue == 0 )
      return maxPreviousLabel;

    int newMaxLabel = maxPreviousLabel + 1;
    bool flagOneLabelNewMax = 0;

    // Relabel
    for( int jj=minValue;jj<=maxValue;++jj )
    {
      for( int ii=0;ii<tracked_images.size();++ii )
      {
        LabelIteratorType it(tracked_images.at(ii),tracked_images.at(ii)->GetLargestPossibleRegion());
        for(it.GoToBegin();!it.IsAtEnd(); ++it)
        {
          if( it.Get() == jj )
          {
            it.Set(newMaxLabel);
            flagOneLabelNewMax = 1;
          }
        }
      }
      if( flagOneLabelNewMax == 1 )
      {
        newMaxLabel++;
      }
      flagOneLabelNewMax = 0;
    }
    newMaxLabel--;

    return newMaxLabel;

  }

} // end of namespace helpers
// 		// Write the output:
// 		for(int t =0; t<tracked_images.size(); t++)
// 		{
// 			writeImage<LabelImageType>(tracked_images[t],trackfnames[t].c_str());
// 		}
