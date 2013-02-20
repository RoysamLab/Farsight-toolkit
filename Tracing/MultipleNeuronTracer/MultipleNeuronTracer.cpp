#include "MultipleNeuronTracer.h"
#include <ctime>

#ifdef _OPENMP
#include "omp.h"
#endif

MultipleNeuronTracer::MultipleNeuronTracer()
{
}

MultipleNeuronTracer::~MultipleNeuronTracer()
{
}

void MultipleNeuronTracer::LoadParameters(const char* parametersFileName,int _argc)
{
  std::map<std::string, std::string> opts;  

  if(_argc == 6){
    std::cout << " _argc: " << _argc << std::endl;
    this->optionsCreate(parametersFileName, opts);
  }

  std::map<std::string,std::string>::iterator mi;

  mi = opts.find("-intensity_threshold"); 
  if(!this->_isCoverageOptimized){
    if(mi!=opts.end())
    { std::istringstream ss((*mi).second); ss>>this->intensity_threshold; 
    }
    else
    { this->intensity_threshold = 0.005; printf("Chose intensity_threshold = 0.005 as default\n");}
  }
  else
    this->intensity_threshold += 0.02;

  if(!this->_isCoverageOptimized){
    mi = opts.find("-contrast_threshold");
    if(mi!=opts.end())
    { std::istringstream ss((*mi).second); ss>>this->contrast_threshold; }
    else
    {	  this->contrast_threshold = 0.0003; printf("Chose contrast_threshold = 0.0003 as default\n"); }
  }
  else
    this->contrast_threshold += 0.0004;

  mi = opts.find("-cost_threshold"); 
  if(mi!=opts.end())
  { std::istringstream ss((*mi).second); ss>>this->cost_threshold; }
  else
  { this->cost_threshold = 700; printf("Chose cost_threshold = 700 as default\n");}

  mi = opts.find("-debris_threshold"); 
  if(mi!=opts.end())
  { std::istringstream ss((*mi).second); ss>>this->debris_threshold; }
  else
  { this->debris_threshold = 0.8; printf("Chose debris_threshold = 0.8 as default\n"); }

  mi = opts.find("-offshoot"); 
  if(mi!=opts.end())
  { std::istringstream ss((*mi).second); ss>>this->offshoot; }
  else
  { this->offshoot = 10; printf("Chose offshoot = 10 as default\n"); }

  mi = opts.find("-device"); 
  if(mi!=opts.end())
  { std::istringstream ss((*mi).second); ss>>this->device; }
  else
  { this->device = 1; printf("Chose device = 0 as default\n"); }

  mi = opts.find("-mu"); 
  if(mi!=opts.end())
  { std::istringstream ss((*mi).second); ss>>this->mu; }
  else
  { this->mu = 100; printf("Chose mu = 100 as default\n"); }

  mi = opts.find("-no_of_iteration"); 
  if(mi!=opts.end())
  { std::istringstream ss((*mi).second); ss>>this->noOfIteration; }
  else
  { this->noOfIteration = 15; printf("Chose noOfIteration = 15 as default\n"); }

  mi = opts.find("-tracing_type");  // 1 for LOG; 2 for GVF - Default - GVF Tracing
  if(mi!=opts.end())
  { std::istringstream ss((*mi).second); ss>>this->tracing_type; }
  else
  { this->tracing_type = 2; printf("Chose tracing_type = 2(GVF Tracing) as default\n"); }

  debug = true;
  std::cout<<"tracing_type="<<this->tracing_type<<std::endl;
  std::cout<<"intensity_threshold="<<this->intensity_threshold<<std::endl;
  std::cout<<"contrast_threshold="<<this->contrast_threshold<<std::endl;
  std::cout<<"cost_threshold="<<this->cost_threshold<<std::endl;
  std::cout<<"debris_threshold="<<this->debris_threshold<<std::endl;
  std::cout<<"offshoot="<<this->offshoot<<std::endl;
  std::cout<<"device="<<this->device<<std::endl;
  std::cout<<"mu="<<this->mu<<std::endl;
  std::cout<<"no_of_iteration="<<this->noOfIteration<<std::endl;
  

}
void MultipleNeuronTracer::LoadParameters_1(const char* parametersFileName,float intensityThreshold,float contrastThreshold,int costThreshold)
{
  std::map<std::string, std::string> opts;  

  std::map<std::string,std::string>::iterator mi;

  this->intensity_threshold = intensityThreshold;
  std::cout<<"intensity_threshold "<<intensity_threshold<<std::endl;	

  this->contrast_threshold = contrastThreshold; 
  std::cout<<"contrast_threshold "<<contrastThreshold<<std::endl;	

  mi = opts.find("-cost_threshold"); 
  if(mi!=opts.end())
  { std::istringstream ss((*mi).second); ss>>this->cost_threshold; }
  else
  { this->cost_threshold = costThreshold; printf("Chose cost_threshold = 700 as default\n");}

  mi = opts.find("-debris_threshold"); 
  if(mi!=opts.end())
  { std::istringstream ss((*mi).second); ss>>this->debris_threshold; }
  else
  { this->debris_threshold = 0.8; printf("Chose debris_threshold = 0.8 as default\n"); }

  mi = opts.find("-offshoot"); 
  if(mi!=opts.end())
  { std::istringstream ss((*mi).second); ss>>this->offshoot; }
  else
  { this->offshoot = 10; printf("Chose offshoot = 10 as default\n"); }

  mi = opts.find("-device"); 
  if(mi!=opts.end())
  { std::istringstream ss((*mi).second); ss>>this->device; }
  else
  { this->device = 1; printf("Chose device = 0 as default\n"); }

  std::cout<<"intensity_threshold="<<this->intensity_threshold<<std::endl;
  std::cout<<"contrast_threshold="<<this->contrast_threshold<<std::endl;
  std::cout<<"cost_threshold="<<this->cost_threshold<<std::endl;
  std::cout<<"debris_threshold="<<this->debris_threshold<<std::endl;
  std::cout<<"offshoot="<<this->offshoot<<std::endl;
  std::cout<<"device="<<this->device<<std::endl;

}
void MultipleNeuronTracer::LoadCurvImage(std::string fname, unsigned int pad) 
{
  std::cout << "Reading input file "<< fname << std::endl;
  ReaderType::GlobalWarningDisplayOff();
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(fname);
  ImageType3D::Pointer image = reader->GetOutput();
  image->Update();

  ////Binarize the image and use as this as the mask instead of the intensity and contrast thresholds
  //MinMaxImageCalculatorType::Pointer minMaxImageCalFilter = MinMaxImageCalculatorType::New();
  //minMaxImageCalFilter->SetImage(image);
  //minMaxImageCalFilter->Compute();	
  //

  //OtsuThresholdImageFilterType::Pointer otsuThresholdImageFilter = OtsuThresholdImageFilterType::New();
  //otsuThresholdImageFilter->SetInput(image);
  //otsuThresholdImageFilter->Update();
  //std::cout << (int)(otsuThresholdImageFilter->GetThreshold()) << std::endl;
  //
  //float lowerThreshold = otsuThresholdImageFilter->GetThreshold();
  //lowerThreshold = lowerThreshold *0.9;

  //BinaryThresholdImageFilterType::Pointer thresholdFilter = BinaryThresholdImageFilterType::New();
  //thresholdFilter->SetInput(image);
  //thresholdFilter->SetLowerThreshold(lowerThreshold);
  //thresholdFilter->SetUpperThreshold(minMaxImageCalFilter->GetMaximum());
  //thresholdFilter->SetInsideValue(255);
  //thresholdFilter->SetOutsideValue(0);
  //thresholdFilter->Update();




  //_MaskedImage = ImageType3D::New();
  //_MaskedImage = thresholdFilter->GetOutput();
  //

  //itk::CastImageFilter< ImageType3D, CharImageType3D>::Pointer caster2 = itk::CastImageFilter< ImageType3D, CharImageType3D>::New();
  //caster2->SetInput(rescaler2->GetOutput());

  //imgIn

  //	float alpha_B, alpha_F, P_I;
  //alpha_B = alpha_F = P_I = 0;	
  ////MinErrorThresholding(imgIn, &alpha_B, &alpha_F, &P_I, R, C, 1, shd,imgOut); 		
  //float alpha_C, P_I2;
  //alpha_C = P_I2 = 0.0;
  //threeLevelMinErrorThresh(imgIn, &alpha_B, &alpha_F, &alpha_C, &P_I, &P_I2, R, C, 1);
  //threeLevelMinErrorThresh()

  std::cout << "Entering LoadCurvImage" << std::endl;

  // add the code for Huang Threshold 
  // get the Masked Image
  LoadCurvImage_1(image, pad);
}

void MultipleNeuronTracer::LoadCurvImage_1(ImageType3D::Pointer &image, unsigned int pad)  
{
  _flagPipeline = false; // By default pipeline off
  _flagOutLog = false;
  ImageType3D::Pointer CurvImage = image;
  _padz = pad;

  RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetOutputMinimum(0.0);
  rescaler->SetOutputMaximum(1.0);
  rescaler->SetInput(CurvImage);

  //Median filter
  std::cout << "Running Median Filter" << std::endl;
  MedianFilterType::Pointer medfilt = MedianFilterType::New();
  // 	medfilt->SetNumberOfThreads(16);
  medfilt->SetInput(rescaler->GetOutput());
  ImageType3D::SizeType rad = { {1, 1, 1} };
  medfilt->SetRadius(rad);
  medfilt->Update();
  CurvImage = medfilt->GetOutput();

  //pad z slices
  std::cout << "pad z slices" << std::endl;
  itk::Size<3> isz = CurvImage->GetBufferedRegion().GetSize();
  itk::Size<3> osz = isz;
  osz[2] += 2*_padz;
  itk::Index<3> indx, ondx;

  _PaddedCurvImage = ImageType3D::New();
  _PaddedCurvImage->SetRegions(osz);
  _PaddedCurvImage->Allocate();
  _PaddedCurvImage->SetSpacing(CurvImage->GetSpacing());

  for(ondx[2] = 0; ondx[2] < osz[2]; ++ondx[2]) 
  {
    indx[2] = (ondx[2] < _padz) ? 0 : ondx[2] - _padz;
    indx[2] = (ondx[2] >= osz[2]-_padz) ? isz[2]-1 : indx[2];
    for(ondx[1] = 0; ondx[1] < osz[1]; ++ondx[1]) 
    {
      indx[1] = ondx[1];
      for(ondx[0] = 0; ondx[0] < osz[0]; ++ondx[0]) 
      {
        indx[0] = ondx[0];
        _PaddedCurvImage->SetPixel(ondx, CurvImage->GetPixel(indx));
      }
    }
  }

  std::cout << "Input file size (after zero padding) is " << _PaddedCurvImage->GetBufferedRegion().GetSize() << std::endl;
  _size = _PaddedCurvImage->GetBufferedRegion().GetSize();
  //CurvImage->Delete();
}


void MultipleNeuronTracer::LoadCurvImage_2(ImageType3D::Pointer &image)
{
  _flagPipeline = false; // By default pipeline off
  _flagOutLog = false;
  ImageType3D::Pointer CurvImage = image;

  unsigned int padz = 0;

  //pad z slices
  std::cout << "pad z slices" << std::endl;
  itk::Size<3> isz = CurvImage->GetBufferedRegion().GetSize();
  itk::Size<3> osz = isz;
  osz[2] += 2*padz;
  itk::Index<3> indx, ondx;

  _PaddedCurvImage = ImageType3D::New();
  _PaddedCurvImage->SetRegions(osz);
  _PaddedCurvImage->Allocate();
  _PaddedCurvImage->SetSpacing(CurvImage->GetSpacing());

  for(ondx[2] = 0; ondx[2] < osz[2]; ++ondx[2]) 
  {
    indx[2] = (ondx[2] < padz) ? 0 : ondx[2] - padz;
    indx[2] = (ondx[2] >= osz[2]-padz) ? isz[2]-1 : indx[2];
    for(ondx[1] = 0; ondx[1] < osz[1]; ++ondx[1]) 
    {
      indx[1] = ondx[1];
      for(ondx[0] = 0; ondx[0] < osz[0]; ++ondx[0]) 
      {
        indx[0] = ondx[0];
        _PaddedCurvImage->SetPixel(ondx, CurvImage->GetPixel(indx));
      }
    }
  }

  std::cout << "Input file size (after zero padding) is " << _PaddedCurvImage->GetBufferedRegion().GetSize() << std::endl;
  _size = _PaddedCurvImage->GetBufferedRegion().GetSize();
}

ObjectnessMeasures_micro::ObjectnessMeasures_micro(){

  this->alpha = 0.5;
  this->beta = 0.5;
  this->gamma = 0.25;

  this->sigma_min = 0.5;
  this->sigma_max = 2.0;
  this->sigma_intervals = 1;
  this->objectness_type = 1;

  this->ballness = 0.0;
  this->plateness = 0.0;
  this->vesselness = 0.0;
  this->noiseness = 0.0;
}

ObjectnessMeasures_micro::ObjectnessMeasures_micro(float alpha, float beta, float gamma){

  this->alpha = alpha;
  this->beta = beta;
  this->gamma = gamma;

  this->sigma_min = 0.5;
  this->sigma_max = 2.0;
  this->sigma_intervals = 1;
  this->objectness_type = 1;

  this->ballness = 0.0;
  this->plateness = 0.0;
  this->vesselness = 0.0;
  this->noiseness = 0.0;
}

ObjectnessMeasures_micro::ObjectnessMeasures_micro(float sigma_min, float sigma_max, float sigma_intervals, int obj_type){

  this->alpha = 0.5;
  this->beta = 0.5;
  this->gamma = 0.25;

  this->sigma_min = sigma_min;
  this->sigma_max = sigma_max;
  this->sigma_intervals = sigma_intervals;
  this->objectness_type = obj_type;

  this->ballness = 0.0;
  this->plateness = 0.0;
  this->vesselness = 0.0;
  this->noiseness = 0.0;
}

void MultipleNeuronTracer::Set_isCoverageOptimized(bool opt_cov){
  this->_isCoverageOptimized = opt_cov;
}

void MultipleNeuronTracer::OptimizeCoverage(std::string coverageFileName, bool writeResult){

  std::cout << std::endl<< "Optimizing feature coverage for the image." << std::endl;

  // Optimizing coverage at single scale	
  float sigma = 5.6569f;
  float sigma_min = sigma;
  float sigma_max = sigma;
  int sigma_intervals = 1;

  double intensity_weight = 0.2; //0.1; //0.4; //0.2;
  double objectness_weight = 1.0 - intensity_weight;

  StatisticsFilterType::Pointer stats_filter = StatisticsFilterType::New();
  stats_filter->SetInput(this->_PaddedCurvImage);
  stats_filter->Update();
  double img_max_val = stats_filter->GetMaximum();


  //ObjectnessMeasures_micro obj_measures(sigma_min, sigma_max, sigma_intervals, 0); // use this for astrocytes
  ObjectnessMeasures_micro obj_measures(sigma_min, sigma_max, sigma_intervals, 1); // use this for microglia
  obj_measures.alpha = 0.5 * img_max_val;
  obj_measures.beta = 0.5 * img_max_val;
  obj_measures.gamma = 0.25 * img_max_val; //0.25 * img_max_val;

  //std::cout << "Max image value: " << img_max_val << " Mean img value: " << img_mean_val << std::endl;

  this->ComputeObjectnessImage(obj_measures);

  //Write out the vesselenss image
  if(writeResult){
    std::string vesselnessPointsFileName = coverageFileName;
    vesselnessPointsFileName.erase(vesselnessPointsFileName.length()-4, vesselnessPointsFileName.length());
    vesselnessPointsFileName.append("_vesselness.mhd");

    //std::cout << vesselnessPointsFileName << std::endl;

    typedef itk::ImageFileWriter<ImageType3D> ImageWriterType;
    ImageWriterType::Pointer image_writer = ImageWriterType::New();
    image_writer->SetFileName(vesselnessPointsFileName);
    image_writer->SetInput(this->ObjectnessImage);
    image_writer->Update();
  }

  StatisticsFilterType::Pointer stats_filter2 = StatisticsFilterType::New();
  stats_filter2->SetInput(this->ObjectnessImage);
  stats_filter2->Update();
  double max_objectness = stats_filter2->GetMaximum();

  MultiplyImageFilter::Pointer img_multiplier = MultiplyImageFilter::New();
  img_multiplier->SetInput(this->ObjectnessImage);
  img_multiplier->SetConstant2(1.0/(max_objectness));
  img_multiplier->Update();

  MultiplyImageFilter::Pointer img_multiplier2 = MultiplyImageFilter::New();
  img_multiplier2->SetInput(this->_PaddedCurvImage);
  img_multiplier2->SetConstant2(1.0/(img_max_val));
  img_multiplier2->Update();

  ImageType3D::Pointer normalized_img = img_multiplier2->GetOutput();

  ImageType3D::Pointer normalized_obj_img = img_multiplier->GetOutput();

  StatisticsFilterType::Pointer stats_filter3 = StatisticsFilterType::New();
  stats_filter3->SetInput(normalized_obj_img);
  stats_filter3->Update();
  double mean_objectness = stats_filter3->GetMean();
  double std_objectness = stats_filter3->GetVariance();

  StatisticsFilterType::Pointer stats_filter4 = StatisticsFilterType::New();
  stats_filter4->SetInput(normalized_img);
  stats_filter4->Update();
  double img_mean_val = stats_filter4->GetMean();
  double img_std_val = stats_filter4->GetVariance();


  MultiplyImageFilter::Pointer img_multiplier3 = MultiplyImageFilter::New();
  img_multiplier3->SetInput1(normalized_img);
  img_multiplier3->SetInput2(normalized_obj_img);
  img_multiplier3->Update();
  ImageType3D::Pointer int_obj_prod_img = img_multiplier3->GetOutput();

  MultiplyImageFilter::Pointer img_multiplier4 = MultiplyImageFilter::New();
  img_multiplier4->SetInput(normalized_img);
  img_multiplier4->SetConstant2(intensity_weight);
  img_multiplier4->Update();
  MultiplyImageFilter::Pointer img_multiplier5 = MultiplyImageFilter::New();
  img_multiplier5->SetInput(normalized_obj_img);
  img_multiplier5->SetConstant2(objectness_weight);
  img_multiplier5->Update();
  AddImageFilter::Pointer img_adder = AddImageFilter::New();
  img_adder->SetInput1(img_multiplier4->GetOutput());
  img_adder->SetInput2(img_multiplier5->GetOutput());
  img_adder->Update();
  ImageType3D::Pointer int_obj_sum_img = img_adder->GetOutput();

  StatisticsFilterType::Pointer stats_filter5 = StatisticsFilterType::New();
  stats_filter5->SetInput(int_obj_sum_img);
  stats_filter5->Update();
  double added_img_mean_val = stats_filter5->GetMean();

  typedef itk::LaplacianRecursiveGaussianImageFilter< ImageType3D , ImageType3D> LoGFilterType;
  LoGFilterType::Pointer gauss = LoGFilterType::New();
  gauss->SetInput( _PaddedCurvImage );
  gauss->SetSigma( sigma );
  gauss->SetNormalizeAcrossScale(false);
  gauss->GetOutput()->Update();

  float tot = 0.0f, num = 0.0f;
  itk::ImageRegionIterator<ImageType3D> ittemp(gauss->GetOutput(), gauss->GetOutput()->GetBufferedRegion());
  float gamma = 1.6f;
  float tnorm = vcl_pow(sigma,gamma);
  for(ittemp.GoToBegin(); !ittemp.IsAtEnd(); ++ittemp){
    float q = ittemp.Get()*tnorm;
    ittemp.Set(-1.0f*q);
    tot += q*q;
    num ++;
  }
  //std::cout << "Scale "<< sigma << " had average Energy: " << tot <<std::endl;


  // set the diagonal terms in neighborhood iterator
  itk::Offset<3>
    xp =  {{2 ,  0 ,   0}},
       xn =  {{-2,  0,    0}},
       yp =  {{0,   2,   0}},
       yn =  {{0,  -2,    0}},
       zp =  {{0,   0,    2}},
       zn =  {{0,   0,   -2}};

  itk::Size<3> rad = {{1,1,1}};
  itk::NeighborhoodIterator<ImageType3D> nit(rad , gauss->GetOutput(), gauss->GetOutput()->GetBufferedRegion());
  itk::ImageRegionIterator<ImageType3D> it(gauss->GetOutput(), gauss->GetOutput()->GetBufferedRegion());

  unsigned int
    xy1 =  17, //{ 1 ,   1 ,  0 },
        xy2 =  9,  //{ -1,  -1 ,  0 },
        xy3 =  15, //{ -1,   1 ,  0 },
        xy4 =  11, //{ 1 ,  -1 ,  0 },

        yz1 =  25, //{ 0 ,   1 ,  1 },
        yz2 =  1,  //{ 0 ,  -1 , -1 },
        yz3 =  19, //{ 0 ,  -1 ,  1 },
        yz4 =  7,  //{ 0 ,   1 , -1 },

        xz1 =  23, //{ 1 ,   0 ,  1 },
        xz2 =  3,  //{-1 ,   0 , -1 },
        xz3 =  21, //{-1 ,   0 ,  1 },
        xz4 =  5;  //{ 1 ,   0 , -1 };

  typedef itk::FixedArray< double, 3 > EigenValuesArrayType;
  typedef itk::Matrix< double, 3, 3 > EigenVectorMatrixType;
  typedef itk::SymmetricSecondRankTensor<double,3> TensorType;

  itk::Size<3> sz = _PaddedCurvImage->GetBufferedRegion().GetSize();
  sz[0] = sz[0] - 3;
  sz[1] = sz[1] - 3; 
  sz[2] = sz[2] - 3;

  itk::Vector<float,3> sp = _PaddedCurvImage->GetSpacing();


  // This parameter will affect the density og LoG points detected at each scale.
  float win_scale = 2;

  long win = long(sigma)/2;
  //long win = win_scale * long(sigma);
  if (win <2) 
    win = 2;

  int opt_iter = 0, max_opt_iter = 10; //20;
  float max_coverage = 0.002; //0.01; //0.0;
  float min_coverage = 0.001; // 0.0;
  float coverage_upper_limit = max_coverage + (0.2*max_coverage); //(0.05*max_coverage);
  float coverage_lower_limit = min_coverage - (0.25*min_coverage); //(0.2*min_coverage); //(0.05*min_coverage);

  std::cout << "Coverage limits: [" << coverage_lower_limit << ", " << coverage_upper_limit << "] " << std::endl;

  float thresh1_step_size = 0.02;
  float thresh2_step_size = 0.0004;

  float thresh1 = 0.01; //0.03; //0.01; //0.05; //0.08; //0.005; //0.03   // 3% of maximum theshold from Lowe 2004
  float thresh2 = 0.0001; //0.0005; //0.0001; //0.0009; //0.015; //0.0003;  //0.001 -0.1 percent of range

  double Mi_coverage4 = 0.0;
  double Mi_coverage5 = 0.0;
  double Mi_coverage6 = 0.0;


  std::ofstream optim_internal_file;

  if(writeResult){
    std::string internalOptimFileName = coverageFileName;
    internalOptimFileName.erase(internalOptimFileName.length()-4, internalOptimFileName.length());
    internalOptimFileName.append("_internal.txt");
    optim_internal_file.open(internalOptimFileName.c_str(), std::ios::out);
  }

  while(opt_iter < max_opt_iter){

    _NDXImage = ImageType3D::New();
    _NDXImage->SetRegions(_PaddedCurvImage->GetBufferedRegion());
    _NDXImage->Allocate();
    _NDXImage->FillBuffer(0.0f);

    it.GoToBegin();
    nit.GoToBegin();

    long ctCnt = 0;
    long rejectCtCnt = 0;

    float total_fg_objectness = 0.0;
    float total_fg_intensity = 0.0;
    float total_fg_obj_int = 0.0;
    float total_fg_wt_obj_int = 0.0;

    while(!nit.IsAtEnd()){

      itk::Index<3> ndx = it.GetIndex();
      if ( (ndx[0] < 2) || (ndx[1] < 2) || (ndx[2] < 2) ||
          (ndx[0] > (unsigned int)sz[0]) || (ndx[1] > (unsigned int)sz[1]) ||
          (ndx[2] > (unsigned int)sz[2]) ){
        ++it;
        ++nit;
        continue;
      }

      float a1 = 0.0;
      for (unsigned int i=0; i < 13; ++i)
        a1 += vnl_math_max(nit.GetPixel(i), nit.GetPixel(26 - i));

      float val = nit.GetPixel(13);

      if ( ((val - a1/13.0f) > thresh2 ) && ( val > thresh1 )){

        TensorType h;
        h[0] = gauss->GetOutput()->GetPixel( ndx + xp ) + gauss->GetOutput()->GetPixel( ndx + xn ) - 2*nit.GetPixel( 13 );
        h[3] = gauss->GetOutput()->GetPixel( ndx + yp ) + gauss->GetOutput()->GetPixel( ndx + yn ) - 2*nit.GetPixel( 13 );
        h[5] = gauss->GetOutput()->GetPixel( ndx + zp ) + gauss->GetOutput()->GetPixel( ndx + zn ) - 2*nit.GetPixel( 13 );
        h[1] = nit.GetPixel(xy1) + nit.GetPixel(xy2) - nit.GetPixel(xy3) - nit.GetPixel(xy4);
        h[2] = nit.GetPixel(xz1) + nit.GetPixel(xz2) - nit.GetPixel(xz3) - nit.GetPixel(xz4);
        h[4] = nit.GetPixel(yz1) + nit.GetPixel(yz2) - nit.GetPixel(yz3) - nit.GetPixel(yz4);

        EigenValuesArrayType ev;
        EigenVectorMatrixType em;
        h.ComputeEigenAnalysis (ev, em);


        unsigned int w;
        if(true){ //(IsSeed(ev, w)){

          float value = vnl_math_abs(ev[0]) + vnl_math_abs(ev[1]) + vnl_math_abs(ev[2]) - vnl_math_abs(ev[w]); // How is this score derived?
          if (RegisterIndex(value, ndx, sz, win)){


            _NDXImage->SetPixel(ndx, value);
            ctCnt++;


            total_fg_objectness += normalized_obj_img->GetPixel(ndx); //this->ObjectnessImage->GetPixel(ndx);
            total_fg_intensity += normalized_img->GetPixel(ndx);  //this->_PaddedCurvImage->GetPixel(ndx);
            //total_fg_obj_int += this->ObjectnessImage->GetPixel(ndx) * this->_PaddedCurvImage->GetPixel(ndx);
            total_fg_obj_int += normalized_obj_img->GetPixel(ndx) * this->_PaddedCurvImage->GetPixel(ndx);
            total_fg_wt_obj_int += int_obj_sum_img->GetPixel(ndx);

          }	
          else
            rejectCtCnt++;
        }
        }
        ++it;
        ++nit;
      }

      double Mi_coverage = (double)ctCnt / img_mean_val;
      double Mi_coverage2 = (double)ctCnt / mean_objectness;
      double Mi_coverage3 = (double)ctCnt / (img_mean_val * mean_objectness);

      ImageType3D::SizeType im_size = this->_PaddedCurvImage->GetBufferedRegion().GetSize();
      double total_objectness = mean_objectness * (im_size[0]*im_size[1]*im_size[2]);
      double total_intensity = img_mean_val * (im_size[0]*im_size[1]*im_size[2]);
      double total_int_obj_sum = added_img_mean_val * (im_size[0]*im_size[1]*im_size[2]);

      if(added_img_mean_val < 0.0001)
        total_int_obj_sum = 0.001 * (im_size[0]*im_size[1]*im_size[2]);


      Mi_coverage4 = total_fg_objectness / total_objectness;
      Mi_coverage5 = total_fg_intensity / total_intensity;
      Mi_coverage6 = total_fg_wt_obj_int / total_int_obj_sum;

      std::cout << "num: " << total_fg_wt_obj_int << " den: " << added_img_mean_val << std::endl;

      std::cout << std::endl;
      std::cout << "Number of CTs at this stage: " << ctCnt <<std::endl;
      std::cout << "Number of CTs rejected by RegisterIndex() are: " << rejectCtCnt << std::endl;

      //std::cout << "Total foreground objectness: " << total_fg_objectness << std::endl;
      //std::cout << "Total foreground intensity: " << total_fg_intensity << std::endl;
      //std::cout << "Total foreground obj*intensity: " << total_fg_obj_int << std::endl;
      //std::cout << "Mi_coverage: " << Mi_coverage << std::endl;
      //std::cout << "Mi_coverage2: " << Mi_coverage2 << std::endl;
      //std::cout << "Mi_coverage3: " << Mi_coverage3 << std::endl;
      //std::cout << "Mi_coverage4: " << Mi_coverage4 << std::endl;
      //std::cout << "Mi_coverage5: " << Mi_coverage5 << std::endl;
      //std::cout << "Mi_coverage6: " << Mi_coverage6 << std::endl;

      std::cout << "Iter: " << opt_iter << " Mi_coverage6: " << Mi_coverage6 << " Intensity threshold: " << thresh1 << ", Contrast threshold: " << thresh2 << std::endl;

      if(writeResult){
        if(optim_internal_file.good())
          optim_internal_file << Mi_coverage6 << '\t' << thresh1 << '\t' << thresh2 << '\t' << ctCnt << std::endl;
        else
          std::cout << "Error writing coverage internal file. " << std::endl;
      }

      opt_iter++;


      if(Mi_coverage6 > coverage_upper_limit){
        thresh1 += thresh1_step_size;
        thresh2 += thresh2_step_size;
      }
      else if(Mi_coverage6 < coverage_upper_limit && Mi_coverage6 > max_coverage){
        thresh1 += thresh1_step_size / 5.0;
        thresh2 += thresh2_step_size / 5.0;
      }
      else if(Mi_coverage6 < coverage_lower_limit){

        if(thresh1 <= 0.01){
          thresh1 -= thresh1_step_size / 5.0;
          thresh2 -= thresh2_step_size / 5.0;
        }
        else if(thresh1 <= 0.0){
          thresh1 = 0.0;
          thresh2 = 0.0;
          break;
        }
        else{
          thresh1 -= thresh1_step_size;
          thresh2 -= thresh2_step_size;
        }
      }
      else
        break;		
    }

    if(writeResult)
      optim_internal_file.close();

    // Setting thresholds based on the optimized coverage
    this->intensity_threshold = thresh1;
    this->contrast_threshold = thresh2;

    // Printing the results
    if(opt_iter == max_opt_iter){
      std::cout << "Coverage might not be correctly optimized. Manual adjustments might be needed. " << std::endl;
      this->_isCoverageOptimized = false;
    }
    else if(thresh1 = 0.0){
      std::cout << "Coverage might not be correctly optimized. Manual adjustments might be needed. " << std::endl;
      this->_isCoverageOptimized = false;
    }
    else{
      std::cout << "Done with optimizing coverage. " <<  std::endl;
      this->_isCoverageOptimized = true;
    }


    if(writeResult){

      std::ofstream coverage_file;
      coverage_file.open(coverageFileName.c_str(), std::ios::out);
      if(coverage_file.good()){
        coverage_file << "intensity threshold: " /*<< Mi_coverage6 << std::endl*/ << this->intensity_threshold + 0.02 << std::endl << "contrast threshold: " << this->contrast_threshold + 0.0004 << std::endl;
        //coverage_file << img_mean_val << std::endl << img_std_val << std::endl << mean_objectness << std::endl << std_objectness << std::endl;

        coverage_file.close();
      }
      else
        std::cout << "Error writing coverage file. " << std::endl;


      std::string LoGPointsFileName = coverageFileName;
      LoGPointsFileName.erase(LoGPointsFileName.length()-4, LoGPointsFileName.length());
      LoGPointsFileName.append("_optimum_LOG_points.tif");

      RescalerType::Pointer rescaler2 = RescalerType::New();
      rescaler2->SetInput(_NDXImage);
      rescaler2->SetOutputMaximum( 255 );
      rescaler2->SetOutputMinimum( 0 );
      rescaler2->Update();
      itk::CastImageFilter< ImageType3D, CharImageType3D>::Pointer caster2 = itk::CastImageFilter< ImageType3D, CharImageType3D>::New();
      caster2->SetInput(rescaler2->GetOutput());

      itk::ImageFileWriter< CharImageType3D >::Pointer LoGwriter2 = itk::ImageFileWriter< CharImageType3D >::New();

      LoGwriter2->SetFileName(LoGPointsFileName.c_str());	
      LoGwriter2->SetInput(caster2->GetOutput());
      LoGwriter2->Update();			
    }
  }

  void MultipleNeuronTracer::ComputeObjectnessImage(ObjectnessMeasures_micro obj_measures){

    //float sigma_min = 2.0f; //0.5f;
    //float sigma_max = 10.0f; //4.0f;
    //int sigma_steps = 5;

    //float alpha = 0.5, beta = 0.5, gamma = 0.25; //5.0;

    //int obj_dim = objectness_type; //1; //0: Blobness, 1: Vesselness, 2: Plateness

    //////MultiScaleHessianFilterType::Pointer multi_scale_Hessian = MultiScaleHessianFilterType::New();
    //////multi_scale_Hessian->SetInput(this->_PaddedCurvImage);
    //////multi_scale_Hessian->SetSigmaMin(obj_measures.sigma_min);
    //////multi_scale_Hessian->SetSigmaMax(obj_measures.sigma_max);
    //////multi_scale_Hessian->SetNumberOfSigmaSteps(obj_measures.sigma_intervals);

    ////////ObjectnessFilterType::Pointer objectness_filter = ObjectnessFilterType::New();
    //////ObjectnessFilterType::Pointer objectness_filter = multi_scale_Hessian->GetHessianToMeasureFilter();
    //////
    //////objectness_filter->SetScaleObjectnessMeasure(false);
    //////objectness_filter->SetBrightObject(true);
    //////objectness_filter->SetAlpha(obj_measures.alpha);
    //////objectness_filter->SetBeta(obj_measures.beta);
    //////objectness_filter->SetGamma(obj_measures.gamma);
    //////objectness_filter->SetObjectDimension(obj_measures.objectness_type);
    //////
    ////////std::cout << obj_measures.alpha << std::endl << obj_measures.beta << std::endl << obj_measures.gamma << std::endl;

    //////multi_scale_Hessian->Update();
    //////
    //////this->ObjectnessImage = multi_scale_Hessian->GetOutput();

    /*typedef itk::ImageFileWriter<ImageType3D> ImageWriterType;
      ImageWriterType::Pointer image_writer = ImageWriterType::New();

      image_writer->SetFileName("C:\\Prathamesh\\Astrocytes\\CoverageExp\\VesselnessImage.mhd");
      image_writer->SetInput(multi_scale_Hessian->GetOutput());
      image_writer->Update();

      ImageWriterType::Pointer image_writer2 = ImageWriterType::New();
      image_writer2->SetFileName("C:\\Prathamesh\\Astrocytes\\CoverageExp\\VesselnessImage_scales.mhd");
      image_writer2->SetInput(multi_scale_Hessian->GetScalesOutput());
      image_writer2->Update();*/
  }


  void MultipleNeuronTracer::RunMask()
  {
    //Binarize the image and use as this as the mask instead of the intensity and contrast thresholds
    MinMaxImageCalculatorType::Pointer minMaxImageCalFilter = MinMaxImageCalculatorType::New();
    minMaxImageCalFilter->SetImage(_PaddedCurvImage);
    minMaxImageCalFilter->Compute();	



    OtsuThresholdImageFilterType::Pointer otsuThresholdImageFilter = OtsuThresholdImageFilterType::New();
    otsuThresholdImageFilter->SetInput(_PaddedCurvImage);
    otsuThresholdImageFilter->Update();
    std::cout << (int)(otsuThresholdImageFilter->GetThreshold()) << std::endl;

    float lowerThreshold = otsuThresholdImageFilter->GetThreshold();
    lowerThreshold = lowerThreshold *0.9;

    BinaryThresholdImageFilterType::Pointer thresholdFilter = BinaryThresholdImageFilterType::New();
    thresholdFilter->SetInput(_PaddedCurvImage);
    thresholdFilter->SetLowerThreshold(lowerThreshold);
    thresholdFilter->SetUpperThreshold(minMaxImageCalFilter->GetMaximum());
    thresholdFilter->SetInsideValue(255);
    thresholdFilter->SetOutsideValue(0);
    thresholdFilter->Update();

    _MaskedImage = thresholdFilter->GetOutput();
  }


  ///////////////////////////////////////////////////////////////////////
  void MultipleNeuronTracer::ReadStartPoints(std::string fname, unsigned int pad) 
  {
    _padz = pad;

    std::string temp, num;
    std::ifstream infile;
    infile.open(fname.c_str());
    if(!infile.good())
    {
      std::cout << "Error reading seed points" << std::endl;
      exit(1);
    }
    size_t x1, x2;
    std::cout << "Reading start points " << std::endl;

    while(!infile.eof()) 
    {
      std::getline(infile,temp);
      if (temp.length() < 1)
        continue;

      std::cout<<temp; // Prints our STRING.
      x1 = temp.find_first_of("0123456789.");
      x2 = temp.find_first_not_of("0123456789.",x1);
      if ((x2 - x1) > 10)
        continue;

      num = temp.substr(x1,x2-x1);
      float x = atof(num.c_str());

      x1 = temp.find_first_of("0123456789.",x2+1);
      x2 = temp.find_first_not_of("0123456789.",x1);
      if ((x2 - x1) > 10)
        continue;

      num = temp.substr(x1,x2-x1);
      float y = atof(num.c_str());

      x1 = temp.find_first_of("0123456789.",x2+1);
      x2 = temp.find_first_not_of("0123456789.",x1);
      if (x2 > temp.length())
        x2 = temp.length();

      if ((x2 - x1) > 10)
        continue;

      num = temp.substr(x1,x2-x1);
      float z = atof(num.c_str());

      itk::Size<3> osz = _size;  //original size padz
      osz[2] = osz[2]-_padz;
      std::cout <<" after conversion " << x <<" "<< y <<" "<< z << std::endl;

      if ( (x>=0.0) && (y>=0.0) && (z>=0.0) )
      {
        itk::Index<3> n;
        n[0] = long(x + 0.5); 
        if (n[0] >= (unsigned int)osz[0]) 
          n[0] = osz[0]-1;
        n[1] = long(y + 0.5);
        if (n[1] >= (unsigned int)osz[1])
          n[1] = osz[1]-1;
        n[2] = long(z + 0.5);
        if (n[2] >= (unsigned int)osz[2])
          n[2] = osz[2]-1;
        _StartPoints.push_back(n);
        // 			std::cout << " is read as " << n << std::endl;
      }
      else
        std::cout << " is discarded (Recommended format XXX YYY ZZZ , Try removing decimal points, add leading zeros in the input text file)" << std::endl;
    }
    infile.close();
  }

  void MultipleNeuronTracer::ReadStartPoints_1(std::vector< itk::Index<3> > somaCentroids, unsigned int pad) 
  {
    _padz = pad;

    std::cout << "Reading start points " << std::endl;
    for(int i=0; i<(int)somaCentroids.size(); ++i)
    {
      float x = (float)somaCentroids.at(i)[0];
      float y = (float)somaCentroids.at(i)[1];
      float z = (float)somaCentroids.at(i)[2];

      // 		std::cout << x <<" "<< y <<" "<< z << std::endl;
      itk::Size<3> osz = _size;  //original size padz
      osz[2] = osz[2]-_padz;

      if ( (x>=0.0) && (y>=0.0) && (z>=0.0) )
      {
        itk::Index<3> n;
        n[0] = long(x + 0.5); 
        if (n[0] >= (unsigned int)osz[0]) 
          n[0] = osz[0]-1;
        n[1] = long(y + 0.5);
        if (n[1] >= (unsigned int)osz[1])
          n[1] = osz[1]-1;
        n[2] = long(z + 0.5);
        if (n[2] >= (unsigned int)osz[2])
          n[2] = osz[2]-1;
        _StartPoints.push_back(n);
        // 			std::cout << " is read as " << n << std::endl;
      }
      else
        std::cout << " is discarded (Recommended format XXX YYY ZZZ , Try removing decimal points, add leading zeros in the input text file)" << std::endl;
    }
  }

  ///////////////////////////////////////////////////////////////////////
  void MultipleNeuronTracer::ReadStartPoints_2(std::string fname, unsigned int pad,float startx,float starty,float startz, float widthx,float widthy,float widthz)
  {

    _padz = pad;

    std::string temp, num;
    std::ifstream infile;
    infile.open(fname.c_str());
    if(!infile.good())
    {
      std::cout << "Error reading seed points" << std::endl;
      exit(1);
    }
    size_t x1, x2;
    //std::cout << "Reading start points " << std::endl;
    std::cout << " Start x" << startx << std::endl;
    std::cout << " Start y" << starty << std::endl;
    std::cout << " Start Z" << startz << std::endl;
    std::cout << " width x" << widthx << std::endl;
    std::cout << " width y" << widthy << std::endl;
    std::cout << " width Z" << widthz << std::endl;
    while(!infile.eof()) 
    {
      std::getline(infile,temp);
      if (temp.length() < 1)
        continue;

      //std::cout<<temp; // Prints our STRING.
      x1 = temp.find_first_of("0123456789.");
      x2 = temp.find_first_not_of("0123456789.",x1);
      if ((x2 - x1) > 10)
        continue;

      num = temp.substr(x1,x2-x1);
      float x = atof(num.c_str());

      x1 = temp.find_first_of("0123456789.",x2+1);
      x2 = temp.find_first_not_of("0123456789.",x1);
      if ((x2 - x1) > 10)
        continue;

      num = temp.substr(x1,x2-x1);
      float y = atof(num.c_str());

      x1 = temp.find_first_of("0123456789.",x2+1);
      x2 = temp.find_first_not_of("0123456789.",x1);
      if (x2 > temp.length())
        x2 = temp.length();

      if ((x2 - x1) > 10)
        continue;

      num = temp.substr(x1,x2-x1);
      float z = atof(num.c_str());

      itk::Size<3> osz = _size;  //original size padz
      osz[2] = osz[2]-_padz;
      //std::cout <<" after conversion " << x <<" "<< y <<" "<< z << std::endl;

      if ( (x>=0.0) && (y>=0.0) && (z>=0.0) ){
        /*	if( (x>=startx && x<=(startx+widthx)) &&
            (y>=starty && x<=(starty+widthy)) &&
            (z>=startz && x<=(startz+widthz))
            )*/
        {
          itk::Index<3> n;
          n[0] = long(x + 0.5); 
          if (n[0] >= (unsigned int)osz[0]) 
            n[0] = osz[0]-1;
          n[1] = long(y + 0.5);
          if (n[1] >= (unsigned int)osz[1])
            n[1] = osz[1]-1;
          n[2] = long(z + 0.5);
          if (n[2] >= (unsigned int)osz[2])
            n[2] = osz[2]-1;
          _StartPoints.push_back(n);
          std::cout << " is read as " << n << std::endl;
        }
      }
      else
        std::cout << " is discarded (Recommended format XXX YYY ZZZ , Try removing decimal points, add leading zeros in the input text file)" << std::endl;
    }
    infile.close();
  }
  ///////////////////////////////////////////////////////////////////////////////////
  void MultipleNeuronTracer::runNDX(void)
  {
    FeatureMain();
  }

  ///////////////////////////////////////////////////////////////////////////////////
  void MultipleNeuronTracer::RunTracing(void)
  {
    if( _flagPipeline == false )
    {
      std::cout<<std::endl<<"FALSE PIPELINE";
      FeatureMain();			//Nice function here that is easy to miss....
    }
	//WriteImage3D(std::string("D:\\Data\\GVF_Traces\\temp\\LOG_seedPoint_Image.mhd"), _NDXImage);
    _CurrentID = 1;

    //set up the connection image and swc image
    _ConnImage = ImageType3D::New();
    _ConnImage->SetRegions(_PaddedCurvImage->GetBufferedRegion());
    _ConnImage->Allocate();
    _ConnImage->FillBuffer(MAXVAL);	//MAXVAL is ... needs to be replaced with std::numeric_limit< float >::max()...

    _SWCImage = SWCImageType3D::New(); //major memory
    _SWCImage->SetRegions(_PaddedCurvImage->GetBufferedRegion());
    _SWCImage->Allocate();
    _SWCImage->FillBuffer(NULL);

    // fill the SWCImage image with start points
    std::vector<IndexType>::iterator startIt;
    int tID = 1;

    clock_t fillSWCImage1_start_time = clock();
    for (startIt = _StartPoints.begin(); startIt != _StartPoints.end(); ++startIt, ++tID)
    {
      itk::Index<3> startIndex = (*startIt);
      startIndex[2] += _padz;													//Convert to padded image index
      SWCNode* start_node = new SWCNode(_CurrentID++, -1, tID, startIndex);	//This is the seed points SWCNode
      _SWCImage->SetPixel(startIndex,start_node);								//Adding all seed points to the SWCImage
      _ConnImage->SetPixel(startIndex,0.0f);									//Set the ConnectedImage to 0.0 at all the seed nodes (remember that the Connected image is all initialized with MAXVAL)... 
      _SWCNodeContainer.push_back(start_node);									//Fill the _SWCNodeContainer with start points
      HeapNode *h = new HeapNode(start_node->ndx, 0.0);						//Heap nodes hold an (index, value) pair
      _PQ.push(h);																//Priority Queue contains the seed nodes now...
    }
    std::cout << "fillSWCImage1 took: " << (clock() - fillSWCImage1_start_time)/(float) CLOCKS_PER_SEC << std::endl;

    clock_t fillSWCImage2_start_time = clock();

    long eCounter = 0, TotalePoints;
    itk::ImageRegionConstIterator<ImageType3D> Nit(_NDXImage, _NDXImage->GetBufferedRegion());
    for (Nit.GoToBegin(); !Nit.IsAtEnd(); ++Nit) 
    {
      if (Nit.Get() > 0)	//Vesselness value is greater than 0
      {
        itk::Index<3> endx = Nit.GetIndex();
        SWCNode* s2 = new SWCNode(0, -1, -1*(++eCounter), endx);	//id = 0, parent_id = -1, tree id = -1 * eCounter, index that this vesselness value is greater than 0
        _SWCImage->SetPixel(endx,s2);								//Adding all critical points where vesselness value is greater than 0 to the SWC image
      }
    }
    std::cout << "fillSWCImage2 took: " << (clock() - fillSWCImage2_start_time)/(float) CLOCKS_PER_SEC << std::endl;

    bool print_out_critical_point_image = false;
    if (print_out_critical_point_image)
    {
      //Make a unsigned char image to print out the critical points image
      typedef itk::Image< unsigned char, 3 > CriticalPointsImageType;
      CriticalPointsImageType::Pointer critical_point_image = CriticalPointsImageType::New();
      critical_point_image->SetRegions(_SWCImage->GetLargestPossibleRegion());
      critical_point_image->Allocate();
      critical_point_image->FillBuffer(0);

      //Iterate through SWCImage and setting critical points to 255 in critical_point_image
      itk::ImageRegionConstIterator< SWCImageType3D > SWCImage_iter(_SWCImage, _SWCImage->GetLargestPossibleRegion());
      SWCImage_iter.GoToBegin();

      while (!SWCImage_iter.IsAtEnd())
      {
        SWCNode* critical_point_node = SWCImage_iter.Get();	
        critical_point_image->SetPixel(critical_point_node->ndx, 255);
        ++SWCImage_iter;
      }

      typedef itk::ImageFileWriter< CriticalPointsImageType > CriticalPointsWriterType;
      CriticalPointsWriterType::Pointer crit_pts_writer = CriticalPointsWriterType::New();
      crit_pts_writer->SetInput(critical_point_image);
      crit_pts_writer->SetFileName("critical_point_image.mhd");
      crit_pts_writer->Update();
    }

    TotalePoints = eCounter;
    std::cout<<"eCounter = "<<eCounter<<std::endl;	//eCounter is just number of nodes that are critical points (but not seed points)
    //std::cout << "No of CTs inserted : " <<  TotalePoints << std::endl;

    //Generating some kind of offset neighborhood... this needs to be done with itkNeighborhoodIterator
    itk::Offset<3> x1 = {{-1, 0 ,0}};
    _off.push_back( x1 );
    x1[0] = 1;					// x1 = {{1, 0, 0}}
    _off.push_back( x1 );
    x1[0] = 0; 
    x1[1] = -1;					// x1 = {{0, -1, 0}}
    _off.push_back( x1 );
    x1[1] = 1;					// x1 = {{0, 1, 0}}
    _off.push_back( x1 );
    x1[1] = 0; 
    x1[2] = -1;					// x1 = {{0, 0, -1}}
    _off.push_back( x1 );
    x1[2] = 1;					// x1 = {{0, 0, 1}}
    _off.push_back( x1 );

    std::vector<OffsetType>::iterator oit;
    bool showMessage = false;
    //std::cout << " Heap size: " << PQ.size() << std::endl;
    float KeyValue;

    clock_t PQ_popping_start_time = clock();

    while(!_PQ.empty())	//For each seed node
    {
      //Take the top HeapNode and remove it from the Priority Queue 
      HeapNode *h = _PQ.top();
      _PQ.pop();

      //Temporarily store the index and value of the node
      itk::Index<3> ndx = h->ndx;
      KeyValue = h->KeyValue;
      delete h;

      //Don't do anything if the heapnode value is larger than the one in the connected image
      if ( KeyValue > _ConnImage->GetPixel(ndx) ) 
        continue;


      if ((eCounter <= 0) || (KeyValue > _CostThreshold) ) 
      {
        if (showMessage == true) 
        {
          std::cout << "NOTE: Exiting the search at cost " << _CostThreshold << " However, " << (100*eCounter)/TotalePoints << "%% of the image is still not covered, change cost if necessary!!\r"<< std::endl;
          //std::cout << "Cleaning Heap size: " << PQ.size() << std::endl;
          //std::cout<<"keyvalue = "<<KeyValue<<std::endl;
          showMessage = false;
        }

        SWCNode* t  = _SWCImage->GetPixel(ndx);
        if ( t != NULL) 
        {
          if (t->TreeID < 0) 
          {
            delete t;
          }
        }
        continue;
      }

      SWCNode* s = _SWCImage->GetPixel(ndx);
      if (s != NULL) 
      {
        if (s->TreeID < 0) 
        {
          std::vector<IndexType> Chain;

          SWCNode* L = TBack(ndx, Chain);

          if ( L  != NULL ) 
          {
            float costFactor = GetCostLocal( L , ndx);

            std::vector<IndexType>::reverse_iterator cit;
            SWCNode* par = L;

            for (cit = Chain.rbegin(); cit != Chain.rend(); ++cit) 
            {
              SWCNode* t = _SWCImage->GetPixel(*cit);
              if (t == NULL) 
              {
                float val = _ConnImage->GetPixel(*cit) * costFactor;
                _ConnImage->SetPixel((*cit),val);
                SWCNode* s = new SWCNode(_CurrentID++, par, L->TreeID, (*cit));
                _SWCImage->SetPixel((*cit),s);
                _SWCNodeContainer.push_back(s);
                par->children.push_back(s);
                par = s;
                HeapNode *h = new HeapNode((*cit), val);
                _PQ.push(h);
              }
              else 
              {
                if (t->TreeID < 0) 
                {
                  delete t;
                  eCounter--;
                  float val = _ConnImage->GetPixel(*cit) * costFactor;
                  _ConnImage->SetPixel((*cit),val);
                  SWCNode* s = new SWCNode(_CurrentID++, par, L->TreeID, (*cit));
                  _SWCImage->SetPixel((*cit),s);
                  _SWCNodeContainer.push_back(s);
                  par->children.push_back(s);
                  //std::cout<<"SWCImage Node @ " << (*cit) << "(" << s->ID << ") with parent " << par->ID << "  Cost: " << val << "  " << (100*eCounter)/TotalePoints << "% Remaining.\r";// << std::endl;
                  par = s;
                  HeapNode *h = new HeapNode((*cit), val);
                  _PQ.push(h);
                }
              }
            }
          } 
        }
      }

      for (oit = _off.begin(); oit < _off.end(); ++oit) 
      {
        itk::Index<3> ndx2 = ndx + (*oit);
        if ( (ndx2[0] < 2) || (ndx2[1] < 2) || (ndx2[2] < 2) || (ndx2[0] >= unsigned(_size[0] - 2)) || (ndx2[1] >= unsigned(_size[1] - 2)) || (ndx2[2] >= unsigned(_size[2] - 2)) )  
          continue;

        if (_SWCImage->GetPixel(ndx2) != NULL) 
        {
          if (_SWCImage->GetPixel(ndx2)->TreeID > 0) 
          {
            continue;			
          }
        }
        PixelType P = 1/(_PaddedCurvImage->GetPixel(ndx2) + 0.001f);  // consider taking inverse here
        PixelType a1, a2, a3;
        ScanNeighbors(a1,a2,a3, ndx2);
        PixelType aa = Update( a1, a2, a3, P );
        if ( _ConnImage->GetPixel(ndx2) > aa )  
        {
          _ConnImage->SetPixel(ndx2, aa);
          HeapNode *h = new HeapNode(ndx2, aa);
          _PQ.push(h);
        }
      }
    }


    std::cout << "PQ popping took: " << (clock() - PQ_popping_start_time)/(float) CLOCKS_PER_SEC << std::endl;

    clock_t Interpolate1_start_time = clock();
    Interpolate(2.0);
    std::cout << "Interpolate1 took: " << (clock() - Interpolate1_start_time)/(float) CLOCKS_PER_SEC << std::endl;

    clock_t Decimate_start_time = clock();
    Decimate();
    std::cout << "Decimate took: " << (clock() - Decimate_start_time)/(float) CLOCKS_PER_SEC << std::endl;

    clock_t Interpolate2_start_time = clock();
    Interpolate(2.0);	
    std::cout << "Interpolate2 took: " << (clock() - Interpolate2_start_time)/(float) CLOCKS_PER_SEC << std::endl;

    // 	clock_t RemoveIntraSomaNodes_start_time = clock();
    // 	RemoveIntraSomaNodes();
    // 	std::cout << "RemoveIntraSomaNodes took: " << (clock() - RemoveIntraSomaNodes_start_time)/(float) CLOCKS_PER_SEC << std::endl;

  }



  /**
    Run Tracing USING GVF seed detection technique
   **/

  void MultipleNeuronTracer::RunGVFTracing(bool preComputedGVF)
  {
    if( _flagPipeline == false )
    {
      std::cout<<std::endl<<"FALSE PIPELINE";
      UpdateNDXImage_GVF(preComputedGVF);			//Nice function here that is easy to miss....
    }
    bool yuWangTest = true;
    //////////if(yuWangTest == true){


	//WriteImage3D(std::string("D:\\Data\\GVF_Traces\\temp\\GVF_seedPoint_Image.mhd"), _NDXImage);

    //////////	//ReaderType::GlobalWarningDisplayOff();
    //////////	//ReaderType::Pointer reader = ReaderType::New();
    //////////	////reader->SetFileName("D:\\Data\\YuWang_Test\\GVF_Seed_Test\\vesselimage.nrrd");
    //////////	////_NDXImage = reader->GetOutput();
    //////////	////_NDXImage->Update();



    //////////	_NDXImage = ImageType3D::New();
    //////////	_NDXImage->SetRegions(_PaddedCurvImage->GetBufferedRegion());
    //////////	_NDXImage->Allocate();
    //////////	_NDXImage->FillBuffer(0.0f);





    //////////	std::string fileName = "D:\\Data\\YuWang_Test\\yan_xu\\SeedPoint.txt";
    //////////	FILE * fp = fopen(fileName.c_str(), "r");
    //////////	char buff[1024];
    //////////	if(fp==NULL)
    //////////	{
    //////////		printf("Couldn't open file %s for parsing\n");
    //////////		//	return false;
    //////////	}
    //////////	double i,j,k,r;
    //////////	while(!feof(fp))
    //////////	{
    //////////		if(fgets(buff,1024,fp)==NULL)
    //////////		{
    //////////			break;
    //////////		}
    //////////		int pc = 0;
    //////////		while(buff[pc]==' '&&(pc<1023))
    //////////		{
    //////////			pc++;
    //////////		}
    //////////		if(buff[pc]=='#') // ignoring comment lines for now. We have to read scale from it! TODO
    //////////		{
    //////////			continue;
    //////////		}
    //////////		sscanf(buff,"%lf %lf %lf",&i,&j,&k);
    //////////		double radius = r;
    //////////		ImageType3D::IndexType pixelIndex;
    //////////		pixelIndex[0] = i;
    //////////		pixelIndex[1] = j;
    //////////		pixelIndex[2] = k;
    //////////		_NDXImage->SetPixel(pixelIndex, 255);
    //////////		WriteImage3D(std::string("seedPoint_Image.mhd"), _NDXImage);
    //////////	}
    //////////}



    //WriteImage3D(std::string("seedPoint_Image.mhd"), _NDXImage);








    _CurrentID = 1;

    //set up the connection image and swc image
    _ConnImage = ImageType3D::New();
    _ConnImage->SetRegions(_PaddedCurvImage->GetBufferedRegion());
    _ConnImage->Allocate();
    _ConnImage->FillBuffer(MAXVAL);	//MAXVAL is ... needs to be replaced with std::numeric_limit< float >::max()...

    _SWCImage = SWCImageType3D::New(); //major memory
    _SWCImage->SetRegions(_PaddedCurvImage->GetBufferedRegion());
    _SWCImage->Allocate();
    _SWCImage->FillBuffer(NULL);

    // fill the SWCImage image with start points
    std::vector<IndexType>::iterator startIt;
    int tID = 1;

    clock_t fillSWCImage1_start_time = clock();
    for (startIt = _StartPoints.begin(); startIt != _StartPoints.end(); ++startIt, ++tID)
    {
      itk::Index<3> startIndex = (*startIt);
      startIndex[2] += _padz;													//Convert to padded image index
      SWCNode* start_node = new SWCNode(_CurrentID++, -1, tID, startIndex);	//This is the seed points SWCNode
      _SWCImage->SetPixel(startIndex,start_node);								//Adding all seed points to the SWCImage
      _ConnImage->SetPixel(startIndex,0.0f);									//Set the ConnectedImage to 0.0 at all the seed nodes (remember that the Connected image is all initialized with MAXVAL)... 
      _SWCNodeContainer.push_back(start_node);									//Fill the _SWCNodeContainer with start points
      HeapNode *h = new HeapNode(start_node->ndx, 0.0);						//Heap nodes hold an (index, value) pair
      _PQ.push(h);																//Priority Queue contains the seed nodes now...
    }
    std::cout << "fillSWCImage1 took: " << (clock() - fillSWCImage1_start_time)/(float) CLOCKS_PER_SEC << std::endl;

    clock_t fillSWCImage2_start_time = clock();

    long eCounter = 0, TotalePoints;

    itk::ImageRegionConstIterator<ImageType3D> Nit(_NDXImage, _NDXImage->GetBufferedRegion());
    for (Nit.GoToBegin(); !Nit.IsAtEnd(); ++Nit) 
    {
      if (Nit.Get() > 0)	//Vesselness value is greater than 0
      {
        itk::Index<3> endx = Nit.GetIndex();
        SWCNode* s2 = new SWCNode(0, -1, -1*(++eCounter), endx);	//id = 0, parent_id = -1, tree id = -1 * eCounter, index that this vesselness value is greater than 0
        _SWCImage->SetPixel(endx,s2);								//Adding all critical points where vesselness value is greater than 0 to the SWC image
      }
    }
    std::cout << "fillSWCImage2 took: " << (clock() - fillSWCImage2_start_time)/(float) CLOCKS_PER_SEC << std::endl;

    bool print_out_critical_point_image = false;
    if (print_out_critical_point_image)
    {
      //Make a unsigned char image to print out the critical points image
      typedef itk::Image< unsigned char, 3 > CriticalPointsImageType;
      CriticalPointsImageType::Pointer critical_point_image = CriticalPointsImageType::New();
      critical_point_image->SetRegions(_SWCImage->GetLargestPossibleRegion());
      critical_point_image->Allocate();
      critical_point_image->FillBuffer(0);

      //Iterate through SWCImage and setting critical points to 255 in critical_point_image
      itk::ImageRegionConstIterator< SWCImageType3D > SWCImage_iter(_SWCImage, _SWCImage->GetLargestPossibleRegion());
      SWCImage_iter.GoToBegin();

      while (!SWCImage_iter.IsAtEnd())
      {
        SWCNode* critical_point_node = SWCImage_iter.Get();	
        critical_point_image->SetPixel(critical_point_node->ndx, 255);
        ++SWCImage_iter;
      }

      typedef itk::ImageFileWriter< CriticalPointsImageType > CriticalPointsWriterType;
      CriticalPointsWriterType::Pointer crit_pts_writer = CriticalPointsWriterType::New();
      crit_pts_writer->SetInput(critical_point_image);
      crit_pts_writer->SetFileName("critical_point_image.mhd");
      crit_pts_writer->Update();
    }

    TotalePoints = eCounter;
    std::cout<<"eCounter = "<<eCounter<<std::endl;	//eCounter is just number of nodes that are critical points (but not seed points)
    //std::cout << "No of CTs inserted : " <<  TotalePoints << std::endl;

    //Generating some kind of offset neighborhood... this needs to be done with itkNeighborhoodIterator
    itk::Offset<3> x1 = {{-1, 0 ,0}};
    _off.push_back( x1 );
    x1[0] = 1;					// x1 = {{1, 0, 0}}
    _off.push_back( x1 );
    x1[0] = 0; 
    x1[1] = -1;					// x1 = {{0, -1, 0}}
    _off.push_back( x1 );
    x1[1] = 1;					// x1 = {{0, 1, 0}}
    _off.push_back( x1 );
    x1[1] = 0; 
    x1[2] = -1;					// x1 = {{0, 0, -1}}
    _off.push_back( x1 );
    x1[2] = 1;					// x1 = {{0, 0, 1}}
    _off.push_back( x1 );

    std::vector<OffsetType>::iterator oit;
    bool showMessage = false;
    //std::cout << " Heap size: " << PQ.size() << std::endl;
    float KeyValue;

    clock_t PQ_popping_start_time = clock();

    while(!_PQ.empty())	//For each seed node
    {
      //Take the top HeapNode and remove it from the Priority Queue 
      HeapNode *h = _PQ.top();
      _PQ.pop();

      //Temporarily store the index and value of the node
      itk::Index<3> ndx = h->ndx;
      KeyValue = h->KeyValue;
      delete h;

      //Don't do anything if the heapnode value is larger than the one in the connected image
      if ( KeyValue > _ConnImage->GetPixel(ndx) ) 
        continue;


      if ((eCounter <= 0) || (KeyValue > _CostThreshold) ) 
      {
        if (showMessage == true) 
        {
          std::cout << "NOTE: Exiting the search at cost " << _CostThreshold << " However, " << (100*eCounter)/TotalePoints << "%% of the image is still not covered, change cost if necessary!!\r"<< std::endl;
          //std::cout << "Cleaning Heap size: " << PQ.size() << std::endl;
          //std::cout<<"keyvalue = "<<KeyValue<<std::endl;
          showMessage = false;
        }

        SWCNode* t  = _SWCImage->GetPixel(ndx);
        if ( t != NULL) 
        {
          if (t->TreeID < 0) 
          {
            delete t;
          }
        }
        continue;
      }

      SWCNode* s = _SWCImage->GetPixel(ndx);
      if (s != NULL) 
      {
        if (s->TreeID < 0) 
        {
          std::vector<IndexType> Chain;

          SWCNode* L = TBack(ndx, Chain);

          if ( L  != NULL ) 
          {
            float costFactor = GetCostLocal( L , ndx);

            std::vector<IndexType>::reverse_iterator cit;
            SWCNode* par = L;

            for (cit = Chain.rbegin(); cit != Chain.rend(); ++cit) 
            {
              SWCNode* t = _SWCImage->GetPixel(*cit);
              if (t == NULL) 
              {
                float val = _ConnImage->GetPixel(*cit) * costFactor;
                _ConnImage->SetPixel((*cit),val);
                SWCNode* s = new SWCNode(_CurrentID++, par, L->TreeID, (*cit));
                _SWCImage->SetPixel((*cit),s);
                _SWCNodeContainer.push_back(s);
                par->children.push_back(s);
                par = s;
                HeapNode *h = new HeapNode((*cit), val);
                _PQ.push(h);
              }
              else 
              {
                if (t->TreeID < 0) 
                {
                  delete t;
                  eCounter--;
                  float val = _ConnImage->GetPixel(*cit) * costFactor;
                  _ConnImage->SetPixel((*cit),val);
                  SWCNode* s = new SWCNode(_CurrentID++, par, L->TreeID, (*cit));
                  _SWCImage->SetPixel((*cit),s);
                  _SWCNodeContainer.push_back(s);
                  par->children.push_back(s);
                  //std::cout<<"SWCImage Node @ " << (*cit) << "(" << s->ID << ") with parent " << par->ID << "  Cost: " << val << "  " << (100*eCounter)/TotalePoints << "% Remaining.\r";// << std::endl;
                  par = s;
                  HeapNode *h = new HeapNode((*cit), val);
                  _PQ.push(h);
                }
              }
            }
          } 
        }
      }

      for (oit = _off.begin(); oit < _off.end(); ++oit) 
      {
        itk::Index<3> ndx2 = ndx + (*oit);
        if ( (ndx2[0] < 2) || (ndx2[1] < 2) || (ndx2[2] < 2) || (ndx2[0] >= unsigned(_size[0] - 2)) || (ndx2[1] >= unsigned(_size[1] - 2)) || (ndx2[2] >= unsigned(_size[2] - 2)) )  
          continue;

        if (_SWCImage->GetPixel(ndx2) != NULL) 
        {
          if (_SWCImage->GetPixel(ndx2)->TreeID > 0) 
          {
            continue;			
          }
        }
        PixelType P = 1/(_PaddedCurvImage->GetPixel(ndx2) + 0.001f);  // consider taking inverse here
        PixelType a1, a2, a3;
        ScanNeighbors(a1,a2,a3, ndx2);
        PixelType aa = Update( a1, a2, a3, P );
        if ( _ConnImage->GetPixel(ndx2) > aa )  
        {
          _ConnImage->SetPixel(ndx2, aa);
          HeapNode *h = new HeapNode(ndx2, aa);
          _PQ.push(h);
        }
      }
    }


    std::cout << "PQ popping took: " << (clock() - PQ_popping_start_time)/(float) CLOCKS_PER_SEC << std::endl;

    clock_t Interpolate1_start_time = clock();
    Interpolate(2.0);
    std::cout << "Interpolate1 took: " << (clock() - Interpolate1_start_time)/(float) CLOCKS_PER_SEC << std::endl;

    clock_t Decimate_start_time = clock();
    Decimate();
    std::cout << "Decimate took: " << (clock() - Decimate_start_time)/(float) CLOCKS_PER_SEC << std::endl;

    clock_t Interpolate2_start_time = clock();
    Interpolate(2.0);	
    std::cout << "Interpolate2 took: " << (clock() - Interpolate2_start_time)/(float) CLOCKS_PER_SEC << std::endl;

    // 	clock_t RemoveIntraSomaNodes_start_time = clock();
    //RemoveIntraSomaNodes();
    // 	std::cout << "RemoveIntraSomaNodes took: " << (clock() - RemoveIntraSomaNodes_start_time)/(float) CLOCKS_PER_SEC << std::endl;

  }



  void MultipleNeuronTracer::UpdateNDXImage_GVF(bool preComputedGVFAndVessel)
  {

	// int num_iteration = 15;
	// float noise_level = 100;
	
	int num_iteration = this->noOfIteration;
	float noise_level = this->mu;
    int smoothing_scale = 1;
    int detection_method = 1;
    int radius = 0.1;
    int iter_num = 20;
    double sigma_min = 0.2;
    double sigma_max = 8;
    int sigma_step = 10;
    bool upSample = false;
    int scaleFactor = 2;
    std::cout<<std::endl<<"upscale the image";
    //WriteImage3D(std::string("before_Upsampling.tif"), _PaddedCurvImage);
    // upsample the InputImage and the somaImage
    if(upSample){
      _PaddedCurvImage = Upsampling(_PaddedCurvImage,scaleFactor);
      //_SomaImage = UpsamplingLabel(_SomaImage,scaleFactor);
      //std::vector<IndexType> tempStartPoints;
      //std::vector<IndexType>::iterator startIt;
      //for (startIt = _StartPoints.begin(); startIt != _StartPoints.end(); ++startIt)
      //{
      //	itk::Index<3> startIndex = (*startIt);
      //	std::cout<<startIndex[0]<<" "<<startIndex[1]<<" "<<startIndex[2]<<std::endl;
      //	startIndex[0] = startIndex[0]*scaleFactor;													
      //	startIndex[1] = startIndex[1]*scaleFactor;
      //	startIndex[2] = startIndex[2]*scaleFactor;
      //	std::cout<<startIndex[0]<<" "<<startIndex[1]<<" "<<startIndex[2]<<std::endl;
      //	//_StartPoints.erase(startIt);
      //	tempStartPoints.push_back(startIndex);
      //}
      //_StartPoints.swap(tempStartPoints);
      //std::cout<<"***************************************************************"<<std::endl;
      ////print the start points
      //std::vector<IndexType>::iterator startIt1;
      //for (startIt1 = _StartPoints.begin(); startIt1!= _StartPoints.end(); ++startIt1)
      //{
      //	itk::Index<3> startIndex = (*startIt1);
      //	std::cout<<startIndex[0]<<" "<<startIndex[1]<<" "<<startIndex[2]<<std::endl;
      //}

      WriteImage3D(std::string("after_Upsampling.tif"), _PaddedCurvImage);
    }







    //std::cout<<"compute Multi Vesselness"<<std::endl;
    //ComputeMultiVesselness(sigma_min, sigma_max, sigma_step);
    //WriteImage3D(std::string("multi_Vesselness_enhancement.tif"), _PaddedCurvImage);

    if(!preComputedGVFAndVessel){
      std::cout<<"compute GVF"<<std::endl;
      this->computeGVF(noise_level,num_iteration,smoothing_scale);

      std::cout<<"compute GVF Vesselness"<<std::endl;
      this->ComputeGVFVesselness();
    }
    //WriteImage3D(std::string("GVF_Vesselness_enhancement.tif"), _IVessel);
    std::cout<<"compute seed Detection"<<std::endl;
    _v_threshold = 0.1;// this value is calculated in the computeGVFVesselness() function by 
    this->SeedDetection(_v_threshold,detection_method,radius);
    std::cout<<"compute seed Adjustment"<<std::endl;
    this->SeedAdjustment(iter_num);



    _NDXImage = ImageType3D::New();
    _NDXImage->SetRegions(_PaddedCurvImage->GetBufferedRegion());
    _NDXImage->Allocate();
    _NDXImage->FillBuffer(0.0f);

    // update the NDXImage with the seedpoints
    int i = 0;
    while( i < SeedPt.NP )
    {
      ImageType3D::IndexType index; 

      index[0] = ceil(SeedPt.Pt[i].x);
      index[1] = ceil(SeedPt.Pt[i].y);
      index[2] = ceil(SeedPt.Pt[i].z);
      _NDXImage->SetPixel(index,255);
      i++;
    }







  }


  void MultipleNeuronTracer::computeGVF_2(ImageType3D::Pointer &image,int noise_level, int num_iteration, int smoothing_scale)
  {
		_PaddedCurvImage = image;
		computeGVF(noise_level, num_iteration, smoothing_scale);
  }
  void MultipleNeuronTracer::ComputeGVFVesselness_2(ImageType3D::Pointer &image)
  {
		_PaddedCurvImage = image;
		ComputeGVFVesselness();
  }







  void MultipleNeuronTracer::ComputeMultiVesselness(double sigma_min, double sigma_max, int sigma_step)
  {

    typedef itk::RescaleIntensityImageFilter< ImageType3D, ImageType3D> RescaleFilterType;
    bool SatoVesselness = false;
    /** NOT USED ???? **/
    if( SatoVesselness )
    {
      //// Declare the type of enhancement filter - use ITK's 3D vesselness (Sato)
      //typedef itk::Hessian3DToVesselnessMeasureImageFilter<float> VesselnessFilterType;

      //// Declare the type of multiscale enhancement filter
      //typedef itk::MultiScaleHessianBasedMeasureImageFilter<ImageType,VesselnessFilterType> MultiScaleEnhancementFilterType;

      //// Instantiate the multiscale filter and set the input image
      //MultiScaleEnhancementFilterType::Pointer multiScaleEnhancementFilter = MultiScaleEnhancementFilterType::New();
      //multiScaleEnhancementFilter->SetInput( _PaddedCurvImage );
      //multiScaleEnhancementFilter->SetSigmaMin( sigma_min );
      //multiScaleEnhancementFilter->SetSigmaMax( sigma_max );
      //multiScaleEnhancementFilter->SetNumberOfSigmaSteps( sigma_step );

      //// Get the vesselness filter and set the parameters
      //VesselnessFilterType* vesselnessFilter = 
      //	multiScaleEnhancementFilter->GetHessianToMeasureFilter();
      //vesselnessFilter->SetAlpha1(0.5);
      //vesselnessFilter->SetAlpha2(2);

      //try
      //{
      //	multiScaleEnhancementFilter->Update();
      //}
      //		                     
      //{
      //	std::cerr << "Exception caught: "<< err << std::endl;
      //}
      //RescaleFilterType::Pointer rescale1 = RescaleFilterType::New();
      //rescale1->SetInput( multiScaleEnhancementFilter->GetOutput() );
      //rescale1->SetOutputMinimum( 0 );
      //rescale1->SetOutputMaximum( 255 );
      //rescale1->Update();
      //_PaddedCurvImage = rescale1->GetOutput();
    }
    else
    {
      typedef itk::MultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter<ImageType3D,ImageType3D> MultiScaleVesselnessFilterType;
      MultiScaleVesselnessFilterType::Pointer MultiScaleVesselnessFilter =
        MultiScaleVesselnessFilterType::New();
      MultiScaleVesselnessFilter->SetInput( _PaddedCurvImage );
      MultiScaleVesselnessFilter->SetSigmaMin( sigma_min );
      MultiScaleVesselnessFilter->SetSigmaMax( sigma_max );
      MultiScaleVesselnessFilter->SetNumberOfSigmaSteps( sigma_step );

      try
      {
        MultiScaleVesselnessFilter->Update();
      }
      catch( itk::ExceptionObject & err )
      {
        std::cerr << "Exception caught: "<< err << std::endl;
      }


      RescaleFilterType::Pointer rescale = RescaleFilterType::New();
      rescale->SetInput( MultiScaleVesselnessFilter->GetOutput() );
      rescale->SetOutputMinimum( 0 );
      rescale->SetOutputMaximum( 1 );
      rescale->Update();
      _PaddedCurvImage = rescale->GetOutput();

    }
  }



  void MultipleNeuronTracer::computeGVF(int noise_level, int num_iteration, int smoothing_scale)
  { 
    
	  std::cout<<"noise_level	"<<noise_level<<std::endl;
	  std::cout<<"num_iteration	"<<num_iteration<<std::endl;
	if( smoothing_scale == 0 )
    {
      typedef itk::GradientImageFilter<ImageType3D, float, float> GradientImageFilterType;
      GradientImageFilterType::Pointer gradientFilter = GradientImageFilterType::New();

      //typedef itk::RescaleIntensityImageFilter< ImageType3D, ImageType3D> RescaleFilterType;
      //RescaleFilterType::Pointer rescale = RescaleFilterType::New();
      //rescale->SetInput( _PaddedCurvImage );
      //rescale->SetOutputMinimum( 0 );
      //rescale->SetOutputMaximum( 1 );
      //rescale->Update();

      //gradientFilter->SetInput(rescale->GetOutput());
      // Not required image is already 0 to 1
      gradientFilter->SetInput( _PaddedCurvImage );
      //gradientFilter->SetInput(I);

      try
      {
        gradientFilter->Update();
      }
      catch( itk::ExceptionObject & err )
      {
        std::cerr << "Exception caught: " << err << std::endl;
      }

      //IG = gradientFilter->GetOutput();
      if( num_iteration == 0 )
      {
        _IGVF = gradientFilter->GetOutput();
      }
      else
      {
        typedef itk::GradientVectorFlowImageFilter<GradientImageType, GradientImageType> GradientVectorFlowFilterType;
        GradientVectorFlowFilterType::Pointer GVFFilter = GradientVectorFlowFilterType::New();

        GVFFilter->SetInput(gradientFilter->GetOutput());
        GVFFilter->SetNoiseLevel(noise_level);
        GVFFilter->SetIterationNum(num_iteration);

        try
        {
          GVFFilter->Update();
        }
        catch( itk::ExceptionObject & err )
        {
          std::cerr << "Exception caught: " << err << std::endl;
        }
        _IGVF = GVFFilter->GetOutput();
      }

    }
    else
    {
      typedef itk::GradientRecursiveGaussianImageFilter<ImageType3D, GradientImageType> GradientImageFilterType;
      GradientImageFilterType::Pointer gradientFilter = GradientImageFilterType::New();

      //typedef itk::RescaleIntensityImageFilter< ImageType3D, ImageType3D> RescaleFilterType;
      //RescaleFilterType::Pointer rescale = RescaleFilterType::New();
      //rescale->SetInput( _PaddedCurvImage );
      //rescale->SetOutputMinimum( 0 );
      //rescale->SetOutputMaximum( 1 );
      //rescale->Update();

      gradientFilter->SetSigma(smoothing_scale);
      //gradientFilter->SetInput(rescale->GetOutput());
      gradientFilter->SetInput( _PaddedCurvImage );

      try
      {
        gradientFilter->Update();
      }
      catch( itk::ExceptionObject & err )
      {
        std::cerr << "Exception caught: " << err << std::endl;
      }

      //IG = gradientFilter->GetOutput();

      typedef itk::GradientVectorFlowImageFilter<GradientImageType, GradientImageType> GradientVectorFlowFilterType;
      GradientVectorFlowFilterType::Pointer GVFFilter = GradientVectorFlowFilterType::New();

      GVFFilter->SetInput(gradientFilter->GetOutput());
      GVFFilter->SetNoiseLevel(noise_level);
      GVFFilter->SetIterationNum(num_iteration);

      try
      {
        GVFFilter->Update();
      }
      catch( itk::ExceptionObject & err )
      {
        std::cerr << "Exception caught: " << err << std::endl;
      }
      _IGVF = GVFFilter->GetOutput();

      /////// Write GVF Image////////
      //typedef itk::ImageFileWriter<GradientImageType> ImageWriterType;
      //ImageWriterType::Pointer image_writer = ImageWriterType::New();
      //image_writer->SetFileName("D:\\Data\\FSData\\GVF_Test\\gvf_image.mhd");
      //image_writer->SetInput(_IGVF);
      //try
      //{
      //	image_writer->Update();
      //}
      //catch( itk::ExceptionObject & err )
      //{
      //	std::cerr << "Exception caught: " << err << std::endl;
      //}
    }


  }



  void MultipleNeuronTracer::ComputeGVFVesselness()
  {


    double FrangiAlpha = 0.5;
    double FrangiBeta = 0.5;
    double FrangiC = 10;
    double A = 2 * pow(FrangiAlpha,2);
    double B = 2 * pow(FrangiBeta,2);
    double C = 2 * pow(FrangiC,2);

	// This should be a duplicator
    typedef itk::CastImageFilter<ImageType3D,ImageType3D> CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(_PaddedCurvImage);
    caster->Update();
    _IVessel = caster->GetOutput();
    //_IVessel = _PaddedCurvImage;

    /*CasterType::Pointer caster1 = CasterType::New();
      caster1->SetInput(I);
      caster1->Update();
      Vx1 = caster1->GetOutput();


      CasterType::Pointer caster2 = CasterType::New();
      caster2->SetInput(I);
      caster2->Update();
      Vy1 = caster2->GetOutput();


      CasterType::Pointer caster3 = CasterType::New();
      caster3->SetInput(I);
      caster3->Update();
      Vz1 = caster3->GetOutput();*/


    typedef itk::RecursiveGaussianImageFilter<
      ImageType3D, ImageType3D > FilterType;
    //typedef itk::GradientImageFilter<ImageType3D, float, float> FilterType;
    typedef ImageType3D::Pointer ProbImagePointer;
    ProbImagePointer Dx = extract_one_component(0,_IGVF);
    ProbImagePointer Dy = extract_one_component(1,_IGVF);
    ProbImagePointer Dz = extract_one_component(2,_IGVF); 

    /* WriteImage3D("DX.tif",Dx);
       WriteImage3D("DX.tif",Dy);
       WriteImage3D("DX.tif",Dz);
       */ 

    /*FilterType::Pointer filter_x = FilterType::New();
      filter_x->SetInput(Dx);
      filter_x->Update();
      Dxx = extract_one_component(0,filter_x->GetOutput());
      Dxy = extract_one_component(1,filter_x->GetOutput());
      Dxz = extract_one_component(2,filter_x->GetOutput());

      FilterType::Pointer filter_y = FilterType::New();
      filter_y->SetInput(Dy);
      filter_y->Update();
      Dyy = extract_one_component(1,filter_y->GetOutput());
      Dyz = extract_one_component(2,filter_y->GetOutput());

      FilterType::Pointer filter_z = FilterType::New();
      filter_z->SetInput(Dz);
      filter_z->Update();
      Dzz = extract_one_component(2,filter_z->GetOutput());*/

    //Dxx
    FilterType::Pointer filter1 = FilterType::New();
    filter1->SetDirection(0);
    filter1->SetOrder(FilterType::FirstOrder);
    filter1->SetInput(Dx);
    filter1->SetSigma(1);
    filter1->Update();
    ProbImagePointer Dxx = filter1->GetOutput();

    //Dxy
    FilterType::Pointer filter2 = FilterType::New();
    filter2->SetDirection(1);
    filter2->SetOrder(FilterType::FirstOrder);
    filter2->SetInput(Dx);
    filter2->SetSigma(1);
    filter2->Update();
    ProbImagePointer Dxy = filter2->GetOutput();

    //Dxz
    FilterType::Pointer filter3 = FilterType::New();
    filter3->SetDirection(2);
    filter3->SetOrder(FilterType::FirstOrder);
    filter3->SetInput(Dx);
    filter3->SetSigma(1);
    filter3->Update();
    ProbImagePointer Dxz = filter3->GetOutput();

    //Dyy
    FilterType::Pointer filter4 = FilterType::New();
    filter4->SetDirection(1);
    filter4->SetOrder(FilterType::FirstOrder);
    filter4->SetInput(Dy);
    filter4->SetSigma(1);
    filter4->Update();
    ProbImagePointer Dyy = filter4->GetOutput();

    //Dyz
    FilterType::Pointer filter5 = FilterType::New();
    filter5->SetDirection(2);
    filter5->SetOrder(FilterType::FirstOrder);
    filter5->SetInput(Dy);
    filter5->SetSigma(1);
    filter5->Update();
    ProbImagePointer Dyz = filter5->GetOutput();

    //Dzz
    FilterType::Pointer filter6 = FilterType::New();
    filter6->SetDirection(2);
    filter6->SetOrder(FilterType::FirstOrder);
    filter6->SetInput(Dz);
    filter6->SetSigma(1);
    filter6->Update();
    ProbImagePointer Dzz = filter6->GetOutput(); 


    typedef itk::ImageSliceIteratorWithIndex< ImageType3D > SliceIteratorType;
    typedef itk::ImageSliceIteratorWithIndex< GradientImageType > GradientSliceIteratorType;
    typedef itk::ImageSliceIteratorWithIndex< ImageType3D > ISliceIteratorType;

    SliceIteratorType inputIt1( Dxx, Dxx->GetRequestedRegion() );
    inputIt1.SetFirstDirection( 0 );
    inputIt1.SetSecondDirection( 1 );
    SliceIteratorType inputIt2( Dxy, Dxy->GetRequestedRegion() );
    inputIt2.SetFirstDirection( 0 );
    inputIt2.SetSecondDirection( 1 );
    SliceIteratorType inputIt3( Dxz, Dxz->GetRequestedRegion() );
    inputIt3.SetFirstDirection( 0 );
    inputIt3.SetSecondDirection( 1 );
    SliceIteratorType inputIt4( Dyy, Dyy->GetRequestedRegion() );
    inputIt4.SetFirstDirection( 0 );
    inputIt4.SetSecondDirection( 1 );
    SliceIteratorType inputIt5( Dyz, Dyz->GetRequestedRegion() );
    inputIt5.SetFirstDirection( 0 );
    inputIt5.SetSecondDirection( 1 );
    SliceIteratorType inputIt6( Dzz, Dzz->GetRequestedRegion() );
    inputIt6.SetFirstDirection( 0 );
    inputIt6.SetSecondDirection( 1 );

    ISliceIteratorType inputIt7(_PaddedCurvImage, _PaddedCurvImage->GetRequestedRegion() );
    inputIt7.SetFirstDirection( 0 );
    inputIt7.SetSecondDirection( 1 );

    /*SliceIteratorType outputIt1(Vx1, Vx1->GetRequestedRegion() );
      outputIt1.SetFirstDirection( 2 );
      outputIt1.SetSecondDirection( 1 );
      SliceIteratorType outputIt2(Vy1, Vy1->GetRequestedRegion() );
      outputIt2.SetFirstDirection( 2 );
      outputIt2.SetSecondDirection( 1 );
      SliceIteratorType outputIt3(Vz1, Vz1->GetRequestedRegion() );
      outputIt3.SetFirstDirection( 2 );
      outputIt3.SetSecondDirection( 1 );*/

    GradientSliceIteratorType inputIt0( _IGVF, _IGVF->GetRequestedRegion() );
    inputIt0.SetFirstDirection( 0 );
    inputIt0.SetSecondDirection( 1 );

    SliceIteratorType outputIt4( _IVessel, _IVessel->GetRequestedRegion() );
    outputIt4.SetFirstDirection( 0 );
    outputIt4.SetSecondDirection( 1 );

    inputIt0.GoToBegin();
    inputIt1.GoToBegin();
    inputIt2.GoToBegin();
    inputIt3.GoToBegin();
    inputIt4.GoToBegin();
    inputIt5.GoToBegin();
    inputIt6.GoToBegin();
    inputIt7.GoToBegin();

    outputIt4.GoToBegin();

    double Ma[3][3];
    double V[3][3];
    double d[3];

    while( !inputIt1.IsAtEnd() )
    {
      while ( !inputIt1.IsAtEndOfSlice() )
      {
        while ( !inputIt1.IsAtEndOfLine() )
        {
          if( inputIt7.Get() == 0 )
            //if( 0 )
          {
            outputIt4.Set(0);
          }
          else
          {
            Ma[0][0] = inputIt1.Get();
            Ma[0][1] = inputIt2.Get();
            Ma[0][2] = inputIt3.Get();
            Ma[1][0] = inputIt2.Get();
            Ma[1][1] = inputIt4.Get();
            Ma[1][2] = inputIt5.Get();
            Ma[2][0] = inputIt3.Get();
            Ma[2][1] = inputIt5.Get();
            Ma[2][2] = inputIt6.Get();
            eigen_decomposition(Ma,V,d);

            double Lambda1 = d[0];
            double Lambda2 = d[1];
            double Lambda3 = d[2];
            if(Lambda2>=0.0 || Lambda3>=0.0)
            {
              outputIt4.Set(0);
            }
            else
            {
              double Ra  = Lambda2 / Lambda3; 
              double Rb  = Lambda1 / vcl_sqrt ( vnl_math_abs( Lambda2 * Lambda3 )); 
              double S  = vcl_sqrt( pow(Lambda1,2) + pow(Lambda2,2) + pow(Lambda3,2) );
              double vesMeasure_1  = 
                ( 1 - vcl_exp(-1.0*(( vnl_math_sqr( Ra ) ) / ( A ))));
              double vesMeasure_2  = 
                vcl_exp ( -1.0 * ((vnl_math_sqr( Rb )) /  ( B )));
              double vesMeasure_3  = 
                ( 1 - vcl_exp( -1.0 * (( vnl_math_sqr( S )) / ( C ))));

              float V_Saliency = vesMeasure_1 * vesMeasure_2 * vesMeasure_3;
              //float V_Saliency = fabs(inputIt0.Get()[0] * (float)V[0][0] + inputIt0.Get()[1] * (float)V[0][1] + inputIt0.Get()[2] * (float)V[0][2]);
              //std::cout<<Lambda1<<","<<Lambda2<<","<<Lambda3<<std::endl;
              //std::cout<<V_Saliency<<std::endl;
              outputIt4.Set(V_Saliency);
            }
          }
          ++inputIt0;
          ++inputIt1;
          ++inputIt2;
          ++inputIt3;
          ++inputIt4;
          ++inputIt5;
          ++inputIt6;
          ++inputIt7;
          ++outputIt4;
        }
        inputIt0.NextLine();
        inputIt1.NextLine();
        inputIt2.NextLine();
        inputIt3.NextLine();
        inputIt4.NextLine();
        inputIt5.NextLine();
        inputIt6.NextLine();
        inputIt7.NextLine();
        outputIt4.NextLine();
      }
      inputIt0.NextSlice();
      inputIt1.NextSlice();
      inputIt2.NextSlice();
      inputIt3.NextSlice();
      inputIt4.NextSlice();
      inputIt5.NextSlice();
      inputIt6.NextSlice();
      inputIt7.NextSlice();
      outputIt4.NextSlice();
    }

    //rescale the vesselness to [0,255]
    typedef itk::RescaleIntensityImageFilter< ImageType3D, ImageType3D> RescaleFilterType;
    RescaleFilterType::Pointer rescale = RescaleFilterType::New();
    rescale->SetInput( _IVessel );
    rescale->SetOutputMinimum( 0 );
    rescale->SetOutputMaximum( 255 );
    rescale->Update();
    _IVessel = rescale->GetOutput();

    //typedef itk::OtsuThresholdImageFilter<
    //  ImageType3D, ImageType3D > TH_FilterType;
    //TH_FilterType::Pointer thresholder = TH_FilterType::New();
    //thresholder->SetInput( _IVessel );
    //thresholder->SetOutsideValue( 1 );
    //thresholder->SetInsideValue( 0 );
    //thresholder->SetNumberOfHistogramBins( 256 );
    //thresholder->Update();

    //v_threshold = thresholder->GetThreshold();

    //std::cout<<" - Selected Otsu Threshold:"<<v_threshold<<std::endl;




  }



  void MultipleNeuronTracer::SeedDetection(float th, int detection_method, int seed_radius)
  {

    SM = _PaddedCurvImage->GetLargestPossibleRegion().GetSize()[0];
    SN = _PaddedCurvImage->GetLargestPossibleRegion().GetSize()[1];
    SZ = _PaddedCurvImage->GetLargestPossibleRegion().GetSize()[2];

    SeedType seed;
    SamplePointer seeds_candidate = SampleType::New();
    seeds_candidate->SetMeasurementVectorSize( 3 );

    //int max_num_seed = 100000;
    //SeedPt.SetN(max_num_seed);

    //int remove_radius = 5;

    int rad = 1;

    typedef itk::CastImageFilter<ImageType3D,LabelImageType3D> CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(_PaddedCurvImage);
    caster->Update();
    LabelImageType3D::Pointer I_Seed = caster->GetOutput();

    typedef itk::NeighborhoodIterator< ImageType3D > NeighborhoodIteratorType;
    typedef itk::NeighborhoodIterator< ImageType3D > NeighborhoodIteratorType1;
    typedef itk::NeighborhoodIterator< GradientImageType > NeighborhoodIteratorType2;

    NeighborhoodIteratorType::RadiusType radius;
    NeighborhoodIteratorType1::RadiusType radius1;
    NeighborhoodIteratorType2::RadiusType radius2;
    radius.Fill(rad);
    radius1.Fill(rad);
    radius2.Fill(rad);

    NeighborhoodIteratorType it( radius, _IVessel, _IVessel->GetRequestedRegion());
    //NeighborhoodIteratorType2 it21( radius2, IGVF, IGVF->GetRequestedRegion());
    //NeighborhoodIteratorType2 it22( radius2, V1, V1->GetRequestedRegion());
    //std::cout<<"SeedDetection-Threshold:"<<th<<std::endl;///
    //std::cout<<"SeedDetection-Seed Radius:"<<seed_radius<<std::endl;///

    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
      /*if( detection_method == 0 ) //skeleton seeds
        {
        if(it1.GetCenterPixel() == 1)
        {
      //SeedPt.AddPt(it.GetIndex()[0],it.GetIndex()[1],it.GetIndex()[2]);
      seed[0] = it.GetIndex()[0];
      seed[1] = it.GetIndex()[1];
      seed[2] = it.GetIndex()[2];
      if( !SeedSparsify( seeds_candidate, seed, seed_radius ) )
      seeds_candidate->PushBack( seed );
      }
      }*/

      if(it.GetCenterPixel() >= th )
        //if( it.GetCenterPixel() >= th )
      {
        //SeedPt.AddPt(it.GetIndex()[0],it.GetIndex()[1],it.GetIndex()[2]);
        seed[0] = it.GetIndex()[0];
        seed[1] = it.GetIndex()[1];
        seed[2] = it.GetIndex()[2];
        if( seed_radius != 0 )
        {
          if( !SeedSparsify( seeds_candidate, seed, seed_radius ) )
            seeds_candidate->PushBack( seed );
        }
        else
        {
          seeds_candidate->PushBack( seed );
        }
      } 

    }

    //SeedPt.SetN( seeds_candidate->Size() );

    /***** Update the code to write it to NDXImage  ***** MURAD ****/
    SeedPt.RemoveAllPts();
    for( unsigned int i = 0; i < seeds_candidate->Size(); i++ )
    {
      SeedPt.AddPt( seeds_candidate->GetMeasurementVector(i)[0], seeds_candidate->GetMeasurementVector(i)[1], seeds_candidate->GetMeasurementVector(i)[2] );	
    }

    visit_label.set_size(SeedPt.GetSize());
    visit_label.fill(0);

    //SeedPt_mg = SeedPt;
    //std::cout<<"Detected Seed Points:"<<SeedPt.NP<<std::endl;///

  }



  bool MultipleNeuronTracer::SeedSparsify(SamplePointer seeds_candidate, SeedType query_point, int radius)
  {
    bool removal = false;
    Point3D temp_pt1, temp_pt2;
    temp_pt1.x = query_point[0];
    temp_pt1.y = query_point[1];
    temp_pt1.z = query_point[2];

    for( unsigned int i = 0; i < seeds_candidate->Size(); i++ )
    {
      SeedType seed;
      seed = seeds_candidate->GetMeasurementVector(i);
      temp_pt2.x = seed[0];
      temp_pt2.y = seed[1];
      temp_pt2.z = seed[2];
      if( temp_pt1.GetDistTo(temp_pt2) < radius )
      {
        removal = true;
        break;
      }
    }

    /*typedef itk::Statistics::KdTreeGenerator< SampleType > TreeGeneratorType;
      TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();

      treeGenerator->SetSample( seeds_candidate );
      treeGenerator->SetBucketSize( 16 );
      treeGenerator->Update();

      typedef TreeGeneratorType::KdTreeType TreeType;
      TreeType::Pointer tree = treeGenerator->GetOutput();

      TreeType::InstanceIdentifierVectorType neighbors;
      tree->Search( query_point, (double)radius, neighbors );

      if( neighbors.size() > 0 )
      {
      removal = true;
      } */
    return removal;
  }




  ProbImagePointer extract_one_component(int index, GradientImagePointer G)
  {
    typedef itk::ImageAdaptor<  GradientImageType, 
            VectorPixelAccessor > ImageAdaptorType;
    ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
    VectorPixelAccessor  accessor;
    accessor.SetIndex( index );
    adaptor->SetPixelAccessor( accessor );
    adaptor->SetImage( G );
    typedef itk::CastImageFilter<ImageAdaptorType,ImageType3D> CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(adaptor);
    caster->Update();
    ProbImagePointer output = caster->GetOutput();

    return output;
  }






  void MultipleNeuronTracer::SeedAdjustment(int iter_num)
  {

    normalizeGVF();

    if( iter_num == 0 )
    {
      //sort seeds by their saliency
      PointList3D New_SeedPt;

      vnl_vector<double> saliency(SeedPt.NP);
      for( int i = 0; i < SeedPt.NP; i++)
      {
        GradientImageType::IndexType index; 
        SeedPt.Pt[i].check_out_of_range_3D(SM,SN,SZ);
        index[0] = ceil(SeedPt.Pt[i].x);
        index[1] = ceil(SeedPt.Pt[i].y);
        index[2] = ceil(SeedPt.Pt[i].z);
        saliency(i) = _IVessel->GetPixel(index);
      }

      for( unsigned int i = 0; i < saliency.size(); i++)
      {
        int index = saliency.arg_max();
        New_SeedPt.AddPt(SeedPt.Pt[index]);
        saliency(index) = -1;
      } 

      SeedPt = New_SeedPt;

      visit_label.set_size(SeedPt.GetSize());
      visit_label.fill(0);
      return;
    }

    //int iter_num = 100;

    typedef itk::VectorLinearInterpolateImageFunction< 
      GradientImageType, float >  GradientInterpolatorType;

    GradientInterpolatorType::Pointer interpolator = GradientInterpolatorType::New();
    interpolator->SetInputImage(_IGVF);


    //move seeds along gradient vector flow
    int j = 0;

    Point3D temp_pt;

    while( j < iter_num )
    {
      for( int i = 0; i < SeedPt.NP; i++ )
      {
        GradientImageType::IndexType index; 
        // std::cout<<SeedPt.Pt[i].x<<","<<SeedPt.Pt[i].y<<","<<SeedPt.Pt[i].z<<std::endl;
        SeedPt.Pt[i].check_out_of_range_3D(SM,SN,SZ);
        index[0] = (SeedPt.Pt[i].x);
        index[1] = (SeedPt.Pt[i].y);
        index[2] = (SeedPt.Pt[i].z);
        //std::cout<<index[0]<<","<<index[1]<<","<<index[2]<<std::endl;
        GradientPixelType gradient = interpolator->EvaluateAtIndex(index);
        //std::cout<<"check point 20"<<std::endl;
        //GradientPixelType gradient = IGVF->GetPixel(index);
        //std::cout<<gradient[0]<<","<<gradient[1]<<","<<gradient[2]<<std::endl;
        //Point3D temp_pt(gradient[0],gradient[1],gradient[2]);
        temp_pt.x = gradient[0];
        temp_pt.y = gradient[1];
        temp_pt.z = gradient[2];
        SeedPt.Pt[i] = SeedPt.Pt[i] + temp_pt;
      }
      j++;
    }

    /*//filter out seeds in the background
      int i = 0;
      while( i < SeedPt.NP )
      {
      ImageType::IndexType index; 

      SeedPt.Pt[i].check_out_of_range_3D(SM,SN,SZ);

      index[0] = ceil(SeedPt.Pt[i].x);
      index[1] = ceil(SeedPt.Pt[i].y);
      index[2] = ceil(SeedPt.Pt[i].z);

      if( VBW->GetPixel(index) == 0)
      {    
      SeedPt.RemovePt(i);
      continue;
      }
      i++;
      }*/


    seed_centroid();

    //sort seeds by their saliency
    PointList3D New_SeedPt;

    vnl_vector<double> saliency(SeedPt.NP);
    unsigned long int number_of_out_of_range=0;///
    //std::cout<<"Number of seeds:"<<SeedPt.NP<<std::endl;///
    for( unsigned long int i = 0; i < SeedPt.NP; i++)///
    {
      GradientImageType::IndexType index; 
      SeedPt.Pt[i].check_out_of_range_3D(SM,SN,SZ);
      if (SeedPt.Pt[i].check_out_of_range_3D(SM,SN,SZ))
        number_of_out_of_range++;
      //std::cout<<"Out of Range"<<std::endl;
      index[0] = ceil(SeedPt.Pt[i].x);
      index[1] = ceil(SeedPt.Pt[i].y);
      index[2] = ceil(SeedPt.Pt[i].z);
      saliency(i) = _IVessel->GetPixel(index);
    }
    //std::cout<<"After for on SeedPt.NP"<<std::endl;
    //std::cout<<"Number of out of range:"<<number_of_out_of_range<<std::endl;

    for( unsigned long int i = 0; i < saliency.size(); i++)////
    {
      unsigned long int index = saliency.arg_max();
      New_SeedPt.AddPt(SeedPt.Pt[index]);
      saliency(index) = -1;
    } 
    //std::cout<<"After for on saliency size"<<std::endl;///

    SeedPt = New_SeedPt;

    visit_label.set_size(SeedPt.GetSize());
    //std::cout<<"After set size"<<std::endl;///

    visit_label.fill(0);

    //SeedPt_mg = SeedPt;
    //std::cout<<"Detected Seed Points:"<<SeedPt.NP<<std::endl;
  } 


  void MultipleNeuronTracer::seed_centroid()
  {
    typedef itk::CastImageFilter<ImageType3D,LabelImageType3D> CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(_PaddedCurvImage);
    caster->Update();
    LabelImageType3D::Pointer I_Seed = caster->GetOutput();

    typedef itk::ImageSliceIteratorWithIndex< LabelImageType3D > SliceIteratorType;
    SliceIteratorType inputIt( I_Seed, I_Seed->GetRequestedRegion() );
    inputIt.SetFirstDirection( 0 );
    inputIt.SetSecondDirection( 1 );
    inputIt.GoToBegin();
    while( !inputIt.IsAtEnd() )
    {
      while ( !inputIt.IsAtEndOfSlice() )
      {
        while ( !inputIt.IsAtEndOfLine() )
        {
          inputIt.Set(0);
          ++inputIt;
        }
        inputIt.NextLine();
      }
      inputIt.NextSlice();
    }

    for( int i = 0; i < SeedPt.GetSize(); i++ )
    {
      LabelImageType3D::IndexType index;
      SeedPt.Pt[i].check_out_of_range_3D(SM,SN,SZ);
      /*for( int j = -1; j <= 1; j++ )
        {
        for( int k = -1; k <= 1; k++ )
        {
        for( int z = -1; z <= 1; z++)
        {
        index[0] = SeedPt.Pt[i].x + j;
        index[1] = SeedPt.Pt[i].y + k;
        index[2] = SeedPt.Pt[i].z + z;

        if( index[0] < 0 || index[1] < 0 || index[2] < 0 || 
        index[0] >= SM || index[1] >= SN || index[2] >= SZ )
        continue;

        I_Seed->SetPixel(index, 1);
        }
        }
        }*/

      index[0] = SeedPt.Pt[i].x;
      index[1] = SeedPt.Pt[i].y;
      index[2] = SeedPt.Pt[i].z;
      I_Seed->SetPixel(index,1);
    }

    // Set up a connected components filter to label the binary objects.
    typedef itk::ConnectedComponentImageFilter< LabelImageType3D, LabelImageTypeCCIF > ConnectedComponentType;//////////////////
    ConnectedComponentType::Pointer connectedComponentFilter = ConnectedComponentType::New();


    try 
    { 
      //std::cout<<"Before connectedComponentFilter->SetInput()"<< std::endl;
      connectedComponentFilter->SetInput( I_Seed );
    } 
    catch( itk::ExceptionObject & err ) 
    { 
      std::cerr << "ExceptionObject caught at connectedComponentFilter GetOutput!" << std::endl; 
      std::cerr << err << std::endl; 
    } 


    // Relabel the components in order of size.
    typedef itk::RelabelComponentImageFilter< LabelImageTypeCCIF, LabelImageTypeCCIF > RelabelType;///////////
    RelabelType::Pointer relabeler = RelabelType::New();

    try 
    { 
      //std::cout<<"Before relaberer->GetOutput()"<< std::endl;
      relabeler->SetInput( connectedComponentFilter->GetOutput() );
    } 
    catch( itk::ExceptionObject & err ) 
    { 
      std::cerr << "ExceptionObject caught at relaberer GetOutput!" << std::endl; 
      std::cerr << err << std::endl; 
    } 



    typedef itk::LabelGeometryImageFilter< LabelImageTypeCCIF > LabelGeometryType;///////////////////
    LabelGeometryType::Pointer labelGeometryFilter = LabelGeometryType::New();

    try 
    { 
      //std::cout<<"Before labelGeometryFilter->GetOutput()"<< std::endl;
      labelGeometryFilter->SetInput( relabeler->GetOutput() ); 
    } 
    catch( itk::ExceptionObject & err ) 
    { 
      std::cerr << "ExceptionObject caught at labelGeometryFilter GetOutput!" << std::endl; 
      std::cerr << err << std::endl; 
    } 





    try 
    { 
      //std::cout<<"Before labelGeometryFilter->Update()"<< std::endl;
      labelGeometryFilter->Update(); 
    } 
    catch( itk::ExceptionObject & err ) 
    { 
      std::cerr << "ExceptionObject caught at labelGeometryFilter GetOutput!" << std::endl; 
      std::cerr << err << std::endl; 
    } 

    //std::cout<<"After labelGeometryFilter->Update()"<< std::endl;


    LabelGeometryType::LabelPointType index_temp;
    int labelValue = labelGeometryFilter->GetNumberOfLabels()-1;

    //std::cout<<"After labelGeometryFilter->GetNumberofLabels()"<< std::endl;

    Point3D temp;
    //PointList3D new_SeedPt(labelValue+1);
    SeedPt.RemoveAllPts();
    for( int i = 1; i <= labelValue; i++)
    {
      index_temp = labelGeometryFilter->GetCentroid(i);
      temp.x = index_temp[0];
      temp.y = index_temp[1];
      temp.z = index_temp[2];
      SeedPt.AddPt(temp);
    }
    std::cout<<"Seed Centroid End"<< std::endl;

  }
  void MultipleNeuronTracer::normalizeGVF()
  {
    typedef itk::ImageSliceIteratorWithIndex< GradientImageType > SliceIteratorType;
    SliceIteratorType inputIt( _IGVF, _IGVF->GetRequestedRegion() );
    inputIt.SetFirstDirection( 0 );
    inputIt.SetSecondDirection( 1 );
    inputIt.GoToBegin();
    while( !inputIt.IsAtEnd() )
    {
      while ( !inputIt.IsAtEndOfSlice() )
      {
        while ( !inputIt.IsAtEndOfLine() )
        {
          if(inputIt.Get().GetNorm() != 0)
          {
            inputIt.Set( inputIt.Get()/( inputIt.Get().GetNorm() + std::numeric_limits<float>::epsilon() ) );
          }
          ++inputIt;
        }
        inputIt.NextLine();
      }
      inputIt.NextSlice();
    }

    /* typedef itk::ImageDuplicator< GradientImageType > DuplicatorType;
       DuplicatorType::Pointer Duplicator = DuplicatorType::New();
       Duplicator->SetInputImage(IGVF);
       Duplicator->Update();
       IGVF_Norm = Duplicator->GetOutput();


       typedef itk::CastImageFilter<ImageType,ProbImageType> CasterType;
       CasterType::Pointer caster = CasterType::New();

       caster->SetInput(I);
       caster->Update();
       GVF_Mag = caster->GetOutput();


       typedef itk::ImageSliceIteratorWithIndex< GradientImageType > SliceIteratorType;
       typedef itk::ImageSliceIteratorWithIndex< ProbImageType > SliceIteratorType1;

       SliceIteratorType inputIt( IGVF, IGVF->GetRequestedRegion() );
       SliceIteratorType outputIt( IGVF_Norm, IGVF_Norm->GetRequestedRegion() );
       SliceIteratorType1 outputIt1( GVF_Mag, GVF_Mag->GetRequestedRegion() );
       inputIt.SetFirstDirection( 0 );
       inputIt.SetSecondDirection( 1 );
       outputIt.SetFirstDirection( 0 );
       outputIt.SetSecondDirection( 1 );
       outputIt1.SetFirstDirection( 0 );
       outputIt1.SetSecondDirection( 1 );

       inputIt.GoToBegin();
       outputIt.GoToBegin();
       outputIt1.GoToBegin();

       while( !inputIt.IsAtEnd() )
       {
       while ( !inputIt.IsAtEndOfSlice() )
       {
       while ( !inputIt.IsAtEndOfLine() )
       {
       if(inputIt.Get().GetNorm() == 0)
       {
       outputIt.Set( inputIt.Get());
       }
       else
       {
       outputIt.Set( inputIt.Get()/inputIt.Get().GetNorm() );
       }


       outputIt1.Set( inputIt.Get().GetNorm() );
       ++inputIt;
       ++outputIt;
       ++outputIt1;
       }
       inputIt.NextLine();
       outputIt.NextLine();
       outputIt1.NextLine();
       }
       inputIt.NextSlice();
       outputIt.NextSlice();
       outputIt1.NextSlice();
       } */
  }


  ProbImagePointer MultipleNeuronTracer::Upsampling(ProbImagePointer inputImage,int upsampleFactor)
  {

    //		int upsampleFactor = 2;
    typedef itk::IdentityTransform<double, 3>	T_Transform;
    typedef itk::BSplineInterpolateImageFunction<ImageType3D, double, double>	T_Interpolator;
    typedef itk::ResampleImageFilter<ImageType3D, ImageType3D>	T_ResampleFilter;


    // Prepare the resampler.

    // Instantiate the transform and specify it should be the id transform.
    T_Transform::Pointer pTransform = T_Transform::New();
    pTransform->SetIdentity();

    // Instantiate the b-spline interpolator and set it as the third order
    // for bicubic.
    T_Interpolator::Pointer pInterpolator = T_Interpolator::New();
    pInterpolator->SetSplineOrder(3);

    // Instantiate the resampler. Wire in the transform and the interpolator.
    T_ResampleFilter::Pointer pResizeFilter = T_ResampleFilter::New();
    pResizeFilter->SetTransform(pTransform);
    pResizeFilter->SetInterpolator(pInterpolator);

    // Specify the input.

    pResizeFilter->SetInput(inputImage);


    const double vfOutputOrigin[3]  = { 0.0, 0.0, 0.0};
    pResizeFilter->SetOutputOrigin(vfOutputOrigin);


    // Fetch original image size.
    ImageType3D::RegionType inputRegion = inputImage->GetLargestPossibleRegion();
    itk::Size<3> inputSize = inputImage->GetBufferedRegion().GetSize();

    // Fetch original image spacing.
    ImageType3D::SpacingType vfInputSpacing = inputImage->GetSpacing();

    double vfOutputSpacing[3];
    vfOutputSpacing[0] = double(vfInputSpacing[0] / upsampleFactor);
    vfOutputSpacing[1] = double(vfInputSpacing[1] / upsampleFactor);
    vfOutputSpacing[2] = double(vfInputSpacing[2] / upsampleFactor);

    // Set the output spacing. If you comment out the following line, the original
    // image will be simply put in the upper left corner of the new image without
    // any scaling.
    pResizeFilter->SetOutputSpacing(vfOutputSpacing);

    // Set the output size as specified on the command line.

    itk::Size<3> vnOutputSize = inputImage->GetBufferedRegion().GetSize();
    vnOutputSize[0] = inputSize[0]*upsampleFactor;
    vnOutputSize[1] = inputSize[1]*upsampleFactor;
    vnOutputSize[2] = inputSize[2]*upsampleFactor;

    pResizeFilter->SetSize(vnOutputSize);

    // Specify the input.

    pResizeFilter->SetInput(inputImage);



    std::cout<<"in upscale before write"<<std::endl;
    /*ImageType3D::Pointer pWriter = ImageType3D::New();
      pWriter->SetFileName("test.tif");
      pWriter->SetInput(pResizeFilter->GetOutput());
      std::cout<<"in upscale after write"<<std::endl;*/
    //WriteImage3D(std::string("test.tif"), pResizeFilter->GetOutput());
    ProbImagePointer outputImage = pResizeFilter->GetOutput();
    return outputImage;	 



  }
  //LabelImagePointer MultipleNeuronTracer::UpsamplingLabel(ProbImagePointer inputImage,int upsampleFactor){
  //	
  ////		int upsampleFactor = 2;
  //	 typedef itk::IdentityTransform<double, 3>	T_Transform;
  //	 typedef itk::BSplineInterpolateImageFunction<LabelImageType3D, double, double>	T_Interpolator;
  //	 typedef itk::ResampleImageFilter<LabelImageType3D, LabelImageType3D>	T_ResampleFilter;
  //
  //
  //	 // Prepare the resampler.
  //
  //	 // Instantiate the transform and specify it should be the id transform.
  //	 T_Transform::Pointer pTransform = T_Transform::New();
  //	 pTransform->SetIdentity();
  //
  //	 // Instantiate the b-spline interpolator and set it as the third order
  //	 // for bicubic.
  //	 T_Interpolator::Pointer pInterpolator = T_Interpolator::New();
  //	 pInterpolator->SetSplineOrder(3);
  //
  //	 // Instantiate the resampler. Wire in the transform and the interpolator.
  //	 T_ResampleFilter::Pointer pResizeFilter = T_ResampleFilter::New();
  //	 pResizeFilter->SetTransform(pTransform);
  //	 pResizeFilter->SetInterpolator(pInterpolator);
  //
  //	 // Specify the input.
  //
  //	 pResizeFilter->SetInput(inputImage);
  //
  //
  //	 const double vfOutputOrigin[3]  = { 0.0, 0.0, 0.0};
  //	 pResizeFilter->SetOutputOrigin(vfOutputOrigin);
  //
  //
  //	 // Fetch original image size.
  //	 ImageType3D::RegionType inputRegion = inputImage->GetLargestPossibleRegion();
  //	itk::Size<3> inputSize = inputImage->GetBufferedRegion().GetSize();
  //
  //	 // Fetch original image spacing.
  //	 ImageType3D::SpacingType vfInputSpacing = inputImage->GetSpacing();
  //
  //	 double vfOutputSpacing[3];
  //	 vfOutputSpacing[0] = double(vfInputSpacing[0] / upsampleFactor);
  //	 vfOutputSpacing[1] = double(vfInputSpacing[1] / upsampleFactor);
  //	 vfOutputSpacing[2] = double(vfInputSpacing[2] / upsampleFactor);
  //
  //	 // Set the output spacing. If you comment out the following line, the original
  //	 // image will be simply put in the upper left corner of the new image without
  //	 // any scaling.
  //	 pResizeFilter->SetOutputSpacing(vfOutputSpacing);
  //
  //	 // Set the output size as specified on the command line.
  //
  //	 itk::Size<3> vnOutputSize = inputImage->GetBufferedRegion().GetSize();
  //	 vnOutputSize[0] = inputSize[0]*upsampleFactor;
  //	 vnOutputSize[1] = inputSize[1]*upsampleFactor;
  //	 vnOutputSize[2] = inputSize[2]*upsampleFactor;
  //
  //	 pResizeFilter->SetSize(vnOutputSize);
  //
  //	 // Specify the input.
  //
  //	 pResizeFilter->SetInput(inputImage);
  //
  //
  //
  //	std::cout<<"in upscale before write"<<std::endl;
  //	/*ImageType3D::Pointer pWriter = ImageType3D::New();
  //	pWriter->SetFileName("test.tif");
  //	pWriter->SetInput(pResizeFilter->GetOutput());
  //	std::cout<<"in upscale after write"<<std::endl;*/
  //	//WriteImage3D(std::string("test.tif"), pResizeFilter->GetOutput());
  //	LabelImagePointer outputImage = pResizeFilter->GetOutput();
  //	return outputImage;	 
  //}
  //

  //}
void eigen_decomposition2D(double A[2][2], double V[2][2], double d[2])
{
  double tmp = sqrt(pow(A[0][0]-A[1][1],2) + 4 * pow(A[0][1],2));
  double v1x = 2 * A[1][0];
  double v1y = A[1][1] - A[0][0] + tmp;
  double mag = sqrt(pow(v1x,2)+pow(v1y,2));

  v1x /= (mag + std::numeric_limits<float>::epsilon());
  v1y /= (mag + std::numeric_limits<float>::epsilon());

  double v2x = -1 * v1y;
  double v2y = v1x;

  double mu1 = 0.5 * (A[0][0] + A[1][1] + tmp);
  double mu2 = 0.5 * (A[0][0] + A[1][1] - tmp);

  if( fabs(mu1) > fabs(mu2) )
  {
    V[0][0] = v2x;
    V[1][0] = v2y;
    V[0][1] = v1x;
    V[1][1] = v1y;
    d[0] = mu2;
    d[1] = mu1;
  }
  else
  {
    V[0][0] = v1x;
    V[1][0] = v1y;
    V[0][1] = v2x;
    V[1][1] = v2y;
    d[0] = mu1;
    d[1] = mu2;
  }

}


void eigen_decomposition(double A[3][3], double V[3][3], double d[3]) {
  int n = 3;
  double e[3];
  double da[3];
  double dt, dat;
  double vet[3];
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      V[i][j] = A[i][j];
    }
  }
  tred2(V, d, e);
  tql2(V, d, e);

  /* Sort the eigen values and vectors by abs eigen value */
  da[0]=absd(d[0]); da[1]=absd(d[1]); da[2]=absd(d[2]);
  if((da[0]>=da[1])&&(da[0]>da[2]))
  {
    dt=d[2];   dat=da[2];    vet[0]=V[0][2];    vet[1]=V[1][2];    vet[2]=V[2][2];
    d[2]=d[0]; da[2]=da[0];  V[0][2] = V[0][0]; V[1][2] = V[1][0]; V[2][2] = V[2][0];
    d[0]=dt;   da[0]=dat;    V[0][0] = vet[0];  V[1][0] = vet[1];  V[2][0] = vet[2]; 
  }
  else if((da[1]>=da[0])&&(da[1]>da[2]))  
  {
    dt=d[2];   dat=da[2];    vet[0]=V[0][2];    vet[1]=V[1][2];    vet[2]=V[2][2];
    d[2]=d[1]; da[2]=da[1];  V[0][2] = V[0][1]; V[1][2] = V[1][1]; V[2][2] = V[2][1];
    d[1]=dt;   da[1]=dat;    V[0][1] = vet[0];  V[1][1] = vet[1];  V[2][1] = vet[2]; 
  }
  if(da[0]>da[1])
  {
    dt=d[1];   dat=da[1];    vet[0]=V[0][1];    vet[1]=V[1][1];    vet[2]=V[2][1];
    d[1]=d[0]; da[1]=da[0];  V[0][1] = V[0][0]; V[1][1] = V[1][0]; V[2][1] = V[2][0];
    d[0]=dt;   da[0]=dat;    V[0][0] = vet[0];  V[1][0] = vet[1];  V[2][0] = vet[2]; 
  }
}

static void tred2(double V[3][3], double d[3], double e[3]) {

  int n = 3;
  /*  This is derived from the Algol procedures tred2 by */
  /*  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for */
  /*  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding */
  /*  Fortran subroutine in EISPACK. */
  int i, j, k;
  double scale;
  double f, g, h;
  double hh;
  for (j = 0; j < n; j++) {d[j] = V[n-1][j]; }

  /* Householder reduction to tridiagonal form. */

  for (i = n-1; i > 0; i--) {
    /* Scale to avoid under/overflow. */
    scale = 0.0;
    h = 0.0;
    for (k = 0; k < i; k++) { scale = scale + fabs(d[k]); }
    if (scale == 0.0) {
      e[i] = d[i-1];
      for (j = 0; j < i; j++) { d[j] = V[i-1][j]; V[i][j] = 0.0;  V[j][i] = 0.0; }
    } else {

      /* Generate Householder vector. */

      for (k = 0; k < i; k++) { d[k] /= scale; h += d[k] * d[k]; }
      f = d[i-1];
      g = sqrt(h);
      if (f > 0) { g = -g; }
      e[i] = scale * g;
      h = h - f * g;
      d[i-1] = f - g;
      for (j = 0; j < i; j++) { e[j] = 0.0; }

      /* Apply similarity transformation to remaining columns. */

      for (j = 0; j < i; j++) {
        f = d[j];
        V[j][i] = f;
        g = e[j] + V[j][j] * f;
        for (k = j+1; k <= i-1; k++) { g += V[k][j] * d[k]; e[k] += V[k][j] * f; }
        e[j] = g;
      }
      f = 0.0;
      for (j = 0; j < i; j++) { e[j] /= h; f += e[j] * d[j]; }
      hh = f / (h + h);
      for (j = 0; j < i; j++) { e[j] -= hh * d[j]; }
      for (j = 0; j < i; j++) {
        f = d[j]; g = e[j];
        for (k = j; k <= i-1; k++) { V[k][j] -= (f * e[k] + g * d[k]); }
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
      }
    }
    d[i] = h;
  }

  /* Accumulate transformations. */

  for (i = 0; i < n-1; i++) {
    V[n-1][i] = V[i][i];
    V[i][i] = 1.0;
    h = d[i+1];
    if (h != 0.0) {
      for (k = 0; k <= i; k++) { d[k] = V[k][i+1] / h;}
      for (j = 0; j <= i; j++) {
        g = 0.0;
        for (k = 0; k <= i; k++) { g += V[k][i+1] * V[k][j]; }
        for (k = 0; k <= i; k++) { V[k][j] -= g * d[k]; }
      }
    }
    for (k = 0; k <= i; k++) { V[k][i+1] = 0.0;}
  }
  for (j = 0; j < n; j++) { d[j] = V[n-1][j]; V[n-1][j] = 0.0; }
  V[n-1][n-1] = 1.0;
  e[0] = 0.0;
}

/* Symmetric tridiagonal QL algorithm. */
static void tql2(double V[3][3], double d[3], double e[3]) {

  int n = 3;
  /*  This is derived from the Algol procedures tql2, by */
  /*  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for */
  /*  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding */
  /*  Fortran subroutine in EISPACK. */

  int i, j, k, l, m;
  double f;
  double tst1;
  double eps;
  int iter;
  double g, p, r;
  double dl1, h, c, c2, c3, el1, s, s2;

  for (i = 1; i < n; i++) { e[i-1] = e[i]; }
  e[n-1] = 0.0;

  f = 0.0;
  tst1 = 0.0;
  eps = pow(2.0, -52.0);
  for (l = 0; l < n; l++) {

    /* Find small subdiagonal element */

    tst1 = MAX_1(tst1, fabs(d[l]) + fabs(e[l]));
    m = l;
    while (m < n) {
      if (fabs(e[m]) <= eps*tst1) { break; }
      m++;
    }

    /* If m == l, d[l] is an eigenvalue, */
    /* otherwise, iterate. */

    if (m > l) {
      iter = 0;
      do {
        iter = iter + 1;  /* (Could check iteration count here.) */
        /* Compute implicit shift */
        g = d[l];
        p = (d[l+1] - g) / (2.0 * e[l]);
        r = hypot2(p, 1.0);
        if (p < 0) { r = -r; }
        d[l] = e[l] / (p + r);
        d[l+1] = e[l] * (p + r);
        dl1 = d[l+1];
        h = g - d[l];
        for (i = l+2; i < n; i++) { d[i] -= h; }
        f = f + h;
        /* Implicit QL transformation. */
        p = d[m]; c = 1.0; c2 = c; c3 = c;
        el1 = e[l+1]; s = 0.0; s2 = 0.0;
        for (i = m-1; i >= l; i--) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = hypot2(p, e[i]);
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * (c * g + s * d[i]);
          /* Accumulate transformation. */
          for (k = 0; k < n; k++) {
            h = V[k][i+1];
            V[k][i+1] = s * V[k][i] + c * h;
            V[k][i] = c * V[k][i] - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;

        /* Check for convergence. */
      } while (fabs(e[l]) > eps*tst1);
    }
    d[l] = d[l] + f;
    e[l] = 0.0;
  }

  /* Sort eigenvalues and corresponding vectors. */
  for (i = 0; i < n-1; i++) {
    k = i;
    p = d[i];
    for (j = i+1; j < n; j++) {
      if (d[j] < p) {
        k = j;
        p = d[j];
      }
    }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (j = 0; j < n; j++) {
        p = V[j][i];
        V[j][i] = V[j][k];
        V[j][k] = p;
      }
    }
  }
}


double MAX_1(double a, double b)
{
  return ((a)>(b)?(a):(b));
}

static double hypot2(double x, double y) { return sqrt(x*x+y*y); }

double absd(double val){ if(val>0){ return val;} else { return -val;} };









void MultipleNeuronTracer::RemoveSoma( LabelImageType3D::Pointer image )
{
  _SomaImage = image;
  RemoveIntraSomaNodes();

}

//////////////////////////////////////////////////////////////////////////////////////
//	INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////

void MultipleNeuronTracer::FeatureMain(void)
{
  time_t FeatureMain_start_time = clock();
  std::cout << std::endl<< "Feature detection 3D" << std::endl;
  _NDXImage = ImageType3D::New();
  _NDXImage->SetRegions(_PaddedCurvImage->GetBufferedRegion());
  _NDXImage->Allocate();
  _NDXImage->FillBuffer(0.0f);

  _NDXImage2 = CharImageType3D::New();///////////////////
  _NDXImage2->SetRegions(_PaddedCurvImage->GetBufferedRegion());//////////////////////
  _NDXImage2->Allocate();/////////////////////////
  _NDXImage2->FillBuffer(0.0f);///////////////////////////


  float sigmas[] =  { 2.0f, 2.8284f, 4.0f, 5.6569f, 8.0f, 11.31f};	//LoG scales
  for (unsigned int i = 0; i < 6; ++i)
  {
    std::cout << "Analysis at " << sigmas[i] << std::endl;
    if( _flagOutLog == false )
    {
      GetFeature( sigmas[i] );			//I guess this is finding all the critical points and throwing their vesselness values into _NDXImage
    }
    else
    {
      GetFeature_2( sigmas[i], i );			//I guess this is finding all the critical points and throwing their vesselness values into _NDXImage
    }
  }

  //itk::RescaleIntensityImageFilter<ImageType3D, ImageType3D>::Pointer rescaler = itk::RescaleIntensityImageFilter<ImageType3D, ImageType3D>::New();
  //rescaler->SetInput(_NDXImage);;
  //rescaler->SetOutputMaximum(255);
  //rescaler->SetOutputMinimum(0);

  //itk::CastImageFilter<ImageType3D, CharImageType3D>::Pointer caster = itk::CastImageFilter<ImageType3D, CharImageType3D>::New();
  //caster->SetInput(rescaler->GetOutput());

  //itk::ImageFileWriter<CharImageType3D>::Pointer writer = itk::ImageFileWriter<CharImageType3D>::New();
  //writer->SetInput(caster->GetOutput());
  //writer->SetFileName("C:\\Data\\Darpa\\TEST_FOR_PIPELINE\\23_2100_4200\\seed_points.tif");
  //writer->Update();

  /*ImageType3D::Pointer temp = ImageType3D::New();
    temp->SetRegions(_PaddedCurvImage->GetBufferedRegion());
    temp->Allocate();

    itk::ImageRegionConstIterator<ImageType3D> Nit(_NDXImage, _NDXImage->GetBufferedRegion());
    itk::ImageRegionIterator<ImageType3D> tit(temp, temp->GetBufferedRegion());
    itk::ImageRegionConstIterator<ImageType3D> Cit(_PaddedCurvImage, _PaddedCurvImage->GetBufferedRegion());
    for (Nit.GoToBegin(), Cit.GoToBegin(), tit.GoToBegin(); !Nit.IsAtEnd(); ++Nit, ++tit, ++Cit)
    {
    if (Nit.Get() > 0) 
    {
    tit.Set(1.0);		
    }
    }*/

  ////saving debris points in Debris_Points.tif
  //CharRescalerType::Pointer rescaler2 = CharRescalerType::New();////
  //rescaler2->SetInput( _NDXImage2 );/////
  //rescaler2->SetOutputMaximum( 255 );//////
  //rescaler2->SetOutputMinimum( 0 );/////
  //rescaler2->Update();/////

  //itk::CastImageFilter< CharImageType3D, CharImageType3D>::Pointer caster2 = itk::CastImageFilter< CharImageType3D, CharImageType3D>::New();///////
  //caster2->SetInput(rescaler2->GetOutput());/////

  //itk::ImageFileWriter< CharImageType3D >::Pointer seedsWriter2 = itk::ImageFileWriter< CharImageType3D >::New();//////
  //seedsWriter2->SetFileName("C:\\Users\\msavelon\\Desktop\\crop0131\\Debris_Points.tif");/////
  //seedsWriter2->SetInput(caster2->GetOutput());////
  //seedsWriter2->Update();/////

}

void MultipleNeuronTracer::GetFeature( float sigma ) 
{
  //std::ofstream out_seeds;////////////
  //out_seeds.open("C:\\Users\\msavelon\\Desktop\\crop0131\\seeds.txt", std::ios::app);//////////////////

  std::cout<<std::endl<<"Get Feature 1";

  clock_t LoG_start_time = clock();
  typedef itk::LaplacianRecursiveGaussianImageFilter< ImageType3D , ImageType3D> GFilterType;
  GFilterType::Pointer gauss = GFilterType::New();
  gauss->SetInput( _PaddedCurvImage );
  gauss->SetSigma( sigma );
  gauss->SetNormalizeAcrossScale(false);
  //ImageType3D::Pointer smoothedCurvImage = gauss->GetOutput();
  gauss->GetOutput()->Update();
  std::cout << "325Laplacian of Gaussian at " << sigma << " took " << (clock() - LoG_start_time)/(float) CLOCKS_PER_SEC << std::endl;


  //itk::RescaleIntensityImageFilter<ImageType3D, ImageType3D>::Pointer rescaler = itk::RescaleIntensityImageFilter<ImageType3D, ImageType3D>::New();
  //rescaler->SetInput(gauss->GetOutput());;
  //rescaler->SetOutputMaximum(255);
  //rescaler->SetOutputMinimum(0);

  //itk::CastImageFilter<ImageType3D, CharImageType3D>::Pointer caster = itk::CastImageFilter<ImageType3D, CharImageType3D>::New();
  //caster->SetInput(rescaler->GetOutput());

  //std::stringstream ss;
  //ss << ceil(sigma);
  //itk::ImageFileWriter<CharImageType3D>::Pointer writer = itk::ImageFileWriter<CharImageType3D>::New();
  //writer->SetInput(caster->GetOutput());
  //writer->SetFileName("C:\\Data\\Darpa\\TEST_FOR_PIPELINE\\23_2100_4200\\LOG_" + ss.str() + ".tif");
  //writer->Update();

  float tot = 0.0f, num = 0.0f;
  itk::ImageRegionIterator<ImageType3D> ittemp(gauss->GetOutput(), gauss->GetOutput()->GetBufferedRegion());
  float gamma = 1.6f;
  float tnorm = vcl_pow(sigma,gamma);

  for(ittemp.GoToBegin(); !ittemp.IsAtEnd(); ++ittemp)
  {
    float q = ittemp.Get()*tnorm;
    ittemp.Set(-1.0f*q);
    tot += q*q;
    num ++;
  }




  //if(debug){
  //	itk::ImageFileWriter<ImageType3D>::Pointer writer = itk::ImageFileWriter<ImageType3D>::New();
  //	writer->SetInput(gauss->GetOutput());
  //	writer->SetFileName("D:\\Data\\FSData\\LOG_TEST_ROYSAM\\LOG_" + ss.str() + ".tif");
  //	writer->Update();
  //}


  //std::cout << "Scale "<< sigma << " had average Energy: " << tot <<std::endl;

  // set the diagonal terms in neighborhood iterator
  itk::Offset<3>
    xp =  {{2 ,  0 ,   0}},
       xn =  {{-2,  0,    0}},
       yp =  {{0,   2,   0}},
       yn =  {{0,  -2,    0}},
       zp =  {{0,   0,    2}},
       zn =  {{0,   0,   -2}};

  itk::Size<3> rad = {{1,1,1}};
  itk::NeighborhoodIterator<ImageType3D> nit(rad , gauss->GetOutput(), gauss->GetOutput()->GetBufferedRegion());
  itk::ImageRegionIterator<ImageType3D> it(gauss->GetOutput(), gauss->GetOutput()->GetBufferedRegion());

  unsigned int
    xy1 =  17, //{ 1 ,   1 ,  0 },
        xy2 =  9,  //{ -1,  -1 ,  0 },
        xy3 =  15, //{ -1,   1 ,  0 },
        xy4 =  11, //{ 1 ,  -1 ,  0 },

        yz1 =  25, //{ 0 ,   1 ,  1 },
        yz2 =  1,  //{ 0 ,  -1 , -1 },
        yz3 =  19, //{ 0 ,  -1 ,  1 },
        yz4 =  7,  //{ 0 ,   1 , -1 },

        xz1 =  23, //{ 1 ,   0 ,  1 },
        xz2 =  3,  //{-1 ,   0 , -1 },
        xz3 =  21, //{-1 ,   0 ,  1 },
        xz4 =  5;  //{ 1 ,   0 , -1 };

  typedef itk::FixedArray< double, 3 > EigenValuesArrayType;
  typedef itk::Matrix< double, 3, 3 > EigenVectorMatrixType;
  typedef itk::SymmetricSecondRankTensor<double,3> TensorType;

  itk::Size<3> sz = _PaddedCurvImage->GetBufferedRegion().GetSize();
  sz[0] = sz[0] - 3;
  sz[1] = sz[1] - 3; 
  sz[2] = sz[2] - 3;

  it.GoToBegin();
  nit.GoToBegin();
  itk::Vector<float,3> sp = _PaddedCurvImage->GetSpacing();

  long win = long(sigma)/2;
  if (win < 2) 
  {
    win = 2;
  }

  //typedef itk::StatisticsImageFilter< ImageType3D > StatisticsImageFilterType;
  //StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New ();
  //statisticsImageFilter->SetInput(_PaddedCurvImage);
  //statisticsImageFilter->Update();
  //double image_mean = statisticsImageFilter->GetMean();
  //double image_stddev = statisticsImageFilter->GetSigma();

  //const float thresh1 = image_mean - (image_stddev/3);   // 3% of maximum theshold from Lowe 2004
  //const float thresh2 = image_mean/45;  // -0.1 percent of range
  const float thresh1 = intensity_threshold;// 0.01;//0.005;   // 3% of maximum theshold from Lowe 2004
  const float thresh2 = contrast_threshold;//0.003;//0.0003;  // -0.1 percent of range

  //std::cout << "Intensity thr: "  << intensity_threshold << " Contrast thr: " << contrast_threshold << std::endl;

  long ctCnt = 0;
  int inrt = 0; // niclas testing
  while(!nit.IsAtEnd()) 
  {
    // 		std::cout<<inrt++<<"\n "<<std::flush;
    itk::Index<3> ndx = it.GetIndex();
    if ( (ndx[0] < 2) || (ndx[1] < 2) || (ndx[2] < 2) ||
        (ndx[0] > (unsigned int)sz[0]) || (ndx[1] > (unsigned int)sz[1]) ||
        (ndx[2] > (unsigned int)sz[2]) )
    {
      ++it;
      ++nit;
      continue;
    }

    float a1 = 0.0;
    for (unsigned int i=0; i < 13; ++i)
    {
      a1 += vnl_math_max(nit.GetPixel(i), nit.GetPixel(26 - i));
    }

    float val = nit.GetPixel(13) ;

    //float mask_value = _MaskedImage->GetPixel(ndx);

    if ( ((val - a1/13.0f) > thresh2 ) && ( val > thresh1 ))  
      //if(mask_value > 0)
    {
      TensorType h;
      h[0] = gauss->GetOutput()->GetPixel( ndx + xp ) + gauss->GetOutput()->GetPixel( ndx + xn ) - 2*nit.GetPixel( 13 );
      h[3] = gauss->GetOutput()->GetPixel( ndx + yp ) + gauss->GetOutput()->GetPixel( ndx + yn ) - 2*nit.GetPixel( 13 );
      h[5] = gauss->GetOutput()->GetPixel( ndx + zp ) + gauss->GetOutput()->GetPixel( ndx + zn ) - 2*nit.GetPixel( 13 );
      h[1] = nit.GetPixel(xy1) + nit.GetPixel(xy2) - nit.GetPixel(xy3) - nit.GetPixel(xy4);
      h[2] = nit.GetPixel(xz1) + nit.GetPixel(xz2) - nit.GetPixel(xz3) - nit.GetPixel(xz4);
      h[4] = nit.GetPixel(yz1) + nit.GetPixel(yz2) - nit.GetPixel(yz3) - nit.GetPixel(yz4);

      EigenValuesArrayType ev;
      EigenVectorMatrixType em;
      h.ComputeEigenAnalysis (ev, em);

      unsigned int w;

      bool IsPlateCheck;
      if (device)
        IsPlateCheck=(IsPlate(ev, w));
      else
        IsPlateCheck=(IsPlate_control(ev, w));

      if (IsPlateCheck) 
      {
        float value = vnl_math_abs(ev[0]) + vnl_math_abs(ev[1]) + vnl_math_abs(ev[2]) - vnl_math_abs(ev[w]);
        if (RegisterIndex(value, ndx, sz, win))	//RegisterIndex returns true if this value is the highest in the neighborhood, otherwise it will return false
        {
          _NDXImage->SetPixel(ndx,value);
          ctCnt++;			//CriTical Counter I guess
          //std::cout<<ctCnt<<" ";
        }
      }

      if (IsDebris(ev, w, val))
      {
        float value = vnl_math_abs(ev[0]) + vnl_math_abs(ev[1]) + vnl_math_abs(ev[2]) - vnl_math_abs(ev[w]);//////
        //if (IsDebris(ev, w, val))
        if (RegisterIndexDebris(value, ndx, sz, win)) 
        {
          _NDXImage2->SetPixel(ndx,255);

          //if(out_seeds.good()){
          //	out_seeds << ndx[0] << " " << ndx[1] << " " << ndx[2] << " " << std::endl;}

          itk::Index<3> global_index;
          global_index[0]=ndx[0]+_indxDice[0];
          global_index[1]=ndx[1]+_indxDice[1];
          global_index[2]=ndx[2]+_indxDice[2];

          DebrisNode *db = new DebrisNode(global_index);//storing GLOBAL INDEX // same as ndx when executing as standalone
          _DebrisNodeContainer.push_back(db);

        }

      }


    }
    ++it;
    ++nit;
  }
  //if(debug){
  //	itk::ImageFileWriter<ImageType3D>::Pointer writer2 = itk::ImageFileWriter<ImageType3D>::New();
  //	writer2->SetInput(_NDXImage);
  //	writer2->SetFileName("D:\\Data\\FSData\\LOG_TEST_ROYSAM\\NDX_IMAGE" + ss.str() + ".tif");
  //	writer2->Update();
  //}
  std::cout <<"asdfNumber of CTs at this stage: " << ctCnt <<std::endl<<std::flush;
  //out_seeds.close();
}


void MultipleNeuronTracer::setDiceSize( itk::Size<3> sizeDice )
{
  _sizeDice[0] = sizeDice[0];
  _sizeDice[1] = sizeDice[1];
  _sizeDice[2] = sizeDice[2];
}

void MultipleNeuronTracer::setDiceIndex( itk::Index<3> indxDice )
{
  _indxDice[0] = indxDice[0];
  _indxDice[1] = indxDice[1];
  _indxDice[2] = indxDice[2];

}

void MultipleNeuronTracer::setLogScale( ImageType3D::Pointer inputImageLoG, int scale )
{
  if( scale == 0 )
    _logScale_1 = inputImageLoG;
  else if( scale == 1 )
    _logScale_2 = inputImageLoG;
  else if( scale == 2 )
    _logScale_3 = inputImageLoG;
  else if( scale == 3 )
    _logScale_4 = inputImageLoG;
  else if( scale == 4 )
    _logScale_5 = inputImageLoG;
  else if( scale == 5 )
    _logScale_6 = inputImageLoG;

}


void MultipleNeuronTracer::GetFeature_2( float sigma, int scale ) 
{
  // // 	clock_t LoG_start_time = clock();
  // 	typedef itk::LaplacianRecursiveGaussianImageFilter< ImageType3D , ImageType3D> GFilterType;
  // 	GFilterType::Pointer gauss = GFilterType::New();
  // 	gauss->SetInput( _PaddedCurvImage );
  // 	gauss->SetSigma( sigma );
  // 	gauss->SetNormalizeAcrossScale(false);
  // 	//ImageType3D::Pointer smoothedCurvImage = gauss->GetOutput();
  // 	gauss->GetOutput()->Update();
  // 	std::cout << "Laplacian of Gaussian at " << sigma << " took " << (clock() - LoG_start_time)/(float) CLOCKS_PER_SEC << std::endl;
  std::cout<<std::endl<<"Get Feature 2";


  ImageType3D::Pointer gauss_3 = ImageType3D::New();
  ImageType3D::PointType originGaussLocal;
  originGaussLocal[0] = 0; 
  originGaussLocal[1] = 0;
  originGaussLocal[2] = 0;
  gauss_3->SetOrigin( originGaussLocal );
  ImageType3D::IndexType startGaussLocal;
  startGaussLocal[0] = 0;
  startGaussLocal[1] = 0;
  startGaussLocal[2] = 0;
  ImageType3D::RegionType regionGaussLocal;
  regionGaussLocal.SetSize ( _sizeDice  );
  regionGaussLocal.SetIndex( startGaussLocal );
  gauss_3->SetRegions( regionGaussLocal );
  gauss_3->Allocate();
  gauss_3->FillBuffer(0);
  // 			SetSpacing(_logScale_1->GetSpacing());
  gauss_3->Update();


  ImageType3D::RegionType region;
  region.SetSize(_sizeDice);
  region.SetIndex(_indxDice);



  itk::ImageRegionIterator<ImageType3D> itGauss_3(gauss_3, gauss_3->GetLargestPossibleRegion());
  int gauss_slice_size = _sizeDice[1] * _sizeDice[0];
  ImageType3D::PixelType * gauss_3_Array = gauss_3->GetBufferPointer();

  itk::Index<3> local_origin = _logScale_1->GetRequestedRegion().GetIndex();
  // 	#pragma omp critical
  std::cout << std::endl << "Start of the Block : " << local_origin << std::endl;
  itk::Index<3> local_offset;
  local_offset[0] = _indxDice[0] - local_origin[0];
  local_offset[1] = _indxDice[1] - local_origin[1];
  local_offset[2] = _indxDice[2] - local_origin[2];

  itk::Size<3> block_size = _logScale_1->GetRequestedRegion().GetSize();
  // 	#pragma omp critical
  std::cout << std::endl << "Size of the Block : " << block_size << std::endl;
  int block_slice_size = block_size[1] * block_size[0];

  ImageType3D::PixelType * logScale_Array;
  if( scale == 0 )
  {
    logScale_Array = _logScale_1->GetBufferPointer();

    // 		itk::ImageRegionIterator<ImageType3D> itPreLoG(_logScale_1, region);
    // 		for(itPreLoG.GoToBegin(); !itPreLoG.IsAtEnd(); ++itPreLoG,++itGauss_3)
    // 		{
    // 			itGauss_3.Set(itPreLoG.Get());
    // 		}
  }
  else if( scale == 1 )
  {
    logScale_Array = _logScale_2->GetBufferPointer();

    // 		itk::ImageRegionIterator<ImageType3D> itPreLoG(_logScale_2, region);
    // 		for(itPreLoG.GoToBegin(); !itPreLoG.IsAtEnd(); ++itPreLoG,++itGauss_3)
    // 		{
    // 			itGauss_3.Set(itPreLoG.Get());
    // 		}
  }
  else if( scale == 2 )
  {
    logScale_Array = _logScale_3->GetBufferPointer();

    // 		itk::ImageRegionIterator<ImageType3D> itPreLoG(_logScale_3, region);
    // 		for(itPreLoG.GoToBegin(); !itPreLoG.IsAtEnd(); ++itPreLoG,++itGauss_3)
    // 		{
    // 			itGauss_3.Set(itPreLoG.Get());
    // 		}
  }
  else if( scale == 3 )
  {
    logScale_Array = _logScale_4->GetBufferPointer();

    // 		itk::ImageRegionIterator<ImageType3D> itPreLoG(_logScale_4, region);
    // 		for(itPreLoG.GoToBegin(); !itPreLoG.IsAtEnd(); ++itPreLoG,++itGauss_3)
    // 		{
    // 			itGauss_3.Set(itPreLoG.Get());
    // 		}
  }
  else if( scale == 4 )
  {
    logScale_Array = _logScale_5->GetBufferPointer();

    // 		itk::ImageRegionIterator<ImageType3D> itPreLoG(_logScale_5, region);
    // 		for(itPreLoG.GoToBegin(); !itPreLoG.IsAtEnd(); ++itPreLoG,++itGauss_3)
    // 		{
    // 			itGauss_3.Set(itPreLoG.Get());
    // 		}
  }
  else if( scale == 5 )
  {
    logScale_Array = _logScale_6->GetBufferPointer();

    // 		itk::ImageRegionIterator<ImageType3D> itPreLoG(_logScale_6, region);
    // 		for(itPreLoG.GoToBegin(); !itPreLoG.IsAtEnd(); ++itPreLoG,++itGauss_3)
    // 		{
    // 			itGauss_3.Set(itPreLoG.Get());
    // 		}
  }

  for(int z=0; z<_sizeDice[2]; ++z)
  {
    for(int y=0; y<_sizeDice[1]; ++y)
    {
      for(int x=0; x<_sizeDice[0]; ++x)
      {
        gauss_3_Array[(z * gauss_slice_size) + (y * _sizeDice[0]) + (x)] = logScale_Array[((z+local_offset[2]) * block_slice_size) + ((y+local_offset[1]) * block_size[0]) + (x+local_offset[0])];
      }
    }
  }

  ImageType3D::Pointer gauss_2 = gauss_3;

  // 	ImageType3D::Pointer gauss_2;
  // 	if( scale == 0 )
  // 	{
  // 		#pragma omp critical
  // 		{
  // 		typedef itk::RegionOfInterestImageFilter< ImageType3D, ImageType3D > ROIFilterType;
  // 		ROIFilterType::Pointer ROIfilter = ROIFilterType::New();
  // 		ROIfilter->SetRegionOfInterest(region);
  // 		ROIfilter->SetInput(_logScale_1);
  // 		ROIfilter->Update();
  // 		typedef itk::ImageDuplicator< ImageType3D > DuplicatorType;
  // 		DuplicatorType::Pointer duplicator = DuplicatorType::New();
  // 		duplicator->SetInputImage(ROIfilter->GetOutput());
  // 		duplicator->Update();
  // 		gauss_2 = duplicator->GetOutput();
  // 		}
  // 	}
  // 	else if( scale == 1 )
  // 	{
  // 		#pragma omp critical
  // 		{
  // 		typedef itk::RegionOfInterestImageFilter< ImageType3D, ImageType3D > ROIFilterType;
  // 		ROIFilterType::Pointer ROIfilter = ROIFilterType::New();
  // 		ROIfilter->SetRegionOfInterest(region);
  // 		ROIfilter->SetInput(_logScale_2);
  // 		ROIfilter->Update();
  // 		typedef itk::ImageDuplicator< ImageType3D > DuplicatorType;
  // 		DuplicatorType::Pointer duplicator = DuplicatorType::New();
  // 		duplicator->SetInputImage(ROIfilter->GetOutput());
  // 		duplicator->Update();
  // 		gauss_2 = duplicator->GetOutput();
  // 		}
  // 	}
  // 	else if( scale == 2 )
  // 	{
  // 		#pragma omp critical
  // 		{
  // 		typedef itk::RegionOfInterestImageFilter< ImageType3D, ImageType3D > ROIFilterType;
  // 		ROIFilterType::Pointer ROIfilter = ROIFilterType::New();
  // 		ROIfilter->SetRegionOfInterest(region);
  // 		ROIfilter->SetInput(_logScale_3);
  // 		ROIfilter->Update();
  // 		typedef itk::ImageDuplicator< ImageType3D > DuplicatorType;
  // 		DuplicatorType::Pointer duplicator = DuplicatorType::New();
  // 		duplicator->SetInputImage(ROIfilter->GetOutput());
  // 		duplicator->Update();
  // 		gauss_2 = duplicator->GetOutput();
  // 		}
  // 	}
  // 	else if( scale == 3 )
  // 	{
  // 		#pragma omp critical
  // 		{
  // 		typedef itk::RegionOfInterestImageFilter< ImageType3D, ImageType3D > ROIFilterType;
  // 		ROIFilterType::Pointer ROIfilter = ROIFilterType::New();
  // 		ROIfilter->SetRegionOfInterest(region);
  // 		ROIfilter->SetInput(_logScale_4);
  // 		ROIfilter->Update();
  // 		typedef itk::ImageDuplicator< ImageType3D > DuplicatorType;
  // 		DuplicatorType::Pointer duplicator = DuplicatorType::New();
  // 		duplicator->SetInputImage(ROIfilter->GetOutput());
  // 		duplicator->Update();
  // 		gauss_2 = duplicator->GetOutput();
  // 		}
  // 	}
  // 	else if( scale == 4 )
  // 	{
  // 		#pragma omp critical
  // 		{
  // 		typedef itk::RegionOfInterestImageFilter< ImageType3D, ImageType3D > ROIFilterType;
  // 		ROIFilterType::Pointer ROIfilter = ROIFilterType::New();
  // 		ROIfilter->SetRegionOfInterest(region);
  // 		ROIfilter->SetInput(_logScale_5);
  // 		ROIfilter->Update();
  // 		typedef itk::ImageDuplicator< ImageType3D > DuplicatorType;
  // 		DuplicatorType::Pointer duplicator = DuplicatorType::New();
  // 		duplicator->SetInputImage(ROIfilter->GetOutput());
  // 		duplicator->Update();
  // 		gauss_2 = duplicator->GetOutput();
  // 		}
  // 	}
  // 	else if( scale == 5 )
  // 	{
  // 		#pragma omp critical
  // 		{
  // 		typedef itk::RegionOfInterestImageFilter< ImageType3D, ImageType3D > ROIFilterType;
  // 		ROIFilterType::Pointer ROIfilter = ROIFilterType::New();
  // 		ROIfilter->SetRegionOfInterest(region);
  // 		ROIfilter->SetInput(_logScale_6);
  // 		ROIfilter->Update();
  // 		typedef itk::ImageDuplicator< ImageType3D > DuplicatorType;
  // 		DuplicatorType::Pointer duplicator = DuplicatorType::New();
  // 		duplicator->SetInputImage(ROIfilter->GetOutput());
  // 		duplicator->Update();
  // 		gauss_2 = duplicator->GetOutput();
  // 		}
  // 	}


  // #pragma omp critical
  // 	ROIfilter->Update();
  // 	typedef itk::ImageDuplicator< ImageType3D > DuplicatorType;
  // 	DuplicatorType::Pointer duplicator = DuplicatorType::New();
  // 	duplicator->SetInputImage(ROIfilter->GetOutput());
  // #pragma omp critical
  // 	duplicator->Update();
  // 	ImageType3D::Pointer gauss_2 = duplicator->GetOutput();



  //itk::RescaleIntensityImageFilter<ImageType3D, ImageType3D>::Pointer rescaler = itk::RescaleIntensityImageFilter<ImageType3D, ImageType3D>::New();
  //rescaler->SetInput(gauss->GetOutput());;
  //rescaler->SetOutputMaximum(255);
  //rescaler->SetOutputMinimum(0);

  //itk::CastImageFilter<ImageType3D, CharImageType3D>::Pointer caster = itk::CastImageFilter<ImageType3D, CharImageType3D>::New();
  //caster->SetInput(rescaler->GetOutput());

  //std::stringstream ss;
  //ss << ceil(sigma);
  //itk::ImageFileWriter<CharImageType3D>::Pointer writer = itk::ImageFileWriter<CharImageType3D>::New();
  //writer->SetInput(caster->GetOutput());
  //writer->SetFileName("C:\\Data\\Darpa\\TEST_FOR_PIPELINE\\23_2100_4200\\LOG_" + ss.str() + ".tif");
  //writer->Update();

  float tot = 0.0f, num = 0.0f;
  itk::ImageRegionIterator<ImageType3D> ittemp(gauss_2, gauss_2->GetBufferedRegion());
  float gamma = 1.6f;
  float tnorm = vcl_pow(sigma,gamma);

  for(ittemp.GoToBegin(); !ittemp.IsAtEnd(); ++ittemp)
  {
    float q = ittemp.Get()*tnorm;
    ittemp.Set(-1.0f*q);
    tot += q*q;
    num ++;
  }
  //std::cout << "Scale "<< sigma << " had average Energy: " << tot <<std::endl;

  // set the diagonal terms in neighborhood iterator
  itk::Offset<3>
    xp =  {{2 ,  0 ,   0}},
       xn =  {{-2,  0,    0}},
       yp =  {{0,   2,   0}},
       yn =  {{0,  -2,    0}},
       zp =  {{0,   0,    2}},
       zn =  {{0,   0,   -2}};

  itk::Size<3> rad = {{1,1,1}};
  itk::NeighborhoodIterator<ImageType3D> nit(rad , gauss_2, gauss_2->GetBufferedRegion());
  itk::ImageRegionIterator<ImageType3D> it(gauss_2, gauss_2->GetBufferedRegion());

  unsigned int
    xy1 =  17, //{ 1 ,   1 ,  0 },
        xy2 =  9,  //{ -1,  -1 ,  0 },
        xy3 =  15, //{ -1,   1 ,  0 },
        xy4 =  11, //{ 1 ,  -1 ,  0 },

        yz1 =  25, //{ 0 ,   1 ,  1 },
        yz2 =  1,  //{ 0 ,  -1 , -1 },
        yz3 =  19, //{ 0 ,  -1 ,  1 },
        yz4 =  7,  //{ 0 ,   1 , -1 },

        xz1 =  23, //{ 1 ,   0 ,  1 },
        xz2 =  3,  //{-1 ,   0 , -1 },
        xz3 =  21, //{-1 ,   0 ,  1 },
        xz4 =  5;  //{ 1 ,   0 , -1 };

  typedef itk::FixedArray< double, 3 > EigenValuesArrayType;
  typedef itk::Matrix< double, 3, 3 > EigenVectorMatrixType;
  typedef itk::SymmetricSecondRankTensor<double,3> TensorType;

  itk::Size<3> sz = _PaddedCurvImage->GetBufferedRegion().GetSize();
  sz[0] = sz[0] - 3;
  sz[1] = sz[1] - 3; 
  sz[2] = sz[2] - 3;

  it.GoToBegin();
  nit.GoToBegin();
  itk::Vector<float,3> sp = _PaddedCurvImage->GetSpacing();

  long win = long(sigma)/2;
  if (win < 2) 
  {
    win = 2;
  }

  //typedef itk::StatisticsImageFilter< ImageType3D > StatisticsImageFilterType;
  //StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New ();
  //statisticsImageFilter->SetInput(_PaddedCurvImage);
  //statisticsImageFilter->Update();
  //double image_mean = statisticsImageFilter->GetMean();
  //double image_stddev = statisticsImageFilter->GetSigma();

  //const float thresh1 = image_mean - (image_stddev/3);   // 3% of maximum theshold from Lowe 2004
  //const float thresh2 = image_mean/45;  // -0.1 percent of range
  const float thresh1 = intensity_threshold;// 0.01;//0.005;   // 3% of maximum theshold from Lowe 2004
  const float thresh2 = contrast_threshold;//0.003;//0.0003;  // -0.1 percent of range

  long ctCnt = 0;
  while(!nit.IsAtEnd()) 
  {
    itk::Index<3> ndx = it.GetIndex();
    if ( (ndx[0] < 2) || (ndx[1] < 2) || (ndx[2] < 2) ||
        (ndx[0] > (unsigned int)sz[0]) || (ndx[1] > (unsigned int)sz[1]) ||
        (ndx[2] > (unsigned int)sz[2]) )
    {
      ++it;
      ++nit;
      continue;
    }

    float a1 = 0.0;
    for (unsigned int i=0; i < 13; ++i)
    {
      a1 += vnl_math_max(nit.GetPixel(i), nit.GetPixel(26 - i));
    }

    float val = nit.GetPixel(13) ;

    if ( ((val - a1/13.0f) > thresh2 ) && ( val > thresh1 ))  
    {
      TensorType h;
      h[0] = gauss_2->GetPixel( ndx + xp ) + gauss_2->GetPixel( ndx + xn ) - 2*nit.GetPixel( 13 );
      h[3] = gauss_2->GetPixel( ndx + yp ) + gauss_2->GetPixel( ndx + yn ) - 2*nit.GetPixel( 13 );
      h[5] = gauss_2->GetPixel( ndx + zp ) + gauss_2->GetPixel( ndx + zn ) - 2*nit.GetPixel( 13 );
      h[1] = nit.GetPixel(xy1) + nit.GetPixel(xy2) - nit.GetPixel(xy3) - nit.GetPixel(xy4);
      h[2] = nit.GetPixel(xz1) + nit.GetPixel(xz2) - nit.GetPixel(xz3) - nit.GetPixel(xz4);
      h[4] = nit.GetPixel(yz1) + nit.GetPixel(yz2) - nit.GetPixel(yz3) - nit.GetPixel(yz4);

      EigenValuesArrayType ev;
      EigenVectorMatrixType em;
      h.ComputeEigenAnalysis (ev, em);

      unsigned int w;

      bool IsPlateCheck;
      if (device)
        IsPlateCheck=(IsPlate(ev, w));
      else
        IsPlateCheck=(IsPlate_control(ev, w));

      if (IsPlateCheck) 
      {
        float value = vnl_math_abs(ev[0]) + vnl_math_abs(ev[1]) + vnl_math_abs(ev[2]) - vnl_math_abs(ev[w]);
        if (RegisterIndex(value, ndx, sz, win))	//RegisterIndex returns true if this value is the highest in the neighborhood, otherwise it will return false
        {
          _NDXImage->SetPixel(ndx,value);
          ctCnt++;			//CriTical Counter I guess
        }
      }
    }
    ++it;
    ++nit;
  }
  std::cout <<"Number of CTs at this stage: " << ctCnt <<std::endl;
}

bool MultipleNeuronTracer::IsPlate(const itk::FixedArray<float, 3> &ev, unsigned int &w)  
{
  float L1, L2, L;
  if ( (ev[0] > ev[1]) && (ev[0] > ev[2]) ) 
  {
    w = 0;
    L = ev[0];
    L1 = ev[1]; 
    L2 = ev[2];
    if (ev[1] > ev[2])
    {
      L1 = ev[2]; 
      L2 = ev[1];		
    }
  }

  else if( (ev[1] > ev[0]) && (ev[1] > ev[2]) ) 
  {
    w = 1;
    L = ev[1];
    L1 = ev[0];
    L2 = ev[2];
    if (ev[0] > ev[2]) 
    {
      L1 = ev[2];
      L2 = ev[0];		
    }
  }

  else  
  {
    w = 2;
    L = ev[2];
    L1 = ev[0];
    L2 = ev[1];
    if (ev[0] > ev[1]) 
    {
      L1 = ev[1];
      L2 = ev[0];		
    }
  }

  if (std::abs(L2)/std::sqrt(std::abs(L1*L))<debris_threshold) //*( (abs(L2)/sqrt(abs(L1*L))<0.25) && (abs(L)+abs(L1)+abs(L2)>0.05) )*/(abs(L2)/sqrt(abs(L1*L))<0.5)//((L - L2) > (L2 - L1) && (L - L2) > vnl_math_abs(L)) 
  {
    return true;
  }

  return false;  /// right now this is turned off (Amit)
}


bool MultipleNeuronTracer::IsPlate_control(const itk::FixedArray<float, 3> &ev, unsigned int &w)  
{
  float L1, L2, L;
  if ( (ev[0] > ev[1]) && (ev[0] > ev[2]) ) 
  {
    w = 0;
    L = ev[0];
    L1 = ev[1]; 
    L2 = ev[2];
    if (ev[1] > ev[2])
    {
      L1 = ev[2]; 
      L2 = ev[1];		
    }
  }

  else if( (ev[1] > ev[0]) && (ev[1] > ev[2]) ) 
  {
    w = 1;
    L = ev[1];
    L1 = ev[0];
    L2 = ev[2];
    if (ev[0] > ev[2]) 
    {
      L1 = ev[2];
      L2 = ev[0];		
    }
  }

  else  
  {
    w = 2;
    L = ev[2];
    L1 = ev[0];
    L2 = ev[1];
    if (ev[0] > ev[1]) 
    {
      L1 = ev[1];
      L2 = ev[0];		
    }
  }

  // 	if /*( (abs(L2)/sqrt(abs(L1*L))<0.25) && (abs(L)+abs(L1)+abs(L2)>0.05) )*//*(abs(L2)/sqrt(abs(L1*L))<0.5)*/((L - L2) > (L2 - L1) && (L - L2) > vnl_math_abs(L)) 

  if ( (std::abs(L2)/std::sqrt(std::abs(L1*L))<0.25) && (std::abs(L)+std::abs(L1)+std::abs(L2)>0.05) )
  {
    return true;
  }

  return true;  /// right now this is turned off (Amit)
}


//Searches in some specified window in _NDXImage around ndx and see if something larger than value is present
//If there is a value in the neighborhood larger than value then return false, else set the _NDXImage at ndx to 0 and return true 
bool MultipleNeuronTracer::IsDebris(const itk::FixedArray<float, 3> &ev, unsigned int &w, float val)  
{
  float L1, L2, L;
  if ( (ev[0] > ev[1]) && (ev[0] > ev[2]) ) 
  {
    w = 0;
    L = ev[0];
    L1 = ev[1]; 
    L2 = ev[2];
    if (ev[1] > ev[2])
    {
      L1 = ev[2]; 
      L2 = ev[1];		
    }
  }

  else if( (ev[1] > ev[0]) && (ev[1] > ev[2]) ) 
  {
    w = 1;
    L = ev[1];
    L1 = ev[0];
    L2 = ev[2];
    if (ev[0] > ev[2]) 
    {
      L1 = ev[2];
      L2 = ev[0];		
    }
  }

  else  
  {
    w = 2;
    L = ev[2];
    L1 = ev[0];
    L2 = ev[1];
    if (ev[0] > ev[1]) 
    {
      L1 = ev[1];
      L2 = ev[0];		
    }
  }

  if (( std::abs(L-L1)/std::abs(L2)<0.7 )&& (std::abs(L)+std::abs(L1)+std::abs(L2)>0.002) && (val>0.006) ) // 0.7 0.005 0.09
  {
    return true;
  }

  return false;//true;  /// right now this is turned off (Amit)
}


bool MultipleNeuronTracer::RegisterIndex(const float value, itk::Index<3> &ndx, itk::Size<3>& sz, long h = 2) 
{
  itk::Index<3> n;
  bool higherPresent = false;
  for (n[0] = ndx[0]-h; n[0] <= ndx[0]+h; ++n[0]) 
  {
    for (n[1] = ndx[1]-h; n[1] <= ndx[1]+h; ++n[1]) 
    {
      for (n[2] = ndx[2]-h; n[2] <= ndx[2]+h; ++n[2]) 
      {
        if ( (n[0] < 2) || (n[1] < 2) || (n[2] < 2) || (n[0] > (unsigned int)sz[0]) ||
            (n[1] > (unsigned int)sz[1]) || (n[2] > (unsigned int)sz[2]) )
        {
          continue;
        }

        float curval = _NDXImage->GetPixel(n);
        if (value > curval) 
        {
          _NDXImage->SetPixel(n,0.0f);	//Why do we set this to 0.0? We overwrite with value later anyways...
        }
        else if (value < curval) 
        {
          higherPresent = true;				
        }
      }
    }
  }
  if (higherPresent == true) 
  {
    return false;	
  }

  return true;
}

bool MultipleNeuronTracer::RegisterIndexDebris(const float value, itk::Index<3> &ndx, itk::Size<3>& sz, long h = 2) 
{
  itk::Index<3> n;
  bool higherPresent = false;
  for (n[0] = ndx[0]-h; n[0] <= ndx[0]+h; ++n[0]) 
  {
    for (n[1] = ndx[1]-h; n[1] <= ndx[1]+h; ++n[1]) 
    {
      for (n[2] = ndx[2]-h; n[2] <= ndx[2]+h; ++n[2]) 
      {
        if ( (n[0] < 2) || (n[1] < 2) || (n[2] < 2) || (n[0] > (unsigned int)sz[0]) ||
            (n[1] > (unsigned int)sz[1]) || (n[2] > (unsigned int)sz[2]) )
        {
          continue;
        }

        float curval = (float)(_NDXImage2->GetPixel(n));
        if (value > curval) 
        {
          _NDXImage2->SetPixel(n,0);	//Why do we set this to 0.0? We overwrite with value later anyways...
        }
        else if (value < curval) 
        {
          higherPresent = true;				
        }
      }
    }
  }
  if (higherPresent == true) 
  {
    return false;	
  }

  return true;
}


/////////////////////////////////////////////////////////////////////////////////////////////

SWCNode* MultipleNeuronTracer::TBack(itk::Index<3> &ndx, std::vector<IndexType>& Chain)  
{

  SWCNode* Label = NULL;
  itk::Index<3> n;
  itk::Vector<float,3> p, x, d, dold;
  for (int i=0; i<3; i++) 
  {
    p[i] = static_cast<PixelType>(ndx[i]);
    dold[i] = 0.0f;
  }
  bool done = false;
  if (_SWCImage->GetPixel(ndx)->TreeID > 0) 
  {
    done = true;
  }
  const float MAXDERV = 10000.0f;

  Chain.push_back(ndx);

  while (done == false) 
  {
    //x
    x = p; 
    x[0]++;
    n.CopyWithRound(x);
    if (n[0] < (unsigned int)_size[0])
    {
      d[0] = _ConnImage->GetPixel(n);
    }
    else
    {
      d[0] = MAXDERV;
    }

    x = p; 
    x[0]--;
    n.CopyWithRound(x);
    if (n[0] >= 0)    
    {
      d[0] -= _ConnImage->GetPixel(n);   
    }
    else 
    {
      d[0] -= MAXDERV; 
    }

    // y
    x = p; 
    x[1]++;
    n.CopyWithRound(x);
    if (n[1] < (unsigned int)_size[1]) 
    {
      d[1] = _ConnImage->GetPixel(n);
    }
    else
    {
      d[1] = MAXDERV;
    }

    x = p; 
    x[1]--;
    n.CopyWithRound(x);
    if (n[1] >= 0)
    {
      d[1] -= _ConnImage->GetPixel(n);
    }
    else
    {
      d[1] -= MAXDERV; 
    }

    // z
    x = p; 
    x[2]++;
    n.CopyWithRound(x);
    if (n[2] < (unsigned int)_size[2]) 
    {
      d[2] = _ConnImage->GetPixel(n); 
    }
    else
    {
      d[2] = MAXDERV;
    }

    x = p; 
    x[2]--;
    n.CopyWithRound(x);
    if (n[2] >= 0)
    {
      d[2] -= _ConnImage->GetPixel(n);
    }
    else
    {
      d[2] -= MAXDERV;
    }

    double norm2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
    if (norm2>0.001)
    {
      d.Normalize();
      d += dold;
      d.Normalize();
      dold = d;
      d *= 0.5;
      p -= d;
    }

    n.CopyWithRound(p);
    Chain.push_back(n);
    //check termination
    SWCNode *t = _SWCImage->GetPixel(n);
    if (t != NULL ) 
    {
      if (t->TreeID > 0) 
      {
        done = true;
        Label = _SWCImage->GetPixel(n);
        break;
      }
    }
    if (Chain.size() > 500) 
    {
      //std::cout << "Tree not found for " << ndx << " in 500 steps, exiting!! " << std::endl;
      Chain.clear();
      Label = NULL;
      done = true;
      break;
    }
  }
  //

  float costFactorLabel = 0.0;

  if (Chain.size()!=0 )
    costFactorLabel = GetCostLocalLabel( Label , ndx);

  if (costFactorLabel>=10.0)
  {	
    Chain.clear();
    done=true;
    return NULL;}
  else 
    return Label;
  //


  //return Label;
}

///////////////////////////////////////////////////////////////////////////////////
float MultipleNeuronTracer::GetCost(SWCNode* s, itk::Index<3>& endx ) 
{
  itk::Index<3> base = endx, ndx = s->ndx;
  float cost = 0.0f, angsum = 0.0f, count = 0.01f;
  itk::Vector<float,3> d1, d2 , gd, gd1;
  d1.Filled(0.0);
  gd1.Fill(0.0);
  bool first = true;

  while (count < 500.0f) 
  {
    float d = (ndx[0] - base[0])*(ndx[0] - base[0]) + (ndx[1] - base[1])*(ndx[1] - base[1]) + (ndx[2] - base[2])*(ndx[2] - base[2]) ;
    if ( vcl_sqrt(d) > 6.0f) 
    {
      d2 = d1;
      d1[0] = float(ndx[0] - base[0]);
      d1[1] = float(ndx[1] - base[1]);
      d1[2] = float(ndx[2] - base[2]);
      d1.Normalize();
      if (first == true) 
      {
        first = false;
        gd1 = d1;
      }
      else 
      {
        PixelType w = dot_product(d1.Get_vnl_vector(),d2.Get_vnl_vector());
        if (w < 0.99f) 
        {
          angsum += vcl_acos(vnl_math_abs(w));
        }
        count ++;
      }
      base = ndx;
    }
    s = s->parent;
    if (s == NULL) 
    {
      break;
    }
    ndx = s->ndx;
  }
  gd[0] = float(ndx[0] - endx[0]);
  gd[1] = float(ndx[1] - endx[1]);
  gd[2] = float(ndx[2] - endx[2]);
  gd.Normalize();
  float allowedTurns = 1.0f;
  if (dot_product(gd.Get_vnl_vector(),gd1.Get_vnl_vector() ) >= 0) 
  {
    cost = (angsum/allowedTurns);
    if ( cost > 1.0) 
    {
      cost = 1.0f;		
    }
  }
  else 
  {
    cost = 1.0f;
  }

  return cost;
}

float MultipleNeuronTracer::GetCostLocal(SWCNode* s, itk::Index<3>& endx ) 
{
  itk::Index<3> base = endx, ndx = s->ndx;
  float cost = 0.0f, count = 0.01f;
  itk::Vector<float,3> d1, d2;
  d2.Filled(0.0);

  d1[0] = float(ndx[0] - base[0]);
  d1[1] = float(ndx[1] - base[1]);
  d1[2] = float(ndx[2] - base[2]);
  d1.Normalize();

  base = ndx;

  while (count < 500.0f) 
  {
    float d = (ndx[0] - base[0])*(ndx[0] - base[0]) + (ndx[1] - base[1])*(ndx[1] - base[1]) + (ndx[2] - base[2])*(ndx[2] - base[2]) ;
    if ( vcl_sqrt(d) > 6.0f) 
    {
      d2[0] = float(ndx[0] - base[0]);
      d2[1] = float(ndx[1] - base[1]);
      d2[2] = float(ndx[2] - base[2]);
      d2.Normalize();

      PixelType w = dot_product(d1.Get_vnl_vector(),d2.Get_vnl_vector());
      if ( w <= 0.0f) 
      {
        cost = 1.0f;
      }
      else if (( w > 0.0f) && (w <= 0.98f)) 
      {
        cost = 1.0 - w;
      }
      else 
      {
        cost = 0.0f;
      }
      break;
    }
    count++;
    s = s->parent;
    if (s == NULL) 
    {
      break;
    }
    ndx = s->ndx;
  }

  return cost;
}



float MultipleNeuronTracer::GetCostLocalLabel(SWCNode* s, itk::Index<3>& endx ) //this functions is used in TBack to prevent nodes to join trees, in cases of abrupt directional changes 
{
  itk::Index<3> base = endx, ndx = s->ndx;
  float cost = 0.0f, count = 0.01f, local_count=0.01f;
  itk::Vector<float,3> d1, d2;
  d2.Filled(0.0);

  d1[0] = float(ndx[0] - base[0]); //  leaf-current
  d1[1] = float(ndx[1] - base[1]);
  d1[2] = float(ndx[2] - base[2]);
  d1.Normalize();

  base = ndx;

  float local_abrupt=0.0;//flag for marking if there is an abrupt directional change as we cross the chain towards close ancestors

  while (count < 500.0f) //500 (5)
  {
    float d = (ndx[0] - base[0])*(ndx[0] - base[0]) + (ndx[1] - base[1])*(ndx[1] - base[1]) + (ndx[2] - base[2])*(ndx[2] - base[2]) ;
    if ( vcl_sqrt(d) > 6.0f) //6.0 (0)
    {
      d2[0] = float(ndx[0] - base[0]); //  ancestor-leaf
      d2[1] = float(ndx[1] - base[1]);
      d2[2] = float(ndx[2] - base[2]);
      d2.Normalize();

      PixelType w = dot_product(d1.Get_vnl_vector(),d2.Get_vnl_vector());

      if ( w <= 0.4f) //0.0 (0.2)
      {
        cost = 1.0f; //1.0 (10.0)
        local_abrupt=1.0;
      }
      else if (( w > 0.4f) && (w <= 0.98f))//0.0 && 0.98  (0.2-0.98)
      {
        cost = 1.0 - w;//-
      }
      else 
      {
        cost = 0.0f;//0
      }
      break;
    }
    count++;
    s = s->parent;
    if (s == NULL) 
    {
      break;
    }
    ndx = s->ndx;
  }

  if (local_abrupt==1.0)
    cost=10.0;

  return cost;


}

void MultipleNeuronTracer::ScanNeighbors( PixelType &a1, PixelType &a2, PixelType &a3, itk::Index<3> &ndx) 
{
  a1 = MAXVAL;
  if(ndx[0] > 0)
  {
    a1 = _ConnImage->GetPixel(ndx + _off.at(0));
  }	
  if (ndx[0] < (unsigned int)_size[0]-1) 
  {
    a1 = vnl_math_min(_ConnImage->GetPixel(ndx + _off.at(1)), a1 );
  }

  a2 = MAXVAL;
  if(ndx[1] > 0)  
  {
    a2 = _ConnImage->GetPixel(ndx + _off.at(2));
  }
  if (ndx[1] < (unsigned int)_size[1]-1) 
  {
    a2 = vnl_math_min(_ConnImage->GetPixel(ndx + _off.at(3)), a2 );
  }

  a3 = MAXVAL;
  if(ndx[2] > 0)  
  {
    a3 = _ConnImage->GetPixel(ndx + _off.at(4));
  }
  if (ndx[2] < (unsigned int)_size[2]-1) 
  {
    a3 = vnl_math_min(_ConnImage->GetPixel(ndx + _off.at(5)), a3 );
  }
}

///////////////////////////////////////////////////////////////////////
PixelType MultipleNeuronTracer::Update( PixelType a1,  PixelType a2,  PixelType a3,  PixelType P )  
{
  if (a1 > a2)  
  {
    PixelType temp = a2;
    a2 = a1;
    a1 = temp;
  }

  if (a2 > a3)  
  {
    PixelType temp = a3;
    a3 = a2;
    a2 = temp;
  }

  if (a1 > a2)  
  {
    PixelType temp = a2;
    a2 = a1;
    a1 = temp;
  }

  PixelType A1 = 0;
  PixelType delta = (a2+a1+a3)*(a2+a1+a3) - 3*(a1*a1 + a2*a2 + a3*a3 - P*P);

  if( delta>=0 ) 
  {
    A1 = ( a2+a1+a3 + vcl_sqrt(delta) ) / 3.0;
  }

  if( A1 <= a3 )  
  {
    delta = (a2+a1)*(a2+a1) - 2*(a1*a1 + a2*a2 - P*P);
    A1 = 0;
    if( delta>=0 )  
    {
      A1 = ( a2+a1 + vcl_sqrt(delta) ) / 2.0;
    }
    if( A1 <= a2 ) 
    {
      A1 = a1 + P;
    }
  }
  return A1;
}

///////////////////////////////////////////////////////////////////////////////////

void MultipleNeuronTracer::Decimate() 
{
  //std::cout << "Decimating the tree of size: " << _SWCNodeContainer.size() << std::endl;
  std::vector<SWCNode*>::iterator sit;
  for (sit = _SWCNodeContainer.begin(); sit != _SWCNodeContainer.end(); ++sit) 
  {
    if((*sit)->children.size() >= 2) 
    {
      (*sit)->IsBranch = true;
      (*sit)->IsActive = true;
    }
    else if((*sit)->children.size() == 0) 
    {
      (*sit)->IsLeaf = true;
      (*sit)->IsActive = true;
    }
    else if ((*sit)->parent == NULL) 
    {
      (*sit)->IsActive = true;
    }
    else 
    {
      (*sit)->IsActive = false;
    }
  }

  //std::cout << "Tree labeled: 1" << std::endl;
  for (sit = _SWCNodeContainer.begin(); sit != _SWCNodeContainer.end(); ++sit) 
  {
    if ((*sit)->IsActive == false) 
    {
      if ((*sit)->parent->IsActive == false)  
      {
        bool chActive = false;
        for (unsigned int i = 0; i < (*sit)->children.size(); ++i) 
        {
          if ((*sit)->children[i]->IsActive == true) 
          {
            chActive = true;					
          }
        }
        if (chActive == false) 
        {
          (*sit)->IsActive = true;				
        }
      }
    }
  }

  const float minOffshootLength = offshoot;//10;
  //std::cout << "Removing offshoots of length less than " << minOffshootLength  << std::endl;

  for (sit = _SWCNodeContainer.begin(); sit != _SWCNodeContainer.end(); ++sit) 
  {
    if ((*sit)->IsLeaf == true) 
    {
      //std::cout << "Leaf at" << (*sit)->ndx << " ID:" << (*sit)->ID << " ParentID:" << (*sit)->PID << std::endl;
      SWCNode* par = (*sit)->parent;
      if (par == NULL) 
      {
        continue;
      }

      if (par->PID == -1) 
      {
        continue;
      }

      itk::Vector<PixelType,3> p1 = (*sit)->pos;
      itk::Vector<PixelType,3> p2 = par->pos;
      itk::Vector<PixelType,3> dp = p1 - p2;
      float d = dp.GetNorm();

      while ( par->IsBranch == false ) 
      {
        p1 = p2;
        par = par->parent;
        if (par == NULL) 
        {
          break;
        }
        p2 = par->pos;
        dp = p1 - p2;
        d += dp.GetNorm();
      }

      if (d < minOffshootLength) 
      {
        SWCNode* n = (*sit);
        while ( n != par ) 
        {
          n->IsActive = false;
          n = n->parent;
        }
        if(par != NULL)
        {
          par->IsBranch = false;
        }
      }
    }
  }

  //std::cout << "Tree labeled: 3" << std::endl;

  std::vector<SWCNode*> NewContainer;
  NewContainer.reserve(_SWCNodeContainer.size());

  long newID = 1;
  itk::Array<long> IDLookUp(_SWCNodeContainer.size());
  IDLookUp.Fill(0);

  for (unsigned int i=0; i < _SWCNodeContainer.size(); ++i) 
  {
    if (_SWCNodeContainer[i]->IsActive == true) 
    {
      IDLookUp[i] = newID++;		
    }
  }
  //std::cout << "Lookup generated: " << std::endl;

  for (sit = _SWCNodeContainer.begin(); sit != _SWCNodeContainer.end(); ++sit) 
  {
    if ((*sit)->IsActive == true) 
    {
      long PID;
      itk::Index<3> ndx;
      for (int i = 0; i < 3; ++i) 
      {
        ndx[i] = long((*sit)->pos[i]);
      }

      SWCNode* par = (*sit)->parent;
      if (par == NULL) 
      {
        PID = -1;
      }
      else 
      {
        while (par->IsActive == false)
        {
          par = par->parent;
          if(par == NULL)
          {
            break;
          }
        }
        if(par == NULL)
        {
          PID = -1;
        }

        else
        {
          PID = IDLookUp(par->ID - 1);
          if(PID < 1 || (unsigned int)PID > NewContainer.size())
          {
            continue;
          }
          par = NewContainer[PID-1];
        }
      }

      SWCNode* s = new SWCNode();
      s->ID = IDLookUp((*sit)->ID - 1);
      s->PID = PID;
      s->IsActive = true;
      s->IsBranch = (*sit)->IsBranch;
      s->IsLeaf = (*sit)->IsLeaf;
      s->ndx = ndx;
      s->parent = par;
      s->pos = (*sit)->pos;
      s->TreeID = (*sit)->TreeID;
      if (par != NULL) 
      {
        par->children.push_back(s);
      }		
      NewContainer.push_back(s);
    }
  }
  //std::cout << "NewContainer created: " << std::endl;

  for (sit = _SWCNodeContainer.begin(); sit != _SWCNodeContainer.end(); ++sit) 
  {
    delete (*sit);
  }

  _SWCNodeContainer = NewContainer;
}

///////////////////////////////////////////////////////////////////////////////////
void MultipleNeuronTracer::Interpolate(float sigma) 
{
  //std::cout << "Interpolating the tree: " << std::endl;
  typedef itk::SmoothingRecursiveGaussianImageFilter< ImageType3D , ImageType3D> GFilterType;
  GFilterType::Pointer gauss = GFilterType::New();
  gauss->SetInput( _PaddedCurvImage );
  gauss->SetSigma( sigma );
  gauss->SetNormalizeAcrossScale(false);
  //ImageType3D::Pointer smoothedCurvImage = gauss->GetOutput();
  gauss->GetOutput()->Update();

  std::vector<SWCNode*>::iterator sit;
  for (sit = _SWCNodeContainer.begin(); sit != _SWCNodeContainer.end(); ++sit) 
  {
    float w,x,y,z;
    if (((*sit)->children.size() > 0) && ((*sit)->parent != NULL)) 
    {
      w = vnl_math_max(gauss->GetOutput()->GetPixel((*sit)->ndx), 0.1f);
      x = w * float((*sit)->ndx[0]);
      y = w * float((*sit)->ndx[1]);
      z = w * float((*sit)->ndx[2]);

      if ((*sit)->parent != NULL) 
      {
        float w1 = vnl_math_max(gauss->GetOutput()->GetPixel((*sit)->parent->ndx), 0.1f);
        w += w1;
        x += (w1 * float((*sit)->parent->ndx[0]));
        y += (w1 * float((*sit)->parent->ndx[1]));
        z += (w1 * float((*sit)->parent->ndx[2]));
      }
      for (unsigned int i = 0; i < (*sit)->children.size() ; ++i) 
      {
        float w1 = vnl_math_max(gauss->GetOutput()->GetPixel((*sit)->children[i]->ndx), 0.1f);
        w += w1;
        x += (w1 * float((*sit)->children[i]->ndx[0]));
        y += (w1 * float((*sit)->children[i]->ndx[1]));
        z += (w1 * float((*sit)->children[i]->ndx[2]));
      }
      (*sit)->pos[0] = x / w;
      (*sit)->pos[1] = y / w;
      (*sit)->pos[2] = z / w;
      (*sit)->ndx.CopyWithRound((*sit)->pos);
    }
    else 
    {
      (*sit)->pos[0] = float((*sit)->ndx[0]);
      (*sit)->pos[1] = float((*sit)->ndx[1]);
      (*sit)->pos[2] = float((*sit)->ndx[2]);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////
void MultipleNeuronTracer::LoadSomaImage(std::string somaFileName)
{
  /*typedef itk::ImageFileReader<CharImageType3D> SomaReaderType;
    SomaReaderType::Pointer somaReader = SomaReaderType::New();
    somaReader->SetFileName(somaFileName);
    _SomaImage = somaReader->GetOutput();
    somaReader->Update();
    */

  // Now reading a labeled image for the somas

  typedef itk::ImageFileReader<LabelImageType3D> SomaReaderType;
  SomaReaderType::Pointer somaReader = SomaReaderType::New();
  somaReader->SetFileName(somaFileName);
  _SomaImage = somaReader->GetOutput();
  somaReader->Update();
}

void MultipleNeuronTracer::RemoveIntraSomaNodes(void)
{
  std::cout << "Removing nodes that fall inside the somas of the Curvelets Image" << std::endl;

  unsigned int originalSize = _SWCNodeContainer.size();
  LabelArrayType somaArray = _SomaImage->GetBufferPointer();
  itk::Size<3> im_size = _SomaImage->GetBufferedRegion().GetSize();
  int slice_size = im_size[0] * im_size[1];
  int row_size = im_size[0];

  //find the root nodes of each tree
  std::cout << "Finding the root nodes of each tree" << std::endl;
  std::map<long, SWCNode*> treeIDToRootMap;
  std::vector<SWCNode*>::iterator sit;
  for (sit = _SWCNodeContainer.begin(); sit != _SWCNodeContainer.end(); ++sit)
  {
    //assume that a node with no parent is a centroid
    if( (*sit)->parent == NULL )
    {
      treeIDToRootMap[(*sit)->TreeID] = (*sit);
    }
  }

  if(treeIDToRootMap.size() != this->_StartPoints.size()){
    std::cout << "Centroids missing!!" << std::endl;

    /*std::ofstream cent_out_1("cent1.txt");
      std::ofstream cent_out_2("cent2.txt");
      for(int i = 0; i < this->StartPoints.size(); i++)
      cent_out_1 << this->StartPoints[i][0] << "," << this->StartPoints[i][1] << "," << this->StartPoints[i][2] << std::endl;
      for(int i = 0; i < treeIDToRootMap.size(); i++)
      cent_out_2 << treeIDToRootMap[i]->ndx[0] << "," << treeIDToRootMap[i]->ndx[1] << "," << treeIDToRootMap[i]->ndx[2] << std::endl;
      cent_out_1.close();
      cent_out_2.close();*/
  }

  itk::Index<3> dummy_index;
  //Removing nodes
  for (sit = _SWCNodeContainer.begin(); sit != _SWCNodeContainer.end();)
  {
    //don't check nodes that are outside the extent of the soma image
    if ( !_SomaImage->GetLargestPossibleRegion().IsInside( (*sit)->ndx ) )
    {
      ++sit;
      continue;
    }

    //don't remove centroid nodes
    if( (*sit)->parent == NULL )
    {
      ++sit;
      continue;
    }

    //remove any other node that falls within a soma
    /*if ( _SomaImage->GetPixel( (*sit)->ndx ) != 0 )
      {
      delete (*sit);
      sit = _SWCNodeContainer.erase(sit);
      }*/

    // Removing nodes only lying in the foreground of the current soma
    itk::Index<3> Node_0 = (*sit)->ndx;
    itk::Index<3> Node_1 = treeIDToRootMap[(*sit)->TreeID]->ndx;
    //if ( somaArray[(slice_size * Node_0[2]) + (row_size * Node_0[1]) + Node_0[0]] != 0 )
    if ( _SomaImage->GetPixel( (*sit)->ndx ) != 0 )
    {
      //if( somaArray[(slice_size * Node_0[2]) + (row_size * Node_0[1]) + Node_0[0]] == somaArray[(slice_size * Node_1[2]) + (row_size * Node_1[1]) + Node_1[0]])
      if( _SomaImage->GetPixel((*sit)->ndx) == _SomaImage->GetPixel(treeIDToRootMap[(*sit)->TreeID]->ndx) )
      {
        for(int i = 0; i < this->_StartPoints.size(); i++)
        {
          if((*sit)->ndx[0] == this->_StartPoints[i][0] && (*sit)->ndx[1] == this->_StartPoints[i][1] && (*sit)->ndx[2] == this->_StartPoints[i][2])
            std::cout << "Centroid " << (*sit)->ndx[0] << ", " << (*sit)->ndx[1] << ", " << (*sit)->ndx[2] << " deleted!! " <<std::endl;
        }

        delete (*sit);
        sit = _SWCNodeContainer.erase(sit);
        //std::cout << "Deleted node. " << std::endl;				
      }	

      else
      {
        SWCNode *parent = (*sit)->parent;
        SWCNode *root = treeIDToRootMap[(*sit)->TreeID];

        if(parent->ndx == root->ndx)
        {
          ++sit;
          continue;
        }					
        if(_SomaImage->GetPixel(parent->ndx) == _SomaImage->GetPixel(root->ndx))
        {
          (*sit)->parent = root;
          (*sit)->PID = root->ID;

          ++sit;
        }
        else
        {
          ++sit;
          continue;
        }		
      }
    }

    //otherwise if its parent lies within a soma reassign it to be a child
    //of the centroid instead.
    else
    {
      SWCNode *parent = (*sit)->parent;
      if ( !_SomaImage->GetLargestPossibleRegion().IsInside( parent->ndx ) )
      {
        ++sit;
        continue;
      }			

      itk::Index<3> Node_2 = parent->ndx;
      //if( somaArray[(slice_size * Node_2[2]) + (row_size * Node_2[1]) + Node_2[0]] != 0)

      if( _SomaImage->GetPixel( parent->ndx ) != 0)
      {					
        if( _SomaImage->GetPixel(parent->ndx) == _SomaImage->GetPixel(treeIDToRootMap[(*sit)->TreeID]->ndx) )
        {
          (*sit)->parent = treeIDToRootMap[(*sit)->TreeID];
          (*sit)->PID = treeIDToRootMap[(*sit)->TreeID]->ID;
        }					
      }
      ++sit;
    }
  }

  size_t newSize = _SWCNodeContainer.size();
  std::cout << "Just removed " << originalSize - newSize
    << " nodes (" << originalSize << " to " << newSize << ")"
    << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////

void MultipleNeuronTracer::WriteMultipleSWCFiles(std::string fname, unsigned int padz) 
{
  // check number of start points to determine number of files to write, with new filename eachtime
  std::cout << "Total " << _SWCNodeContainer.size() << " nodes..." <<std::endl;
  std::vector<SWCNode*>::iterator sit;
  float SCALE = 1.0f;

  for (unsigned int i = 0; i < _StartPoints.size(); ++i) 
  {
    std::stringstream ss;
    ss << "_" << i+1 << ".swc";
    std::string fname1 = fname;
    fname1.replace(fname.length()-4,8,ss.str());
    std::cout << "Writing SWCImage file " << fname1 << " \n Tree ID " << i+1 <<std::endl;
    std::ofstream ofile(fname1.c_str());
    //ofile << "#Neuron Tracing Code 3D, RPI" << std::endl;
    //ofile << "#author: AM" << std::endl;

    //make the LookUp table
    std::map<long, long> NodeIDToSWCIDMap;
    long ID = 1;
    long rootID = 1;
    for (sit = _SWCNodeContainer.begin(); sit != _SWCNodeContainer.end(); ++sit) 
    {
      if ((*sit)->TreeID == i+1) 
        NodeIDToSWCIDMap[(*sit)->ID] = ID++;			
    }
    std::cout << ID << " Nodes found  ";

    //create the SWCImage file
    for (sit = _SWCNodeContainer.begin(); sit != _SWCNodeContainer.end(); ++sit) 
    {
      if ((*sit)->TreeID == i+1) 
      {
        long id = NodeIDToSWCIDMap[(*sit)->ID];
        long pid = -1;
        long type = 3;
        if ((*sit)->PID > 0) 
          pid = NodeIDToSWCIDMap[(*sit)->PID];

        if(pid == -1) 
        {
          type = 1;
          rootID = NodeIDToSWCIDMap[(*sit)->ID];
        }

        //hack for when your parent was deleted but you didn't get assigned as
        //a child of the root
        if(pid == 0)
          pid = rootID;

        //get radius estimate for this node
        float radius = getRadius((*sit)->pos);

        ofile << id << " " << type << " " << SCALE*(*sit)->pos[0] << " "
          << SCALE*(*sit)->pos[1] << " " << SCALE*(*sit)->pos[2]-padz
          << " " <<  " " << radius << " " << pid << std::endl;
      }
    }
    ofile.close();
    std::cout << " file written. " << std::endl;
  }

  for (sit = _SWCNodeContainer.begin(); sit != _SWCNodeContainer.end(); ++sit) 
    delete (*sit);

  std::cout << " done! " << std::endl;
}

float MultipleNeuronTracer::getRadius(itk::Vector<float,3>& pos) 
{
  float r = 2.0f;
  itk::Vector<float,3> m1, m2, m;
  itk::Index<3> ndx;

  for (int iter = 0; iter < 20; ++iter) 
  {
    for( int i = 0; i<3; i++) 
    {
      m1[i] = pos[i] - vnl_math_max(2.0f*r, 5.0f);
      m2[i] = pos[i] + vnl_math_max(2.0f*r, 5.0f);
    }

    std::vector<float> c;
    c.reserve(4*4*int(r*r));
    float i1 = 0.0f, i2 = 0.0f, i1s = 0.0f, i2s = 0.0f;
    for (m[2] = m1[2]; m[2] <= m2[2]; m[2]++) 
    {
      for (m[1] = m1[1]; m[1] <= m2[1]; m[1]++) 
      {
        for (m[0] = m1[0]; m[0] <= m2[0]; m[0]++) 
        {
          ndx.CopyWithRound(m);
          itk::Vector<float,3> mm = pos - m;
          float d = mm.GetNorm();
          if (_PaddedCurvImage->GetBufferedRegion().IsInside(ndx)) 
          {
            float val = _PaddedCurvImage->GetPixel(ndx);
            if (d < r) 
            {
              i1 += val;
              ++i1s;
            }
            else 
            {
              i2 += val;
              ++i2s;
            }

            if (vnl_math_abs(d - r) < 0.7f) 
              c.push_back(val);						
          }
        }
      }
    }
    i1 /= i1s;
    i2 /= i2s;
    float dr = 0.0f;
    for (std::vector<float>::iterator it = c.begin(); it < c.end(); ++it) 
      dr += vnl_math_abs((*it) - i1) - vnl_math_abs((*it) - i2);

    dr *= 1.0f / float(c.size()); //rate
    dr = vnl_math_max(dr , -1.0f);
    dr = vnl_math_min(dr , 1.0f);
    r -= dr;
    r = vnl_math_max(r , 1.0f) ; 
  }
  return r;
}

///////////////////////////////////////////////////////////////////////////////////
void MultipleNeuronTracer::WriteSWCFile(std::string fname, unsigned int padz) 
{
  std::vector<SWCNode*>::iterator sit;
  std::cout << "Writing SWCImage file " << fname << " with " << _SWCNodeContainer.size() << " nodes...";
  std::ofstream ofile(fname.c_str());
  //ofile << "#Neuron Tracing Code 3D, RPI" << std::endl;
  //ofile << "#author: AM" << std::endl;
  for (sit = _SWCNodeContainer.begin(); sit != _SWCNodeContainer.end(); ++sit) 
  {
    //get radius estimate for this node
    float radius = getRadius((*sit)->pos);
    ofile << (*sit)->ID << " 3 " << (*sit)->pos[0] << " " << (*sit)->pos[1] << " "
      << (*sit)->pos[2]-padz << " " << radius <<" " << (*sit)->PID << std::endl;
    delete (*sit);
  }
  ofile.close();
  std::cout << " done! " << std::endl;
}


vtkSmartPointer< vtkTable > MultipleNeuronTracer::GetSWCTable(unsigned int padz) 
{
  vtkSmartPointer< vtkTable > SWCTable = vtkSmartPointer< vtkTable >::New();
  SWCTable->Initialize();

  vtkSmartPointer< vtkDoubleArray > column = vtkSmartPointer< vtkDoubleArray >::New();
  column->SetName("A");
	SWCTable->AddColumn(column);
	column = vtkSmartPointer< vtkDoubleArray >::New();
	column->SetName("B");
	SWCTable->AddColumn(column);
	column = vtkSmartPointer< vtkDoubleArray >::New();
	column->SetName("C");
	SWCTable->AddColumn(column);
	column = vtkSmartPointer< vtkDoubleArray >::New();
	column->SetName("D");
	SWCTable->AddColumn(column);
	column = vtkSmartPointer< vtkDoubleArray >::New();
	column->SetName("E");
	SWCTable->AddColumn(column);
	column = vtkSmartPointer< vtkDoubleArray >::New();
	column->SetName("F");
	SWCTable->AddColumn(column);
	column = vtkSmartPointer< vtkDoubleArray >::New();
	column->SetName("G");
	SWCTable->AddColumn(column);

	std::vector<SWCNode*>::iterator sit;
	for (sit = _SWCNodeContainer.begin(); sit != _SWCNodeContainer.end(); ++sit) 
	{
		//get radius estimate for this node
		float radius = getRadius((*sit)->pos);
		vtkSmartPointer< vtkVariantArray > row = vtkSmartPointer< vtkVariantArray >::New();
		row->InsertNextValue(vtkVariant((*sit)->ID));
		row->InsertNextValue(vtkVariant(3));
		row->InsertNextValue(vtkVariant((*sit)->pos[0]));
		row->InsertNextValue(vtkVariant((*sit)->pos[1]));
		row->InsertNextValue(vtkVariant((*sit)->pos[2]-padz));
		row->InsertNextValue(vtkVariant(radius));
		row->InsertNextValue(vtkVariant((*sit)->PID));
		SWCTable->InsertNextRow(row);
		delete (*sit);
	}

	return SWCTable;	
}

///////////////////////////////////////////////////////////////////////
vtkSmartPointer< vtkTable > MultipleNeuronTracer::GetDebrisTable(unsigned int padz) 
{
	vtkSmartPointer< vtkTable > DebrisTable = vtkSmartPointer< vtkTable >::New();
	DebrisTable->Initialize();

	vtkSmartPointer< vtkDoubleArray > column = vtkSmartPointer< vtkDoubleArray >::New();
	column->SetName("A");
	DebrisTable->AddColumn(column);
	column = vtkSmartPointer< vtkDoubleArray >::New();
	column->SetName("B");
	DebrisTable->AddColumn(column);
	column = vtkSmartPointer< vtkDoubleArray >::New();
	column->SetName("C");
	DebrisTable->AddColumn(column);


	std::vector<DebrisNode*>::iterator dit;
	for (dit = _DebrisNodeContainer.begin(); dit != _DebrisNodeContainer.end(); ++dit) 
	{
		vtkSmartPointer< vtkVariantArray > row = vtkSmartPointer< vtkVariantArray >::New();
		row->InsertNextValue(vtkVariant((*dit)->ndx[0]));
		row->InsertNextValue(vtkVariant((*dit)->ndx[1]));
		row->InsertNextValue(vtkVariant((*dit)->ndx[2]));

		DebrisTable->InsertNextRow(row);
		delete (*dit);
	}

	return DebrisTable;	
}


void MultipleNeuronTracer::GenerateTestImage(void) 
{
	_PaddedCurvImage = ImageType3D::New();
	_size[0] = 20; 
	_size[1] = 20; 
	_size[2] = 20;
	_PaddedCurvImage->SetRegions(_size);
	_PaddedCurvImage->Allocate();
	_PaddedCurvImage->FillBuffer(0.0);

	itk::Vector<float,3> dir; 
	dir.Fill(1.0f); 
	dir.Normalize();
	itk::Vector<float,3> acc, pos;
	pos.Fill(3.0f);
	itk::Index<3> ndx;
	ndx.CopyWithRound(pos);

	_PaddedCurvImage->SetPixel(ndx, 1.0f);

	for (int i=0; i<15; i++) 
	{
		float val = float(rand()%100) / 100.0f;
		for (int j = 0;j<3;j++) 
			acc[j] = (float(rand()%100) / 100.0f) - 0.5f;
		
		dir += acc*0.5;
		dir.Normalize();

		pos += dir;
		ndx.CopyWithRound(pos);
		_PaddedCurvImage->SetPixel(ndx,val);
	}

	WriteImage3D(std::string("GeneratedImage.mhd"), _PaddedCurvImage);
}

void MultipleNeuronTracer::WriteImage3D(std::string fname, MultipleNeuronTracer::ImageType3D::Pointer image)  
{
	std::cout << "Writing output file "<< fname << std::endl;
	typedef itk::ImageFileWriter<ImageType3D> WriterType;
	WriterType::GlobalWarningDisplayOff();
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(fname);
	writer->SetInput(image);
	writer->Update();
}

///////////////////////////////////////////////////////////////////////////////
void MultipleNeuronTracer::BlackOut(itk::Index<3> &stndx)
{
	for (long z = -3; z <=3 ; ++z) 
	{
		for (long y = -5; y <=5 ; ++y) 
		{
			for (long x = -5; x <=5 ; ++x) 
			{
				itk::Offset<3> _off = { {x,y,z} };
				itk::Index<3> n = stndx + _off;
				if ( (n[0] < 0) || (n[1] < 0) || (n[2] < 0) ||
					(n[0] >= (unsigned int)_size[0]) || (n[1] >= (unsigned int)_size[1]) ||
					(n[2] >= (unsigned int)_size[2]) )  
				{
						continue;
				}
				_PaddedCurvImage->SetPixel(n,1.0f);
				_NDXImage->SetPixel(n,0);
			}
		}
	}
}

int MultipleNeuronTracer::optionsCreate(const char* optfile, std::map<std::string,std::string>& options)
{
	options.clear();
	ifstream fin(optfile); assert(fin.good());
	std::string name;  fin>>name;
	while(fin.good()) {
		char cont[100];	 fin.getline(cont, 99);
		options[name] = std::string(cont);
		fin>>name;
	}
	fin.close();
	return 0;
}

void MultipleNeuronTracer::threeLevelMinErrorThresh(unsigned char* im, float* Alpha1, float* Alpha2, float* Alpha3, float* P_I1, float* P_I2, size_t r, size_t c, size_t z)
{
	//create a normalized image histogram
	float Hst[256];
	for(int i=0; i<256; i++)
		Hst[i] = 0.0;

	for(size_t i=0; i<r*c*z; i++)
	{
		int v = (int) im[i];
		Hst[v]++;
	}

	for(int i=0; i<256; i++)
		Hst[i] /= (r*c*z);


	//The three-level min error thresholding algorithm
	float P0, U0, P1, U1, P2, U2, U, J, min_J;
	min_J = 1000000.0;
	// Try this: we need to define a penalty term that depends on the number of parameters
	//The penalty term is given as 0.5*k*ln(n)
	//where k is the number of parameters of the model and n is the number of samples
	//In this case, k=6 and n=256
	double PenTerm3 =  sqrt(6.0)*log(256.0);
	for(int i=0; i<254; i++)//to set the first threshold
	{
		//compute the current parameters of the first component
		P0 = U0 = 0.0;		
		for(int l=0; l<=i; l++)
		{
			P0+=Hst[l];
			U0+=(l+1)*Hst[l];
		}
		U0 /= P0;

		for(int j=i+1; j<255; j++)//to set the second threshold
		{
			//compute the current parameters of the second component
			P1 = U1 = 0.0;		
			for(int l=i+1; l<=j; l++)
			{
				P1+=Hst[l];
				U1+=(l+1)*Hst[l];
			}
			U1 /= P1;

			//compute the current parameters of the third component
			P2 = U2 = 0.0;		
			for(int l=j+1; l<=255; l++)
			{
				P2+=Hst[l];
				U2+=(l+1)*Hst[l];
			}
			U2 /= P2;

			//compute the overall mean
			U = P0*U0 + P1*U1 + P2*U2;

			//Compute the current value of the error criterion function
			J =  U - (P0*(log(P0)+U0*log(U0))+ P1*(log(P1)+U1*log(U1)) + P2*(log(P2)+U2*log(U2)));
			//Add the penalty term
			J +=PenTerm3;

			if(J<min_J)
			{
				min_J = J;
				Alpha1[0] = U0;
				P_I1[0] = P0;
				Alpha2[0] = U1;
				P_I2[0] = P1;
				Alpha3[0] = U2;				
			}
		}
	}

	//try this: see if using two components is better
	//The penalty term is given as sqrt(k)*ln(n)	
	//In this case, k=4 and n=256
	double PenTerm2 =  2*log(256.0);
	for(int i=0; i<254; i++)//to set the first threshold
	{
		//compute the current parameters of the first component
		P0 = U0 = 0.0;		
		for(int l=0; l<=i; l++)
		{
			P0+=Hst[l];
			U0+=(l+1)*Hst[l];
		}
		U0 /= P0;

		for(int j=i+1; j<255; j++)//to set the second threshold
		{
			//compute the current parameters of the second component
			P1 = U1 = 0.0;		
			for(int l=j; l<=255; l++)
			{
				P1+=Hst[l];
				U1+=(l+1)*Hst[l];
			}
			U1 /= P1;

			//compute the overall mean
			U = P0*U0 + P1*U1;

			//Compute the current value of the error criterion function
			J =  U - (P0*(log(P0)+U0*log(U0))+ P1*(log(P1)+U1*log(U1)));
			//Add the penalty term
			J +=PenTerm2;
			if(J<min_J)
			{
				min_J = J;
				Alpha1[0] = U0;
				P_I1[0] = P0;
				Alpha2[0] = U1;
				P_I2[0] = P1;
				Alpha3[0] = -1; //Just a negative number to let the program knows that two levels will be used		
			}
		}
	}

}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	SWCImage NODE and HEAP NODE
//
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
SWCNode::SWCNode()
{
	this->ID = -1;
	this->PID = -1;
	this->TreeID = -1;
	this->IsLeaf = false;
	this->IsBranch = false;
	this->parent = NULL;
	this->children.reserve(2);
}

SWCNode::SWCNode(long id, long parent_id, long tree_id, itk::Index<3> index)
{
	this->ID = id;
	this->PID = parent_id;
	this->TreeID = tree_id;
	this->ndx = index;
	this->IsLeaf = false;
	this->IsBranch = false;
	this->parent = NULL;
	this->children.reserve(2);
}

SWCNode::SWCNode(long id, SWCNode * parent, long tree_id, itk::Index<3> index)
{
	this->ID = id;
	this->PID = parent->ID;
	this->TreeID = tree_id;
	this->ndx = index;
	this->IsLeaf = false;
	this->IsBranch = false;
	this->parent = parent;
	this->children.reserve(2);
}

HeapNode::HeapNode(itk::Index<3> n1, PixelType d)
{
	ndx = n1;
	KeyValue = d;
}

DebrisNode::DebrisNode(itk::Index<3> dndx)
{
	this->ndx=dndx;
}







