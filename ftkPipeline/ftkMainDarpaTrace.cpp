#include "ftkMainDarpaTrace.h"

void ftkMainDarpaTrace::readParameters( std::string segmentParams )
{
  _segmentParams = segmentParams;

  std::map< std::string, std::string > options;
  std::map< std::string, std::string >::iterator iter;

  options.clear();
  ifstream fin(segmentParams.c_str());
  assert(fin.good());
  std::string name;
  fin>>name;

  while(fin.good()) {
    char cont[1000];
    fin.getline(cont, 999);
    options[name] = std::string(cont);
    // 		std::cout << std::endl << name << " " << cont;
    fin>>name;
  }
  fin.close();

  // 	iter = options.find("-xSize");
  // 	if(iter!=options.end())
  // 	{ std::istringstream ss((*iter).second); ss >> _xSize;}
  // 	else
  // 	{ _xSize = 1024; printf("Choose _xSize = 1024 as default\n");}
  //
  // 	iter = options.find("-ySize");
  // 	if(iter!=options.end())
  // 	{ std::istringstream ss((*iter).second); ss >> _ySize;}
  // 	else
  // 	{ _ySize = 1024; printf("Choose _ySize = 1024 as default\n");}
  //
  // 	iter = options.find("-zSize");
  // 	if(iter!=options.end())
  // 	{ std::istringstream ss((*iter).second); ss >> _zSize;}
  // 	else
  // 	{ _zSize = 1024; printf("Choose _zSize = 1024 as default\n");}

  iter = options.find("-xTile");
  if(iter!=options.end())
  { std::istringstream ss((*iter).second); ss >> _xTile;}
  else
  { _xTile = 1024; printf("Choose xTile = 1024 as default\n");}

  iter = options.find("-yTile");
  if(iter!=options.end())
  { std::istringstream ss((*iter).second); ss >> _yTile;}
  else
  { _yTile = 1024; printf("Choose yTile = 1024 as default\n");}

  iter = options.find("-zTile");
  if(iter!=options.end())
  { std::istringstream ss((*iter).second); ss >> _zTile;}
  else
  { _zTile = 1024; printf("Choose _zTile = 1024 as default\n");}

  iter = options.find("-num_threads");
  if(iter!=options.end())
  { std::istringstream ss((*iter).second); ss >> _num_threads;}
  else
  { _num_threads = 80; printf("Choose _num_threads = 80 as default\n");}

  iter = options.find("-Cy5_Image");
  if(iter!=options.end())
  { std::istringstream ss((*iter).second); ss >> _Cy5_Image;}
  else
  { _Cy5_Image.clear(); printf("Choose _Cy5_Image = NULL as default\n");}

  iter = options.find("-TRI_Image");
  if(iter!=options.end())
  { std::istringstream ss((*iter).second); ss >> _TRI_Image;}
  else
  { _TRI_Image.clear(); printf("Choose _TRI_Image = NULL as default\n");}

  iter = options.find("-GFP_Image");
  if(iter!=options.end())
  { std::istringstream ss((*iter).second); ss >> _GFP_Image;}
  else
  { _GFP_Image.clear(); printf("Choose _GFP_Image = NULL as default\n");}

  iter = options.find("-DAP_Image");
  if(iter!=options.end())
  { std::istringstream ss((*iter).second); ss >> _DAP_Image;}
  else
  { _DAP_Image.clear(); printf("Choose _DAP_Image = NULL as default\n");}

  iter = options.find("-Soma_Centroids");
  if(iter!=options.end())
  { std::istringstream ss((*iter).second); ss >> _Soma_Centroids;}
  else
  { _Soma_Centroids.clear(); printf("Choose _Soma_Centroids = NULL as default\n");}

  iter = options.find("-Soma_Montage");
  if(iter!=options.end())
  { std::istringstream ss((*iter).second); ss >> _Soma_Montage;}
  else
  { _Soma_Montage.clear(); printf("Choose _Soma_Montage = NULL as default\n");}

  iter = options.find("-isSmall");
  if(iter!=options.end())
  { std::istringstream ss((*iter).second); ss >> _isSmall;}
  else
  { _isSmall = 1; printf("Choose _isSmall = 1 as default\n");}

  iter = options.find("-traceParams");
  if(iter!=options.end())
  { std::istringstream ss((*iter).second); ss >> _traceParams;}
  else
  { _traceParams.clear(); printf("Choose _traceParams = NULL as default\n");}

  iter = options.find("-outPath");
  if(iter!=options.end())
  { std::istringstream ss((*iter).second); ss >> _outPath;}
  else
  { _outPath.clear(); printf("Choose _outPath = NULL as default\n");}

  iter = options.find("-outPathDebug");
  if(iter!=options.end())
  { std::istringstream ss((*iter).second); ss >> _outPathDebug;}
  else
  { _outPathDebug.clear(); printf("Choose _outPathDebug = NULL as default\n");}

  iter = options.find("-outPathDebugLevel2");
  if(iter!=options.end())
  { std::istringstream ss((*iter).second); ss >> _outPathDebugLevel2;}
  else
  { _outPathDebugLevel2.clear(); printf("Choose _outPathDebugLevel2 = NULL as default\n");}

  iter = options.find("-outPathTemp");
  if(iter!=options.end())
  { std::istringstream ss((*iter).second); ss >> _outPathTemp;}
  else
  { _outPathTemp.clear(); printf("Choose _outPathTemp = NULL as default\n");}

  iter = options.find("-optimizeCoverage?");
  if(iter!=options.end())
  { std::istringstream ss((*iter).second); ss >> _optimizeCoverage;}
  else
  { _optimizeCoverage = false; printf("Choose _optimizeCoverage? = NULL as default\n");}

  iter = options.find("-overridedefaultsTraceParams");
  if(iter!=options.end())
  { std::istringstream ss((*iter).second); ss >> _overridedefaultsTraceParams;}
  else
  { _overridedefaultsTraceParams = "NO"; printf("Choose _overridedefaultsTraceParams = NO as default\n");}

  _Cy5_ImageNRRD = _Cy5_Image+".nrrd";
  _TRI_ImageNRRD = _TRI_Image+".nrrd";
  _GFP_ImageNRRD = _GFP_Image+".nrrd";
  _DAP_ImageNRRD = _DAP_Image+".nrrd";
  _Soma_MontageNRRD = _Soma_Montage+".nrrd";


  // Print Parameters
  std::cout << std::endl << "This are the parameters";
  std::cout << std::endl << "_xSize: " << _xSize;
  std::cout << std::endl << "_ySize: " << _ySize;
  std::cout << std::endl << "_zSize: " << _zSize;
  std::cout << std::endl << "_xTile: " << _xTile;
  std::cout << std::endl << "_yTile: " << _yTile;
  std::cout << std::endl << "_zTile: " << _zTile;
  std::cout << std::endl << "_num_threads: " << _num_threads;
  std::cout << std::endl << "_Cy5_Image: " << _Cy5_Image;
  std::cout << std::endl << "_TRI_Image: " << _TRI_Image;
  std::cout << std::endl << "_GFP_Image: " << _GFP_Image;
  std::cout << std::endl << "_DAP_Image: " << _DAP_Image;
  std::cout << std::endl << "_Soma_Centroids: " << _Soma_Centroids;
  std::cout << std::endl << "_Soma_Montage: " << _Soma_Montage;
  std::cout << std::endl << "_isSmall: " << _isSmall;
  std::cout << std::endl << "_traceParams: " << _traceParams;
  std::cout << std::endl << "_outPath: " << _outPath;
  std::cout << std::endl << "_outPathDebug: " << _outPathDebug;
  std::cout << std::endl << "_outPathTemp: " << _outPathTemp;
  std::cout<< std::endl <<"_overridedefaultsTraceParams "<<_overridedefaultsTraceParams;
  std::cout << std::endl << "_optimizeCoverage? " << _optimizeCoverage;
}

void ftkMainDarpaTrace::runPreprocesing()
{
  rawImageType_flo::Pointer MontageGFP_Image = readImage< rawImageType_flo >(_GFP_ImageNRRD.c_str());


  //RescaleFilterType_floToflo::Pointer rescaler = RescaleFilterType_floToflo::New();
  //rescaler->SetOutputMinimum(0.0);
  //rescaler->SetOutputMaximum(1.0);
  //rescaler->SetInput( MontageGFP_Image );

  //std::cout<<std::endl<<"Rescaler"<<std::flush;

  //MedianFilterType_floToflo::Pointer medfilt = MedianFilterType_floToflo::New();
  //medfilt->SetInput( rescaler->GetOutput() );
  //rawImageType_flo::SizeType rad = { {1, 1, 1} };
  //medfilt->SetRadius(rad);
  //medfilt->Update();

  //std::cout<<std::endl<<"MedFilter"<<std::flush;

  _GFP_ImagePREPMNT = _outPathTemp+"/GFP_MNT_PRE.nrrd";
  //writeImage< rawImageType_flo >( medfilt->GetOutput(), _GFP_ImagePREPMNT.c_str());

  _xSize = MontageGFP_Image->GetLargestPossibleRegion().GetSize()[0];
  _ySize = MontageGFP_Image->GetLargestPossibleRegion().GetSize()[1];
  _zSize = MontageGFP_Image->GetLargestPossibleRegion().GetSize()[2];

  std::cout << std::endl << "ACAASD: " << MontageGFP_Image->GetLargestPossibleRegion().GetSize();
}

void ftkMainDarpaTrace::computeTileGVFAndVesselness()
{
  std::cout<<"In Pre-Compute GVF and Vesselness"<<std::endl;
  int num_iteration = 15;
  int smoothing_scale = 1;
  _GVF_ImagePREMNT = _outPathTemp+"/GFP_MNT_PRE_GVF_";
  _Vesselness_ImagePREMNT = _outPathTemp+"/GFP_MNT_PRE_Vesselness_";
  std::stringstream str_bigTile;
  computeSplitConst();
  for( unsigned int bigTile = 0; bigTile<_numDivisionsInRowCEN ; ++bigTile )
  {
    std::cout<<std::endl<<"bigTile: "<<bigTile;
    str_bigTile << bigTile;
    std::cout << std::endl << "_initialBigTileLOG: " << _initialBigTileLOG[bigTile][0] <<", "<<_initialBigTileLOG[bigTile][1] <<", "<<_initialBigTileLOG[bigTile][2];
    std::cout << std::endl << "_sizeOfBigTilesLOG: " << _sizeOfBigTilesLOG[bigTile][0] <<", "<<_sizeOfBigTilesLOG[bigTile][1] <<", "<<_sizeOfBigTilesLOG[bigTile][2];

    itk::Index<3> initialBigIndexLOG;
    itk::Size<3> sizeOfTheRegionLOG;

    initialBigIndexLOG[0] = _initialBigTileLOG[bigTile][0];
    initialBigIndexLOG[1] = _initialBigTileLOG[bigTile][1];
    initialBigIndexLOG[2] = _initialBigTileLOG[bigTile][2];

    sizeOfTheRegionLOG[0] = _sizeOfBigTilesLOG[bigTile][0];
    sizeOfTheRegionLOG[1] = _sizeOfBigTilesLOG[bigTile][1]+_sizeOfBigTilesLOG[bigTile+1][1];
    sizeOfTheRegionLOG[2] = _sizeOfBigTilesLOG[bigTile][2];

    std::cout << std::endl << "sizeOfTheRegionLOG: " << sizeOfTheRegionLOG[0] <<", "<<sizeOfTheRegionLOG[1] <<", "<<sizeOfTheRegionLOG[2];

    // 		std::cout << std::endl << initialBigIndexLOG;
    // 		std::cout << std::endl << "SIZE: " << sizeOfTheRegionLOG[1];

    rawImageType_flo::RegionType desiredRegionBigTileLOG;
    desiredRegionBigTileLOG.SetSize(sizeOfTheRegionLOG);
    desiredRegionBigTileLOG.SetIndex(initialBigIndexLOG);

    itk::Index<3> initialBigIndexCEN;
    itk::Size<3> sizeOfTheRegionCEN;

    initialBigIndexCEN[0] = _initialBigTileCEN[bigTile][0];
    initialBigIndexCEN[1] = _initialBigTileCEN[bigTile][1];
    initialBigIndexCEN[2] = _initialBigTileCEN[bigTile][2];
    sizeOfTheRegionCEN[0] = _sizeOfBigTilesCEN[bigTile][0];
    sizeOfTheRegionCEN[1] = _sizeOfBigTilesCEN[bigTile][1];
    sizeOfTheRegionCEN[2] = _sizeOfBigTilesCEN[bigTile][2];
    std::cout << std::endl << "_sizeOfBigTilesCEN: " << _sizeOfBigTilesCEN[bigTile][0] <<", "<<_sizeOfBigTilesCEN[bigTile][1] <<", "<<_sizeOfBigTilesCEN[bigTile][2];
    std::cout << std::endl << desiredRegionBigTileLOG;
    _img_traceDesiredRegion = readImageRegion< rawImageType_flo >( _GFP_ImagePREPMNT.c_str(), desiredRegionBigTileLOG );
    std::cout<<"done Reading***********************************************************"<<std::endl;

    //Logic to compute the GVF and Vesselness
    MultipleNeuronTracer * MNT = new MultipleNeuronTracer();
    //this->computeGVF(100,num_iteration,smoothing_scale);
    //MNT->LoadCurvImage_2(_img_traceDesiredRegion);
    double origin_zero[3];
    origin_zero[0] = 0.0;
    origin_zero[1] = 0.0;
    origin_zero[2] = 0.0;
    _img_traceDesiredRegion->SetOrigin(origin_zero);
    std::cout<<"compute GVF*******************"<<std::endl;
    MNT->computeGVF_2(_img_traceDesiredRegion,100,num_iteration,smoothing_scale);
    std::string gvfPath = _GVF_ImagePREMNT+str_bigTile.str()+".nrrd";
    GradientImageType::Pointer gvfImage = MNT->getGVFImage();
    double origin_indx[3];
    origin_indx[0] = initialBigIndexLOG[0];
    origin_indx[1] = initialBigIndexLOG[1];
    origin_indx[2] = initialBigIndexLOG[2];
    gvfImage->SetOrigin(origin_indx);
    gvfImage->Update();
    writeImage< GradientImageType >( gvfImage, gvfPath.c_str());

    std::cout<<"compute Vesselness******************"<<std::endl;
    MNT->ComputeGVFVesselness_2(_img_traceDesiredRegion);
    ImageType3D::Pointer vesselImage = MNT->getVessleness();
    vesselImage->SetOrigin(origin_indx);
    vesselImage->Update();
    std::string vesselPath = _Vesselness_ImagePREMNT+str_bigTile.str()+".nrrd";
    writeImage< rawImageType_flo >( vesselImage, vesselPath.c_str());
    delete MNT;
  }



}
void ftkMainDarpaTrace::runTracing()
{
  std::cout << std::endl << "LETS TRACE";

  if( _isSmall == 1 )
  {
    itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1); // This one can not be changed
    //itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1); // This one can chenga
    int num_threads = 1;
#ifdef _OPENMP
    omp_set_nested(1);
#if _OPENMP >= 200805L
    omp_set_max_active_levels(2);
#endif
    omp_set_num_threads(num_threads);
#endif
  }

  int counterCentro = 0;
  std::vector< itk::Index<3> > centroid_list = getCentroidList();
  std::cout << "Number of cells to be traced : " << centroid_list.size() << "\n";

  std::string SWCFilename = _outPath + "/OnlySWC.xml";
  std::ofstream outSWCFile;
  outSWCFile.open(SWCFilename.c_str());
  outSWCFile << "<?xml\tversion=\"1.0\"\t?>\n";
  outSWCFile << "<Source>\n\n";

  computeSplitConst();

  // CAN NOT BE RUN PARALLEL
  for( unsigned int bigTile = 0; bigTile<_numDivisionsInRowCEN ; ++bigTile )
  {
    std::cout<<std::endl<<"bigTile: "<<bigTile;
    std::cout << std::endl << "_initialBigTileLOG: " << _initialBigTileLOG[bigTile][0] <<", "<<_initialBigTileLOG[bigTile][1] <<", "<<_initialBigTileLOG[bigTile][2];
    std::cout << std::endl << "_sizeOfBigTilesLOG: " << _sizeOfBigTilesLOG[bigTile][0] <<", "<<_sizeOfBigTilesLOG[bigTile][1] <<", "<<_sizeOfBigTilesLOG[bigTile][2];

    // 			stringstream out34;
    // 			out34<<bigTile;
    // 			string srr = out34.str();
    // 			std::string SWCFilenameDivided = _outPath + "/TracesAndSomasDivided/OnlySWC_" + srr +".xml";
    // 			std::ofstream outfileDivided;
    // 			outfileDivided.open(SWCFilenameDivided.c_str());
    // 			outfileDivided << "<?xml\tversion=\"1.0\"\t?>\n";
    // 			outfileDivided << "<Source>\n\n";



    std::stringstream str_bigTile;
    str_bigTile << bigTile;
    itk::Index<3> initialBigIndexLOG;
    itk::Size<3> sizeOfTheRegionLOG;

    initialBigIndexLOG[0] = _initialBigTileLOG[bigTile][0];
    initialBigIndexLOG[1] = _initialBigTileLOG[bigTile][1];
    initialBigIndexLOG[2] = _initialBigTileLOG[bigTile][2];

    sizeOfTheRegionLOG[0] = _sizeOfBigTilesLOG[bigTile][0];
    sizeOfTheRegionLOG[1] = _sizeOfBigTilesLOG[bigTile][1]+_sizeOfBigTilesLOG[bigTile+1][1];
    sizeOfTheRegionLOG[2] = _sizeOfBigTilesLOG[bigTile][2];

    std::cout << std::endl << "sizeOfTheRegionLOG: " << sizeOfTheRegionLOG[0] <<", "<<sizeOfTheRegionLOG[1] <<", "<<sizeOfTheRegionLOG[2];

    // 		std::cout << std::endl << initialBigIndexLOG;
    // 		std::cout << std::endl << "SIZE: " << sizeOfTheRegionLOG[1];

    itk::Index<3> initialBigIndexLOG_WithoutIndex;
    initialBigIndexLOG_WithoutIndex[0] = 0;
    initialBigIndexLOG_WithoutIndex[1] = 0;
    initialBigIndexLOG_WithoutIndex[2] = 0;

    rawImageType_flo::RegionType desiredRegionBigTileLOG;
    rawImageType_flo::RegionType desiredRegionBigTileLOG_WithoutIndex;
    GradientImageType::RegionType desiredRegionBigTileLOGGradient;

    desiredRegionBigTileLOG.SetSize(sizeOfTheRegionLOG);
    desiredRegionBigTileLOG.SetIndex(initialBigIndexLOG);

    desiredRegionBigTileLOGGradient.SetSize(sizeOfTheRegionLOG);
    desiredRegionBigTileLOGGradient.SetIndex(initialBigIndexLOG);

    desiredRegionBigTileLOG_WithoutIndex.SetSize(sizeOfTheRegionLOG);
    desiredRegionBigTileLOG_WithoutIndex.SetIndex(initialBigIndexLOG_WithoutIndex);

    itk::Index<3> initialBigIndexCEN;
    itk::Size<3> sizeOfTheRegionCEN;

    initialBigIndexCEN[0] = _initialBigTileCEN[bigTile][0];
    initialBigIndexCEN[1] = _initialBigTileCEN[bigTile][1];
    initialBigIndexCEN[2] = _initialBigTileCEN[bigTile][2];

    sizeOfTheRegionCEN[0] = _sizeOfBigTilesCEN[bigTile][0];
    sizeOfTheRegionCEN[1] = _sizeOfBigTilesCEN[bigTile][1];
    sizeOfTheRegionCEN[2] = _sizeOfBigTilesCEN[bigTile][2];

    std::cout << std::endl << "_sizeOfBigTilesCEN: " << _sizeOfBigTilesCEN[bigTile][0] <<", "<<_sizeOfBigTilesCEN[bigTile][1] <<", "<<_sizeOfBigTilesCEN[bigTile][2];
    std::cout << std::endl << desiredRegionBigTileLOG;

    // 			std::vector< rawImageType_flo::Pointer > LoGDesiredRegion;
    // 			LoGDesiredRegion.resize(6);

    // 			rawImageType_uint::Pointer _somaMontageDesiredRegion;
    // 			rawImageType_flo::Pointer _img_traceDesiredRegion;

    // 			std::string tempFileName_42 = _GFP_ImagePREPMNT
    if( _isSmall == 1 )
    {
      itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1); // This one can not be changed
      //itk::MultiThreader::SetGlobalMaximumNumberOfThreads(80); // This one can chenga
      int num_threads = 1;
#ifdef _OPENMP
      omp_set_nested(1);
#if _OPENMP >= 200805L
      omp_set_max_active_levels(2);
#endif
      omp_set_num_threads(num_threads);
#endif
    }
    _img_traceDesiredRegion = readImageRegion< rawImageType_flo >( _GFP_ImagePREPMNT.c_str(), desiredRegionBigTileLOG );
    _somaMontageDesiredRegion = readImageRegion< rawImageType_uint >( _Soma_MontageNRRD.c_str(), desiredRegionBigTileLOG );


    std::string vesselPath = _Vesselness_ImagePREMNT+str_bigTile.str()+".nrrd";
    std::string gvfPath = _GVF_ImagePREMNT+str_bigTile.str()+".nrrd";

    // Compute GVF and Vesselnes
    std::cout << std::endl << "Now the GVF and Vessel will be computed " << std::flush;
    int num_iteration = 15;
    int smoothing_scale = 1;
    itk::MultiThreader::SetGlobalDefaultNumberOfThreads(80); // This one can chenga

    MultipleNeuronTracer * MNT = new MultipleNeuronTracer();
    MNT->computeGVF_2(_img_traceDesiredRegion,100,num_iteration,smoothing_scale);
    _img_GVFDesiredRegion = MNT->getGVFImage();
    std::cout << std::endl << "GVF DONE" << std::flush;

    MNT->ComputeGVFVesselness_2(_img_traceDesiredRegion);
    _img_VesselDesiredRegion = MNT->getVessleness();

    delete MNT;
    itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1); // This one can chenga


    //Load the pre-computed the GVF and Vesselness
    //_img_GVFDesiredRegion = readImage< GradientImageType >( gvfPath.c_str());
    //_img_VesselDesiredRegion = readImage< rawImageType_flo >( vesselPath.c_str() );


#pragma omp parallel for num_threads(_num_threads) schedule(dynamic, 1)
    for(  long long a=0; a<centroid_list.size(); ++a )
      if( _isSmall == 1 )
      {
        itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1); // This one can not be changed
        //itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1); // This one can chenga
        int num_threads = 1;
#ifdef _OPENMP
        omp_set_nested(1);
#if _OPENMP >= 200805L
        omp_set_max_active_levels(2);
#endif
        omp_set_num_threads(num_threads);
#endif
      }



#pragma omp parallel for num_threads(_num_threads) schedule(dynamic, 1)
#if _OPENMP >= 200805L
    for( unsigned long long a=0; a<centroid_list.size(); ++a )
#else
      for( long long a=0; a<centroid_list.size(); ++a )	//OpenMP2.5 requires signed intergral type
#endif
      {
        int x, y, z;
        std::stringstream ssx, ssy, ssz;

        x = centroid_list[a][0];
        y = centroid_list[a][1];
        z = centroid_list[a][2];

        unsigned long long yMin = initialBigIndexCEN[1];
        unsigned long long yMax = initialBigIndexCEN[1] + sizeOfTheRegionCEN[1];

        if(!( (yMin <= y) && (y < yMax) ) )
        {
          continue;
        }

        ssx << x; ssy << y; ssz << z;

#pragma omp critical
        {
          counterCentro++;
          std::cout<<std::endl<<"\t\t\t\t asdfasdf ----->>>>> " << counterCentro << ", of " << centroid_list.size();
          std::cout<<", x: "<<centroid_list[a][0]<<", y: "<<centroid_list[a][1]<<", z: "<<centroid_list[a][2]<<std::endl;
        }

        // 			std::cout << std::endl << "x: " << x << ", y: " << y << ", z: " << z;
        // 			std::cout << std::endl << "xTile: " << _xTile << ", yTile: " << _yTile << ", zTile: " << _zTile;

        std::stringstream ssx_off, ssy_off, ssz_off;
        std::stringstream ssx_offBig, ssy_offBig, ssz_offBig;
        ssx_offBig << 0;
        ssy_offBig << 0;
        ssz_offBig << 0;
        if(x >= _xTile/2)
          ssx_off << x - _xTile/2;
        else
          ssx_off << 0;
        if(y >= _yTile/2)
          ssy_off << y - _yTile/2;
        else
          ssy_off << 0;
        if(z >= _zTile/2)
          ssz_off << z - _zTile/2;
        else
          ssz_off << 0;

        // 			std::cout << std::endl << "_initialBigTileLOG: " << _initialBigTileLOG[bigTile][1];

        int x_local = x;
        int y_local = y - _initialBigTileLOG[bigTile][1];
        int z_local = z;

        // 			std::cout << std::endl << "x_local: " << x_local << ", y_local: " << y_local << ", z_local: " << z_local;
        // 			std::cout << std::endl << "x: " << x << ", y: " << y << ", z: " << z;

        if(x_local >= _xTile/2)
          ssx_offBig << x_local - _xTile/2;
        else
          ssx_off << 0;
        if(y_local >= _yTile/2)
          ssy_offBig << y_local - _yTile/2;
        else
          ssy_off << 0;
        if(z_local >= _zTile/2)
          ssz_offBig << z_local - _zTile/2;
        else
          ssz_off << 0;

        //########    CROP THE DESIRED DICE FROM THE GFP AND SOMA MONTAGES   ########


        //########    FETCH ALL CENTROIDS THAT FALL WITHIN THE DICE    ########

        std::vector< itk::Index<3> > soma_Table = getSomaTable(centroid_list, x, y, z );

        // 			//########    RUN TRACING    ########
        MultipleNeuronTracer * MNT = new MultipleNeuronTracer();

        //
        // 			Automatic parameter estimation
        //
        //MNT->LoadCurvImage_1(img_trace, 0);
        //std::cout << std::endl << "LAREGION ES: " << _img_traceDesiredRegion;
        rawImageType_flo::Pointer img_trace;
        GradientImageType::Pointer img_gvf;
        rawImageType_flo::Pointer img_vessel;

#pragma omp critical
        {
          img_trace = cropImages< rawImageType_flo >( _img_traceDesiredRegion, x, y, z);
          img_gvf = cropImages< GradientImageType >( _img_GVFDesiredRegion, x, y, z);;
          img_vessel = cropImages< rawImageType_flo >( _img_VesselDesiredRegion, x, y, z);;
          MNT->LoadCurvImage_2(img_trace);
          MNT->setGVFImage(img_gvf);
          MNT->setVesselImage(img_vessel);
        }

        //MNT->RunMask();
        /*MNT->LoadParameters_1(_traceParams.c_str(),5);*/
        float calc_intensity_threshold = 0;
        float calc_contrast_threshold = 0;
        if(_overridedefaultsTraceParams == "YES")
        {
          std::vector<float> features = this->computeFeatures(img_trace);
          calc_intensity_threshold = getCalcThreshold(features,"intensity");
          calc_contrast_threshold = getCalcThreshold(features,"contrast");
          // For some images the threshold goes to negative in that case use the once that is specified in the
          // option_mnt
          if(calc_intensity_threshold < 0 || calc_contrast_threshold < 0 )
          {
            MNT->LoadParameters(_traceParams.c_str(),5);
          }
          else
          {
            MNT->LoadParameters_1(_traceParams.c_str(),calc_intensity_threshold,calc_contrast_threshold,350);
          }
        }else
        {
          if(_optimizeCoverage == 1){
            std::string coverageFileName = _outPathDebug + "/Coverage_" + ssx.str() + "_" + ssy.str() + "_" + ssz.str() + ".txt";

#pragma omp critical
            {
              MNT->OptimizeCoverage(coverageFileName, true);
            }
          }

          MNT->LoadParameters(_traceParams.c_str(), 6);
        }


        MNT->ReadStartPoints_1(soma_Table, 0);
        // 				MNT->SetCostThreshold(1000);
        MNT->SetCostThreshold(MNT->cost_threshold);

        // 	// 			MNT->LoadSomaImage_1(img_soma_yan);
        bool flagLog = false;
        bool flagPreComputedGVFVessel = true;
        MNT->setFlagOutLog(flagLog);
        MNT->RunGVFTracing(flagPreComputedGVFVessel);

#pragma omp critical
        {
          rawImageType_uint::Pointer img_soma = cropImages< rawImageType_uint >( _somaMontageDesiredRegion, x, y, z);
          MNT->RemoveSoma( img_soma );
        }
        //
        x = std::min(_xTile/2, x);
        y = std::min(_yTile/2, y);
        z = std::min(_zTile/2, z);
        //
        vtkSmartPointer< vtkTable > swcTable = MNT->GetSWCTable(0);
        delete MNT;
        std::string swcFilename = _outPath + "/Trace_" + ssx.str() + "_" + ssy.str() + "_" + ssz.str() + "_ANT.swc";
        WriteCenterTrace(swcTable, x, y, z, swcFilename);
        //
#pragma omp critical
        {
          outSWCFile << "\t<File\tFileName=\"Trace_" << ssx.str() << "_" << ssy.str() << "_" << ssz.str() << "_ANT.swc\"\tType=\"Trace\"\ttX=\"" << ssx_off.str() << "\"\ttY=\"" << ssy_off.str() << "\"\ttZ=\"" << ssz_off.str() << "\"/>\n";
        }

        // // 				vtkSmartPointer< vtkTable > swcTable = MNT->GetSWCTable(0);
        // 				std::string swcFilenameDivided = _outPath + "/TracesAndSomasDivided/Trace_BigTile_" + srr + "_" + ssx.str() + "_" + ssy.str() + "_" + ssz.str() + "_ANT.swc";
        // 				WriteCenterTrace(swcTable, x, y, z, swcFilenameDivided);
        // //
        // 				#pragma omp critical
        // 				{
        // 					outfileDivided << "\t<File\tFileName=\"Trace_BigTile_" << srr << "_" << ssx.str() << "_" << ssy.str() << "_" << ssz.str() << "_ANT.swc\"\tType=\"Trace\"\ttX=\"" << ssx_offBig.str() << "\"\ttY=\"" << ssy_offBig.str() << "\"\ttZ=\"" << ssz_offBig.str() << "\"/>\n";
        // 				}
      }
  }
}








void ftkMainDarpaTrace::computeSplitConst( )
{
  // 	std::cout << std::endl << "ACA VAMOS";
  // 	std::cout << std::endl << "PP: " << _xSize << " " << _xTile;
  // 	std::cout << std::endl << "PP: " << _ySize << " " << _yTile;
  // 	std::cout << std::endl << "PP: " << _zSize << " " << _zTile;
  /////////////////////////////////////////////////////////////////////////////////////////////
  int numDivisionsInRowLOG = (int)floor((double)_ySize/(2*_yTile+10)); // Minimum value 3, otherwise does not make sense
  if( numDivisionsInRowLOG < 3 )
    numDivisionsInRowLOG = 2;
  _numDivisionsInRowCEN = numDivisionsInRowLOG - 1;

  std::cout << std::endl << "O: " << numDivisionsInRowLOG << " " << _numDivisionsInRowCEN;

  _sizeOfBigTilesLOG.resize(numDivisionsInRowLOG);
  _sizeOfBigTilesCEN.resize(_numDivisionsInRowCEN);
  for( unsigned int jj=0; jj<numDivisionsInRowLOG; ++jj )
  {
    _sizeOfBigTilesLOG[jj][0] = _xSize;
    _sizeOfBigTilesLOG[jj][1] = (unsigned long long)floor((double)_ySize/(double)numDivisionsInRowLOG);
    _sizeOfBigTilesLOG[jj][2] = _zSize;
  }
  std::cout << std::endl << "O: " << _sizeOfBigTilesLOG[numDivisionsInRowLOG-1][1];
  _sizeOfBigTilesLOG[numDivisionsInRowLOG-1][1] = _ySize - (numDivisionsInRowLOG-1)*_sizeOfBigTilesLOG[numDivisionsInRowLOG-1][1];
  std::cout << std::endl << "O: " << _sizeOfBigTilesLOG[numDivisionsInRowLOG-1][1];
  if( _sizeOfBigTilesLOG[0][1] <2*_yTile)
  {
    std::cout<<std::endl<<"MISTAKE_1";
  }
  _sizeOfBigTilesCEN[0][0] = _xSize;
  _sizeOfBigTilesCEN[0][1] = _sizeOfBigTilesLOG[0][1]+(unsigned long long)floor((double)_sizeOfBigTilesLOG[0][1]/2);
  _sizeOfBigTilesCEN[0][2] = _zSize;
  for( unsigned int jj=1; jj<_numDivisionsInRowCEN; ++jj )
  {
    _sizeOfBigTilesCEN[jj][0] = _xSize;
    _sizeOfBigTilesCEN[jj][1] = _sizeOfBigTilesLOG[jj][1];
    _sizeOfBigTilesCEN[jj][2] = _zSize;
  }
  if( _numDivisionsInRowCEN == 1 )
  {
    _sizeOfBigTilesCEN[_numDivisionsInRowCEN-1][1] = _ySize;
  }
  else
  {
    _sizeOfBigTilesCEN[_numDivisionsInRowCEN-1][1] = _ySize -_sizeOfBigTilesCEN[0][1] -  (_numDivisionsInRowCEN-2)*_sizeOfBigTilesLOG[numDivisionsInRowLOG-2][1];
  }


  unsigned long long total = 0;
  for( unsigned int jj=0; jj<numDivisionsInRowLOG; ++jj )
  {
    std::cout<<std::endl<<"Row Size Big: "<<_sizeOfBigTilesLOG[jj][1];
    total = total + _sizeOfBigTilesLOG[jj][1];
  }
  if( total != _ySize )
  {
    std::cout<<std::endl<<"MISTAKE_2";
  }
  std::cout<<std::endl<<"\t\tTotalLOG: "<<total<<", Ysize: "<<_ySize;

  unsigned long long total2 = 0;
  for( unsigned int jj=0; jj<_numDivisionsInRowCEN; ++jj )
  {
    std::cout<<std::endl<<"Col Size Big"<<_sizeOfBigTilesCEN[jj][1];
    total2 = total2 + _sizeOfBigTilesCEN[jj][1];
  }
  if( total2 != _ySize )
  {
    std::cout<<std::endl<<"MISTAKE_3";
  }
  std::cout<<std::endl<<"\t\tTotalCEN: "<<total2<<", Ysize: "<<_ySize;

  _initialBigTileLOG.resize(numDivisionsInRowLOG);
  _initialBigTileLOG[0][0] = 0;
  _initialBigTileLOG[0][1] = 0;
  _initialBigTileLOG[0][2] = 0;
  for( unsigned int jj=1; jj<numDivisionsInRowLOG; ++jj )
  {
    _initialBigTileLOG[jj][0] = 0;
    _initialBigTileLOG[jj][1] = _sizeOfBigTilesLOG[jj-1][1] + _initialBigTileLOG[jj-1][1];
    _initialBigTileLOG[jj][2] = 0;
    std::cout << std::endl << "\t->IniLog:\t" <<_initialBigTileLOG[jj][1];
  }

  _initialBigTileCEN.resize(_numDivisionsInRowCEN);
  _initialBigTileCEN[0][0] = 0;
  _initialBigTileCEN[0][1] = 0;
  _initialBigTileCEN[0][2] = 0;
  for( unsigned int jj=1; jj<_numDivisionsInRowCEN; ++jj )
  {
    _initialBigTileCEN[jj][0] = 0;
    _initialBigTileCEN[jj][1] = _sizeOfBigTilesCEN[jj-1][1] + _initialBigTileCEN[jj-1][1];
    _initialBigTileCEN[jj][2] = 0;
    std::cout << std::endl << "\t->IniCen:\t" <<_initialBigTileCEN[jj][1];
  }
  std::cout << std::endl << "-----------------END";
}

std::vector< itk::Index<3> > ftkMainDarpaTrace::getCentroidList()
{
  std::cout << std::endl << _Soma_Centroids;
  std::cout << std::endl << _Soma_Centroids;
  std::cout << std::endl << _Soma_Centroids;

  vtkSmartPointer<vtkTable> somaCentroidsTable = ftk::LoadTable(_Soma_Centroids);

  std::vector< itk::Index<3> > centroid_list;
  for(int r=0; r<(int)somaCentroidsTable->GetNumberOfRows(); ++r)
  {
    int cx = somaCentroidsTable->GetValue(r, 0).ToInt();
    int cy = somaCentroidsTable->GetValue(r, 1).ToInt();
    int cz = somaCentroidsTable->GetValue(r, 2).ToInt();

    itk::Index<3> cen;
    cen[0] = cx; cen[1] = cy; cen[2] = cz;
    centroid_list.push_back(cen);
  }
  // 	std::cout << std::endl << "-----------------END";
  // 	std::cout << std::endl << "-----------------END";
  return centroid_list;
}

std::vector< itk::Index<3> > ftkMainDarpaTrace::getSomaTable( std::vector< itk::Index<3> > centroid_list, int x, int y, int z )
{
  itk::Index<3> centroid;
  centroid[0] = ((x - _xTile/2)>0) ? _xTile/2:x;
  centroid[1] = ((y - _yTile/2)>0) ? _yTile/2:y;
  centroid[2] = ((z - _zTile/2)>0) ? _zTile/2:z;

  itk::Index<3> start;
  start[0] = ((x - _xTile/2)>0) ? (x - _xTile/2):0;
  start[1] = ((y - _yTile/2)>0) ? (y - _yTile/2):0;
  start[2] = ((z - _zTile/2)>0) ? (z - _zTile/2):0;

  itk::Size<3> size;
  size[0] = ((x+_xTile/2)<_xSize) ? _xTile : (_xTile/2+_xSize-x-1);
  size[1] = ((y+_yTile/2)<_ySize) ? _yTile : (_yTile/2+_ySize-y-1);
  size[2] = ((z+_zTile/2)<_zSize) ? _zTile : (_zTile/2+_zSize-z-1);

  std::vector< itk::Index<3> > soma_Table;
  for(int ctr =0; ctr<centroid_list.size() ; ++ctr)
  {
    itk::Index<3> cen = centroid_list[ctr];
    if( (cen[0]>=start[0]) && (cen[0]<(start[0]+size[0])) && (cen[1]>=start[1]) && (cen[1]<(start[1]+size[1])) && (cen[2]>=start[2]) && (cen[2]<(start[2]+size[2])) )
    {
      itk::Index<3> centroid2;
      centroid2[0] = centroid[0] + cen[0] - x;
      centroid2[1] = centroid[1] + cen[1] - y;
      centroid2[2] = centroid[2] + cen[2] - z;
      soma_Table.push_back(centroid2);
    }
  }
  return soma_Table;
}


void ftkMainDarpaTrace::WriteCenterTrace(vtkSmartPointer< vtkTable > swcNodes, int x, int y, int z, std::string filename)
{
  std::cout << "Writing SWCImage file " << filename << " with " << swcNodes->GetNumberOfRows() << " nodes...";

  std::vector<int> soma_ids;
  std::vector<int> del_ids;

  for(int r=0; r<(int)swcNodes->GetNumberOfRows(); ++r)
  {
    if(swcNodes->GetValue(r,6).ToInt() == -1)
      soma_ids.push_back(swcNodes->GetValue(r,0).ToInt());
    else
      break;
  }

  for(int i=0; i<soma_ids.size(); ++i)
  {
    if( (swcNodes->GetValue(i,2).ToInt() != x) || (swcNodes->GetValue(i,3).ToInt() != y) || (swcNodes->GetValue(i,4).ToInt() != z) )//if( ((int)swcNodes[i][2] + x > 861) && ((int)swcNodes[i][2] + x < 3861) && ((int)swcNodes[i][3] + y < 8610) )
    {
      del_ids.push_back(soma_ids[i]);
    }
  }

  for(int r=0 ; r<(int)swcNodes->GetNumberOfRows(); ++r)
  {
    if(r+1 > (int)swcNodes->GetNumberOfRows()) break;
    std::vector<int>::iterator posn1 = find(del_ids.begin(), del_ids.end(), swcNodes->GetValue(r,6).ToInt());
    if(posn1 != del_ids.end())
    {
      del_ids.push_back(swcNodes->GetValue(r,0).ToInt());
      swcNodes->RemoveRow(r);
      --r;
    }
  }

  for(int i=0; i<soma_ids.size(); ++i)
  {
    if(i+1 > (int)swcNodes->GetNumberOfRows()) break;
    std::vector<int>::iterator posn1 = find(del_ids.begin(), del_ids.end(), swcNodes->GetValue(i,0).ToInt());
    if(posn1 != del_ids.end())
    {
      swcNodes->RemoveRow(i);
      --i;
    }
  }

  std::ofstream outfile(filename.c_str());

  for (int row = 0; row < (int)swcNodes->GetNumberOfRows(); ++row)
  {
    for (int col = 0; col < (int)swcNodes->GetNumberOfColumns(); ++col)
    {
      outfile << swcNodes->GetValue(row,col) << " ";
    }
    outfile << "\n";
  }
  outfile.close();
}
/**
  Calculate L Measures
 **/
void ftkMainDarpaTrace::calcLMeasures(int argc, char *argv[])
{


  std::string LMeasure_FileName = _outPath + "/LMeasures.xls";
  std::string tracesSCreenShot = _outPath + "A_trace_screenShot.tif";
  //QApplication app(argc, argv);
  //ImageFeatureThreshold *featureThreshold = new ImageFeatureThreshold();
  ////std::string inputFilename = "C:\\Data\\mnt_data\\InputImage.tif";
  //ImageViewer *imageView = new ImageViewer();
  //imageView->setImageFeatureThreshold(featureThreshold);
  //imageView->show();
  //imageView->showImage();
  //imageView->showTracesXML();
  //imageView->saveScreenShot(tracesSCreenShot);
  //imageView->saveTraceFeatures(LMeasure_FileName);
  //app.closeAllWindows();
  //app.exec();
  //app.quit();

}

float ftkMainDarpaTrace::getCalcThreshold(std::vector<float> &features, std::string type){
  std::cout<<"in getCalcThreshold"<<std::endl;
  float threshold = 0.00;
  if(type == "intensity"){
    //Const = 0
    //		var	PReg
    //		2	-0.2046
    //		4	25.9873
    //		6	-60.3184
    //		11	-640.2467
    //		12	363.7265
    //		13	987.3842
    //		14	-406.3363
    threshold = features[1]*(-0.2046)
      +features[3]*(25.9873)
      +features[5]*(-60.3184)
      +features[10]*(-640.2467)
      +features[11]*(363.7265)
      +features[12]*(987.3842)
      +features[13]*(-406.3363);

  }else if(type=="contrast"){
    //Const = 0
    //var	PReg
    //1		-0.0109
    //4		-12.4757
    //6		70.7617
    //8		-188.6227
    //9		-13.0697
    //10	330.6868
    //12	-397.2717
    //14	243.9354
    threshold = features[0]*(-0.0109)
      +features[3]*(-12.4757)
      +features[5]*(70.7617)
      +features[7]*(-188.6227)
      +features[8]*(-13.0697)
      +features[9]*(330.6868)
      +features[11]*(-397.2717)
      +features[13]*(243.9354);
  }else{
    threshold = 0.003;//Return default
  }
  return threshold;
}
std::vector<float> ftkMainDarpaTrace::computeFeatures(rawImageType_flo::Pointer &image){

  std::vector<float> features;
  float sigmas[] =  { 2.0f, 2.8284f, 4.0f, 5.6569f, 8.0f, 11.31f };	//LoG scales
  typedef float PixelType;
  typedef itk::Image< PixelType, 3 >  ImageType3D;
  typedef itk::LaplacianRecursiveGaussianImageFilter< ImageType3D , ImageType3D> GFilterType;
  typedef itk::StatisticsImageFilter<ImageType3D> StatisticsImageFilterType;


  //RescalerType::Pointer rescaler = RescalerType::New();
  //rescaler->SetOutputMinimum(0.0);
  //rescaler->SetOutputMaximum(1.0);
  //rescaler->SetInput(_PaddedCurvImage);

  ////Median filter
  //std::cout << "Running Median Filter" << std::endl;
  //MedianFilterType::Pointer medfilt = MedianFilterType::New();
  //medfilt->SetNumberOfThreads(16);
  //medfilt->SetInput(rescaler->GetOutput());

  //ImageType3D::SizeType rad = { {1, 1, 1} };
  //medfilt->SetRadius(rad);
  //medfilt->Update();
  // Image is rescaled
  ImageType3D::Pointer _PaddedCurvImage = image;

  StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New ();
  statisticsImageFilter->SetInput(_PaddedCurvImage);
  statisticsImageFilter->Update();

  features.push_back(statisticsImageFilter->GetMean());
  features.push_back(statisticsImageFilter->GetSigma());

  for (unsigned int i = 0; i < 6; ++i)
  {

    GFilterType::Pointer gauss = GFilterType::New();
    gauss->SetInput( _PaddedCurvImage );
    gauss->SetSigma( sigmas[i] );
    gauss->SetNormalizeAcrossScale(false);
    gauss->GetOutput()->Update();

    StatisticsImageFilterType::Pointer LOGstatisticsImageFilter = StatisticsImageFilterType::New ();
    LOGstatisticsImageFilter->SetInput(gauss->GetOutput());
    LOGstatisticsImageFilter->Update();
    //std::cout << "Mean: " << LOGstatisticsImageFilter->GetMean() << std::endl;
    //std::cout << "Std.: " << LOGstatisticsImageFilter->GetSigma() << std::endl;
    features.push_back(LOGstatisticsImageFilter->GetMean());
    features.push_back(LOGstatisticsImageFilter->GetSigma());
  }
  //std::cout<<"********************END computing features*********************"<<std::endl;

  return features;
}
