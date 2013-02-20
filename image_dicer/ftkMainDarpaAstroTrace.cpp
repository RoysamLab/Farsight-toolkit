
#include "ftkMainDarpaAstroTrace.h"

void ftkMainDarpaAstroTrace::readParameters( std::string astroTraceParams )
{
	_astroTraceParams = astroTraceParams;
	
	std::map< std::string, std::string > options;
	std::map< std::string, std::string >::iterator iter;
	
	options.clear();
	ifstream fin(astroTraceParams.c_str()); 
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

	iter = options.find("-xTileBor"); 
	if(iter!=options.end())
	{ std::istringstream ss((*iter).second); ss >> _xTileBor;}
	else
	{ _xTile = 1024; printf("Choose xTileBor = 1024 as default\n");}

	iter = options.find("-yTileBor"); 
	if(iter!=options.end())
	{ std::istringstream ss((*iter).second); ss >> _yTileBor;}
	else
	{ _xTile = 1024; printf("Choose yTileBor = 1024 as default\n");}

	iter = options.find("-zTileBor"); 
	if(iter!=options.end())
	{ std::istringstream ss((*iter).second); ss >> _zTileBor;}
	else
	{ _xTile = 1024; printf("Choose zTileBor = 1024 as default\n");}
	
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

	iter = options.find("-Dist_Map_Image"); 
	if(iter!=options.end())
	{ std::istringstream ss((*iter).second); ss >> _Dist_Map_Image;}
	else
	{ _Dist_Map_Image.clear(); printf("Choose _Dist_Map_Image = NULL as default\n");}
	
	//iter = options.find("-Label_Centroids"); 
	//if(iter!=options.end())
	//{ std::istringstream ss((*iter).second); ss >> _Label_Centroids;}
	//else
	//{ _Soma_Centroids.clear(); printf("Choose _Soma_Centroids = NULL as default\n");}
	
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
	
	iter = options.find("-astroTraceParams"); 
	if(iter!=options.end())
	{ std::istringstream ss((*iter).second); ss >> _astroTraceParams;}
	else
	{ _astroTraceParams.clear(); printf("Choose _astroTraceParams = NULL as default\n");}
	
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
	
	
	_Cy5_ImageNRRD = _Cy5_Image+".nrrd";
	_TRI_ImageNRRD = _TRI_Image+".nrrd";
	_GFP_ImageNRRD = _GFP_Image+".nrrd";
	_DAP_ImageNRRD = _DAP_Image+".nrrd";
	_Dist_Map_ImageNRRD = _Dist_Map_Image+".nrrd";
	_Soma_MontageNRRD = _Soma_Montage+".nrrd";
	
	
	// Print Parameters
	std::cout << std::endl << "This are the parameters";
	std::cout << std::endl << "_xSize: " << _xSize;
	std::cout << std::endl << "_ySize: " << _ySize;
	std::cout << std::endl << "_zSize: " << _zSize;
	std::cout << std::endl << "_xTile: " << _xTile;
	std::cout << std::endl << "_yTile: " << _yTile;
	std::cout << std::endl << "_zTile: " << _zTile;
	std::cout << std::endl << "_xTileBor: " << _xTileBor;
	std::cout << std::endl << "_yTileBor: " << _yTileBor;
	std::cout << std::endl << "_zTileBor: " << _zTileBor;
	std::cout << std::endl << "_num_threads: " << _num_threads;
	std::cout << std::endl << "_Cy5_Image: " << _Cy5_Image;
	std::cout << std::endl << "_TRI_Image: " << _TRI_Image;
	std::cout << std::endl << "_GFP_Image: " << _GFP_Image;
	std::cout << std::endl << "_DAP_Image: " << _DAP_Image;
	std::cout << std::endl << "_Dist_Map_Image: " << _Dist_Map_Image;
	std::cout << std::endl << "_Soma_Centroids: " << _Soma_Centroids;
	std::cout << std::endl << "_Soma_Montage: " << _Soma_Montage;
	std::cout << std::endl << "_isSmall: " << _isSmall;
	std::cout << std::endl << "_astroTraceParams: " << _astroTraceParams;
	std::cout << std::endl << "_outPath: " << _outPath;
	std::cout << std::endl << "_outPathDebug: " << _outPathDebug;
	std::cout << std::endl << "_outPathTemp: " << _outPathTemp;
}

void ftkMainDarpaAstroTrace::runPreprocesing()
{
	rawImageType_flo::Pointer MontageTRI_Image = readImage< rawImageType_flo >(_TRI_ImageNRRD.c_str());

	
	RescaleFilterType_floToflo::Pointer rescaler = RescaleFilterType_floToflo::New();
	rescaler->SetOutputMinimum(0.0);
	rescaler->SetOutputMaximum(1.0);
	rescaler->SetInput( MontageTRI_Image );

	std::cout<<std::endl<<"Rescaler"<<std::flush;

	MedianFilterType_floToflo::Pointer medfilt = MedianFilterType_floToflo::New();
	medfilt->SetInput( rescaler->GetOutput() );
	rawImageType_flo::SizeType rad = { {1, 1, 1} };
	medfilt->SetRadius(rad);
	medfilt->Update();
	
	std::cout<<std::endl<<"MedFilter"<<std::flush;
	
	_TRI_ImagePREPMNT = _outPathTemp+"/TRI_MNT_PRE.nrrd";
	writeImage< rawImageType_flo >( medfilt->GetOutput(), _TRI_ImagePREPMNT.c_str());
	
	_xSize = MontageTRI_Image->GetLargestPossibleRegion().GetSize()[0];
	_ySize = MontageTRI_Image->GetLargestPossibleRegion().GetSize()[1];
	_zSize = MontageTRI_Image->GetLargestPossibleRegion().GetSize()[2];
	
	std::cout << std::endl << "ACAASD: " << MontageTRI_Image->GetLargestPossibleRegion().GetSize();
}


void ftkMainDarpaAstroTrace::runSplitting()
{
	//if( !_Cy5_Image.empty() )
	//{
	//	rawImageType_8bit::Pointer ImageMontage_Cy5 = readImage< rawImageType_8bit >(_Cy5_ImageNRRD.c_str());
	//	_ImageMontage_Cy5Size = ImageMontage_Cy5->GetLargestPossibleRegion().GetSize();
	//	splitStore< rawImageType_8bit >( ImageMontage_Cy5, _Cy5_Image );
	//}
	if( !_TRI_Image.empty() )
	{
		rawImageType_16bit::Pointer ImageMontage_TRI = readImage< rawImageType_16bit >(_TRI_ImageNRRD.c_str());
		_ImageMontage_TRISize = ImageMontage_TRI->GetLargestPossibleRegion().GetSize();
		splitStore< rawImageType_16bit >( ImageMontage_TRI, _TRI_Image );
	}
	//if( !_GFP_Image.empty() )
	//{
	//	rawImageType_8bit::Pointer ImageMontage_GFP = readImage< rawImageType_8bit >(_GFP_ImageNRRD.c_str());
	//	_ImageMontage_GFPSize = ImageMontage_GFP->GetLargestPossibleRegion().GetSize();
	//	splitStore< rawImageType_8bit >( ImageMontage_GFP, _GFP_Image );
	//}
	//if( !_DAP_Image.empty() )
	//{
	//	rawImageType_8bit::Pointer ImageMontage_DAP = readImage< rawImageType_8bit >(_DAP_ImageNRRD.c_str());
	//	_ImageMontage_DAPSize = ImageMontage_DAP->GetLargestPossibleRegion().GetSize();
	//	splitStore< rawImageType_8bit >( ImageMontage_DAP, _DAP_Image );
	//}
	if( !_Dist_Map_Image.empty() )
	{
		rawImageType_16bit::Pointer ImageMontage_Dist_Map = readImage< rawImageType_16bit >(_Dist_Map_ImageNRRD.c_str());
		_ImageMontage_Dist_MapSize = ImageMontage_Dist_Map->GetLargestPossibleRegion().GetSize();
		splitStore< rawImageType_16bit >( ImageMontage_Dist_Map, _Dist_Map_Image );
	}
	if( !_Dist_Map_Image.empty() && !_TRI_Image.empty() )
	{
		if( (_ImageMontage_TRISize[0]!= _ImageMontage_Dist_MapSize[0]) || (_ImageMontage_TRISize[1]!= _ImageMontage_Dist_MapSize[1]) || (_ImageMontage_GFPSize[2]!= _ImageMontage_DAPSize[2]) ) 
			std::cout << std::endl << "IMAGE GFP AND DAPI ARA OF DIFFERENTE TAMANO, VOY A FALLAR";
	}
	
// 	// Compare if they are equal
// 	int flagEqual = 1;
// 	for( int i=0; i<3; ++i )
// 	{
// 		if( (_ImageMontage_Cy5Size[i] != _ImageMontage_Cy5Size[i]) || (_ImageMontage_TRISize[i] != _ImageMontage_TRISize[i]) || (_ImageMontage_GFPSize[i] != _ImageMontage_GFPSize[i]) || (_ImageMontage_DAPSize[i] != _ImageMontage_DAPSize[i]) )
// 			flagEqual = 0;
// 	}
// 	if( flagEqual == 0 )
// 		std::cout << std::endl << "THE IMAGE SIZES ARE NOT THE SAME";

	// COMPAR QUE AL MENOS DAPI OR GFP EXISTS
}


void ftkMainDarpaAstroTrace::runInterestPoints(  )
{
	
// 	if( _isSmall == 1 )
// 	{
// 		itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1); // This one can not be changed
// 		itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1); // This one can chenga
// 		omp_set_nested(1);
// 		omp_set_max_active_levels(2);
// 		int num_threads = 1;
// 		omp_set_num_threads(num_threads);
// 	}
	
	rawImageType_8bit::RegionType ImageMontageRegion;
	itk::Size<3> ImageMontageSize;
	// Assume TRI or dapi always exist
	if( !_TRI_Image.empty() )
	{
		rawImageType_8bit::Pointer ImageMontage_TRI = readImage< rawImageType_8bit >(_TRI_ImageNRRD.c_str());
		ImageMontageRegion = ImageMontage_TRI->GetLargestPossibleRegion();
		ImageMontageSize = ImageMontageRegion.GetSize();
		computeSplitConst( ImageMontage_TRI );
	}
	if( !_DAP_Image.empty() )
	{
		rawImageType_8bit::Pointer ImageMontage_DAP = readImage< rawImageType_8bit >(_DAP_ImageNRRD.c_str());
		ImageMontageRegion = ImageMontage_DAP->GetLargestPossibleRegion();
		ImageMontageSize = ImageMontageRegion.GetSize();
		computeSplitConst( ImageMontage_DAP );
	}
	
	if( _kx*_ky*_kz < _num_threads )
	{
		itk::MultiThreader::SetGlobalDefaultNumberOfThreads(int(_num_threads/_kx*_ky*_kz)); // This one can not be changed
		itk::MultiThreader::SetGlobalMaximumNumberOfThreads(int(_num_threads/_kx*_ky*_kz)); // This one can chenga
	}
	else
	{
		itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1); // This one can not be changed
		itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1); // This one can chenga
	}
		
	omp_set_nested(1);
	omp_set_max_active_levels(2);
	int num_threads = 1;
	omp_set_num_threads(num_threads);
	
	std::cout << std::endl << "fir MAX NUMBER OF CORES: " << itk::MultiThreader::GetGlobalMaximumNumberOfThreads() << ", fir DEFAULT NUMBER OF CORES: " << itk::MultiThreader::GetGlobalDefaultNumberOfThreads();
	

	int contadorSegment = 0;
	int contadorImageDoneSegment = 0;
	int contadorStich = 0;
	unsigned int maxValue = 0;
	unsigned int maxValueOld = 0;
	int flagFirstStich = 1;
	
#pragma omp parallel for collapse(3) num_threads(_num_threads) schedule(dynamic, 1)
	for(int xco = 0; xco < _kx; xco++)
	{
		for(int yco = 0; yco < _ky; yco++)
		{
			for(int zco = 0; zco < _kz; zco++)
			{
				#pragma omp critical
				{
					++contadorSegment;
					std::cout<<std::endl<< "\t\t--->>> ImageSegment " << contadorSegment << " of " << _kx*_ky*_kz;
				}
				
				rawImageType_8bit::RegionType regionLocal_inside = ComputeLocalRegionSegment( ImageMontageSize, xco, yco, zco ); // The inside
				rawImageType_uint::RegionType regionMontage_inside = ComputeGlobalRegionSegment( ImageMontageSize, xco, yco, zco ); // The inside
				
				rawImageType_8bit::RegionType regionLocal_all = ComputeLocalRegionSplit( ImageMontageSize, xco, yco, zco );
				rawImageType_8bit::RegionType regionMontage_all = ComputeGlobalRegionSplit( ImageMontageSize, xco, yco, zco );

				// Store local image
				std::stringstream out_x;
				std::stringstream out_y;
				std::stringstream out_z;
				out_x<<xco;
				out_y<<yco;
				out_z<<zco;
				std::string xStr = out_x.str();
				std::string yStr = out_y.str();
				std::string zStr = out_z.str();
				
// 				std::cout << std::endl << "REGION_LOC: " << regionLocal_inside;
// 				std::cout << std::endl << "REGION_MON: " << regionMontage_inside;
				
				//rawImageType_8bit::Pointer imageLocalCy5;
				//rawImageType_8bit::Pointer imageLocalTRI;
				//rawImageType_8bit::Pointer imageLocalGFP;
				//rawImageType_8bit::Pointer imageLocalDAP;
				rawImageType_flo::Pointer imageLocalTRI;
				rawImageType_16bit::Pointer imageLocalDist_Map;

				std::vector< rawImageType_flo::Pointer > Images_Tiles;
				Images_Tiles.resize(1);

				std::vector< rawImageType_16bit::Pointer > Dist_Map_Tiles;
				Dist_Map_Tiles.resize(1);

				std::vector< rawImageType_16bit::Pointer > Label_Tiles;
				Label_Tiles.resize(1);
				
				std::vector< vtkSmartPointer< vtkTable > > Table_Tiles;
				Table_Tiles.resize(1);

				std::vector< vtkSmartPointer< vtkTable > > feature_Vector_Tables;
				feature_Vector_Tables.resize(1);
				
				std::vector< rawImageType_16bit::Pointer > ID_Images;
				ID_Images.resize(1);

				std::vector< std::map< unsigned int, itk::Index<3> > > Centroids_Tiles;
				Centroids_Tiles.resize(1);
				
				// Reading part
				#pragma omp critical
				{
				if( !_TRI_Image.empty( ) /*&& !_DAP_Image.empty()*/ )
				{
					int foundTRI = _TRI_Image.find_last_of("/\\");
					std::string _TRI_ImageNoPath = _TRI_Image.substr(foundTRI+1);
					// !!! this has to be nrrd is just to test
					std::string tempTRI = _outPathTemp+"/"+_TRI_ImageNoPath+"_"+xStr+"_"+yStr+"_"+zStr+".nrrd";
					imageLocalTRI = readImage< rawImageType_flo >(tempTRI.c_str());
					
					//int foundDAP = _DAP_Image.find_last_of("/\\");
					//std::string _DAP_ImageNoPath = _DAP_Image.substr(foundDAP+1);
					//// !!! this has to be nrrd is just to test
					//std::string tempDAP = _outPathTemp+"/"+_DAP_ImageNoPath+"_"+xStr+"_"+yStr+"_"+zStr+".nrrd";
					//imageLocalDAP = readImage< rawImageType_8bit >(tempDAP.c_str());
				}
				if( !_Dist_Map_Image.empty( ) )
				{
					int foundDist_Map = _Dist_Map_Image.find_last_of("/\\");
					std::string _Dist_Map_ImageNoPath = _Dist_Map_Image.substr(foundDist_Map+1);
					// !!! this has to be nrrd is just to test
					std::string tempDist_Map = _outPathTemp+"/"+_Dist_Map_ImageNoPath+"_"+xStr+"_"+yStr+"_"+zStr+".nrrd";
					imageLocalDist_Map = readImage< rawImageType_16bit >(tempDist_Map.c_str());
				}
				}
				//Images_Tiles[0] = imageLocalCy5;
				//Images_Tiles[1] = imageLocalTRI;
				//Images_Tiles[2] = imageLocalGFP;
				//Images_Tiles[3] = imageLocalDAP;
				Images_Tiles[0] = imageLocalTRI;

				Dist_Map_Tiles[0] = imageLocalDist_Map;
								
// 				RunSegmentation(regionLocal_inside, Images_Tiles, Label_Tiles, Table_Tiles, Centroids_Tiles);
				
				//Label_Tiles[0] = RunNuclearSegmentation( Images_Tiles[3] );
				//Table_Tiles[0] = ComputeFeaturesAndAssociations( Images_Tiles, Label_Tiles );
				//Centroids_Tiles[0] = GetLabelToCentroidMap(Table_Tiles[0]);

				AstroTracer * AT = new AstroTracer();
				AT->LoadCurvImage(Images_Tiles[0], 0);
				//AT->LoadSomaImage(std::string(argv[2]));
				std::string coverageFileName = _outPathTemp+"/coverage_"+xStr+"_"+yStr+"_"+zStr+".nrrd";
				AT->OptimizeCoverage(coverageFileName, true);
				AT->LoadParameters(_astroTraceParams.c_str());	
				AT->SetScaleRange(4, 4); //(2, 5); //(2, 2)
				AT->CallFeatureMainExternal();
				AT->Set_DistanceMapImage(imageLocalDist_Map);
				std::string featureVectorFileName = _outPathTemp+"/featureVector_"+xStr+"_"+yStr+"_"+zStr+".nrrd";
				std::string IDImageFileName = _outPathTemp+"/IDImage_"+xStr+"_"+yStr+"_"+zStr+".nrrd";
				AT->ComputeAstroFeaturesPipeline(featureVectorFileName, IDImageFileName, 0, regionLocal_inside, feature_Vector_Tables, ID_Images, true);
						

				RemoveLabelNearBorder(regionLocal_inside, Label_Tiles, Table_Tiles, Centroids_Tiles );
				
					
					
// 					// TEST SAVE SEGMENTATION WITH LABELS REMOVED
				#pragma omp critical
				{
					std::string tempTABLERE = _outPathTemp+"/segTable_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REMO.txt";
					ftk::SaveTable(tempTABLERE, Table_Tiles[0]);
					std::string tempLABELRE = _outPathTemp+"/segLabel_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REMO.nrrd";
					writeImage<rawImageType_16bit>(Label_Tiles[0],tempLABELRE.c_str());

	// 				std::string tempLABELRE = _outPathDebugLevel2+"/segLabel_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REMO.tif";
					ftkMainDarpa objftkMainDarpa;
					objftkMainDarpa.projectImage<rawImageType_16bit, rawImageType_16bit>( Label_Tiles[0], tempLABELRE, _outPathDebugLevel2, "ORG_RES_BIN", "TIFF" );
				}
				#pragma omp critical
				{
					contadorImageDoneSegment++;
					std::cout<<std::endl<< "\t\t--->>> ImageDoneSegment " << contadorImageDoneSegment << " of " << _kx*_ky*_kz;
				}
			}
		}
	}
	k
}


//void ftkMainDarpaAstroTrace::runTracing()
//{
//	std::cout << std::endl << "LETS TRACE";
//	
//	if( _isSmall == 1 )
//	{
//		itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1); // This one can not be changed
//		itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1); // This one can chenga
//		omp_set_nested(1); // For the crop
//		omp_set_max_active_levels(2); // For the crop
//		int num_threads = 1;
//		omp_set_num_threads(num_threads);
//	}
//
//	int counterCentro = 0;
//	std::vector< itk::Index<3> > centroid_list = getCentroidList();
//	std::cout << "Number of cells to be traced : " << centroid_list.size() << "\n";
//
//	std::string SWCFilename = _outPath + "/OnlySWC.xml";
//	std::ofstream outSWCFile;
//	outSWCFile.open(SWCFilename.c_str());
//	outSWCFile << "<?xml\tversion=\"1.0\"\t?>\n";
//	outSWCFile << "<Source>\n\n";
//		
//	computeSplitConst();
//		
//	// CAN NOT BE RUN PARALLEL
//	for( unsigned int bigTile = 0; bigTile<_numDivisionsInRowCEN ; ++bigTile )
//	{
//			std::cout<<std::endl<<"bigTile: "<<bigTile;
//			std::cout << std::endl << "_initialBigTileLOG: " << _initialBigTileLOG[bigTile][0] <<", "<<_initialBigTileLOG[bigTile][1] <<", "<<_initialBigTileLOG[bigTile][2];
//			std::cout << std::endl << "_sizeOfBigTilesLOG: " << _sizeOfBigTilesLOG[bigTile][0] <<", "<<_sizeOfBigTilesLOG[bigTile][1] <<", "<<_sizeOfBigTilesLOG[bigTile][2];
//		
//// 			stringstream out34;
//// 			out34<<bigTile;
//// 			string srr = out34.str();
//// 			std::string SWCFilenameDivided = _outPath + "/TracesAndSomasDivided/OnlySWC_" + srr +".xml";
//// 			std::ofstream outfileDivided;
//// 			outfileDivided.open(SWCFilenameDivided.c_str());
//// 			outfileDivided << "<?xml\tversion=\"1.0\"\t?>\n";
//// 			outfileDivided << "<Source>\n\n";
//
//		itk::Index<3> initialBigIndexLOG;
//		itk::Size<3> sizeOfTheRegionLOG;
//		
//		initialBigIndexLOG[0] = _initialBigTileLOG[bigTile][0];
//		initialBigIndexLOG[1] = _initialBigTileLOG[bigTile][1];
//		initialBigIndexLOG[2] = _initialBigTileLOG[bigTile][2];
//		
//		sizeOfTheRegionLOG[0] = _sizeOfBigTilesLOG[bigTile][0];
//		sizeOfTheRegionLOG[1] = _sizeOfBigTilesLOG[bigTile][1]+_sizeOfBigTilesLOG[bigTile+1][1];
//		sizeOfTheRegionLOG[2] = _sizeOfBigTilesLOG[bigTile][2];
//		
//		std::cout << std::endl << "sizeOfTheRegionLOG: " << sizeOfTheRegionLOG[0] <<", "<<sizeOfTheRegionLOG[1] <<", "<<sizeOfTheRegionLOG[2];
//		
//// 		std::cout << std::endl << initialBigIndexLOG;
//// 		std::cout << std::endl << "SIZE: " << sizeOfTheRegionLOG[1];
//		
//		rawImageType_flo::RegionType desiredRegionBigTileLOG;
//		desiredRegionBigTileLOG.SetSize(sizeOfTheRegionLOG);
//		desiredRegionBigTileLOG.SetIndex(initialBigIndexLOG);
//		
//		itk::Index<3> initialBigIndexCEN;
//		itk::Size<3> sizeOfTheRegionCEN;
//		
//		initialBigIndexCEN[0] = _initialBigTileCEN[bigTile][0];
//		initialBigIndexCEN[1] = _initialBigTileCEN[bigTile][1];
//		initialBigIndexCEN[2] = _initialBigTileCEN[bigTile][2];
//		
//		sizeOfTheRegionCEN[0] = _sizeOfBigTilesCEN[bigTile][0];
//		sizeOfTheRegionCEN[1] = _sizeOfBigTilesCEN[bigTile][1];
//		sizeOfTheRegionCEN[2] = _sizeOfBigTilesCEN[bigTile][2];
//		
//		std::cout << std::endl << "_sizeOfBigTilesCEN: " << _sizeOfBigTilesCEN[bigTile][0] <<", "<<_sizeOfBigTilesCEN[bigTile][1] <<", "<<_sizeOfBigTilesCEN[bigTile][2];
//		std::cout << std::endl << desiredRegionBigTileLOG;
//		
//// 			std::vector< rawImageType_flo::Pointer > LoGDesiredRegion;
//// 			LoGDesiredRegion.resize(6);
//		
//// 			rawImageType_uint::Pointer _somaMontageDesiredRegion;
//// 			rawImageType_flo::Pointer _img_traceDesiredRegion;
//		
//// 			std::string tempFileName_42 = _GFP_ImagePREPMNT
//		_img_traceDesiredRegion = readImageRegion< rawImageType_flo >( _TRI_ImagePREPMNT.c_str(), desiredRegionBigTileLOG );
//		_somaMontageDesiredRegion = readImageRegion< rawImageType_uint >( _Soma_MontageNRRD.c_str(), desiredRegionBigTileLOG );
//
//	#pragma omp parallel for num_threads(_num_threads) schedule(dynamic, 1)
//		for( unsigned long long a=0; a<centroid_list.size(); ++a )
//		{
//			int x, y, z;
//			std::stringstream ssx, ssy, ssz;
//
//			x = centroid_list[a][0];
//			y = centroid_list[a][1];
//			z = centroid_list[a][2];
//			
//			unsigned long long yMin = initialBigIndexCEN[1];
//			unsigned long long yMax = initialBigIndexCEN[1] + sizeOfTheRegionCEN[1];
//
//			if(!( (yMin <= y) && (y < yMax) ) )
//			{
//				continue;
//			}
//			
//			ssx << x; ssy << y; ssz << z;
//			
//	#pragma omp critical
//			{
//				counterCentro++;
//				std::cout<<std::endl<<"\t\t\t\t asdfasdf ----->>>>> " << counterCentro << ", of " << centroid_list.size();
//				std::cout<<", x: "<<centroid_list[a][0]<<", y: "<<centroid_list[a][1]<<", z: "<<centroid_list[a][2]<<std::endl;
//			}
//
//// 			std::cout << std::endl << "x: " << x << ", y: " << y << ", z: " << z;
//// 			std::cout << std::endl << "xTile: " << _xTile << ", yTile: " << _yTile << ", zTile: " << _zTile;
//
//			std::stringstream ssx_off, ssy_off, ssz_off;		
//			std::stringstream ssx_offBig, ssy_offBig, ssz_offBig;
//			ssx_offBig << 0;
//			ssy_offBig << 0;
//			ssz_offBig << 0;
//			if(x >= _xTile/2)
//				ssx_off << x - _xTile/2;
//			else 
//				ssx_off << 0;
//			if(y >= _yTile/2)
//				ssy_off << y - _yTile/2;
//			else 
//				ssy_off << 0;
//			if(z >= _zTile/2)
//				ssz_off << z - _zTile/2;
//			else 
//				ssz_off << 0;
//
//// 			std::cout << std::endl << "_initialBigTileLOG: " << _initialBigTileLOG[bigTile][1];
//			
//			int x_local = x;
//			int y_local = y - _initialBigTileLOG[bigTile][1];
//			int z_local = z;
//			
//// 			std::cout << std::endl << "x_local: " << x_local << ", y_local: " << y_local << ", z_local: " << z_local;
//// 			std::cout << std::endl << "x: " << x << ", y: " << y << ", z: " << z;
//			
//			if(x_local >= _xTile/2)
//				ssx_offBig << x_local - _xTile/2;
//			else 
//				ssx_off << 0;
//			if(y_local >= _yTile/2)
//				ssy_offBig << y_local - _yTile/2;
//			else 
//				ssy_off << 0;
//			if(z_local >= _zTile/2)
//				ssz_offBig << z_local - _zTile/2;
//			else 
//				ssz_off << 0;
//
//			//########    CROP THE DESIRED DICE FROM THE GFP AND SOMA MONTAGES   ########
//
//			
//			//########    FETCH ALL CENTROIDS THAT FALL WITHIN THE DICE    ########
//
//			std::vector< itk::Index<3> > soma_Table = getSomaTable(centroid_list, x, y, z );
//			
//// 			//########    RUN TRACING    ########
//			MultipleNeuronTracer * MNT = new MultipleNeuronTracer();
//			MNT->LoadParameters(_astroTraceParams.c_str(),5);
//// 				
//// 				MNT->LoadCurvImage_1(img_trace, 0);
//			#pragma omp critical
//			{
//// 				std::cout << std::endl << "LAREGION ES: " << _img_traceDesiredRegion;
//				rawImageType_flo::Pointer img_trace = cropImages< rawImageType_flo >( _img_traceDesiredRegion, x, y, z);
//				MNT->LoadCurvImage_2(img_trace);
//			}
//			MNT->ReadStartPoints_1(soma_Table, 0);
//// 				MNT->SetCostThreshold(1000);
//			MNT->SetCostThreshold(MNT->cost_threshold);
//			
//// 	// 			MNT->LoadSomaImage_1(img_soma_yan);
//			bool flagLog = false;
//			MNT->setFlagOutLog(flagLog);
//			MNT->RunTracing();
//			
//			#pragma omp critical
//			{
//				rawImageType_uint::Pointer img_soma = cropImages< rawImageType_uint >( _somaMontageDesiredRegion, x, y, z);
//				MNT->RemoveSoma( img_soma );
//			}
//// 
//			x = min(_xTile/2, x);
//			y = min(_yTile/2, y);
//			z = min(_zTile/2, z);
////
//			vtkSmartPointer< vtkTable > swcTable = MNT->GetSWCTable(0);
//			delete MNT;
//			std::string swcFilename = _outPath + "/Trace_" + ssx.str() + "_" + ssy.str() + "_" + ssz.str() + "_ANT.swc";
//			WriteCenterTrace(swcTable, x, y, z, swcFilename);
//// 				
//			#pragma omp critical
//			{
//				outSWCFile << "\t<File\tFileName=\"Trace_" << ssx.str() << "_" << ssy.str() << "_" << ssz.str() << "_ANT.swc\"\tType=\"Trace\"\ttX=\"" << ssx_off.str() << "\"\ttY=\"" << ssy_off.str() << "\"\ttZ=\"" << ssz_off.str() << "\"/>\n";
//			}
//			
//// // 				vtkSmartPointer< vtkTable > swcTable = MNT->GetSWCTable(0);
//// 				std::string swcFilenameDivided = _outPath + "/TracesAndSomasDivided/Trace_BigTile_" + srr + "_" + ssx.str() + "_" + ssy.str() + "_" + ssz.str() + "_ANT.swc";
//// 				WriteCenterTrace(swcTable, x, y, z, swcFilenameDivided);
//// // 				
//// 				#pragma omp critical
//// 				{
//// 					outfileDivided << "\t<File\tFileName=\"Trace_BigTile_" << srr << "_" << ssx.str() << "_" << ssy.str() << "_" << ssz.str() << "_ANT.swc\"\tType=\"Trace\"\ttX=\"" << ssx_offBig.str() << "\"\ttY=\"" << ssy_offBig.str() << "\"\ttZ=\"" << ssz_offBig.str() << "\"/>\n";
//// 				}
//		}
//	}
//}

	
	
	
	



void ftkMainDarpaAstroTrace::computeSplitConst( rawImageType_8bit::Pointer ImageMontage )
{
	itk::Size<3> ImageMontageSize = ImageMontage->GetLargestPossibleRegion().GetSize();

	_kx = ImageMontageSize[0] /(_xTile-_xTileBor);
	_ky = ImageMontageSize[1] /(_yTile-_yTileBor);
	_kz = ImageMontageSize[2] /(_zTile-_zTileBor);
// 	std::cout << std::endl << _kx << " " << _ky << " " << _kz;

	int remx = ImageMontageSize[0] % (_xTile-_xTileBor);
	int remy = ImageMontageSize[1] % (_yTile-_yTileBor);
	int remz = ImageMontageSize[2] % (_zTile-_zTileBor);
// 	std::cout << std::endl << remx << " " << remy << " " << remz;

	itk::Size<3> size_test;
	size_test[0] =  MINNIC((_kx)*(_xTile-_xTileBor)+_xTile-1,ImageMontageSize[0]-1) -  (_kx) * (_xTile-_xTileBor) +1;
	size_test[1] =  MINNIC((_ky)*(_yTile-_yTileBor)+_yTile-1,ImageMontageSize[1]-1) -  (_ky) * (_yTile-_yTileBor) +1;
	size_test[2] =  MINNIC((_kz)*(_zTile-_zTileBor)+_zTile-1,ImageMontageSize[2]-1) -  (_kz) * (_zTile-_zTileBor) +1;
// 	std::cout << std::endl << size_test[0] << " " << size_test[1] << " " << size_test[2];

	if( size_test[0] > _xTileBor  )
	{
		if ( remx > 0 )
			_kx ++;
	}
	if( size_test[1] > _yTileBor  )
	{
		if ( remy > 0 )
			_ky ++;
	}
	if( size_test[2] > _zTileBor  )
	{
		if ( remz > 0 )
			_kz ++;
	}
	
	std::cout << std::endl << size_test[0] << " " << size_test[1] << " " << size_test[2];
	std::cout << std::endl << _kx << " " << _ky << " " << _kz;
	std::cout << std::endl << remx << " " << remy << " " << remz;
}

rawImageType_8bit::RegionType ftkMainDarpaAstroTrace::ComputeLocalRegionSegment( itk::Size<3> ImageMontageSize, int xco, int yco, int zco )
{
	rawImageType_8bit::SizeType size;
	size[0] =  MINNIC((xco)*(_xTile-_xTileBor)+_xTile-1,ImageMontageSize[0]-1) -  xco * (_xTile-_xTileBor) +1;
	size[1] =  MINNIC((yco)*(_yTile-_yTileBor)+_yTile-1,ImageMontageSize[1]-1) -  yco * (_yTile-_yTileBor) +1;
	size[2] =  MINNIC((zco)*(_zTile-_zTileBor)+_zTile-1,ImageMontageSize[2]-1) -  zco * (_zTile-_zTileBor) +1;
	
// 				std::cout << std::endl << "INDICE: " << xco << " " << yco << " " << zco;
// 				std::cout << std::endl << "Tamano Entrada: " << size[0] << " " << size[1] << " " << size[2];

	if( size[0] <= _xTileBor/2 )
		std::cout<<std::endl<<"IT SHOULD NEVER ENTER HERE";
	if( size[1] <= _yTileBor/2 )
		std::cout<<std::endl<<"IT SHOULD NEVER ENTER HERE";
	if( size[2] <= _zTileBor/2 )
		std::cout<<std::endl<<"IT SHOULD NEVER ENTER HERE";

	if((xco != 0) )
		size[0] = size[0] - _xTileBor/2;
	if(xco != _kx-1)
		size[0] = size[0] - _xTileBor/2;
	
	if((yco != 0) )
		size[1] = size[1] - _yTileBor/2;
	if(yco != _ky-1)
		size[1] = size[1] - _yTileBor/2;
	
	if((zco != 0) )
		size[2] = size[2] - _zTileBor/2;
	if(zco != _kz-1)
		size[2] = size[2] - _zTileBor/2;
// 
// 				std::cout << std::endl << "Tamano Salida: " << size[0] << " " << size[1] << " " << size[2];
// 
	rawImageType_8bit::IndexType indexLocal_1;
// 				rawImageType_8bit::IndexType indexMontage_1;
// 					
	rawImageType_8bit::RegionType regionLocal_inside;
// 				rawImageType_8bit::RegionType regionMontage_inside;
// 			
	indexLocal_1[0] = xco *(_xTile-_xTileBor);
// 				indexMontage_1[0] = xco *(_xTile-_xTileBor);
	indexLocal_1[1] = yco *(_yTile-_yTileBor);
// 				indexMontage_1[1] = yco *(_yTile-_yTileBor);
	indexLocal_1[2] = zco *(_zTile-_zTileBor);
// 				indexMontage_1[2] = zco *(_zTile-_zTileBor);
// 					
	if((xco != 0) )
		indexLocal_1[0] = _xTileBor/2;
	if((yco != 0) )
		indexLocal_1[1] = _yTileBor/2;
	if((zco != 0) )
		indexLocal_1[2] = _zTileBor/2;
// 					
	regionLocal_inside.SetIndex(indexLocal_1);
	regionLocal_inside.SetSize(size);
	
	return regionLocal_inside;
// 					
// 				if(xco!=0)
// 					indexMontage_1[0] = xco *(_xTile-_xTileBor)+_xTileBor/2;
// 				if(yco!=0)
// 					indexMontage_1[1] = yco *(_yTile-_yTileBor)+_yTileBor/2;
// 				if(zco!=0)
// 					indexMontage_1[2] = zco *(_zTile-_zTileBor)+_zTileBor/2;
// // 					
// 				regionMontage_inside.SetIndex(indexMontage_1);
// 				regionMontage_inside.SetSize(size);
}

rawImageType_uint::RegionType ftkMainDarpaAstroTrace::ComputeGlobalRegionSegment( itk::Size<3> ImageMontageSize, int xco, int yco, int zco )
{
	rawImageType_8bit::SizeType size;
	size[0] =  MINNIC((xco)*(_xTile-_xTileBor)+_xTile-1,ImageMontageSize[0]-1) -  xco * (_xTile-_xTileBor) +1;
	size[1] =  MINNIC((yco)*(_yTile-_yTileBor)+_yTile-1,ImageMontageSize[1]-1) -  yco * (_yTile-_yTileBor) +1;
	size[2] =  MINNIC((zco)*(_zTile-_zTileBor)+_zTile-1,ImageMontageSize[2]-1) -  zco * (_zTile-_zTileBor) +1;
	
// 				std::cout << std::endl << "INDICE: " << xco << " " << yco << " " << zco;
// 				std::cout << std::endl << "Tamano Entrada: " << size[0] << " " << size[1] << " " << size[2];

	if( size[0] <= _xTileBor/2 )
		std::cout<<std::endl<<"IT SHOULD NEVER ENTER HERE";
	if( size[1] <= _yTileBor/2 )
		std::cout<<std::endl<<"IT SHOULD NEVER ENTER HERE";
	if( size[2] <= _zTileBor/2 )
		std::cout<<std::endl<<"IT SHOULD NEVER ENTER HERE";

	if((xco != 0) )
		size[0] = size[0] - _xTileBor/2;
	if(xco != _kx-1)
		size[0] = size[0] - _xTileBor/2;
	
	if((yco != 0) )
		size[1] = size[1] - _yTileBor/2;
	if(yco != _ky-1)
		size[1] = size[1] - _yTileBor/2;
	
	if((zco != 0) )
		size[2] = size[2] - _zTileBor/2;
	if(zco != _kz-1)
		size[2] = size[2] - _zTileBor/2;
// 
// 				std::cout << std::endl << "Tamano Salida: " << size[0] << " " << size[1] << " " << size[2];
// 
// 				rawImageType_8bit::IndexType indexLocal_1;
	rawImageType_8bit::IndexType indexMontage_1;
// 					
// 				rawImageType_8bit::RegionType regionLocal_inside;
	rawImageType_8bit::RegionType regionMontage_inside;
// 			
// 				indexLocal_1[0] = xco *(_xTile-_xTileBor);
	indexMontage_1[0] = xco *(_xTile-_xTileBor);
// 				indexLocal_1[1] = yco *(_yTile-_yTileBor);
	indexMontage_1[1] = yco *(_yTile-_yTileBor);
// 				indexLocal_1[2] = zco *(_zTile-_zTileBor);
	indexMontage_1[2] = zco *(_zTile-_zTileBor);
// 					
// 				if((xco != 0) )
// 					indexLocal_1[0] = _xTileBor/2;
// 				if((yco != 0) )
// 					indexLocal_1[1] = _yTileBor/2;
// 				if((zco != 0) )
// 					indexLocal_1[2] = _zTileBor/2;
// // 					
// 				regionLocal_inside.SetIndex(indexLocal_1);
// 				regionLocal_inside.SetSize(size);
	
// 				return regionLocal_inside;
// 					
	if(xco!=0)
		indexMontage_1[0] = xco *(_xTile-_xTileBor)+_xTileBor/2;
	if(yco!=0)
		indexMontage_1[1] = yco *(_yTile-_yTileBor)+_yTileBor/2;
	if(zco!=0)
		indexMontage_1[2] = zco *(_zTile-_zTileBor)+_zTileBor/2;
// 					
	regionMontage_inside.SetIndex(indexMontage_1);
	regionMontage_inside.SetSize(size);

	return regionMontage_inside;
}

rawImageType_8bit::RegionType ftkMainDarpaAstroTrace::ComputeLocalRegionSplit( itk::Size<3> ImageMontageSize, int xco, int yco, int zco )
{
	rawImageType_8bit::SizeType size;
	size[0] =  MINNIC((xco)*(_xTile-_xTileBor)+_xTile-1,ImageMontageSize[0]-1) -  xco * (_xTile-_xTileBor) +1;
	size[1] =  MINNIC((yco)*(_yTile-_yTileBor)+_yTile-1,ImageMontageSize[1]-1) -  yco * (_yTile-_yTileBor) +1;
	size[2] =  MINNIC((zco)*(_zTile-_zTileBor)+_zTile-1,ImageMontageSize[2]-1) -  zco * (_zTile-_zTileBor) +1;

	if( size[0] <= _xTileBor/2 )
		std::cout<<std::endl<<"IT SHOULD NEVER ENTER HERE";
	if( size[1] <= _yTileBor/2 )
		std::cout<<std::endl<<"IT SHOULD NEVER ENTER HERE";
	if( size[2] <= _zTileBor/2 )
		std::cout<<std::endl<<"IT SHOULD NEVER ENTER HERE";

	rawImageType_8bit::IndexType index;
	index.Fill(0);
	rawImageType_8bit::RegionType region;
	region.SetIndex(index);
	region.SetSize(size);
	return region;
}

rawImageType_8bit::RegionType ftkMainDarpaAstroTrace::ComputeGlobalRegionSplit( itk::Size<3> ImageMontageSize, int xco, int yco, int zco )
{
	rawImageType_8bit::SizeType size;
	size[0] =  MINNIC((xco)*(_xTile-_xTileBor)+_xTile-1,ImageMontageSize[0]-1) -  xco * (_xTile-_xTileBor) +1;
	size[1] =  MINNIC((yco)*(_yTile-_yTileBor)+_yTile-1,ImageMontageSize[1]-1) -  yco * (_yTile-_yTileBor) +1;
	size[2] =  MINNIC((zco)*(_zTile-_zTileBor)+_zTile-1,ImageMontageSize[2]-1) -  zco * (_zTile-_zTileBor) +1;

	if( size[0] <= _xTileBor/2 )
		std::cout<<std::endl<<"IT SHOULD NEVER ENTER HERE";
	if( size[1] <= _yTileBor/2 )
		std::cout<<std::endl<<"IT SHOULD NEVER ENTER HERE";
	if( size[2] <= _zTileBor/2 )
		std::cout<<std::endl<<"IT SHOULD NEVER ENTER HERE";
	
	rawImageType_8bit::IndexType index;
	rawImageType_8bit::RegionType regionMontage;
	index[0] = xco *(_xTile-_xTileBor);
	index[1] = yco *(_yTile-_yTileBor);
	index[2] = zco *(_zTile-_zTileBor);
	regionMontage.SetIndex(index);
	regionMontage.SetSize(size);
	
	return regionMontage;
}


//
//std::vector< itk::Index<3> > ftkMainDarpaAstroTrace::getCentroidList()
//{	
//	std::cout << std::endl << _Soma_Centroids;
//	std::cout << std::endl << _Soma_Centroids;
//	std::cout << std::endl << _Soma_Centroids;
//	
//	vtkSmartPointer<vtkTable> somaCentroidsTable = ftk::LoadTable(_Soma_Centroids);
//
//	std::vector< itk::Index<3> > centroid_list;
//	for(int r=0; r<(int)somaCentroidsTable->GetNumberOfRows(); ++r)
//	{
//		int cx = somaCentroidsTable->GetValue(r, 0).ToInt();
//		int cy = somaCentroidsTable->GetValue(r, 1).ToInt();
//		int cz = somaCentroidsTable->GetValue(r, 2).ToInt();
//		
//		itk::Index<3> cen;
//		cen[0] = cx; cen[1] = cy; cen[2] = cz; 
//		centroid_list.push_back(cen);
//	}
//// 	std::cout << std::endl << "-----------------END";
//// 	std::cout << std::endl << "-----------------END";
//	return centroid_list;
//}
//
//std::vector< itk::Index<3> > ftkMainDarpaAstroTrace::getSomaTable( std::vector< itk::Index<3> > centroid_list, int x, int y, int z )
//{
//	itk::Index<3> centroid;
//	centroid[0] = ((x - _xTile/2)>0) ? _xTile/2:x; 
//	centroid[1] = ((y - _yTile/2)>0) ? _yTile/2:y;
//	centroid[2] = ((z - _zTile/2)>0) ? _zTile/2:z;
//	
//	itk::Index<3> start;
//	start[0] = ((x - _xTile/2)>0) ? (x - _xTile/2):0;
//	start[1] = ((y - _yTile/2)>0) ? (y - _yTile/2):0;
//	start[2] = ((z - _zTile/2)>0) ? (z - _zTile/2):0;
//
//	itk::Size<3> size;
//	size[0] = ((x+_xTile/2)<_xSize) ? _xTile : (_xTile/2+_xSize-x-1); 
//	size[1] = ((y+_yTile/2)<_ySize) ? _yTile : (_yTile/2+_ySize-y-1);
//	size[2] = ((z+_zTile/2)<_zSize) ? _zTile : (_zTile/2+_zSize-z-1);
//	
//	std::vector< itk::Index<3> > soma_Table;      
//	for(int ctr =0; ctr<centroid_list.size() ; ++ctr)
//	{
//		itk::Index<3> cen = centroid_list[ctr];
//		if( (cen[0]>=start[0]) && (cen[0]<(start[0]+size[0])) && (cen[1]>=start[1]) && (cen[1]<(start[1]+size[1])) && (cen[2]>=start[2]) && (cen[2]<(start[2]+size[2])) )
//		{
//			itk::Index<3> centroid2;
//			centroid2[0] = centroid[0] + cen[0] - x;
//			centroid2[1] = centroid[1] + cen[1] - y;
//			centroid2[2] = centroid[2] + cen[2] - z;
//			soma_Table.push_back(centroid2);
//		}
//	}
//	return soma_Table;
//}
//
//
//void ftkMainDarpaAstroTrace::WriteCenterTrace(vtkSmartPointer< vtkTable > swcNodes, int x, int y, int z, std::string filename)
//{
//	std::cout << "Writing SWCImage file " << filename << " with " << swcNodes->GetNumberOfRows() << " nodes...";
//
//	std::vector<int> soma_ids;
//	std::vector<int> del_ids;
//
//	for(int r=0; r<(int)swcNodes->GetNumberOfRows(); ++r)
//	{
//		if(swcNodes->GetValue(r,6).ToInt() == -1)
//			soma_ids.push_back(swcNodes->GetValue(r,0).ToInt());
//		else
//			break;
//	}
//
//	for(int i=0; i<soma_ids.size(); ++i)
//	{
//		if( (swcNodes->GetValue(i,2).ToInt() != x) || (swcNodes->GetValue(i,3).ToInt() != y) || (swcNodes->GetValue(i,4).ToInt() != z) )//if( ((int)swcNodes[i][2] + x > 861) && ((int)swcNodes[i][2] + x < 3861) && ((int)swcNodes[i][3] + y < 8610) )
//		{
//			del_ids.push_back(soma_ids[i]);
//		}
//	}
//
//	for(int r=0 ; r<(int)swcNodes->GetNumberOfRows(); ++r)
//	{
//		if(r+1 > (int)swcNodes->GetNumberOfRows()) break;
//		std::vector<int>::iterator posn1 = find(del_ids.begin(), del_ids.end(), swcNodes->GetValue(r,6).ToInt());
//		if(posn1 != del_ids.end())
//		{
//			del_ids.push_back(swcNodes->GetValue(r,0).ToInt());
//			swcNodes->RemoveRow(r);
//			--r;
//		}
//	}
//
//	for(int i=0; i<soma_ids.size(); ++i)
//	{
//		if(i+1 > (int)swcNodes->GetNumberOfRows()) break;
//		std::vector<int>::iterator posn1 = find(del_ids.begin(), del_ids.end(), swcNodes->GetValue(i,0).ToInt());
//		if(posn1 != del_ids.end())
//		{
//			swcNodes->RemoveRow(i);
//			--i;
//		}
//	}
//
//	std::ofstream outfile(filename.c_str());
//
//	for (int row = 0; row < (int)swcNodes->GetNumberOfRows(); ++row) 
//	{
//		for (int col = 0; col < (int)swcNodes->GetNumberOfColumns(); ++col) 
//		{
//			outfile << swcNodes->GetValue(row,col) << " ";
//		}
//		outfile << "\n";
//	}
//	outfile.close();
//}