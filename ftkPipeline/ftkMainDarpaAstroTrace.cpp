
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

	iter = options.find("-Label_Image"); 
	if(iter!=options.end())
	{ std::istringstream ss((*iter).second); ss >> _Label_Image;}
	else
	{ _Label_Image.clear(); printf("Choose _Label_Image = NULL as default\n");}

	iter = options.find("-Nuclei_Table"); 
	if(iter!=options.end())
	{ std::istringstream ss((*iter).second); ss >> _Nuclei_Table;}
	else
	{ _Nuclei_Table.clear(); printf("Choose _Nuclei_Table = NULL as default\n");}

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

	iter = options.find("-roots_model_AL"); 
	if(iter!=options.end())
	{ std::istringstream ss((*iter).second); ss >> _roots_model_AL;}
	else
	{ _roots_model_AL.clear(); printf("Choose _roots_model_AL = NULL as default\n");}

	iter = options.find("-final_classification_model"); 
	if(iter!=options.end())
	{ std::istringstream ss((*iter).second); ss >> _final_classification_model;}
	else
	{ _final_classification_model.clear(); printf("Choose _final_classification_model = NULL as default\n");}
	
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
	_Label_ImageNRRD = _Label_Image+".nrrd";
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
	std::cout << std::endl << "_Label_Image: " << _Label_Image;
	std::cout << std::endl << "_Nuclei_Table: " << _Nuclei_Table;
	std::cout << std::endl << "_Dist_Map_Image: " << _Dist_Map_Image;
	std::cout << std::endl << "_Soma_Centroids: " << _Soma_Centroids;
	std::cout << std::endl << "_Soma_Montage: " << _Soma_Montage;
	std::cout << std::endl << "_isSmall: " << _isSmall;
	std::cout << std::endl << "_astroTraceParams: " << _astroTraceParams;
	std::cout << std::endl << "_roots_model_AL: " << _roots_model_AL;
	std::cout << std::endl << "_final_classification_model: " << _final_classification_model;
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
	if( !_Label_Image.empty() )
	{
		rawImageType_uint::Pointer ImageMontage_Label = readImage< rawImageType_uint >(_Label_ImageNRRD.c_str());
		_ImageMontage_LabelSize = ImageMontage_Label->GetLargestPossibleRegion().GetSize();
		splitStore< rawImageType_uint >( ImageMontage_Label, _Label_Image );
	}

	if( !_Dist_Map_Image.empty() )
	{
		rawImageType_flo::Pointer ImageMontage_Dist_Map = readImage< rawImageType_flo >(_Dist_Map_ImageNRRD.c_str());
		_ImageMontage_Dist_MapSize = ImageMontage_Dist_Map->GetLargestPossibleRegion().GetSize();
		splitStore< rawImageType_flo >( ImageMontage_Dist_Map, _Dist_Map_Image );
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
    
    int num_threads = 1;
#ifdef _OPENMP
	omp_set_nested(1);
#if _OPENMP >= 200805L
	omp_set_max_active_levels(2);
#endif
	omp_set_num_threads(num_threads);
#endif
	std::cout << std::endl << "fir MAX NUMBER OF CORES: " << itk::MultiThreader::GetGlobalMaximumNumberOfThreads() << ", fir DEFAULT NUMBER OF CORES: " << itk::MultiThreader::GetGlobalDefaultNumberOfThreads();
	
	vtkSmartPointer< vtkTable > AllNucleiTable = ftk::LoadTable( _Nuclei_Table );

	int contadorSegment = 0;
	int contadorImageDoneSegment = 0;
	int contadorStich = 0;
	unsigned int maxValue = 0;
	unsigned int maxValueOld = 0;
	int flagFirstStich = 1;
	
#if _OPENMP < 200805L
	#pragma omp parallel
#else
	#pragma omp parallel for collapse(3) num_threads(_num_threads) schedule(dynamic, 1)
#endif
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
				rawImageType_uint::Pointer imageLocalLabel;
				rawImageType_flo::Pointer imageLocalTRI;
				rawImageType_flo::Pointer imageLocalDist_Map;
				vtkSmartPointer< vtkTable > nucleiTable;

				std::vector< rawImageType_flo::Pointer > Images_Tiles;
				Images_Tiles.resize(1);

				std::vector< rawImageType_flo::Pointer > Dist_Map_Tiles;
				Dist_Map_Tiles.resize(1);

				std::vector< rawImageType_uint::Pointer > Label_Tiles;
				Label_Tiles.resize(1);
				
				std::vector< vtkSmartPointer< vtkTable > > Table_Tiles;
				Table_Tiles.resize(1);

				std::vector< vtkSmartPointer< vtkTable > > feature_Vector_Tables, root_Vector_Tables;
				feature_Vector_Tables.resize(1);
				root_Vector_Tables.resize(1);
				
				std::vector< rawImageType_16bit::Pointer > ID_Images, root_Images;
				ID_Images.resize(1);
				root_Images.resize(1);

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
					imageLocalDist_Map = readImage< rawImageType_flo >(tempDist_Map.c_str());
				}
				if( !_Label_Image.empty( ) )
				{
					int foundLabel = _Label_Image.find_last_of("/\\");
					std::string _Label_ImageNoPath = _Label_Image.substr(foundLabel+1);
					// !!! this has to be nrrd is just to test
					std::string tempLabel = _outPathTemp+"/"+_Label_ImageNoPath+"_"+xStr+"_"+yStr+"_"+zStr+".nrrd";
					imageLocalLabel = readImage< rawImageType_uint >(tempLabel.c_str());
				}	
				nucleiTable = vtkSmartPointer< vtkTable >::New(); 
				nucleiTable->Initialize();
				for(int col=0; col<(int)AllNucleiTable->GetNumberOfColumns(); col++)
				{
					vtkSmartPointer< vtkDoubleArray > column = vtkSmartPointer< vtkDoubleArray >::New();
					column->SetName(AllNucleiTable->GetColumnName(col));
					nucleiTable->AddColumn(column);
				}
				for(int row=0; row<(int)AllNucleiTable->GetNumberOfRows(); row++)
				{
					itk::Index<3> global_centroid;
					global_centroid[0] = AllNucleiTable->GetValue(row,1).ToUnsignedInt();
					global_centroid[1] = AllNucleiTable->GetValue(row,2).ToUnsignedInt();
					global_centroid[2] = AllNucleiTable->GetValue(row,3).ToUnsignedInt();
					if(regionMontage_inside.IsInside(global_centroid))
					{
						AllNucleiTable->SetValue(row, 1, global_centroid[0] - regionMontage_all.GetIndex()[0]);
						AllNucleiTable->SetValue(row, 2, global_centroid[1] - regionMontage_all.GetIndex()[1]);
						AllNucleiTable->SetValue(row, 3, global_centroid[2] - regionMontage_all.GetIndex()[2]);
						nucleiTable->InsertNextRow(AllNucleiTable->GetRow(row));
						AllNucleiTable->RemoveRow(row);
						row--;
					}
				}
				std::string nucleiTableFileName = _outPathTemp+"/InsideNucleiTable_"+xStr+"_"+yStr+"_"+zStr+".txt";
				ftk::SaveTable(nucleiTableFileName, nucleiTable);
				}


				//Images_Tiles[0] = imageLocalCy5;
				//Images_Tiles[1] = imageLocalTRI;
				//Images_Tiles[2] = imageLocalGFP;
				//Images_Tiles[3] = imageLocalDAP;
				Images_Tiles[0] = imageLocalTRI;
				Label_Tiles[0] = imageLocalLabel;
				Table_Tiles[0] = nucleiTable;
				Dist_Map_Tiles[0] = imageLocalDist_Map;
								
// 				RunSegmentation(regionLocal_inside, Images_Tiles, Label_Tiles, Table_Tiles, Centroids_Tiles);
				
				//Label_Tiles[0] = RunNuclearSegmentation( Images_Tiles[3] );
				//Table_Tiles[0] = ComputeFeaturesAndAssociations( Images_Tiles, Label_Tiles );
				//Centroids_Tiles[0] = GetLabelToCentroidMap(Table_Tiles[0]);

				AstroTracer * AT = new AstroTracer();
				AT->LoadCurvImage(Images_Tiles[0], 0);
				//AT->LoadSomaImage(std::string(argv[2]));
				std::string coverageFileName = _outPathTemp+"/coverage_"+xStr+"_"+yStr+"_"+zStr+".txt";
				AT->OptimizeCoverage(coverageFileName, true);
				AT->LoadParameters(_astroTraceParams.c_str());	
				AT->SetScaleRange(4, 4); //(2, 5); //(2, 2)
				AT->CallFeatureMainExternal();
				AT->Set_DistanceMapImage(imageLocalDist_Map);
				std::string featureVectorFileName = _outPathTemp+"/featureVector_"+xStr+"_"+yStr+"_"+zStr+".txt";
				std::string IDImageFileName = _outPathTemp+"/IDImage_"+xStr+"_"+yStr+"_"+zStr+".nrrd";
				AT->ComputeAstroFeaturesPipeline(featureVectorFileName, IDImageFileName, 0, regionLocal_inside, feature_Vector_Tables, ID_Images, true);
				std::string rootVectorFileName = _outPathTemp+"/rootVector_"+xStr+"_"+yStr+"_"+zStr+".txt";
				std::string rootsImageFileName = _outPathTemp+"/rootsImage_"+xStr+"_"+yStr+"_"+zStr+".nrrd";
				AT->Classification_Roots(root_Vector_Tables, root_Images, _roots_model_AL, rootVectorFileName, rootsImageFileName, true);
				//Centroids_Tiles[0] = GetLabelToCentroidMap(root_Vector_Tables[0]);
				delete AT;
				//RemoveLabelNearBorder(regionLocal_inside, root_Images, root_Vector_Tables, Centroids_Tiles );

				//typedef itk::LabelGeometryImageFilter< rawImageType_uint, rawImageType_flo > LabelGeometryType;
				//LabelGeometryType::Pointer labelGeometryFilter = LabelGeometryType::New();		//Dirk's Filter;
				//labelGeometryFilter->SetInput( imageLocalLabel );
				//labelGeometryFilter->SetIntensityInput( imageLocalTRI );
				////SET ADVANCED (OPTIONAL) ITEMS FOR THIS FILTER:
				//labelGeometryFilter->CalculatePixelIndicesOff();
				//labelGeometryFilter->CalculateOrientedBoundingBoxOff();
				//labelGeometryFilter->CalculateOrientedLabelRegionsOff();
				//labelGeometryFilter->CalculateOrientedIntensityRegionsOff();
				////UPDATE THE FILTER	
				//try
				//{
				//	labelGeometryFilter->Update();
				//}
				//catch (itk::ExceptionObject & e) 
				//{
				//	std::cerr << "Exception in Dirk's Label Geometry Filter: " << e << std::endl;
				//	return false;
				//}				
				//std::vector< LabelGeometryType::LabelPixelType > labels = labelGeometryFilter->GetLabels();
				
				
				// TEST SAVE SEGMENTATION WITH LABELS REMOVED
				if(!root_Vector_Tables[0])
					continue;
				#pragma omp critical
				{
					std::string tempTABLERE = _outPathTemp+"/rootVector_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REMO.txt";
					ftk::SaveTable(tempTABLERE, root_Vector_Tables[0]);
					std::string tempLABELRE = _outPathTemp+"/rootsImage_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REMO.nrrd";
					writeImage<rawImageType_16bit>(root_Images[0],tempLABELRE.c_str());

	// 				std::string tempLABELRE = _outPathDebugLevel2+"/segLabel_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REMO.tif";
					ftkMainDarpa objftkMainDarpa;
					objftkMainDarpa.projectImage<rawImageType_16bit, rawImageType_16bit>( root_Images[0], tempLABELRE, _outPathDebugLevel2, "ORG_BIN", "TIFF" );
				}
				#pragma omp critical
				{
					contadorImageDoneSegment++;
					std::cout<<std::endl<< "\t\t--->>> ImageDoneSegment " << contadorImageDoneSegment << " of " << _kx*_ky*_kz;
				}
			}
		}
	}
}



void ftkMainDarpaAstroTrace::runStitchRoots(  )
{
	std::cout << std::endl << "HERE HERE_2";
	std::cout << std::endl << "HERE HERE_2";
	
	if( _isSmall == 1 )
	{
		itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1); // This one can not be changed
		itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1); // This one can chenga
	
		int num_threads = 1;
	#ifdef _OPENMP
		omp_set_nested(1);
	#if _OPENMP >= 200805L
		omp_set_max_active_levels(2);
	#endif
		omp_set_num_threads(num_threads);
    	#endif
	}
	
	rawImageType_8bit::RegionType ImageMontageRegion;
	itk::Size<3> ImageMontageSize;
	// Assume GFP or dapi always exist
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
	
	std::cout << std::endl << "HERE HERE";
	std::cout << std::endl << "HERE HERE";
	
	rawImageType_uint::Pointer imageRootsMontage = rawImageType_uint::New();
	rawImageType_uint::PointType originz;
	originz[0] = 0; 
	originz[1] = 0;
	originz[2] = 0;
	imageRootsMontage->SetOrigin( originz );
	rawImageType_uint::IndexType indexStich;
	indexStich.Fill(0);
	rawImageType_uint::RegionType regionz;
	regionz.SetSize ( ImageMontageSize  );
	regionz.SetIndex( indexStich );
	imageRootsMontage->SetRegions( regionz );
	
	try
	{
		imageRootsMontage->Allocate();
		imageRootsMontage->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}
	imageRootsMontage->FillBuffer(0);
	
	// GLOBAL TABLE
	vtkSmartPointer< vtkTable > tableRootsMontage;
	
	
	int contadorSegment = 0;
	int contadorStich = 0;
	unsigned int maxValue = 0;
	unsigned int maxValueOld = 0;
	int flagFirstStich = 1;
	
// #pragma omp parallel for collapse(3) num_threads(_num_threads) schedule(dynamic, 1)
	for(int xco = 0; xco < _kx; xco++)
	{
		for(int yco = 0; yco < _ky; yco++)
		{
			for(int zco = 0; zco < _kz; zco++)
			{
// 				#pragma omp critical
// 				{
				++contadorSegment;
				std::cout<<std::endl<< "\t\t--->>> ImageStich " << contadorSegment << " of " << _kx*_ky*_kz;
// 				}
				
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
				
// 				rawImageType_8bit::Pointer imageLocalCy5;
// 				rawImageType_8bit::Pointer imageLocalTRI;
// 				rawImageType_8bit::Pointer imageLocalGFP;
// 				rawImageType_8bit::Pointer imageLocalDAP;

				//std::vector< rawImageType_8bit::Pointer > Images_Tiles;
				//Images_Tiles.resize(4);

				
				std::vector< rawImageType_16bit::Pointer > root_Images;
				root_Images.resize(1);
				
				std::vector< vtkSmartPointer< vtkTable > > root_Vector_Tables;
				root_Vector_Tables.resize(1);
				
				//std::vector< std::map< unsigned int, itk::Index<3> > > Centroids_Tiles;
				//Centroids_Tiles.resize(1);
				
// 				// Reading part
// 				if( !_GFP_Image.empty( ) && !_DAP_Image.empty() )
// 				{
// 					int foundGFP = _GFP_Image.find_last_of("/\\");
// 					std::string _GFP_ImageNoPath = _GFP_Image.substr(foundGFP+1);
// 					// !!! this has to be nrrd is just to test
// 					std::string tempGFP = _outPathTemp+"/"+_GFP_ImageNoPath+"_"+xStr+"_"+yStr+"_"+zStr+".nrrd";
// 					imageLocalGFP = readImage< rawImageType_8bit >(tempGFP.c_str());
// 					
// 					int foundDAP = _DAP_Image.find_last_of("/\\");
// 					std::string _DAP_ImageNoPath = _DAP_Image.substr(foundDAP+1);
// 					// !!! this has to be nrrd is just to test
// 					std::string tempDAP = _outPathTemp+"/"+_DAP_ImageNoPath+"_"+xStr+"_"+yStr+"_"+zStr+".nrrd";
// 					imageLocalDAP = readImage< rawImageType_8bit >(tempDAP.c_str());
// 				}
// 				Images_Tiles[0] = imageLocalCy5;
// 				Images_Tiles[1] = imageLocalTRI;
// 				Images_Tiles[2] = imageLocalGFP;
// 				Images_Tiles[3] = imageLocalDAP;
// 				
// // 				RunSegmentation(regionLocal_inside, Images_Tiles, Label_Tiles, Table_Tiles, Centroids_Tiles);
// 				
// 				Label_Tiles[0] = RunNuclearSegmentation( Images_Tiles[3] );
// 				Table_Tiles[0] = ComputeFeaturesAndAssociations( Images_Tiles, Label_Tiles );
// 				Centroids_Tiles[0] = GetLabelToCentroidMap(Table_Tiles[0]);
// 				
// // 					// TEST PUT A NUMBER IN THE REGION OF INTERES
// // 					IteratorType_16bit iterLocal16_2 = IteratorType_16bit(Label_Tiles[0],regionLocal_inside); iterLocal16_2.GoToBegin();
// // 					for(;!iterLocal16_2.IsAtEnd();++iterLocal16_2)
// // 						iterLocal16_2.Set(50000);
// // 					// TEST SAVE SEGMENTATION
// // 					std::string tempLABEL = _outPathDebugLevel2+"/segLabel_"+"_"+xStr+"_"+yStr+"_"+zStr+".tif";
// // 					std::string tempTABLE = _outPathDebugLevel2+"/segTable_"+"_"+xStr+"_"+yStr+"_"+zStr+".txt";
// // 					writeImage<rawImageType_16bit>(Label_Tiles[0],tempLABEL.c_str());
// // 					ftk::SaveTable(tempTABLE, Table_Tiles[0]);
// 				
// 				RemoveLabelNearBorder(regionLocal_inside, Label_Tiles, Table_Tiles, Centroids_Tiles );
				
					
					
// 					// TEST SAVE SEGMENTATION WITH LABELS REMOVED
				std::string tempTABLERE = _outPathTemp+"/rootVector_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REMO.txt";
				root_Vector_Tables[0] = ftk::LoadTable(tempTABLERE);
				if(!root_Vector_Tables[0])
					continue;
				
				std::string tempLABELRE = _outPathTemp+"/rootsImage_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REMO.nrrd";
				root_Images[0] = readImage<rawImageType_16bit>(tempLABELRE.c_str());

// 				std::string tempLABELRE = _outPathDebugLevel2+"/segLabel_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REMO.tif";
// 				ftkMainDarpa objftkMainDarpa;
// 				objftkMainDarpa.projectImage<rawImageType_16bit, rawImageType_16bit>( root_Images[0], tempLABELRE, _outPathDebugLevel2, "ORG_RES_BIN" );

// /*					// TEST PUT A NUMBER IN THE REGION OF INTERES
// 					IteratorType_16bit iterLocal16_3 = IteratorType_16bit(root_Images[0],regionLocal_inside); iterLocal16_3.GoToBegin();
// 					for(;!iterLocal16_3.IsAtEnd();++iterLocal16_3)
// 						iterLocal16_3.Set(50000);
// 					std::string tempLABELREGIO = _outPathDebugLevel2+"/segLabel_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REGION.nrrd";
// 					writeImage<rawImageType_16bit>(root_Images[0],tempLABELREGIO.c_str());*/
				
// 				StichResults( imageRootsMontage, root_Images, root_Vector_Tables, Centroids_Tiles );

// 				#pragma omp critical
// 				{
				++contadorStich;
				std::cout<<std::endl<< "\t\t--->>> ImageStich " << contadorStich << " of " << _kx*_ky*_kz;
				
				if( flagFirstStich == 1 )
				{
					flagFirstStich = 0;
					tableRootsMontage = vtkSmartPointer<vtkTable>::New();
					tableRootsMontage->Initialize();
					for(int c=0; c<(int)root_Vector_Tables[0]->GetNumberOfColumns(); ++c)
					{
						vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
						column->SetName( root_Vector_Tables[0]->GetColumnName(c) );
						tableRootsMontage->AddColumn(column);
					}
				}
				
				maxValueOld = maxValue;
				if((unsigned long long)root_Vector_Tables[0]->GetNumberOfRows() != 0)
				{
					std::cout << std::endl << "\t\tpiu11 The number of row: " << (int)root_Vector_Tables[0]->GetNumberOfRows() << " in " << contadorStich;
					for(int r=0; r<(int)root_Vector_Tables[0]->GetNumberOfRows(); ++r)
					{
						// Can be done here
						vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
						for(int c=0; c<(int)root_Vector_Tables[0]->GetNumberOfColumns(); ++c)
						{
							if(c == 0)
								model_data1->InsertNextValue(vtkVariant(root_Vector_Tables[0]->GetValue(r,c).ToUnsignedInt() + maxValue));
							else if(c == 1)
								model_data1->InsertNextValue(vtkVariant(root_Vector_Tables[0]->GetValue(r,c).ToInt() + regionMontage_all.GetIndex()[0]));
							else if(c == 2)
								model_data1->InsertNextValue(vtkVariant(root_Vector_Tables[0]->GetValue(r,c).ToInt() + regionMontage_all.GetIndex()[1]));
							else if(c == 3)
								model_data1->InsertNextValue(vtkVariant(root_Vector_Tables[0]->GetValue(r,c).ToInt() + regionMontage_all.GetIndex()[2]));
							else
								model_data1->InsertNextValue(root_Vector_Tables[0]->GetValue(r,c));
						}
						tableRootsMontage->InsertNextRow(model_data1);
					}
					std::cout << std::endl << "\t\tpiu11 The maxvaule before: " << maxValue << " ";
					maxValue = tableRootsMontage->GetValue((int)tableRootsMontage->GetNumberOfRows()-1, 0).ToUnsignedInt();
					std::cout << "after: " << maxValue << " in " << contadorStich;
				}
			
				IteratorType_uint iterMontage_1 = IteratorType_uint(imageRootsMontage,regionMontage_all);
				iterMontage_1.GoToBegin();
				
// 				IteratorType_8bit iterLocal_1 = IteratorType_8bit(imageLocalDAP,regionLocal_inside);
// 				iterLocal_1.GoToBegin();
				
				IteratorType_16bit iterLocal16_1 = IteratorType_16bit(root_Images[0],regionLocal_all);
				iterLocal16_1.GoToBegin();
				
				for(;!iterMontage_1.IsAtEnd();++iterMontage_1)
				{
// 					iterMontage_1.Set(iterLocal_1.Get());
// 					++iterLocal8_1;
					if( iterLocal16_1.Get() != 0 )
						iterMontage_1.Set(iterLocal16_1.Get() + maxValueOld);
					++iterLocal16_1;
				}
				
				std::cout<<std::endl<< "\t\t--->>> ImageStichDone " << contadorStich << " of " << _kx*_ky*_kz;
// 				}
			}
		}
	}

	std::string tempAllRoots_Montage = _outPathTemp+"/Roots_Image.nrrd";
// 	std::cout << std::endl << "SOMA STORED " << temp9;
	writeImage<rawImageType_uint>(imageRootsMontage,tempAllRoots_Montage.c_str());
	
	std::string tempAllRoots_Table = _outPathTemp+"/Roots_Table.txt";
// 	std::cout << std::endl << "SOMA STORED " << temp9;
	ftk::SaveTable(tempAllRoots_Table, tableRootsMontage);
}


void ftkMainDarpaAstroTrace::computeRootFeaturesForNuclei(  )
{
	if( _isSmall == 1 )
	{
		itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1); // This one can not be changed
		itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1); // This one can chenga
		int num_threads = 1;
	#ifdef _OPENMP
		omp_set_nested(1);
	#if _OPENMP >= 200805L
		omp_set_max_active_levels(2);
	#endif
		omp_set_num_threads(num_threads);
    	#endif
	}
	
	rawImageType_8bit::RegionType ImageMontageRegion;
	itk::Size<3> ImageMontageSize;
	// Assume GFP or dapi always exist
	if( !_TRI_Image.empty() )
	{
		rawImageType_8bit::Pointer ImageMontage_TRI = readImage< rawImageType_8bit >(_TRI_ImageNRRD.c_str());
		ImageMontageRegion = ImageMontage_TRI->GetLargestPossibleRegion();
		ImageMontageSize = ImageMontageRegion.GetSize();
		computeSplitConst( ImageMontage_TRI );
	}

	std::string tempAllRoots_Table = _outPathTemp+"/Roots_Table.txt";
	vtkSmartPointer< vtkTable > AllRootsTable = ftk::LoadTable( tempAllRoots_Table );

#if _OPENMP < 200805L
	#pragma omp parallel
#else
	#pragma omp parallel for collapse(3) num_threads(_num_threads) schedule(dynamic, 1)
#endif	
	for(int xco = 0; xco < _kx; xco++)
	{
		for(int yco = 0; yco < _ky; yco++)
		{
			for(int zco = 0; zco < _kz; zco++)
			{
				//#pragma omp critical
				//{
				//	++contadorSegment;
				//	std::cout<<std::endl<< "\t\t--->>> ImageSegment " << contadorSegment << " of " << _kx*_ky*_kz;
				//}
				
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

				rawImageType_uint::Pointer imageLocalLabel;
				rawImageType_flo::Pointer imageLocalTRI;
				rawImageType_flo::Pointer imageLocalDist_Map;
				vtkSmartPointer< vtkTable > nucleiTable, rootsTable;

				std::vector< rawImageType_flo::Pointer > Images_Tiles;
				Images_Tiles.resize(1);

				std::vector< rawImageType_flo::Pointer > Dist_Map_Tiles;
				Dist_Map_Tiles.resize(1);

				std::vector< rawImageType_uint::Pointer > Label_Tiles;
				Label_Tiles.resize(1);
				
				std::vector< vtkSmartPointer< vtkTable > > Table_Tiles;
				Table_Tiles.resize(1);

				std::vector< vtkSmartPointer< vtkTable > > feature_Vector_Tables, root_Vector_Tables;
				feature_Vector_Tables.resize(1);
				root_Vector_Tables.resize(1);
				
				std::vector< rawImageType_16bit::Pointer > ID_Images, root_Images;
				ID_Images.resize(1);
				root_Images.resize(1);

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
					imageLocalDist_Map = readImage< rawImageType_flo >(tempDist_Map.c_str());
				}
				if( !_Label_Image.empty( ) )
				{
					int foundLabel = _Label_Image.find_last_of("/\\");
					std::string _Label_ImageNoPath = _Label_Image.substr(foundLabel+1);
					// !!! this has to be nrrd is just to test
					std::string tempLabel = _outPathTemp+"/"+_Label_ImageNoPath+"_"+xStr+"_"+yStr+"_"+zStr+".nrrd";
					imageLocalLabel = readImage< rawImageType_uint >(tempLabel.c_str());
				}
				rootsTable = vtkSmartPointer< vtkTable >::New(); 
				rootsTable->Initialize();
				for(int col=0; col<(int)AllRootsTable->GetNumberOfColumns(); col++)
				{
					vtkSmartPointer< vtkDoubleArray > column = vtkSmartPointer< vtkDoubleArray >::New();
					column->SetName(AllRootsTable->GetColumnName(col));
					rootsTable->AddColumn(column);
				}
				for(int row=0; row<(int)AllRootsTable->GetNumberOfRows(); row++)
				{
					itk::Index<3> global_root;
					global_root[0] = AllRootsTable->GetValue(row,1).ToUnsignedInt();
					global_root[1] = AllRootsTable->GetValue(row,2).ToUnsignedInt();
					global_root[2] = AllRootsTable->GetValue(row,3).ToUnsignedInt();
					if(regionMontage_all.IsInside(global_root))
					{
						AllRootsTable->SetValue(row, 1, global_root[0] - regionMontage_all.GetIndex()[0]);
						AllRootsTable->SetValue(row, 2, global_root[1] - regionMontage_all.GetIndex()[1]);
						AllRootsTable->SetValue(row, 3, global_root[2] - regionMontage_all.GetIndex()[2]);
						rootsTable->InsertNextRow(AllRootsTable->GetRow(row));
						AllRootsTable->RemoveRow(row);
						row--;
					}
				}
				std::string testRootsTableFileName = _outPathTemp+"/test_roots_table_"+xStr+"_"+yStr+"_"+zStr+".txt";
				ftk::SaveTable(testRootsTableFileName, rootsTable);
				}
				std::string nucleiTableFileName = _outPathTemp+"/InsideNucleiTable_"+xStr+"_"+yStr+"_"+zStr+".txt";
				nucleiTable = ftk::LoadTable(nucleiTableFileName);
				//Images_Tiles[0] = imageLocalCy5;
				//Images_Tiles[1] = imageLocalTRI;
				//Images_Tiles[2] = imageLocalGFP;
				//Images_Tiles[3] = imageLocalDAP;
				Images_Tiles[0] = imageLocalTRI;
				Label_Tiles[0] = imageLocalLabel;
				Table_Tiles[0] = nucleiTable;
				Dist_Map_Tiles[0] = imageLocalDist_Map;
				root_Vector_Tables[0] = rootsTable;
								
// 				RunSegmentation(regionLocal_inside, Images_Tiles, Label_Tiles, Table_Tiles, Centroids_Tiles);
				
				//Label_Tiles[0] = RunNuclearSegmentation( Images_Tiles[3] );
				//Table_Tiles[0] = ComputeFeaturesAndAssociations( Images_Tiles, Label_Tiles );
				//Centroids_Tiles[0] = GetLabelToCentroidMap(Table_Tiles[0]);

				AstroTracer * AT = new AstroTracer();
				AT->LoadCurvImage(Images_Tiles[0], 0);
				AT->LoadParameters(_astroTraceParams.c_str());	
				AT->SetScaleRange(4, 4); //(2, 5); //(2, 2)
				AT->Set_DistanceMapImage(imageLocalDist_Map);
				AT->ReadRootPointsPipeline(root_Vector_Tables);
				AT->ReadNucleiFeaturesPipeline(Table_Tiles);
				std::string AppendedNucleiTableFileName = _outPathTemp+"/AppendedNucleiTable_"+xStr+"_"+yStr+"_"+zStr+".txt";
				AT->ComputeFeaturesFromCandidateRootsPipeline(regionLocal_inside, Table_Tiles, _final_classification_model, AppendedNucleiTableFileName);
				//AT->WriteNucleiFeatures(AppendedNucleiTableFileName);
				delete AT;
				//RemoveLabelNearBorder(regionLocal_inside, root_Images, root_Vector_Tables, Centroids_Tiles );

			}
		}
	}

	vtkSmartPointer< vtkTable > AllNucleiTable = vtkSmartPointer< vtkTable >::New();
	AllNucleiTable->Initialize();

	for(int xco = 0; xco < _kx; xco++)
	{
		for(int yco = 0; yco < _ky; yco++)
		{
			for(int zco = 0; zco < _kz; zco++)
			{
				//#pragma omp critical
				//{
				//	++contadorSegment;
				//	std::cout<<std::endl<< "\t\t--->>> ImageSegment " << contadorSegment << " of " << _kx*_ky*_kz;
				//}
				
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

				std::string AppendedNucleiTableFileName = _outPathTemp+"/AppendedNucleiTable_"+xStr+"_"+yStr+"_"+zStr+".txt";
				vtkSmartPointer< vtkTable > nucleiTable = ftk::LoadTable(AppendedNucleiTableFileName);

				if( (xco==0) && (yco==0) && (zco==0) )
				{
					for(int col=0; col<(int)nucleiTable->GetNumberOfColumns(); col++)
					{
						vtkSmartPointer< vtkDoubleArray > column = vtkSmartPointer< vtkDoubleArray >::New();
						column->SetName(nucleiTable->GetColumnName(col));
						AllNucleiTable->AddColumn(column);
					}
				}
				for(int row=0; row<(int)nucleiTable->GetNumberOfRows(); row++)
				{
					itk::Index<3> local_centroid;
					local_centroid[0] = nucleiTable->GetValue(row,1).ToUnsignedInt();
					local_centroid[1] = nucleiTable->GetValue(row,2).ToUnsignedInt();
					local_centroid[2] = nucleiTable->GetValue(row,3).ToUnsignedInt();
					nucleiTable->SetValue(row, 1, local_centroid[0] + regionMontage_all.GetIndex()[0]);
					nucleiTable->SetValue(row, 2, local_centroid[1] + regionMontage_all.GetIndex()[1]);
					nucleiTable->SetValue(row, 3, local_centroid[2] + regionMontage_all.GetIndex()[2]);
					AllNucleiTable->InsertNextRow(nucleiTable->GetRow(row));
				}				
			}
		}
	}

	//vtkSmartPointer<vtkTable> active_model_table = ftk::LoadTable(_final_classification_model);

	//// to generate the Active Learning Matrix
	//vnl_matrix<double> act_learn_matrix;
	//act_learn_matrix.set_size((int)active_model_table->GetNumberOfColumns() , (int)active_model_table->GetNumberOfRows() - 2);
	//for(int row = 2; row<(int)active_model_table->GetNumberOfRows(); ++row)
	//{
	//	for(int col=0; col<(int)active_model_table->GetNumberOfColumns(); ++col)
	//	{
	//		act_learn_matrix.put(col, row-2, active_model_table->GetValue(row,col).ToDouble());
	//	}
	//}

	////to generate the std_deviation and the mean vectors
	//vnl_vector<double> std_dev_vec, mean_vec; 
	//std_dev_vec.set_size((int)active_model_table->GetNumberOfColumns() - 1);
	//mean_vec.set_size((int)active_model_table->GetNumberOfColumns() - 1);
	//for(int col=1; col<(int)active_model_table->GetNumberOfColumns(); ++col)
	//{
	//	std_dev_vec.put(col-1, active_model_table->GetValue(0,col).ToDouble());
	//	mean_vec.put(col-1, active_model_table->GetValue(1,col).ToDouble());
	//}
	//active_model_table->RemoveRow(0);
	//active_model_table->RemoveRow(0);
	//active_model_table->RemoveColumn(0);

	//std::string classification_name = "multi_class";
	//double confidence_thresh = 1/(int)active_model_table->GetNumberOfRows();

	//MCLR* mclr = new MCLR();
	////Number of features and classes needed in "add_bias" fuction of MCLR
	//mclr->Set_Number_Of_Classes((int)active_model_table->GetNumberOfRows());
	//mclr->Set_Number_Of_Features((int)active_model_table->GetNumberOfColumns());

	////setup the test data
	//vtkSmartPointer<vtkTable> test_table  = vtkSmartPointer<vtkTable>::New();
	//test_table->Initialize();
	//for(int col=0; col<(int)active_model_table->GetNumberOfColumns(); ++col)
	//{
	//	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	//	column->SetName(active_model_table->GetColumnName(col));
	//	column->SetNumberOfValues(AllNucleiTable->GetNumberOfRows());
	//	test_table->AddColumn(column);	
	//}
	//for(int row = 0; row < (int)AllNucleiTable->GetNumberOfRows(); ++row)
	//{	
	//	vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
	//	for(int c=0; c<(int)test_table->GetNumberOfColumns();++c)
	//		model_data1->InsertNextValue(AllNucleiTable->GetValueByName(row,test_table->GetColumnName(c)));
	//	test_table->InsertNextRow(model_data1);
	//}

	//////// Final Data  to classify from the model
	//vnl_matrix<double> data_classify;
	//data_classify =  mclr->Normalize_Feature_Matrix_w(mclr->tableToMatrix_w(test_table), std_dev_vec, mean_vec);
	//data_classify = data_classify.transpose();

	//vnl_matrix<double> currprob;
	//currprob = mclr->Test_Current_Model_w(data_classify, act_learn_matrix);

	//std::string prediction_col_name = "prediction_active_" + classification_name;
	//std::string confidence_col_name = "confidence_" + classification_name;

	//// Add the Prediction Column 
	//std::vector< std::string > prediction_column_names = ftk::GetColumsWithString(prediction_col_name, AllNucleiTable );
	//if(prediction_column_names.size() == 0)
	//{
	//	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	//	column->SetName(prediction_col_name.c_str());
	//	column->SetNumberOfValues( AllNucleiTable->GetNumberOfRows() );
	//	AllNucleiTable->AddColumn(column);
	//}

	//// Add the confidence column
	//std::vector< std::string > confidence_column_names = ftk::GetColumsWithString(confidence_col_name, AllNucleiTable );
	//if(confidence_column_names.size() == 0)
	//{
	//	vtkSmartPointer<vtkDoubleArray> column_confidence = vtkSmartPointer<vtkDoubleArray>::New();
	//	column_confidence->SetName(confidence_col_name.c_str());
	//	column_confidence->SetNumberOfValues( AllNucleiTable->GetNumberOfRows() );
	//	AllNucleiTable->AddColumn(column_confidence);
	//}

	//for(int row = 0; row<(int)AllNucleiTable->GetNumberOfRows(); ++row)  
	//{
	//	vnl_vector<double> curr_col = currprob.get_column(row);
	//	AllNucleiTable->SetValueByName(row, confidence_col_name.c_str(), vtkVariant(curr_col(curr_col.arg_max())));
	//	if(curr_col(curr_col.arg_max()) > confidence_thresh) 
	//	{
	//		AllNucleiTable->SetValueByName(row, prediction_col_name.c_str(), vtkVariant(curr_col.arg_max()+1));						
	//	}
	//	else
	//	{
	//		AllNucleiTable->SetValueByName(row, prediction_col_name.c_str(), vtkVariant(0));
	//	}
	//}
	//
	//delete mclr;
	//	

	std::string tempMulti_Class_Nuclei_Table_Table = _outPathTemp+"/Multi_Class_Nuclei_Table.txt";
	ftk::SaveTable(tempMulti_Class_Nuclei_Table_Table, AllNucleiTable);

	rawImageType_uint::Pointer imageLabelMontage = readImage< rawImageType_uint >(_Label_ImageNRRD.c_str());

	// Soma Montage
	rawImageType_uint::Pointer astrocyteSomaMontage = rawImageType_uint::New();
	rawImageType_uint::RegionType uint_ImageMontageRegion = imageLabelMontage->GetBufferedRegion();
	itk::Size<3> im_size = uint_ImageMontageRegion.GetSize();
	rawImageType_uint::RegionType region;
	region.SetSize( uint_ImageMontageRegion.GetSize() );
	region.SetIndex( uint_ImageMontageRegion.GetIndex() );
	astrocyteSomaMontage->SetRegions( region );
	astrocyteSomaMontage->Allocate();
	astrocyteSomaMontage->FillBuffer(0);
	try
	{
		astrocyteSomaMontage->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}
	rawImageType_uint::PixelType * astrocyteSomaArray = astrocyteSomaMontage->GetBufferPointer();
	rawImageType_uint::PixelType * imageLabelArray = imageLabelMontage->GetBufferPointer();

	unsigned long long sizeXY = im_size[1] * im_size[0];
	unsigned long long sizeX = im_size[0];
	std::map<unsigned int, int> classMap;
	for(int row=0; row<(int)AllNucleiTable->GetNumberOfRows(); ++row)
	{
		classMap[AllNucleiTable->GetValue(row,0).ToUnsignedInt()] = AllNucleiTable->GetValueByName(row, "prediction_active_multi_class").ToInt();
	}

#if _OPENMP < 200805L
	#pragma omp parallel
#else
	#pragma omp parallel for collapse(3)
#endif	
	for(int i=0; i<im_size[2]; ++i)
	{
		for(int j=0; j<im_size[1]; ++j)
		{
			for(int k=0; k<im_size[0]; ++k)
			{
				unsigned long long offset = (i*sizeXY)+(j*sizeX)+k;
				if( imageLabelArray[offset] != 0 )
					if(classMap[imageLabelArray[offset]] == 1)
						astrocyteSomaArray[offset] = imageLabelArray[offset];
			}
		}
	}	
	
	std::string nameAstrocyteSomaMontage = _outPathTemp+"/astrocyte_soma.nrrd";
	writeImage< rawImageType_uint >(astrocyteSomaMontage,nameAstrocyteSomaMontage.c_str());

	rawImageType_uint::Pointer microgliaSomaMontage = rawImageType_uint::New();
	rawImageType_uint::RegionType region1;
	region1.SetSize( uint_ImageMontageRegion.GetSize() );
	region1.SetIndex( uint_ImageMontageRegion.GetIndex() );
	microgliaSomaMontage->SetRegions( region1 );
	microgliaSomaMontage->Allocate();
	microgliaSomaMontage->FillBuffer(0);
	try
	{
		microgliaSomaMontage->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}
	rawImageType_uint::PixelType * microgliaSomaArray = microgliaSomaMontage->GetBufferPointer();
	//rawImageType_uint::PixelType * imageLabelArray = imageLabelMontage->GetBufferPointer();

	//unsigned long long sizeXY = im_size[1] * im_size[0];
	//unsigned long long sizeX = im_size[0];
	//std::map<unsigned int, int> classMap;
	//for(int row=0; row<(int)AllNucleiTable->GetNumberOfRows(); ++row)
	//{
	//	classMap[AllNucleiTable->GetValue(row,0).ToUnsignedInt()] = AllNucleiTable->GetValueByName(row, "prediction_active_multi_class").ToInt();
	//}

#if _OPENMP < 200805L
	#pragma omp parallel
#else
	#pragma omp parallel for collapse(3)
#endif
	for(int i=0; i<im_size[2]; ++i)
	{
		for(int j=0; j<im_size[1]; ++j)
		{
			for(int k=0; k<im_size[0]; ++k)
			{
				unsigned long long offset = (i*sizeXY)+(j*sizeX)+k;
				if( imageLabelArray[offset] != 0 )
					if(classMap[imageLabelArray[offset]] == 2)
						microgliaSomaArray[offset] = imageLabelArray[offset];
			}
		}
	}	
	
	std::string nameMicrogliaSomaMontage = _outPathTemp+"/microglia_soma.nrrd";
	writeImage< rawImageType_uint >(microgliaSomaMontage,nameMicrogliaSomaMontage.c_str());

}

	
	


	
void ftkMainDarpaAstroTrace::RemoveLabelNearBorder(rawImageType_8bit::RegionType regionLocal_inside,std::vector< rawImageType_16bit::Pointer >& Label_Tiles, std::vector< vtkSmartPointer< vtkTable > >& Table_Tiles, std::vector< std::map< unsigned int, itk::Index<3> > >& Centroids_Tiles )
{
	rawImageType_16bit::PixelType * Label_TilesArray = Label_Tiles[0]->GetBufferPointer();
	itk::Size<3> tileSize = Label_Tiles[0]->GetLargestPossibleRegion().GetSize();
	unsigned long long tileSizeXY = tileSize[1] * tileSize[0];
	unsigned long long tileSizeX = tileSize[0];	
	for(unsigned long long z=0; z<tileSize[2]; ++z)
	{
		for(unsigned long long y=0; y<tileSize[1]; ++y)
		{
			for(unsigned long long x=0; x<tileSize[0]; ++x)
			{
				unsigned short value = Label_TilesArray[(tileSizeXY*z) + (tileSizeX*y) + (x)];
				if(value == 0) 
					continue;
				if( !regionLocal_inside.IsInside( Centroids_Tiles[0][value] ) )
					Label_TilesArray[(tileSizeXY*z) + (tileSizeX*y) + (x)] = 0;
			}
		}
	}
	
	vtkSmartPointer< vtkTable > Table_Tiles_out = vtkSmartPointer<vtkTable>::New();
	Table_Tiles_out->Initialize();
	for(int c=0; c<(int)Table_Tiles[0]->GetNumberOfColumns(); ++c)
	{
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( Table_Tiles[0]->GetColumnName(c) );
		Table_Tiles_out->AddColumn(column);
	}
	for(int r=0; r<(int)Table_Tiles[0]->GetNumberOfRows(); ++r)
	{
		rawImageType_8bit::IndexType indexCentroid = Centroids_Tiles[0][ Table_Tiles[0]->GetValue(r, 0).ToUnsignedInt() ];
		if( !regionLocal_inside.IsInside( indexCentroid ) )
			continue;
		
		vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
		for(int c=0; c<(int)Table_Tiles[0]->GetNumberOfColumns(); ++c)
		{
			model_data1->InsertNextValue(Table_Tiles[0]->GetValue(r,c));
		}
		Table_Tiles_out->InsertNextRow(model_data1);
	}
	// Any kind of memory leakage ???
	Table_Tiles[0] = Table_Tiles_out;
}


std::map< unsigned int, itk::Index<3> > ftkMainDarpaAstroTrace::GetLabelToCentroidMap( vtkSmartPointer< vtkTable > table)
{
	//std::cout << std::endl << (int)table->GetNumberOfRows();
	//std::cout << std::endl << (int)table->GetNumberOfRows();
	
	std::map< unsigned int, itk::Index<3> > centroidMap;
	for (int row=0; row<(int)table->GetNumberOfRows(); ++row)
	{
		unsigned int id = table->GetValue(row, 0).ToUnsignedInt();
		itk::Index<3> indx;
		indx[0] = table->GetValue(row, 1).ToUnsignedInt();
		indx[1] = table->GetValue(row, 2).ToUnsignedInt();
		indx[2] = table->GetValue(row, 3).ToUnsignedInt();
		centroidMap[id] = indx;
	}

	return centroidMap;
}
	
	



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
