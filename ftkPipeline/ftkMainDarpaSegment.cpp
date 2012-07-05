
#include "ftkMainDarpaSegment.h"

void ftkMainDarpaSegment::readParameters( std::string segmentParams )
{
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
	{ _xTileBor = 1024; printf("Choose xTileBor = 1024 as default\n");}
	
	iter = options.find("-yTileBor"); 
	if(iter!=options.end())
	{ std::istringstream ss((*iter).second); ss >> _yTileBor;}
	else
	{ _yTileBor = 1024; printf("Choose yTileBor = 1024 as default\n");}
	
	iter = options.find("-zTileBor"); 
	if(iter!=options.end())
	{ std::istringstream ss((*iter).second); ss >> _zTileBor;}
	else
	{ _zTileBor = 1024; printf("Choose _zTileBor = 1024 as default\n");}
	
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
	
	iter = options.find("-isSmall"); 
	if(iter!=options.end())
	{ std::istringstream ss((*iter).second); ss >> _isSmall;}
	else
	{ _isSmall = 1; printf("Choose _isSmall = 1 as default\n");}
	
	iter = options.find("-segParams"); 
	if(iter!=options.end())
	{ std::istringstream ss((*iter).second); ss >> _segParams;}
	else
	{ _segParams.clear(); printf("Choose _segParams = NULL as default\n");}
	
	iter = options.find("-projectDefinition"); 
	if(iter!=options.end())
	{ std::istringstream ss((*iter).second); ss >> _projectDefinition;}
	else
	{ _projectDefinition.clear(); printf("Choose _projectDefinition = NULL as default\n");}
	
	iter = options.find("-optionsMNT"); 
	if(iter!=options.end())
	{ std::istringstream ss((*iter).second); ss >> _optionsMNT;}
	else
	{ _optionsMNT.clear(); printf("Choose optionsMNT = NULL as default\n");}
	
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
	
	iter = options.find("-outPathData"); 
	if(iter!=options.end())
	{ std::istringstream ss((*iter).second); ss >> _outPathData;}
	else
	{ _outPathData.clear(); printf("Choose _outPathData = NULL as default\n");}
	
	
	
	_Cy5_ImageNRRD = _Cy5_Image+".nrrd";
	_TRI_ImageNRRD = _TRI_Image+".nrrd";
	_GFP_ImageNRRD = _GFP_Image+".nrrd";
	_DAP_ImageNRRD = _DAP_Image+".nrrd";
	
	
	// Print Parameters
	std::cout << std::endl << "This are the parameters";
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
	std::cout << std::endl << "_isSmall: " << _isSmall;
	std::cout << std::endl << "_segParams: " << _segParams;
	std::cout << std::endl << "_projectDefinition: " << _projectDefinition;
	std::cout << std::endl << "_outPath: " << _outPath;
	std::cout << std::endl << "_outPathDebug: " << _outPathDebug;
	std::cout << std::endl << "_outPathTemp: " << _outPathTemp;
	std::cout << std::endl << "_outPathData: " << _outPathData;
	
}

void ftkMainDarpaSegment::runSpliting()
{
	if( !_Cy5_Image.empty() )
	{
		rawImageType_8bit::Pointer ImageMontage_Cy5 = readImage< rawImageType_8bit >(_Cy5_ImageNRRD.c_str());
		_ImageMontage_Cy5Size = ImageMontage_Cy5->GetLargestPossibleRegion().GetSize();
		splitStore( ImageMontage_Cy5, _Cy5_Image );
	}
	if( !_TRI_Image.empty() )
	{
		rawImageType_8bit::Pointer ImageMontage_TRI = readImage< rawImageType_8bit >(_TRI_ImageNRRD.c_str());
		_ImageMontage_TRISize = ImageMontage_TRI->GetLargestPossibleRegion().GetSize();
		splitStore( ImageMontage_TRI, _TRI_Image );
	}
	if( !_GFP_Image.empty() )
	{
		rawImageType_8bit::Pointer ImageMontage_GFP = readImage< rawImageType_8bit >(_GFP_ImageNRRD.c_str());
		_ImageMontage_GFPSize = ImageMontage_GFP->GetLargestPossibleRegion().GetSize();
		splitStore( ImageMontage_GFP, _GFP_Image );
	}
	if( !_DAP_Image.empty() )
	{
		rawImageType_8bit::Pointer ImageMontage_DAP = readImage< rawImageType_8bit >(_DAP_ImageNRRD.c_str());
		_ImageMontage_DAPSize = ImageMontage_DAP->GetLargestPossibleRegion().GetSize();
		splitStore( ImageMontage_DAP, _DAP_Image );
	}
	if( !_DAP_Image.empty() && !_GFP_Image.empty() )
	{
		if( (_ImageMontage_GFPSize[0]!= _ImageMontage_DAPSize[0]) || (_ImageMontage_GFPSize[1]!= _ImageMontage_DAPSize[1]) || (_ImageMontage_GFPSize[2]!= _ImageMontage_DAPSize[2]) ) 
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




void ftkMainDarpaSegment::splitStore( rawImageType_8bit::Pointer ImageMontage, std::string nameInput )
{
	int found = nameInput.find_last_of("/\\");
	std::string nameInputNoPath = nameInput.substr(found+1);
	
	computeSplitConst( ImageMontage );
	itk::Size<3> ImageMontageSize = ImageMontage->GetLargestPossibleRegion().GetSize();
	
	int contadorSplit = 0;
#pragma omp parallel for collapse(3) //TEST
	for(int xco = 0; xco < _kx; xco++)
	{
		for(int yco = 0; yco < _ky; yco++)
		{
			for(int zco = 0; zco < _kz; zco++)
			{
// 				#pragma omp critical
// 				{

// 				}
				
				rawImageType_8bit::RegionType regionLocal_all = ComputeLocalRegionSplit( ImageMontageSize, xco, yco, zco );
				rawImageType_8bit::RegionType regionMontage_all = ComputeGlobalRegionSplit( ImageMontageSize, xco, yco, zco );
				
				rawImageType_8bit::Pointer imageLocal = rawImageType_8bit::New();
				imageLocal->SetRegions(regionLocal_all);
				imageLocal->Allocate();
				if(imageLocal->GetBufferPointer()==NULL)
					printf("Couldn't allocate memory - 4 .. going to crash now\n");

				typedef itk::ImageRegionIterator<rawImageType_8bit> IteratorType;
				IteratorType iterMontage(ImageMontage,regionMontage_all);
				IteratorType iterLocal(imageLocal,regionLocal_all);

				iterMontage.GoToBegin();
				iterLocal.GoToBegin();
				for(;!iterMontage.IsAtEnd();++iterMontage,++iterLocal)
				{
					iterLocal.Set(iterMontage.Get());
				}
				
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
				// !!! this has to be nrrd is just to test
				#pragma omp critical
				{
					++contadorSplit;
					std::cout<<std::endl<< "ImageSplit " << contadorSplit << " of " << _kx*_ky*_kz;
					std::string temp9 = _outPathTemp+"/"+nameInputNoPath+"_"+xStr+"_"+yStr+"_"+zStr+".nrrd";
					writeImage< rawImageType_8bit >(imageLocal,temp9.c_str());
				}
			}
		}
	}
}

void ftkMainDarpaSegment::runSegment(  )
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
	// Assume GFP or dapi always exist
	if( !_GFP_Image.empty() )
	{
		rawImageType_8bit::Pointer ImageMontage_GFP = readImage< rawImageType_8bit >(_GFP_ImageNRRD.c_str());
		ImageMontageRegion = ImageMontage_GFP->GetLargestPossibleRegion();
		ImageMontageSize = ImageMontageRegion.GetSize();
		computeSplitConst( ImageMontage_GFP );
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
	
// 	rawImageType_uint::Pointer imageLabelMontage = rawImageType_uint::New();
// 	rawImageType_uint::PointType originz;
// 	originz[0] = 0; 
// 	originz[1] = 0;
// 	originz[2] = 0;
// 	imageLabelMontage->SetOrigin( originz );
// 	rawImageType_uint::IndexType indexStich;
// 	indexStich.Fill(0);
// 	rawImageType_uint::RegionType regionz;
// 	regionz.SetSize ( ImageMontageSize  );
// 	regionz.SetIndex( indexStich );
// 	imageLabelMontage->SetRegions( regionz );
// 	imageLabelMontage->Allocate();
// 	imageLabelMontage->FillBuffer(0);
// 	try
// 	{
// 		imageLabelMontage->Allocate();
// 	}
// 	catch(itk::ExceptionObject &err)
// 	{
// 		std::cerr << "ExceptionObject caught!" <<std::endl;
// 		std::cerr << err << std::endl;
// 	}
	
	// GLOBAL TABLE
// 	vtkSmartPointer< vtkTable > tableLabelMontage;
// 	= vtkSmartPointer<vtkTable>::New();
// 	tableLabelMontage->Initialize();
// 	for(int c=0; c<(int)Table_Tiles[0]->GetNumberOfColumns(); ++c)
// 	{
// 		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
// 		column->SetName( Table_Tiles[0]->GetColumnName(c) );
// 		tableLabelMontage->AddColumn(column);
// 	}
	
	
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
				
				rawImageType_8bit::Pointer imageLocalCy5;
				rawImageType_8bit::Pointer imageLocalTRI;
				rawImageType_8bit::Pointer imageLocalGFP;
				rawImageType_8bit::Pointer imageLocalDAP;

				std::vector< rawImageType_8bit::Pointer > Images_Tiles;
				Images_Tiles.resize(4);

				
				std::vector< rawImageType_16bit::Pointer > Label_Tiles;
				Label_Tiles.resize(1);
				
				std::vector< vtkSmartPointer< vtkTable > > Table_Tiles;
				Table_Tiles.resize(1);
				
				std::vector< std::map< unsigned int, itk::Index<3> > > Centroids_Tiles;
				Centroids_Tiles.resize(1);
				
				// Reading part
				#pragma omp critical
				{
				if( !_GFP_Image.empty( ) && !_DAP_Image.empty() )
				{
					int foundGFP = _GFP_Image.find_last_of("/\\");
					std::string _GFP_ImageNoPath = _GFP_Image.substr(foundGFP+1);
					// !!! this has to be nrrd is just to test
					std::string tempGFP = _outPathTemp+"/"+_GFP_ImageNoPath+"_"+xStr+"_"+yStr+"_"+zStr+".nrrd";
					imageLocalGFP = readImage< rawImageType_8bit >(tempGFP.c_str());
					
					int foundDAP = _DAP_Image.find_last_of("/\\");
					std::string _DAP_ImageNoPath = _DAP_Image.substr(foundDAP+1);
					// !!! this has to be nrrd is just to test
					std::string tempDAP = _outPathTemp+"/"+_DAP_ImageNoPath+"_"+xStr+"_"+yStr+"_"+zStr+".nrrd";
					imageLocalDAP = readImage< rawImageType_8bit >(tempDAP.c_str());
				}
				}
				Images_Tiles[0] = imageLocalCy5;
				Images_Tiles[1] = imageLocalTRI;
				Images_Tiles[2] = imageLocalGFP;
				Images_Tiles[3] = imageLocalDAP;
				
// 				RunSegmentation(regionLocal_inside, Images_Tiles, Label_Tiles, Table_Tiles, Centroids_Tiles);
				
				Label_Tiles[0] = RunNuclearSegmentation( Images_Tiles[3] );
				Table_Tiles[0] = ComputeFeaturesAndAssociations( Images_Tiles, Label_Tiles );
				Centroids_Tiles[0] = GetLabelToCentroidMap(Table_Tiles[0]);
				
// 					// TEST PUT A NUMBER IN THE REGION OF INTERES
// 					IteratorType_16bit iterLocal16_2 = IteratorType_16bit(Label_Tiles[0],regionLocal_inside); iterLocal16_2.GoToBegin();
// 					for(;!iterLocal16_2.IsAtEnd();++iterLocal16_2)
// 						iterLocal16_2.Set(50000);
// 					// TEST SAVE SEGMENTATION
// 					std::string tempLABEL = _outPathDebugLevel2+"/segLabel_"+"_"+xStr+"_"+yStr+"_"+zStr+".tif";
// 					std::string tempTABLE = _outPathDebugLevel2+"/segTable_"+"_"+xStr+"_"+yStr+"_"+zStr+".txt";
// 					writeImage<rawImageType_16bit>(Label_Tiles[0],tempLABEL.c_str());
// 					ftk::SaveTable(tempTABLE, Table_Tiles[0]);
				
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
// /*					// TEST PUT A NUMBER IN THE REGION OF INTERES
// 					IteratorType_16bit iterLocal16_3 = IteratorType_16bit(Label_Tiles[0],regionLocal_inside); iterLocal16_3.GoToBegin();
// 					for(;!iterLocal16_3.IsAtEnd();++iterLocal16_3)
// 						iterLocal16_3.Set(50000);
// 					std::string tempLABELREGIO = _outPathDebugLevel2+"/segLabel_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REGION.nrrd";
// 					writeImage<rawImageType_16bit>(Label_Tiles[0],tempLABELREGIO.c_str());*/
				
// 				StichResults( imageLabelMontage, Label_Tiles, Table_Tiles, Centroids_Tiles );

// 				#pragma omp critical
// 				{
// 					++contadorStich;
// 					std::cout<<std::endl<< "\t\t--->>> ImageStich " << contadorStich << " of " << _kx*_ky*_kz;
// 					
// 					if( flagFirstStich == 1 )
// 					{
// 						flagFirstStich = 0;
// 						tableLabelMontage = vtkSmartPointer<vtkTable>::New();
// 						tableLabelMontage->Initialize();
// 						for(int c=0; c<(int)Table_Tiles[0]->GetNumberOfColumns(); ++c)
// 						{
// 							vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
// 							column->SetName( Table_Tiles[0]->GetColumnName(c) );
// 							tableLabelMontage->AddColumn(column);
// 						}
// 					}
// 					
// 					maxValueOld = maxValue;
// 					if((unsigned long long)Table_Tiles[0]->GetNumberOfRows() != 0)
// 					{
// 						std::cout << std::endl << "\t\tpiu11 The number of row: " << (int)Table_Tiles[0]->GetNumberOfRows() << " in " << contadorStich;
// 						for(int r=0; r<(int)Table_Tiles[0]->GetNumberOfRows(); ++r)
// 						{
// 							vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
// 							for(int c=0; c<(int)Table_Tiles[0]->GetNumberOfColumns(); ++c)
// 							{
// 								if(c == 0)
// 									model_data1->InsertNextValue(vtkVariant(Table_Tiles[0]->GetValue(r,c).ToUnsignedInt() + maxValue));
// 								else if(c == 1)
// 									model_data1->InsertNextValue(vtkVariant(Table_Tiles[0]->GetValue(r,c).ToInt() + regionMontage_all.GetIndex()[0]));
// 								else if(c == 2)
// 									model_data1->InsertNextValue(vtkVariant(Table_Tiles[0]->GetValue(r,c).ToInt() + regionMontage_all.GetIndex()[1]));
// 								else if(c == 3)
// 									model_data1->InsertNextValue(vtkVariant(Table_Tiles[0]->GetValue(r,c).ToInt() + regionMontage_all.GetIndex()[2]));
// 								else
// 									model_data1->InsertNextValue(Table_Tiles[0]->GetValue(r,c));
// 							}
// 							tableLabelMontage->InsertNextRow(model_data1);
// 						}
// 						std::cout << std::endl << "\t\tpiu11 The maxvaule before: " << maxValue << " ";
// 						maxValue = tableLabelMontage->GetValue((int)tableLabelMontage->GetNumberOfRows()-1, 0).ToUnsignedInt();
// 						std::cout << "after: " << maxValue << " in " << contadorStich;
// 					}
// 				
// 					IteratorType_uint iterMontage_1 = IteratorType_uint(imageLabelMontage,regionMontage_all);
// 					iterMontage_1.GoToBegin();
// 					
// 	// 				IteratorType_8bit iterLocal_1 = IteratorType_8bit(imageLocalDAP,regionLocal_inside);
// 	// 				iterLocal_1.GoToBegin();
// 					
// 					IteratorType_16bit iterLocal16_1 = IteratorType_16bit(Label_Tiles[0],regionLocal_all);
// 					iterLocal16_1.GoToBegin();
// 					
// 					for(;!iterMontage_1.IsAtEnd();++iterMontage_1)
// 					{
// 	// 					iterMontage_1.Set(iterLocal_1.Get());
// 	// 					++iterLocal8_1;
// 						if( iterLocal16_1.Get() != 0 )
// 							iterMontage_1.Set(iterLocal16_1.Get() + maxValueOld);
// 						++iterLocal16_1;
// 					}
// 					
// 					std::cout<<std::endl<< "\t\t--->>> ImageStichDone " << contadorStich << " of " << _kx*_ky*_kz;
// 				}
				#pragma omp critical
				{
					contadorImageDoneSegment++;
					std::cout<<std::endl<< "\t\t--->>> ImageDoneSegment " << contadorImageDoneSegment << " of " << _kx*_ky*_kz;
				}
			}
		}
	}
	
// 	std::string temp9 = _GFP_Image+"_label.nrrd";
// // 	std::cout << std::endl << "SOMA STORED " << temp9;
// 	writeImage<rawImageType_uint>(imageLabelMontage,temp9.c_str());
// 	
// 	std::string tempTABLEGLO = _GFP_Image+"_label_table.txt";
// // 	std::cout << std::endl << "SOMA STORED " << temp9;
// 	ftk::SaveTable(tempTABLEGLO, tableLabelMontage);
// 	
// 	// Soma Montage
// 	rawImageType_uint::Pointer imageSomaMontage = rawImageType_uint::New();
// 	itk::Size<3> im_size = imageLabelMontage->GetBufferedRegion().GetSize();
// 	rawImageType_uint::RegionType region;
// 	region.SetSize( ImageMontageRegion.GetSize() );
// 	region.SetIndex( ImageMontageRegion.GetIndex() );
// 	imageSomaMontage->SetRegions( region );
// 	imageSomaMontage->Allocate();
// 	imageSomaMontage->FillBuffer(0);
// 	try
// 	{
// 		imageSomaMontage->Update();
// 	}
// 	catch(itk::ExceptionObject &err)
// 	{
// 		std::cerr << "ExceptionObject caught!" <<std::endl;
// 		std::cerr << err << std::endl;
// 	}
// 	rawImageType_uint::PixelType * imageSomaArray = imageSomaMontage->GetBufferPointer();
// 	rawImageType_uint::PixelType * imageLabelArray = imageLabelMontage->GetBufferPointer();
// 
// 	unsigned long long sizeXY = ImageMontageSize[1] * ImageMontageSize[0];
// 	unsigned long long sizeX = ImageMontageSize[0];
// 	std::map<unsigned int, int> classMap;
// 	for(int row=0; row<(int)tableLabelMontage->GetNumberOfRows(); ++row)
// 	{
// 		classMap[tableLabelMontage->GetValue(row,0).ToUnsignedInt()] = tableLabelMontage->GetValueByName(row, "prediction_active_mg").ToInt();
// 	}
// 
// // 	#pragma omp parallel for collapse(3)
// 	for(int i=0; i<ImageMontageSize[2]; ++i)
// 	{
// 		for(int j=0; j<ImageMontageSize[1]; ++j)
// 		{
// 			for(int k=0; k<ImageMontageSize[0]; ++k)
// 			{
// 				unsigned long long offset = (i*sizeXY)+(j*sizeX)+k;
// 				if(classMap[imageLabelArray[offset]] == 1)
// 					imageSomaArray[offset] = imageLabelArray[offset];
// 			}
// 		}
// 	}	
// 	
// 	std::string nameSomaMontage = _GFP_Image+"_soma.nrrd";
// 	writeImage< rawImageType_uint >(imageSomaMontage,nameSomaMontage.c_str());
// 
// 	vtkSmartPointer<vtkTable> somaCentroidsTable = vtkSmartPointer<vtkTable>::New();
// 	somaCentroidsTable->Initialize();
// 
// 	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
// 	column->SetName( "somaCen_x" );
// 	somaCentroidsTable->AddColumn(column);
// 	column = vtkSmartPointer<vtkDoubleArray>::New();
// 	column->SetName( "somaCen_y" );
// 	somaCentroidsTable->AddColumn(column);
// 	column = vtkSmartPointer<vtkDoubleArray>::New();
// 	column->SetName( "somaCen_z" );
// 	somaCentroidsTable->AddColumn(column);
// 
// 	for(int row = 0; row < (int)tableLabelMontage->GetNumberOfRows(); ++row)
// 	{
// 		if(tableLabelMontage->GetValueByName(row, "prediction_active_mg").ToInt() == 1)
// 		{
// 			vtkSmartPointer<vtkVariantArray> centroid_data = vtkSmartPointer<vtkVariantArray>::New();
// 			centroid_data->InsertNextValue(tableLabelMontage->GetValue(row,1));
// 			centroid_data->InsertNextValue(tableLabelMontage->GetValue(row,2));
// 			centroid_data->InsertNextValue(tableLabelMontage->GetValue(row,3));
// 			somaCentroidsTable->InsertNextRow(centroid_data);			
// 		}
// // 		std::cout << std::endl << tableLabelMontage->GetValueByName(row, "prediction_active_mg").ToInt()
// 	}
// 
// 	for(int row = 0; row < (int)tableLabelMontage->GetNumberOfRows(); ++row)
// 	{
// 		if(tableLabelMontage->GetValueByName(row, "prediction_active_mg").ToInt() != 1)
// 		{
// 			tableLabelMontage->RemoveRow(row);
// 			--row;
// 		}
// 	}
// 	ftk::SaveTable(_GFP_Image + "_soma_table.txt", tableLabelMontage);
// 	ftk::SaveTable(_GFP_Image + "_soma_table_centroids.txt", somaCentroidsTable);
	
// 	// Save RGB image
// 
// 	RGBImageType_8bit::Pointer imageMontageRGBSoma = RGBImageType_8bit::New();
// 	RGBImageType_8bit::IndexType indexMontageRGB;
// 	indexMontageRGB.Fill(0);
// 	RGBImageType_8bit::RegionType regionMontageRGB;
// 	regionMontageRGB.SetIndex(indexMontageRGB);
// 	regionMontageRGB.SetSize(ImageMontageSize);
// 	imageMontageRGBSoma->SetRegions(regionMontageRGB);
// 	try
// 	{
// 		imageMontageRGBSoma->Allocate();
// 	}
// 	catch(itk::ExceptionObject &err)
// 	{
// 		std::cerr << "ExceptionObject caught!" <<std::endl;
// 		std::cerr << err << std::endl;
// 	}
// 	
// 	if( !_GFP_Image.empty() )
// 	{
// 		rawImageType_8bit::Pointer ImageMontage_GFP = readImage< rawImageType_8bit >(_GFP_ImageNRRD.c_str());
// 
// 		IteratorType_8bit itImageMontage_GFP(ImageMontage_GFP, ImageMontage_GFP->GetRequestedRegion()); itImageMontage_GFP.GoToBegin();
// 		IteratorType_uint itimageSomaMontage(imageSomaMontage, imageSomaMontage->GetRequestedRegion()); itimageSomaMontage.GoToBegin();
// 		RGBIteratorType_8bit itImageMontageRGBSoma(imageMontageRGBSoma, imageMontageRGBSoma->GetRequestedRegion()); itImageMontageRGBSoma.GoToBegin();
// 		
// 		for (itImageMontageRGBSoma.GoToBegin(); !itImageMontageRGBSoma.IsAtEnd(); ++itImageMontageRGBSoma, ++itImageMontage_GFP, ++itimageSomaMontage) 
// 		{
// 			RGBPixelType_8bit newPixel;
// 			
// 			newPixel.SetRed( 10/*itImageMontage_GFP.Get()*/ );
// 			newPixel.SetGreen(0);
// 			if( itimageSomaMontage.Get() != 0 )
// 				newPixel.SetBlue(255);
// 			else
// 				newPixel.SetBlue(0);
// 			itImageMontageRGBSoma.Set(newPixel);
// // 			std::cout << std::endl << itImageMontageRGBSoma.Get()[0];
// 		}
// 		
// 		std::string tempRGB = _outPathDebug+"/aRGB_Soma_GFP_"+".nrrd";
// 		writeImage<RGBImageType_8bit>(imageMontageRGBSoma,tempRGB.c_str());
// 	}
	
	
// 	// RGB z
// 	RGBImageType_8bit::Pointer imageMontageRGBSomaZ = RGBImageType_8bit::New();
// 	RGBImageType_8bit::IndexType indexMontageRGB;
// 	indexMontageRGB.Fill(0);
// 	RGBImageType_8bit::RegionType regionMontageRGB;
// 	regionMontageRGB.SetIndex(indexMontageRGB);
// 	regionMontageRGB.SetSize(ImageMontageSize);
// 	imageMontageRGBSoma->SetRegions(regionMontageRGB);
// 	try
// 	{
// 		imageMontageRGBSoma->Allocate();
// 	}
// 	catch(itk::ExceptionObject &err)
// 	{
// 		std::cerr << "ExceptionObject caught!" <<std::endl;
// 		std::cerr << err << std::endl;
// 	}
	
}

void ftkMainDarpaSegment::runStich(  )
{
	std::cout << std::endl << "HERE HERE_2";
	std::cout << std::endl << "HERE HERE_2";
	
	if( _isSmall == 1 )
	{
		itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1); // This one can not be changed
		itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1); // This one can chenga
		omp_set_nested(1);
		omp_set_max_active_levels(2);
		int num_threads = 1;
		omp_set_num_threads(num_threads);
	}
	
	rawImageType_8bit::RegionType ImageMontageRegion;
	itk::Size<3> ImageMontageSize;
	// Assume GFP or dapi always exist
	if( !_GFP_Image.empty() )
	{
		rawImageType_8bit::Pointer ImageMontage_GFP = readImage< rawImageType_8bit >(_GFP_ImageNRRD.c_str());
		ImageMontageRegion = ImageMontage_GFP->GetLargestPossibleRegion();
		ImageMontageSize = ImageMontageRegion.GetSize();
		computeSplitConst( ImageMontage_GFP );
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
	
	rawImageType_uint::Pointer imageLabelMontage = rawImageType_uint::New();
	rawImageType_uint::PointType originz;
	originz[0] = 0; 
	originz[1] = 0;
	originz[2] = 0;
	imageLabelMontage->SetOrigin( originz );
	rawImageType_uint::IndexType indexStich;
	indexStich.Fill(0);
	rawImageType_uint::RegionType regionz;
	regionz.SetSize ( ImageMontageSize  );
	regionz.SetIndex( indexStich );
	imageLabelMontage->SetRegions( regionz );
	
	try
	{
		imageLabelMontage->Allocate();
		imageLabelMontage->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}
	imageLabelMontage->FillBuffer(0);
	
	// GLOBAL TABLE
	vtkSmartPointer< vtkTable > tableLabelMontage;
	
	
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

				std::vector< rawImageType_8bit::Pointer > Images_Tiles;
				Images_Tiles.resize(4);

				
				std::vector< rawImageType_16bit::Pointer > Label_Tiles;
				Label_Tiles.resize(1);
				
				std::vector< vtkSmartPointer< vtkTable > > Table_Tiles;
				Table_Tiles.resize(1);
				
				std::vector< std::map< unsigned int, itk::Index<3> > > Centroids_Tiles;
				Centroids_Tiles.resize(1);
				
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
				std::string tempTABLERE = _outPathTemp+"/segTable_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REMO.txt";
				Table_Tiles[0] = ftk::LoadTable(tempTABLERE);
				
				std::string tempLABELRE = _outPathTemp+"/segLabel_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REMO.nrrd";
				Label_Tiles[0] = readImage<rawImageType_16bit>(tempLABELRE.c_str());

// 				std::string tempLABELRE = _outPathDebugLevel2+"/segLabel_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REMO.tif";
// 				ftkMainDarpa objftkMainDarpa;
// 				objftkMainDarpa.projectImage<rawImageType_16bit, rawImageType_16bit>( Label_Tiles[0], tempLABELRE, _outPathDebugLevel2, "ORG_RES_BIN" );

// /*					// TEST PUT A NUMBER IN THE REGION OF INTERES
// 					IteratorType_16bit iterLocal16_3 = IteratorType_16bit(Label_Tiles[0],regionLocal_inside); iterLocal16_3.GoToBegin();
// 					for(;!iterLocal16_3.IsAtEnd();++iterLocal16_3)
// 						iterLocal16_3.Set(50000);
// 					std::string tempLABELREGIO = _outPathDebugLevel2+"/segLabel_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REGION.nrrd";
// 					writeImage<rawImageType_16bit>(Label_Tiles[0],tempLABELREGIO.c_str());*/
				
// 				StichResults( imageLabelMontage, Label_Tiles, Table_Tiles, Centroids_Tiles );

// 				#pragma omp critical
// 				{
				++contadorStich;
				std::cout<<std::endl<< "\t\t--->>> ImageStich " << contadorStich << " of " << _kx*_ky*_kz;
				
				if( flagFirstStich == 1 )
				{
					flagFirstStich = 0;
					tableLabelMontage = vtkSmartPointer<vtkTable>::New();
					tableLabelMontage->Initialize();
					for(int c=0; c<(int)Table_Tiles[0]->GetNumberOfColumns(); ++c)
					{
						vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
						column->SetName( Table_Tiles[0]->GetColumnName(c) );
						tableLabelMontage->AddColumn(column);
					}
				}
				
				maxValueOld = maxValue;
				if((unsigned long long)Table_Tiles[0]->GetNumberOfRows() != 0)
				{
					std::cout << std::endl << "\t\tpiu11 The number of row: " << (int)Table_Tiles[0]->GetNumberOfRows() << " in " << contadorStich;
					for(int r=0; r<(int)Table_Tiles[0]->GetNumberOfRows(); ++r)
					{
						// Can be done here
						vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
						for(int c=0; c<(int)Table_Tiles[0]->GetNumberOfColumns(); ++c)
						{
							if(c == 0)
								model_data1->InsertNextValue(vtkVariant(Table_Tiles[0]->GetValue(r,c).ToUnsignedInt() + maxValue));
							else if(c == 1)
								model_data1->InsertNextValue(vtkVariant(Table_Tiles[0]->GetValue(r,c).ToInt() + regionMontage_all.GetIndex()[0]));
							else if(c == 2)
								model_data1->InsertNextValue(vtkVariant(Table_Tiles[0]->GetValue(r,c).ToInt() + regionMontage_all.GetIndex()[1]));
							else if(c == 3)
								model_data1->InsertNextValue(vtkVariant(Table_Tiles[0]->GetValue(r,c).ToInt() + regionMontage_all.GetIndex()[2]));
							else
								model_data1->InsertNextValue(Table_Tiles[0]->GetValue(r,c));
						}
						tableLabelMontage->InsertNextRow(model_data1);
					}
					std::cout << std::endl << "\t\tpiu11 The maxvaule before: " << maxValue << " ";
					maxValue = tableLabelMontage->GetValue((int)tableLabelMontage->GetNumberOfRows()-1, 0).ToUnsignedInt();
					std::cout << "after: " << maxValue << " in " << contadorStich;
				}
			
				IteratorType_uint iterMontage_1 = IteratorType_uint(imageLabelMontage,regionMontage_all);
				iterMontage_1.GoToBegin();
				
// 				IteratorType_8bit iterLocal_1 = IteratorType_8bit(imageLocalDAP,regionLocal_inside);
// 				iterLocal_1.GoToBegin();
				
				IteratorType_16bit iterLocal16_1 = IteratorType_16bit(Label_Tiles[0],regionLocal_all);
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
	
	std::string temp9 = _outPathData+"/label.nrrd";
// 	std::cout << std::endl << "SOMA STORED " << temp9;
	writeImage<rawImageType_uint>(imageLabelMontage,temp9.c_str());
	
	std::string tempTABLEGLO = _outPathData+"/label_table.txt";
// 	std::cout << std::endl << "SOMA STORED " << temp9;
	ftk::SaveTable(tempTABLEGLO, tableLabelMontage);
	
	// Soma Montage
	rawImageType_uint::Pointer imageSomaMontage = rawImageType_uint::New();
	itk::Size<3> im_size = imageLabelMontage->GetBufferedRegion().GetSize();
	rawImageType_uint::RegionType region;
	region.SetSize( ImageMontageRegion.GetSize() );
	region.SetIndex( ImageMontageRegion.GetIndex() );
	imageSomaMontage->SetRegions( region );
	imageSomaMontage->Allocate();
	imageSomaMontage->FillBuffer(0);
	try
	{
		imageSomaMontage->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}
	rawImageType_uint::PixelType * imageSomaArray = imageSomaMontage->GetBufferPointer();
	rawImageType_uint::PixelType * imageLabelArray = imageLabelMontage->GetBufferPointer();

	unsigned long long sizeXY = ImageMontageSize[1] * ImageMontageSize[0];
	unsigned long long sizeX = ImageMontageSize[0];
	std::map<unsigned int, int> classMap;
	for(int row=0; row<(int)tableLabelMontage->GetNumberOfRows(); ++row)
	{
		classMap[tableLabelMontage->GetValue(row,0).ToUnsignedInt()] = tableLabelMontage->GetValueByName(row, "prediction_active_mg").ToInt();
	}

	#pragma omp parallel for collapse(3)
	for(int i=0; i<ImageMontageSize[2]; ++i)
	{
		for(int j=0; j<ImageMontageSize[1]; ++j)
		{
			for(int k=0; k<ImageMontageSize[0]; ++k)
			{
				unsigned long long offset = (i*sizeXY)+(j*sizeX)+k;
				if( imageLabelArray[offset] != 0 )
					if(classMap[imageLabelArray[offset]] == 1)
						imageSomaArray[offset] = imageLabelArray[offset];
			}
		}
	}	
	
	std::string nameSomaMontage = _outPathData+"/soma.nrrd";
	writeImage< rawImageType_uint >(imageSomaMontage,nameSomaMontage.c_str());

	vtkSmartPointer<vtkTable> somaCentroidsTable = vtkSmartPointer<vtkTable>::New();
	somaCentroidsTable->Initialize();

	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "somaCen_x" );
	somaCentroidsTable->AddColumn(column);
	column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "somaCen_y" );
	somaCentroidsTable->AddColumn(column);
	column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "somaCen_z" );
	somaCentroidsTable->AddColumn(column);

	for(int row = 0; row < (int)tableLabelMontage->GetNumberOfRows(); ++row)
	{
		if(tableLabelMontage->GetValueByName(row, "prediction_active_mg").ToInt() == 1)
		{
			vtkSmartPointer<vtkVariantArray> centroid_data = vtkSmartPointer<vtkVariantArray>::New();
			centroid_data->InsertNextValue(tableLabelMontage->GetValue(row,1));
			centroid_data->InsertNextValue(tableLabelMontage->GetValue(row,2));
			centroid_data->InsertNextValue(tableLabelMontage->GetValue(row,3));
			somaCentroidsTable->InsertNextRow(centroid_data);			
		}
// 		std::cout << std::endl << tableLabelMontage->GetValueByName(row, "prediction_active_mg").ToInt()
	}

	for(int row = 0; row < (int)tableLabelMontage->GetNumberOfRows(); ++row)
	{
		if(tableLabelMontage->GetValueByName(row, "prediction_active_mg").ToInt() != 1)
		{
			tableLabelMontage->RemoveRow(row);
			--row;
		}
	}
	ftk::SaveTable(_outPathData + "/soma_table.txt", tableLabelMontage);
	ftk::SaveTable(_outPathData + "/soma_table_centroids.txt", somaCentroidsTable);
	
// 	// Save RGB image
// 
// 	RGBImageType_8bit::Pointer imageMontageRGBSoma = RGBImageType_8bit::New();
// 	RGBImageType_8bit::IndexType indexMontageRGB;
// 	indexMontageRGB.Fill(0);
// 	RGBImageType_8bit::RegionType regionMontageRGB;
// 	regionMontageRGB.SetIndex(indexMontageRGB);
// 	regionMontageRGB.SetSize(ImageMontageSize);
// 	imageMontageRGBSoma->SetRegions(regionMontageRGB);
// 	try
// 	{
// 		imageMontageRGBSoma->Allocate();
// 	}
// 	catch(itk::ExceptionObject &err)
// 	{
// 		std::cerr << "ExceptionObject caught!" <<std::endl;
// 		std::cerr << err << std::endl;
// 	}
// 	
// 	if( !_GFP_Image.empty() )
// 	{
// 		rawImageType_8bit::Pointer ImageMontage_GFP = readImage< rawImageType_8bit >(_GFP_ImageNRRD.c_str());
// 
// 		IteratorType_8bit itImageMontage_GFP(ImageMontage_GFP, ImageMontage_GFP->GetRequestedRegion()); itImageMontage_GFP.GoToBegin();
// 		IteratorType_uint itimageSomaMontage(imageSomaMontage, imageSomaMontage->GetRequestedRegion()); itimageSomaMontage.GoToBegin();
// 		RGBIteratorType_8bit itImageMontageRGBSoma(imageMontageRGBSoma, imageMontageRGBSoma->GetRequestedRegion()); itImageMontageRGBSoma.GoToBegin();
// 		
// 		for (itImageMontageRGBSoma.GoToBegin(); !itImageMontageRGBSoma.IsAtEnd(); ++itImageMontageRGBSoma, ++itImageMontage_GFP, ++itimageSomaMontage) 
// 		{
// 			RGBPixelType_8bit newPixel;
// 			
// 			newPixel.SetRed( 10/*itImageMontage_GFP.Get()*/ );
// 			newPixel.SetGreen(0);
// 			if( itimageSomaMontage.Get() != 0 )
// 				newPixel.SetBlue(255);
// 			else
// 				newPixel.SetBlue(0);
// 			itImageMontageRGBSoma.Set(newPixel);
// // 			std::cout << std::endl << itImageMontageRGBSoma.Get()[0];
// 		}
// 		
// 		std::string tempRGB = _outPathDebug+"/aRGB_Soma_GFP_"+".nrrd";
// 		writeImage<RGBImageType_8bit>(imageMontageRGBSoma,tempRGB.c_str());
// 	}
	
	
// 	// RGB z
// 	RGBImageType_8bit::Pointer imageMontageRGBSomaZ = RGBImageType_8bit::New();
// 	RGBImageType_8bit::IndexType indexMontageRGB;
// 	indexMontageRGB.Fill(0);
// 	RGBImageType_8bit::RegionType regionMontageRGB;
// 	regionMontageRGB.SetIndex(indexMontageRGB);
// 	regionMontageRGB.SetSize(ImageMontageSize);
// 	imageMontageRGBSoma->SetRegions(regionMontageRGB);
// 	try
// 	{
// 		imageMontageRGBSoma->Allocate();
// 	}
// 	catch(itk::ExceptionObject &err)
// 	{
// 		std::cerr << "ExceptionObject caught!" <<std::endl;
// 		std::cerr << err << std::endl;
// 	}
	
}


void ftkMainDarpaSegment::computeSplitConst( rawImageType_8bit::Pointer ImageMontage )
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

void ftkMainDarpaSegment::RunSegmentation(rawImageType_8bit::RegionType regionLocal_inside, std::vector< rawImageType_8bit::Pointer >& Images_Tiles, std::vector< rawImageType_16bit::Pointer >& Label_Tiles, std::vector< vtkSmartPointer< vtkTable > >& Table_Tiles, std::vector< std::map< unsigned int, itk::Index<3> > >& Centroids_Tiles)
{
	rawImageType_8bit::Pointer imageLocalCy5 = Images_Tiles[0];
	rawImageType_8bit::Pointer imageLocalTRI = Images_Tiles[1];
	rawImageType_8bit::Pointer imageLocalGFP = Images_Tiles[2];
	rawImageType_8bit::Pointer imageLocalDAP = Images_Tiles[3];
	
	Label_Tiles[0] = RunNuclearSegmentation( imageLocalDAP );
	
//  	Table_Tiles[0] = ComputeFeaturesAndAssociations( Images_Tiles, Label_Tiles );
// 	*myCentroids = GetLabelToCentroidMap(myTable[col]);
	
	

	
	
}


rawImageType_16bit::Pointer ftkMainDarpaSegment::RunNuclearSegmentation(rawImageType_8bit::Pointer inpImg)
{
	rawImageType_8bit::PointType origin;
	origin[0] = 0; 
	origin[1] = 0;
	origin[2] = 0;
	inpImg->SetOrigin( origin );

	itk::SizeValueType size[3];
	size[0] = inpImg->GetLargestPossibleRegion().GetSize()[0];
	size[1] = inpImg->GetLargestPossibleRegion().GetSize()[1];
	size[2] = inpImg->GetLargestPossibleRegion().GetSize()[2];

	unsigned char *in_Image;
	#pragma omp critical
	{
		in_Image = (unsigned char *) malloc (size[0]*size[1]*size[2]);
		if( in_Image == NULL )
		{
			std::cerr<<"Nucleus Seg failed because malloc failed\n";
// 			return NULL;
		}
	}

	typedef itk::ImageRegionIterator< rawImageType_8bit > IteratorType;
	IteratorType pix_buf( inpImg, inpImg->GetRequestedRegion() );
	itk::SizeValueType ind=0;
	for ( pix_buf.GoToBegin(); !pix_buf.IsAtEnd(); ++pix_buf, ++ind )
	{
		in_Image[ind]=(pix_buf.Get());
	}
	yousef_nucleus_seg *NucleusSeg = new yousef_nucleus_seg();
	// Can crash here
	#pragma omp critical
	{
		NucleusSeg->readParametersFromFile(_segParams.c_str());
	}
	NucleusSeg->setDataImage(in_Image,size[0],size[1],size[2],"test.mhd");
	
// 	#pragma omp critical
// 	{
		NucleusSeg->runBinarization();
// 	}
	std::cout<<std::endl << "\t\t Done binarization:" <<std::flush;

//  	#pragma omp critical
	NucleusSeg->runSeedDetection();
	std::cout<<std::endl << "\t\t\t Done Seed Detect: " <<std::flush;



	NucleusSeg->runClustering();
	std::cout<<std::endl << "\t\t\t\t Done Clus: ";
	
	unsigned short *output_img;

	output_img=NucleusSeg->getClustImage();

	rawImageType_16bit::Pointer image = rawImageType_16bit::New();
	rawImageType_16bit::PointType origin1;
	origin1[0] = 0; //Is this OK?
	origin1[1] = 0;
	origin1[2] = 0;
	image->SetOrigin( origin1 );
	rawImageType_16bit::IndexType start1;
	start1[0] = 0;
	start1[1] = 0;
	start1[2] = 0;
	rawImageType_16bit::SizeType size1;
	size1[0] = size[0];
	size1[1] = size[1];
	size1[2] = size[2];
	rawImageType_16bit::RegionType region;
	region.SetSize ( size1  );
	region.SetIndex( start1 );
	image->SetRegions( region );
	image->Allocate();
	image->FillBuffer(0);
	image->Update();

	typedef itk::ImageRegionIteratorWithIndex< rawImageType_16bit > IteratorType2;
	IteratorType2 iterator1(image,image->GetRequestedRegion());
	for(itk::SizeValueType i=0; i<(size[0]*size[1]*size[2]); ++i)
	{
		//unsigned short val = (unsigned short)output_img[i];
		iterator1.Set(output_img[i]);
		++iterator1;
	}

	delete NucleusSeg;
	free(in_Image);
	
	std::cout<<std::endl << "\t\t\t\t HERE2: Done SEG: ";
	
	return image;
}


vtkSmartPointer< vtkTable > ftkMainDarpaSegment::ComputeFeaturesAndAssociations(std::vector< rawImageType_8bit::Pointer >& Images_Tiles, std::vector< rawImageType_16bit::Pointer >& Label_Tiles )
{
	ftk::ProjectDefinition projectDef;
	projectDef.Load( _projectDefinition.c_str() );

	vtkSmartPointer<vtkTable> table = NULL;
	ftk::Image::Pointer nucSegRes = NULL;
	ftk::ProjectProcessor * myProc = new ftk::ProjectProcessor();

	itk::SizeValueType tileSize[3];
	tileSize[0] = Images_Tiles[3]->GetLargestPossibleRegion().GetSize()[0];
	tileSize[1] = Images_Tiles[3]->GetLargestPossibleRegion().GetSize()[1];
	tileSize[2] = Images_Tiles[3]->GetLargestPossibleRegion().GetSize()[2];
	
	std::cout<<std::endl<<"\t\t TRTRTR: "<<_projectDefinition<<" "<<tileSize[0]<<" "<<tileSize[1]<<" "<<tileSize[2];

	ftk::Image::Pointer sourceImages = ftk::Image::New();
	std::vector< unsigned char > color;
	color.push_back(255); color.push_back(0); color.push_back(0);
	sourceImages->AppendChannelFromData3D( Images_Tiles[3]->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), tileSize[0], tileSize[1], tileSize[2], "dapi", color, true);
	color[0] = 0; color[1] = 255; color[2] = 0;
	sourceImages->AppendChannelFromData3D( Images_Tiles[2]->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), tileSize[0], tileSize[1], tileSize[2], "gfp", color, true);
// 	color[0] = 255; color[1] = 255; color[2] = 0;
// 	sourceImages->AppendChannelFromData3D( cy5Img->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), tileSize[0], tileSize[1], tileSize[2], "cy5", color, true);

	ftk::Image::Pointer labelImage = ftk::Image::New();
	color[0] = 255; color[1] = 255; color[2] = 255;
	labelImage->AppendChannelFromData3D( Label_Tiles[0]->GetBufferPointer(), itk::ImageIOBase::USHORT, sizeof(unsigned short), tileSize[0], tileSize[1], tileSize[2], "dapi", color, true);

	myProc->SetInputImage(sourceImages);
	myProc->SetDefinition(&projectDef);
	myProc->SetOutputImage(labelImage);
	myProc->Initialize();

	while(!myProc->DoneProcessing())
		myProc->ProcessNext();

	nucSegRes = myProc->GetOutputImage();
	table = myProc->GetTable();

	//projectDef.Write(fileName);	
	delete myProc;

	return table;
}

std::map< unsigned int, itk::Index<3> > ftkMainDarpaSegment::GetLabelToCentroidMap( vtkSmartPointer< vtkTable > table)
{
	std::cout << std::endl << (int)table->GetNumberOfRows();
	std::cout << std::endl << (int)table->GetNumberOfRows();
	
	std::map< unsigned int, itk::Index<3> > centroidMap;
	for (int row=0; row<(int)table->GetNumberOfRows(); ++row)
	{
		unsigned int id = table->GetValue(row, 0).ToUnsignedInt();
		itk::Index<3> indx;
		indx[0] = table->GetValueByName(row, "centroid_x").ToUnsignedInt();
		indx[1] = table->GetValueByName(row, "centroid_y").ToUnsignedInt();
		indx[2] = table->GetValueByName(row, "centroid_z").ToUnsignedInt();
		centroidMap[id] = indx;
	}

	return centroidMap;
}


void ftkMainDarpaSegment::RemoveLabelNearBorder(rawImageType_8bit::RegionType regionLocal_inside,std::vector< rawImageType_16bit::Pointer >& Label_Tiles, std::vector< vtkSmartPointer< vtkTable > >& Table_Tiles, std::vector< std::map< unsigned int, itk::Index<3> > >& Centroids_Tiles )
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


rawImageType_8bit::RegionType ftkMainDarpaSegment::ComputeLocalRegionSegment( itk::Size<3> ImageMontageSize, int xco, int yco, int zco )
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

rawImageType_uint::RegionType ftkMainDarpaSegment::ComputeGlobalRegionSegment( itk::Size<3> ImageMontageSize, int xco, int yco, int zco )
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

rawImageType_8bit::RegionType ftkMainDarpaSegment::ComputeLocalRegionSplit( itk::Size<3> ImageMontageSize, int xco, int yco, int zco )
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

rawImageType_8bit::RegionType ftkMainDarpaSegment::ComputeGlobalRegionSplit( itk::Size<3> ImageMontageSize, int xco, int yco, int zco )
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

// void ftkMainDarpaSegment::runTracing( )
// {
// 	MultipleNeuronTracer * MNT = new MultipleNeuronTracer()
// 	MNT->LoadParameters(_optionsMNT,5);
// 	
// // 				MNT->setNDXImage(NDXImage);
// 
// // 				MNT->setDiceSize(size);
// // 				MNT->setDiceIndex(start);
// // 				
// // 				MNT->LoadCurvImage_1(img_trace, 0);
// 	MNT->LoadCurvImage_2(img_trace);
// 	MNT->ReadStartPoints_1(soma_Table, 0);
// // 				MNT->SetCostThreshold(1000);
// 	MNT->SetCostThreshold(MNT->cost_threshold);
// 	MNT->LoadSomaImage_1(img_soma);
// 	
// 	bool flagLog = false;
// 	MNT->setFlagOutLog(flagLog);
// 	MNT->RunTracing();
// // 
// 	x = min(tileSizeX/2, x);
// 	y = min(tileSizeY/2, y);
// 	if(size_gfp_montage[2] > tileSizeZ)
// 		z = min(tileSizeZ/2, z);
// // 
// //
// 	vtkSmartPointer< vtkTable > swcTable = MNT->GetSWCTable(0);
// 	std::string swcFilename = filePath + "/TracesAndSomas/Trace_" + ssx.str() + "_" + ssy.str() + "_" + ssz.str() + "_ANT.swc";
// 	WriteCenterTrace(swcTable, x, y, z, swcFilename);
// // 				
// 	#pragma omp critical
// 	{
// 		outfile << "\t<File\tFileName=\"Trace_" << ssx.str() << "_" << ssy.str() << "_" << ssz.str() << "_ANT.swc\"\tType=\"Trace\"\ttX=\"" << ssx_off.str() << "\"\ttY=\"" << ssy_off.str() << "\"\ttZ=\"" << ssz_off.str() << "\"/>\n";
// 	}
// 	
// // 				vtkSmartPointer< vtkTable > swcTable = MNT->GetSWCTable(0);
// 	std::string swcFilenameDivided = filePath + "/TracesAndSomasDivided/Trace_BigTile_" + srr + "_" + ssx.str() + "_" + ssy.str() + "_" + ssz.str() + "_ANT.swc";
// 	WriteCenterTrace(swcTable, x, y, z, swcFilenameDivided);
// }
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 

void ftkMainDarpaSegment::runStichOneAtTheTime(  )
{
	std::cout << std::endl << "HERE HERE_2";
	std::cout << std::endl << "HERE HERE_2";
	
	if( _isSmall == 1 )
	{
		itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1); // This one can not be changed
		itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1); // This one can chenga
		omp_set_nested(1);
		omp_set_max_active_levels(2);
		int num_threads = 1;
		omp_set_num_threads(num_threads);
	}
	
	rawImageType_8bit::RegionType ImageMontageRegion;
	itk::Size<3> ImageMontageSize;
	// Assume GFP or dapi always exist
	if( !_GFP_Image.empty() )
	{
		rawImageType_8bit::Pointer ImageMontage_GFP = readImage< rawImageType_8bit >(_GFP_ImageNRRD.c_str());
		ImageMontageRegion = ImageMontage_GFP->GetLargestPossibleRegion();
		ImageMontageSize = ImageMontageRegion.GetSize();
		computeSplitConst( ImageMontage_GFP );
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
	
	rawImageType_uint::Pointer imageLabelMontage = rawImageType_uint::New();
	rawImageType_uint::PointType originz;
	originz[0] = 0; 
	originz[1] = 0;
	originz[2] = 0;
	imageLabelMontage->SetOrigin( originz );
	rawImageType_uint::IndexType indexStich;
	indexStich.Fill(0);
	rawImageType_uint::RegionType regionz;
	regionz.SetSize ( ImageMontageSize  );
	regionz.SetIndex( indexStich );
	imageLabelMontage->SetRegions( regionz );
	imageLabelMontage->Allocate();
	imageLabelMontage->FillBuffer(0);
	try
	{
		imageLabelMontage->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}
	
	rawImageType_uint::Pointer imageLabelMontageSomaNew = rawImageType_uint::New();
	rawImageType_uint::PointType originzSoma;
	originzSoma[0] = 0; 
	originzSoma[1] = 0;
	originzSoma[2] = 0;
	imageLabelMontage->SetOrigin( originzSoma );
	rawImageType_uint::IndexType indexStichSoma;
	indexStichSoma.Fill(0);
	rawImageType_uint::RegionType regionzSoma;
	regionzSoma.SetSize ( ImageMontageSize  );
	regionzSoma.SetIndex( indexStichSoma );
	imageLabelMontageSomaNew->SetRegions( regionzSoma );
	imageLabelMontageSomaNew->Allocate();
	imageLabelMontageSomaNew->FillBuffer(0);
	try
	{
		imageLabelMontageSomaNew->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}
	
	// GLOBAL TABLE
	vtkSmartPointer< vtkTable > tableLabelMontage;
	vtkSmartPointer< vtkTable > tableLabelMontageSomaNew;
	
	
	int contadorSegment = 0;
	int contadorStich = 0;
	unsigned int maxValue = 0;
	unsigned int maxValueOld = 0;
	int flagFirstStich = 1;
	int firstSomaStich = 1;
	
// #pragma omp parallel for collapse(3) num_threads(_num_threads) schedule(dynamic, 1)
	for(int xco = 0; xco < _kx; xco++)
	{
		for(int yco = 0; yco < _ky; yco++)
		{
			for(int zco = 0; zco < _kz; zco++)
			{
// 				#pragma omp critical
// 				{
// 				++contadorSegment;
// 				std::cout<<std::endl<< "\t\t--->>> ImageStich " << contadorSegment << " of " << _kx*_ky*_kz;
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

				std::vector< rawImageType_8bit::Pointer > Images_Tiles;
				Images_Tiles.resize(4);

				
				std::vector< rawImageType_16bit::Pointer > Label_Tiles;
				Label_Tiles.resize(1);
				
				std::vector< vtkSmartPointer< vtkTable > > Table_Tiles;
				Table_Tiles.resize(1);
				
				std::vector< std::map< unsigned int, itk::Index<3> > > Centroids_Tiles;
				Centroids_Tiles.resize(1);
				
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
				std::string tempTABLERE = _outPathTemp+"/segTable_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REMO.txt";
				Table_Tiles[0] = ftk::LoadTable(tempTABLERE);
				
				std::string tempLABELRE = _outPathTemp+"/segLabel_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REMO.nrrd";
				Label_Tiles[0] = readImage<rawImageType_16bit>(tempLABELRE.c_str());
				
				
	rawImageType_16bit::Pointer imageSomaLocal = rawImageType_16bit::New();
	rawImageType_16bit::IndexType indexSomaLocal;
	indexSomaLocal.Fill(0);
	rawImageType_16bit::RegionType regionzSomaLocal;
	regionzSomaLocal.SetSize ( regionLocal_all.GetSize()  );
	regionzSomaLocal.SetIndex( indexSomaLocal );
	imageSomaLocal->SetRegions( regionzSomaLocal );
	imageSomaLocal->Allocate();
	imageSomaLocal->FillBuffer(0);
	try
	{
		imageSomaLocal->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}
				
// 				rawDuplicatorType_16::Pointer somaDuplicator = rawDuplicatorType_16::New();
// 				somaDuplicator->SetInputImage(Label_Tiles[0]);
// 				somaDuplicator->Update();
// 				rawImageType_16bit::Pointer imageSomaLocal = somaDuplicator->GetOutput();

// 				std::string tempLABELRE = _outPathDebugLevel2+"/segLabel_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REMO.tif";
// 				ftkMainDarpa objftkMainDarpa;
// 				objftkMainDarpa.projectImage<rawImageType_16bit, rawImageType_16bit>( Label_Tiles[0], tempLABELRE, _outPathDebugLevel2, "ORG_RES_BIN" );

// /*					// TEST PUT A NUMBER IN THE REGION OF INTERES
// 					IteratorType_16bit iterLocal16_3 = IteratorType_16bit(Label_Tiles[0],regionLocal_inside); iterLocal16_3.GoToBegin();
// 					for(;!iterLocal16_3.IsAtEnd();++iterLocal16_3)
// 						iterLocal16_3.Set(50000);
// 					std::string tempLABELREGIO = _outPathDebugLevel2+"/segLabel_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REGION.nrrd";
// 					writeImage<rawImageType_16bit>(Label_Tiles[0],tempLABELREGIO.c_str());*/
				
// 				StichResults( imageLabelMontage, Label_Tiles, Table_Tiles, Centroids_Tiles );

// 				#pragma omp critical
// 				{
				++contadorStich;
				std::cout<<std::endl<< "\t\t--->>> ImageStich " << contadorStich << " of " << _kx*_ky*_kz;
				
				if( flagFirstStich == 1 )
				{
					flagFirstStich = 0;
					tableLabelMontage = vtkSmartPointer<vtkTable>::New();
					tableLabelMontage->Initialize();
					for(int c=0; c<(int)Table_Tiles[0]->GetNumberOfColumns(); ++c)
					{
						vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
						column->SetName( Table_Tiles[0]->GetColumnName(c) );
						tableLabelMontage->AddColumn(column);
					}

				}
				vtkSmartPointer< vtkTable > tableLabelSoma = vtkSmartPointer<vtkTable>::New();
				tableLabelSoma->Initialize();
				for(int c=0; c<(int)Table_Tiles[0]->GetNumberOfColumns(); ++c)
				{
					vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
					column->SetName( Table_Tiles[0]->GetColumnName(c) );
					tableLabelSoma->AddColumn(column);
				}
				std::map<unsigned int, int> classMap;
				
				maxValueOld = maxValue;
				if((unsigned long long)Table_Tiles[0]->GetNumberOfRows() != 0)
				{
					std::cout << std::endl << "\t\tpiu11 The number of row: " << (int)Table_Tiles[0]->GetNumberOfRows() << " in " << contadorStich;
					for(int r=0; r<(int)Table_Tiles[0]->GetNumberOfRows(); ++r)
					{
						// Can be done here
						if( Table_Tiles[0]->GetValue(r,0).ToUnsignedInt() + maxValue == 1526 )
							std::cout << std::endl << "LAVELA: " << maxValue << " " << "_" << xStr << "_" << yStr << "_" << zStr << " " << r;
						vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
						for(int c=0; c<(int)Table_Tiles[0]->GetNumberOfColumns(); ++c)
						{
							if(c == 0)
								model_data1->InsertNextValue(vtkVariant(Table_Tiles[0]->GetValue(r,c).ToUnsignedInt() + maxValue));
							else if(c == 1)
								model_data1->InsertNextValue(vtkVariant(Table_Tiles[0]->GetValue(r,c).ToInt() + regionMontage_all.GetIndex()[0]));
							else if(c == 2)
								model_data1->InsertNextValue(vtkVariant(Table_Tiles[0]->GetValue(r,c).ToInt() + regionMontage_all.GetIndex()[1]));
							else if(c == 3)
								model_data1->InsertNextValue(vtkVariant(Table_Tiles[0]->GetValue(r,c).ToInt() + regionMontage_all.GetIndex()[2]));
							else
								model_data1->InsertNextValue(Table_Tiles[0]->GetValue(r,c));
						}
						tableLabelMontage->InsertNextRow(model_data1);
						
						if( Table_Tiles[0]->GetValueByName(r, "prediction_active_mg").ToInt() == 1 )
						{
							vtkSmartPointer<vtkVariantArray> model_data_soma = vtkSmartPointer<vtkVariantArray>::New();
							for(int c=0; c<(int)Table_Tiles[0]->GetNumberOfColumns(); ++c)
							{
								if(c == 0)
									model_data_soma->InsertNextValue(vtkVariant(Table_Tiles[0]->GetValue(r,c).ToUnsignedInt() + maxValue));
								else if(c == 1)
									model_data_soma->InsertNextValue(vtkVariant(Table_Tiles[0]->GetValue(r,c).ToInt()/* + regionMontage_all.GetIndex()[0]*/));
								else if(c == 2)
									model_data_soma->InsertNextValue(vtkVariant(Table_Tiles[0]->GetValue(r,c).ToInt()/* + regionMontage_all.GetIndex()[1]*/));
								else if(c == 3)
									model_data_soma->InsertNextValue(vtkVariant(Table_Tiles[0]->GetValue(r,c).ToInt() /*+ regionMontage_all.GetIndex()[2]*/));
								else
									model_data_soma->InsertNextValue(Table_Tiles[0]->GetValue(r,c));
							}
							tableLabelSoma->InsertNextRow(model_data_soma);
						}
						
						classMap[Table_Tiles[0]->GetValue(r,0).ToUnsignedInt()] = Table_Tiles[0]->GetValueByName(r, "prediction_active_mg").ToInt();
	
					}
					std::cout << std::endl << "\t\tpiu11 The maxvaule before: " << maxValue << " ";
					maxValue = tableLabelMontage->GetValue((int)tableLabelMontage->GetNumberOfRows()-1, 0).ToUnsignedInt();
					std::cout << "after: " << maxValue << " in " << contadorStich;
				}
			
				IteratorType_uint iterMontage_1 = IteratorType_uint(imageLabelMontage,regionMontage_all);
				iterMontage_1.GoToBegin();
				
// 				IteratorType_8bit iterLocal_1 = IteratorType_8bit(imageLocalDAP,regionLocal_inside);
// 				iterLocal_1.GoToBegin();
				
				IteratorType_16bit iterLocal16_1 = IteratorType_16bit(Label_Tiles[0],regionLocal_all);
				iterLocal16_1.GoToBegin();
				
				IteratorType_16bit iterLocal16_soma = IteratorType_16bit(imageSomaLocal,regionLocal_all);
				iterLocal16_soma.GoToBegin();
				
				
				for(;!iterMontage_1.IsAtEnd();++iterMontage_1)
				{
// 					iterLocal16_soma.Set(0);
					if( iterLocal16_1.Get() != 0 )
						iterMontage_1.Set(iterLocal16_1.Get() + maxValueOld);
					if( iterLocal16_1.Get() != 0 )
						if( classMap[iterLocal16_1.Get()]==1) // If soma copy
							iterLocal16_soma.Set(iterLocal16_1.Get() + maxValueOld);
					
					++iterLocal16_1;
					++iterLocal16_soma;
				}
				std::string tempTABLERE_SOMA = _outPathTemp+"/segTable_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REMO_SOMA.txt";
				ftk::SaveTable(tempTABLERE_SOMA, tableLabelSoma);
				
				std::string tempLABELRE_SOMA = _outPathTemp+"/segLabel_"+"_"+xStr+"_"+yStr+"_"+zStr+"_REMO_SOMA.nrrd";
				writeImage<rawImageType_16bit>(imageSomaLocal,tempLABELRE_SOMA.c_str());
				
				
				IteratorType_uint iterMontageSomaNew = IteratorType_uint(imageLabelMontageSomaNew,regionMontage_all);
				iterMontageSomaNew.GoToBegin();
				
				IteratorType_16bit iterimageSomaLocal = IteratorType_16bit(imageSomaLocal,regionLocal_all);
				iterimageSomaLocal.GoToBegin();
				
				for(;!iterMontageSomaNew.IsAtEnd();++iterMontageSomaNew)
				{
					if( iterimageSomaLocal.Get()!= 0)
						iterMontageSomaNew.Set( iterimageSomaLocal.Get() );
					++iterimageSomaLocal;
				}
				
				if( firstSomaStich == 1 )
				{
					firstSomaStich = 0;
					tableLabelMontageSomaNew = vtkSmartPointer<vtkTable>::New();
					tableLabelMontageSomaNew->Initialize();
					for(int c=0; c<(int)tableLabelSoma->GetNumberOfColumns(); ++c)
					{
						vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
						column->SetName( tableLabelSoma->GetColumnName(c) );
						tableLabelMontageSomaNew->AddColumn(column);
					}
				}
				if((unsigned long long)tableLabelSoma->GetNumberOfRows() != 0)
				{
					std::cout << std::endl << "\t\tpiu11 The number of row: " << (int)tableLabelSoma->GetNumberOfRows() << " in " << contadorStich;
					for(int r=0; r<(int)tableLabelSoma->GetNumberOfRows(); ++r)
					{
						// Can be done here
						vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
						for(int c=0; c<(int)Table_Tiles[0]->GetNumberOfColumns(); ++c)
						{
							if(c == 0)
								model_data1->InsertNextValue(vtkVariant(tableLabelSoma->GetValue(r,c).ToUnsignedInt()));
							else if(c == 1)
								model_data1->InsertNextValue(vtkVariant(tableLabelSoma->GetValue(r,c).ToInt() + regionMontage_all.GetIndex()[0]));
							else if(c == 2)
								model_data1->InsertNextValue(vtkVariant(tableLabelSoma->GetValue(r,c).ToInt() + regionMontage_all.GetIndex()[1]));
							else if(c == 3)
								model_data1->InsertNextValue(vtkVariant(tableLabelSoma->GetValue(r,c).ToInt() + regionMontage_all.GetIndex()[2]));
							else
								model_data1->InsertNextValue(tableLabelSoma->GetValue(r,c));
						}
						tableLabelMontageSomaNew->InsertNextRow(model_data1);
					}
				}
				
				
				
				std::cout<<std::endl<< "\t\t--->>> ImageStichDone " << contadorStich << " of " << _kx*_ky*_kz;
			}
		}
	}
	
	std::string temp9 = _outPathData+"/label.nrrd";
// 	std::cout << std::endl << "SOMA STORED " << temp9;
	writeImage<rawImageType_uint>(imageLabelMontage,temp9.c_str());
	
	std::string tempTABLEGLO = _outPathData+"/label_table.txt";
// 	std::cout << std::endl << "SOMA STORED " << temp9;
	ftk::SaveTable(tempTABLEGLO, tableLabelMontage);
	
	// Soma Montage
	rawImageType_uint::Pointer imageSomaMontage = rawImageType_uint::New();
	itk::Size<3> im_size = imageLabelMontage->GetBufferedRegion().GetSize();
	rawImageType_uint::RegionType region;
	region.SetSize( ImageMontageRegion.GetSize() );
	region.SetIndex( ImageMontageRegion.GetIndex() );
	imageSomaMontage->SetRegions( region );
	imageSomaMontage->Allocate();
	imageSomaMontage->FillBuffer(0);
	try
	{
		imageSomaMontage->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}
	rawImageType_uint::PixelType * imageSomaArray = imageSomaMontage->GetBufferPointer();
	rawImageType_uint::PixelType * imageLabelArray = imageLabelMontage->GetBufferPointer();

	unsigned long long sizeXY = ImageMontageSize[1] * ImageMontageSize[0];
	unsigned long long sizeX = ImageMontageSize[0];
	std::map<unsigned int, int> classMap;
	for(int row=0; row<(int)tableLabelMontage->GetNumberOfRows(); ++row)
	{
		classMap[tableLabelMontage->GetValue(row,0).ToUnsignedInt()] = tableLabelMontage->GetValueByName(row, "prediction_active_mg").ToInt();
	}

// 	#pragma omp parallel for collapse(3)
	for(int i=0; i<ImageMontageSize[2]; ++i)
	{
		for(int j=0; j<ImageMontageSize[1]; ++j)
		{
			for(int k=0; k<ImageMontageSize[0]; ++k)
			{
				unsigned long long offset = (i*sizeXY)+(j*sizeX)+k;
				if( imageLabelArray[offset] != 0 )
					if(classMap[imageLabelArray[offset]] == 1)
						imageSomaArray[offset] = imageLabelArray[offset];
			}
		}
	}	
	
	std::string nameSomaMontage = _outPathData+"/soma.nrrd";
	writeImage< rawImageType_uint >(imageSomaMontage,nameSomaMontage.c_str());
	
	std::string nameSomaMontageNew = _outPathData+"/soma_new.nrrd";
	writeImage< rawImageType_uint >(imageLabelMontageSomaNew,nameSomaMontageNew.c_str());
	
	IteratorType_uint iterimageSomaMontage = IteratorType_uint(imageSomaMontage,imageSomaMontage->GetRequestedRegion());
	iterimageSomaMontage.GoToBegin();
	IteratorType_uint iterMontageSomaNew = IteratorType_uint(imageLabelMontageSomaNew,imageLabelMontageSomaNew->GetRequestedRegion());
	iterMontageSomaNew.GoToBegin();

	int saameSoma = 1;
	for(;!iterimageSomaMontage.IsAtEnd();++iterimageSomaMontage, ++iterMontageSomaNew)
	{
		if( iterimageSomaMontage.Get() != iterMontageSomaNew.Get() )
			saameSoma = 0;
	}
	if( saameSoma ==0 )
		std::cout << std::endl << "NOS JODIMOS";
	else
		std::cout << std::endl << "Slimos bien, continuesmo";

	vtkSmartPointer<vtkTable> somaCentroidsTable = vtkSmartPointer<vtkTable>::New();
	somaCentroidsTable->Initialize();

	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "somaCen_x" );
	somaCentroidsTable->AddColumn(column);
	column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "somaCen_y" );
	somaCentroidsTable->AddColumn(column);
	column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName( "somaCen_z" );
	somaCentroidsTable->AddColumn(column);

	for(int row = 0; row < (int)tableLabelMontage->GetNumberOfRows(); ++row)
	{
		if(tableLabelMontage->GetValueByName(row, "prediction_active_mg").ToInt() == 1)
		{
			vtkSmartPointer<vtkVariantArray> centroid_data = vtkSmartPointer<vtkVariantArray>::New();
			centroid_data->InsertNextValue(tableLabelMontage->GetValue(row,1));
			centroid_data->InsertNextValue(tableLabelMontage->GetValue(row,2));
			centroid_data->InsertNextValue(tableLabelMontage->GetValue(row,3));
			somaCentroidsTable->InsertNextRow(centroid_data);			
		}
// 		std::cout << std::endl << tableLabelMontage->GetValueByName(row, "prediction_active_mg").ToInt()
	}

	for(int row = 0; row < (int)tableLabelMontage->GetNumberOfRows(); ++row)
	{
		if(tableLabelMontage->GetValueByName(row, "prediction_active_mg").ToInt() != 1)
		{
			tableLabelMontage->RemoveRow(row);
			--row;
		}
	}
	ftk::SaveTable(_outPathData + "/soma_table.txt", tableLabelMontage);
	ftk::SaveTable(_outPathData + "/soma_table_centroids.txt", somaCentroidsTable);
	
	
// 	vtkSmartPointer<vtkTable> somaCentroidsTableNew = vtkSmartPointer<vtkTable>::New();
// 	somaCentroidsTableNew->Initialize();
// 
// 	vtkSmartPointer<vtkDoubleArray> columnSoma = vtkSmartPointer<vtkDoubleArray>::New();
// 	columnSoma->SetName( "somaCen_x" );
// 	somaCentroidsTable->AddColumn(columnSoma);
// 	columnSoma = vtkSmartPointer<vtkDoubleArray>::New();
// 	columnSoma->SetName( "somaCen_y" );
// 	somaCentroidsTable->AddColumn(columnSoma);
// 	columnSoma = vtkSmartPointer<vtkDoubleArray>::New();
// 	columnSoma->SetName( "somaCen_z" );
// 	somaCentroidsTableNew->AddColumn(columnSoma);
// 
// 	for(int row = 0; row < (int)tableLabelMontageSomaNew->GetNumberOfRows(); ++row)
// 	{
// 		if(tableLabelMontageSomaNew->GetValueByName(row, "prediction_active_mg").ToInt() == 1)
// 		{
// 			vtkSmartPointer<vtkVariantArray> centroid_data = vtkSmartPointer<vtkVariantArray>::New();
// 			centroid_data->InsertNextValue(tableLabelMontageSomaNew->GetValue(row,1));
// 			centroid_data->InsertNextValue(tableLabelMontageSomaNew->GetValue(row,2));
// 			centroid_data->InsertNextValue(tableLabelMontageSomaNew->GetValue(row,3));
// 			somaCentroidsTableNew->InsertNextRow(centroid_data);			
// 		}
// // 		std::cout << std::endl << tableLabelMontage->GetValueByName(row, "prediction_active_mg").ToInt()
// 	}
// 	
// 	ftk::SaveTable(_outPathData + "/soma_table_new.txt", tableLabelMontageSomaNew);
// 	ftk::SaveTable(_outPathData + "/soma_table_centroids_new.txt", somaCentroidsTableNew);
	
	

}