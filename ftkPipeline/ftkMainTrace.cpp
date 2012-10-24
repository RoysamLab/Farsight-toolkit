
#include "ftkMainDarpaTrace.h"

void ftkMainDarpaTrace::readParameters( std::string segmentParams )
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
	
	iter = options.find("-xSize"); 
	if(iter!=options.end())
	{ std::istringstream ss((*iter).second); ss >> _xSize;}
	else
	{ _xSize = 1024; printf("Choose _xSize = 1024 as default\n");}
	
	iter = options.find("-ySize"); 
	if(iter!=options.end())
	{ std::istringstream ss((*iter).second); ss >> _ySize;}
	else
	{ _ySize = 1024; printf("Choose _ySize = 1024 as default\n");}
	
	iter = options.find("-zSize"); 
	if(iter!=options.end())
	{ std::istringstream ss((*iter).second); ss >> _zSize;}
	else
	{ _zSize = 1024; printf("Choose _zSize = 1024 as default\n");}
	
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
	
// 	iter = options.find("-xTileBor"); 
// 	if(iter!=options.end())
// 	{ std::istringstream ss((*iter).second); ss >> _xTileBor;}
// 	else
// 	{ _xTileBor = 1024; printf("Choose xTileBor = 1024 as default\n");}
// 	
// 	iter = options.find("-yTileBor"); 
// 	if(iter!=options.end())
// 	{ std::istringstream ss((*iter).second); ss >> _yTileBor;}
// 	else
// 	{ _yTileBor = 1024; printf("Choose yTileBor = 1024 as default\n");}
// 	
// 	iter = options.find("-zTileBor"); 
// 	if(iter!=options.end())
// 	{ std::istringstream ss((*iter).second); ss >> _zTileBor;}
// 	else
// 	{ _zTileBor = 1024; printf("Choose _zTileBor = 1024 as default\n");}
	
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
	
	iter = options.find("-traceParams"); 
	if(iter!=options.end())
	{ std::istringstream ss((*iter).second); ss >> _traceParams;}
	else
	{ _traceParams.clear(); printf("Choose _traceParams = NULL as default\n");}
	
// 	iter = options.find("-projectDefinition"); 
// 	if(iter!=options.end())
// 	{ std::istringstream ss((*iter).second); ss >> _projectDefinition;}
// 	else
// 	{ _projectDefinition.clear(); printf("Choose _projectDefinition = NULL as default\n");}
	
// 	iter = options.find("-optionsMNT"); 
// 	if(iter!=options.end())
// 	{ std::istringstream ss((*iter).second); ss >> _optionsMNT;}
// 	else
// 	{ _optionsMNT.clear(); printf("Choose optionsMNT = NULL as default\n");}
	
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
	
	
	// Print Parameters
	std::cout << std::endl << "This are the parameters";
	std::cout << std::endl << "_xSize: " << _xSize;
	std::cout << std::endl << "_ySize: " << _ySize;
	std::cout << std::endl << "_zSize: " << _zSize;
	std::cout << std::endl << "_xTile: " << _xTile;
	std::cout << std::endl << "_yTile: " << _yTile;
	std::cout << std::endl << "_zTile: " << _zTile;
// 	std::cout << std::endl << "_xTileBor: " << _xTileBor;
// 	std::cout << std::endl << "_yTileBor: " << _yTileBor;
// 	std::cout << std::endl << "_zTileBor: " << _zTileBor;
	std::cout << std::endl << "_num_threads: " << _num_threads;
	std::cout << std::endl << "_Cy5_Image: " << _Cy5_Image;
	std::cout << std::endl << "_TRI_Image: " << _TRI_Image;
	std::cout << std::endl << "_GFP_Image: " << _GFP_Image;
	std::cout << std::endl << "_DAP_Image: " << _DAP_Image;
	std::cout << std::endl << "_Soma_Centroids: " << _Soma_Centroids;
	std::cout << std::endl << "_Soma_Montage: " << _Soma_Montage;
	std::cout << std::endl << "_isSmall: " << _isSmall;
	std::cout << std::endl << "_traceParams: " << _traceParams;
// 	std::cout << std::endl << "_projectDefinition: " << _projectDefinition;
	std::cout << std::endl << "_outPath: " << _outPath;
	std::cout << std::endl << "_outPathDebug: " << _outPathDebug;
	std::cout << std::endl << "_outPathTemp: " << _outPathTemp;
}

void ftkMainDarpaTrace::runPreprocesing()
{
	Data_flo::Pointer MontageGFP_Image = readImage< Data_flo >(_GFP_ImageNRRD);
	
	RescaleFilterType_floToflo::Pointer rescaler = RescaleFilterType_floToflo::New();
	rescaler->SetOutputMinimum(0.0);
	rescaler->SetOutputMaximum(1.0);
	rescaler->SetInput( MontageGFP_Image );

	std::cout<<std::endl<<"Rescaler"<<std::flush;

	MedianFilterType_floToflo::Pointer medfilt = MedianFilterType_floToflo::New();
	medfilt->SetInput( rescaler->GetOutput() );
	Data_flo::SizeType rad = { {1, 1, 1} };
	medfilt->SetRadius(rad);
	medfilt->Update();
	
	std::cout<<std::endl<<"MedFilter"<<std::flush;
	
	_GFP_ImagePREPMNT = _outPathTemp+"/GFP_MNT_PRE.nrrd";
	std::cout<<_GFP_ImagePREPMNT<<std::endl;
	writeImage< Data_flo >( medfilt->GetOutput(), _GFP_ImagePREPMNT.c_str());
}


void ftkMainDarpaTrace::runTracing()
{
	if( _isSmall == 1 )
	{
		itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1); // This one can not be changed
		itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1); // This one can chenga
		omp_set_nested(0);
		omp_set_max_active_levels(1);
		int num_threads = 1;
		omp_set_num_threads(num_threads);
	}

	int counterCentro = 0;
	std::vector< itk::Index<3> > centroid_list = getCentroidList();
	std::cout << "Number of cells to be traced : " << centroid_list.size() << "\n";

	std::string SWCFilename = _outPath + "/OnlySWC.xml";
	std::ofstream outSWCFile;
	outSWCFile.open(SWCFilename.c_str());
	outSWCFile << "<?xml\tversion=\"1.0\"\t?>\n";
	outSWCFile << "<Source>\n\n";
		

		
		// CAN NOT BE RUN PARALLEL
		for( unsigned int bigTile = 0; bigTile<_numDivisionsInRowCEN ; ++bigTile )
		{
// 			std::cout<<std::endl<<"ITER";
			
// 			stringstream out34;
// 			out34<<bigTile;
// 			string srr = out34.str();
// 			std::string SWCFilenameDivided = _outPath + "/TracesAndSomasDivided/OnlySWC_" + srr +".xml";
// 			std::ofstream outfileDivided;
// 			outfileDivided.open(SWCFilenameDivided.c_str());
// 			outfileDivided << "<?xml\tversion=\"1.0\"\t?>\n";
// 			outfileDivided << "<Source>\n\n";

			itk::Index<3> initialBigIndexLOG;
			itk::Size<3> sizeOfTheRegionLOG;
			
			initialBigIndexLOG[0] = _initialBigTileLOG[bigTile][0];
			initialBigIndexLOG[1] = _initialBigTileLOG[bigTile][1];
			initialBigIndexLOG[2] = _initialBigTileLOG[bigTile][2];
			
			sizeOfTheRegionLOG[0] = _sizeOfBigTilesLOG[bigTile][0];
			sizeOfTheRegionLOG[1] = _sizeOfBigTilesLOG[bigTile][1]+_sizeOfBigTilesLOG[bigTile+1][1];
			sizeOfTheRegionLOG[2] = _sizeOfBigTilesLOG[bigTile][2];
			
			Data_flo::RegionType desiredRegionBigTileLOG;
			desiredRegionBigTileLOG.SetSize(sizeOfTheRegionLOG);
			desiredRegionBigTileLOG.SetIndex(initialBigIndexLOG);
			
			itk::Index<3> initialBigIndexCEN;
			itk::Size<3> sizeOfTheRegionCEN;
			
			initialBigIndexCEN[0] = initialBigTileCEN[bigTile][0];
			initialBigIndexCEN[1] = initialBigTileCEN[bigTile][1];
			initialBigIndexCEN[2] = initialBigTileCEN[bigTile][2];
			
			sizeOfTheRegionCEN[0] = sizeOfBigTilesCEN[bigTile][0];
			sizeOfTheRegionCEN[1] = sizeOfBigTilesCEN[bigTile][1];
			sizeOfTheRegionCEN[2] = sizeOfBigTilesCEN[bigTile][2];
			
// 			std::vector< Data_flo::Pointer > LoGDesiredRegion;
// 			LoGDesiredRegion.resize(6);
			
// 			rawImageType_uint::Pointer _somaMontageDesiredRegion;
// 			Data_flo::Pointer _img_traceDesiredRegion;
			
// 			std::string tempFileName_42 = _GFP_ImagePREPMNT
			_img_traceDesiredRegion = readImageRegion< Data_flo >( _GFP_ImagePREPMNT.c_str(), desiredRegionBigTileLOG );
			_somaMontageDesiredRegion = readImageRegion< rawImageType_uint >( _Soma_Montage.c_str(), desiredRegionBigTileLOG );

		
		#pragma omp parallel for num_threads(_num_threads) schedule(dynamic, 1)
			for( unsigned long long a=0; a<centroid_list.size(); ++a )
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
				
		#pragma omp critical
				{
					counterCentro++;
					std::cout<<std::endl<<"\t\t\t\t asdfasdf ----->>>>> " << counterCentro << ", of " << centroid_list.size();
					std::cout<<" "<<centroid_list[a][0]<<" "<<centroid_list[a][1]<<" "<<centroid_list[a][2];
				}

				ssx << x; ssy << y; ssz << z;

				std::cout << "Tracing Dice " << ssx.str() << "_" << ssy.str() << "_" << ssz.str() << "\n";

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
// 				if(_zSize <= _zTile)
// 					ssz_off << 0;
// 				else
// 				{	
// 					if(z >= _zTile/2)
// 						ssz_off << z - _zTile/2;
// 					else 
// 						ssz_off << 0;
// 				}
				
				int x_local = x;
				int y_local = y - _initialBigTileLOG[bigTile][1];
				int z_local = z;
				
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
// 				if(_zSize <= _zTile)
// 					ssz_offBig << 0;
// 				else
// 				{	
// 					if(z_local >= _zTile/2)
// 						ssz_offBig << z_local - _zTile/2;
// 					else 
// 						ssz_offBig << 0;
// 				}

				//########    CROP THE DESIRED DICE FROM THE GFP AND SOMA MONTAGES   ########

				
				std::cout<<std::cout<<"end of the criticals";
				
				//########    FETCH ALL CENTROIDS THAT FALL WITHIN THE DICE    ########

// 				if(_zSize <= _zTile)
// 					centroid[2] = z;
// 				else
// 					centroid[2] = ((z - _zTile/2)>0) ? _zTile/2:z;

				std::vector< itk::Index<3> > soma_Table = getSomaTable(centroid_list, start2, size2, x, y, z );
				
// 				//########    RUN TRACING    ########
				MultipleNeuronTracer * MNT = new MultipleNeuronTracer();
				MNT->LoadParameters(_traceParams,5);
// 				
// 				MNT->LoadCurvImage_1(img_trace, 0);
				MNT->LoadCurvImage_2(img_trace);
				MNT->ReadStartPoints_1(soma_Table, 0);
// 				MNT->SetCostThreshold(1000);
				MNT->SetCostThreshold(MNT->cost_threshold);
				MNT->LoadSomaImage_1(img_soma);
				
// 	// 			MNT->LoadSomaImage_1(img_soma_yan);
				bool flagLog = false;
				MNT->setFlagOutLog(flagLog);
				MNT->RunTracing();
// 
				x = min(_xTile/2, x);
				y = min(_yTile/2, y);
				z = min(_zTile/2, z);
//
				vtkSmartPointer< vtkTable > swcTable = MNT->GetSWCTable(0);
				std::string swcFilename = _outPath + "/TracesAndSomas/Trace_" + ssx.str() + "_" + ssy.str() + "_" + ssz.str() + "_ANT.swc";
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
	
	
	
	



void ftkMainDarpaTrace::computeSplitConst( rawImageType_8bit::Pointer ImageMontage )
{
		/////////////////////////////////////////////////////////////////////////////////////////////
		int numDivisionsInRowLOG = (int)floor((double)_ySize/(2*_yTile+10)); // Minimum value 3, otherwise does not make sense
		if( numDivisionsInRowLOG < 3 )
			numDivisionsInRowLOG = 2;
		_numDivisionsInRowCEN = numDivisionsInRowLOG - 1;
		
		_sizeOfBigTilesLOG.resize(numDivisionsInRowLOG);
		std::vector< itk::Size<3> > sizeOfBigTilesCEN;
		sizeOfBigTilesCEN.resize(numDivisionsInRowLOG-1);
		for( unsigned int jj=0; jj<numDivisionsInRowLOG; ++jj )
		{
			_sizeOfBigTilesLOG[jj][0] = _xSize;
			_sizeOfBigTilesLOG[jj][1] = (unsigned long long)floor((double)_ySize/(double)numDivisionsInRowLOG);
			_sizeOfBigTilesLOG[jj][2] = _zSize;
		}
		_sizeOfBigTilesLOG[numDivisionsInRowLOG-1][1] = _ySize - (numDivisionsInRowLOG-1)*_sizeOfBigTilesLOG[numDivisionsInRowLOG-1][1];
		if( _sizeOfBigTilesLOG[0][1] <2*_yTile)
		{
			std::cout<<std::endl<<"MISTAKE_1";
		}
		sizeOfBigTilesCEN[0][0] = _xSize;
		sizeOfBigTilesCEN[0][1] = _sizeOfBigTilesLOG[0][1]+(unsigned long long)floor((double)_sizeOfBigTilesLOG[0][1]/2);
		sizeOfBigTilesCEN[0][2] = _zSize;
		for( unsigned int jj=1; jj<_numDivisionsInRowCEN; ++jj )
		{
			sizeOfBigTilesCEN[jj][0] = _xSize;
			sizeOfBigTilesCEN[jj][1] = _sizeOfBigTilesLOG[jj][1];
			sizeOfBigTilesCEN[jj][2] = _zSize;
		}
		if( _numDivisionsInRowCEN == 1 )
		{
			sizeOfBigTilesCEN[_numDivisionsInRowCEN-1][1] = _ySize;
		}
		else
		{
			sizeOfBigTilesCEN[_numDivisionsInRowCEN-1][1] = _ySize -sizeOfBigTilesCEN[0][1] -  (_numDivisionsInRowCEN-2)*_sizeOfBigTilesLOG[numDivisionsInRowLOG-2][1];
		}
		
		
		unsigned long long total = 0;
		for( unsigned int jj=0; jj<numDivisionsInRowLOG; ++jj )
		{
			std::cout<<std::endl<<"Row Size Big"<<_sizeOfBigTilesLOG[jj][1];
			total = total + _sizeOfBigTilesLOG[jj][1];
		}
		if( total != _ySize )
		{
			std::cout<<std::endl<<"MISTAKE_2";
		}
		std::cout<<std::endl<<"\t\t"<<total<<" "<<_ySize<<" "<<_xSize;
		
		unsigned long long total2 = 0;
		for( unsigned int jj=0; jj<_numDivisionsInRowCEN; ++jj )
		{
			std::cout<<std::endl<<"Col Size Big"<<sizeOfBigTilesCEN[jj][1];
			total2 = total2 + sizeOfBigTilesCEN[jj][1];
		}
		if( total2 != _ySize )
		{
			std::cout<<std::endl<<"MISTAKE_3";
		}
		std::cout<<std::endl<<"\t\t"<<total2<<" "<<_ySize<<" "<<_xSize;
		
		_initialBigTileLOG.resize(numDivisionsInRowLOG);
		_initialBigTileLOG[0][0] = 0;
		_initialBigTileLOG[0][1] = 0;
		_initialBigTileLOG[0][2] = 0;
		for( unsigned int jj=1; jj<numDivisionsInRowLOG; ++jj )
		{
			_initialBigTileLOG[jj][0] = 0;
			_initialBigTileLOG[jj][1] = _sizeOfBigTilesLOG[jj-1][1] + _initialBigTileLOG[jj-1][1];
			_initialBigTileLOG[jj][2] = 0;
		}
		
		std::vector< itk::Index<3> > initialBigTileCEN;
		initialBigTileCEN.resize(_numDivisionsInRowCEN);
		initialBigTileCEN[0][0] = 0;
		initialBigTileCEN[0][1] = 0;
		initialBigTileCEN[0][2] = 0;
		for( unsigned int jj=1; jj<_numDivisionsInRowCEN; ++jj )
		{
			initialBigTileCEN[jj][0] = 0;
			initialBigTileCEN[jj][1] = sizeOfBigTilesCEN[jj-1][1] + initialBigTileCEN[jj-1][1];
			initialBigTileCEN[jj][2] = 0;
		}
	
}

std::vector< itk::Index<3> > ftkMainDarpaTrace::getCentroidList()
{
	std::vector< itk::Index<3> > centroid_list;
	
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
	return centroid_list;
}

std::vector< itk::Index<3> > ftkMainDarpaTrace::getSomaTable( std::vector< itk::Index<3> > centroid_list, itk::Index<3> start2, rawImageType_8bit::SizeType size2, int x, int y, int z )
{
	itk::Index<3> centroid;
	centroid[0] = ((x - _xTile/2)>0) ? _xTile/2:x; 
	centroid[1] = ((y - _yTile/2)>0) ? _yTile/2:y;
	centroid[2] = ((z - _zTile/2)>0) ? _zTile/2:z;
	
	std::vector< itk::Index<3> > soma_Table;      
	for(int ctr =0; ctr<centroid_list.size() ; ++ctr)
	{
		itk::Index<3> cen = centroid_list[ctr];
		if( (cen[0]>=start2[0]) && (cen[0]<(start2[0]+size2[0])) && (cen[1]>=start2[1]) && (cen[1]<(start2[1]+size2[1])) && (cen[2]>=start2[2]) && (cen[2]<(start2[2]+size2[2])) )
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

std::vector< itk::Index<3> > ftkMainDarpaTrace::cropImages( std::vector< rawImageType_uint::Pointer > img_soma, std::vector< Data_flo::Pointer > img_trace)
{
	rawImageType_uint::IndexType start;
	start[0] = ((x - _xTile/2)>0) ? (x - _xTile/2):0;
	start[1] = ((y - _yTile/2)>0) ? (y - _yTile/2):0;
	start[2] = ((z - _zTile/2)>0) ? (z - _zTile/2):0;
// 				if(_zSize <= _zTile)
// 					start[2] = 0;
// 				else
// 					start[2] = ((z - _zTile/2)>0) ? (z - _zTile/2):0;

	rawImageType_8bit::IndexType start2;
	start2[0] = ((x - _xTile/2)>0) ? (x - _xTile/2):0;
	start2[1] = ((y - _yTile/2)>0) ? (y - _yTile/2):0;
	start2[2] = ((z - _zTile/2)>0) ? (z - _zTile/2):0;
// 				if(_zSize <= _zTile)
// 					start2[2] = 0;
// 				else
// 					start2[2] = ((z - _zTile/2)>0) ? (z - _zTile/2):0;

	rawImageType_uint::SizeType size;
	size[0] = ((x+_xTile/2)<_xSize) ? _xTile : (_xTile/2+_xSize-x-1); 
	size[1] = ((y+_yTile/2)<_ySize) ? _yTile : (_yTile/2+_ySize-y-1);
	size[2] = ((z+_zTile/2)<_zSize) ? _zTile : (_zTile/2+_zSize-z-1);
// 				if(_zSize <= _zTile)
// 					size[2] = _zSize;
// 				else
// 					size[2] = ((z+_zTile/2)<_zSize) ? _zTile : (_zTile/2+_zSize-z-1);

	rawImageType_8bit::SizeType size2;
	size2[0] = ((x+_xTile/2)<_xSize) ? _xTile : (_xTile/2+_xSize-x-1);
	size2[1] = ((y+_yTile/2)<_ySize) ? _yTile : (_yTile/2+_ySize-y-1);
	size2[2] = ((z+_zTile/2)<_zSize) ? _zTile : (_zTile/2+_zSize-z-1);
// 				if(_zSize <= _zTile)
// 					size2[2] = _zSize;
// 				else
// 					size2[2] = ((z+_zTile/2)<_zSize) ? _zTile : (_zTile/2+_zSize-z-1);

	rawImageType_uint::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);

// 	rawImageType_8bit::RegionType desiredRegion2;
// 	desiredRegion2.SetSize(size2);
// 	desiredRegion2.SetIndex(start2);
// 	
// 	rawImageType_16bit::RegionType desiredRegion2_16bits;
// 	desiredRegion2_16bits.SetSize(size2);
// 	desiredRegion2_16bits.SetIndex(start2);
	
	Data_flo::RegionType desiredRegion2_float;
	desiredRegion2_float.SetSize(size2);
	desiredRegion2_float.SetIndex(start2);
	

	ROIFilterType_uint::Pointer ROIfilter3 = ROIFilterType_uint::New();
	ROIfilter3->SetRegionOfInterest(desiredRegion);
	ROIfilter3->SetInput(_somaMontageDesiredRegion);
#pragma omp critical
	ROIfilter3->Update();
	DuplicatorType_uint::Pointer LabelDuplicator = DuplicatorType_uint::New();
	LabelDuplicator->SetInputImage(ROIfilter3->GetOutput());
#pragma omp critical
	LabelDuplicator->Update();
	img_soma[0] = LabelDuplicator->GetOutput();
// 
	ROIFilterType_flo::Pointer ROIfilter2 = ROIFilterType_flo::New();
	ROIfilter2->SetRegionOfInterest(desiredRegion2_float);
	ROIfilter2->SetInput(_img_traceDesiredRegion);
#pragma omp critical
	ROIfilter2->Update();
	DuplicatorType_flo::Pointer floatDuplicator = DuplicatorType_flo::New();
	floatDuplicator->SetInputImage(ROIfilter2->GetOutput());
#pragma omp critical
	floatDuplicator->Update();

	img_trace[0] = floatDuplicator->GetOutput();
}
