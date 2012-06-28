// make -j80;./exe/darpa_tracer_w_seg ../../data/0103/DAPI_montage.mhd ../../data/0103/GFP_montage.mhd ../../data/0103/Cy5_montage_BS.mhd ../../data/0103/Seg_Params.ini ../../data/0103/ProjectDefinition.xml
// make -j80;((time ./exe/darpa_tracer_w_seg ../../data/0103/DAPI_montage.mhd ../../data/0103/GFP_montage.mhd ../../data/0103/Cy5_montage_BS.mhd ../../data/0103/Seg_Params.ini ../../data/0103/ProjectDefinition.xml) 2>&1) >> salida.log
// set args ../../data/0120/montage_kt01445_w410DAPIdsu_BS.nrrd ../../data/0120/montage_kt01445_w311GFPdsu_BS_CV.mhd ../../data/0120/montage_kt01445_w113Cy5dsu_BS.nrrd ../../data/0120/Parameters/Seg_Params.ini ../../data/0120/Parameters/ProjectDefinition.xml

// scp * nrey@farsight-05.ee.uh.edu:/data/nicolas/farsight_updated/bin/exe/


#include "itkImageFileReader.h"
#include "itkImage.h"
#include "itkIntTypes.h"
#include "boost/tokenizer.hpp"
#include <fstream>
#include "vul/vul_file.h"
#include "iostream"
#include "itkRegionOfInterestImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "../Tracing/MultipleNeuronTracer/MultipleNeuronTracer.h"
#include "../NuclearSegmentation/exe/SomaExtraction.h"
#include "../NuclearSegmentation/Nuclear_Association/VolumeOfInterest.h"

#include "vtkTable.h"
#include "ftkUtils.h"
#include "ftkImage.h"
#include "itkMultiThreader.h"
#include "itkRescaleIntensityImageFilter.h"

#include "../NuclearSegmentation/yousef_core/yousef_seg.h"
#include "../ftkFeatures/ftkLabelImageToFeatures.h"
#include "../NuclearSegmentation/NucleusEditor/ftkProjectProcessor.h"
#include <iostream>
#include <sstream>
#include <time.h>
#include <limits.h>
#include <math.h>
#include <algorithm>

#include <iostream>

#ifdef _OPENMP
#include "omp.h"
#endif

// typedef itk::Image<unsigned char,  3> nucImageType;
// typedef itk::Image<unsigned int, 3> LabelType;
// typedef MultipleNeuronTracer::ImageType3D gfpImageType;


template <typename T>
int writeImage(typename T::Pointer im, const char* filename)
{
	printf("Writing %s ... \n",filename);
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
	std::cout<<" done";
	return EXIT_SUCCESS;
}

template <typename T>
typename T::Pointer readImage(const char* filename)
{
	printf("Reading %s ... \n",filename);
	typedef typename itk::ImageFileReader<T> ReaderType;

	typename ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(filename);
	try
	{
		reader->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
// 		return EXIT_FAILURE;
	}
	std::cout<<" done";
	return reader->GetOutput();
}


template <typename T>
typename T::Pointer readImageRegion(const char* filename, typename T::RegionType region)
{
	printf("Reading %s ... \n",filename);
	typedef typename itk::ImageFileReader<T> ReaderType;

	typename ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(filename);
// 	try
// 	{
// 		reader->Update();
// 	}
// 	catch(itk::ExceptionObject &err)
// 	{
// 		std::cerr << "ExceptionObject caught!" <<std::endl;
// 		std::cerr << err << std::endl;
// // 		return EXIT_FAILURE;
// 	}
	typedef typename itk::ExtractImageFilter< T, T > ROIFilterType;
	typename ROIFilterType::Pointer ROIfilter = ROIFilterType::New();
	ROIfilter->SetExtractionRegion(region);
	ROIfilter->SetInput( reader->GetOutput() );
	ROIfilter->SetDirectionCollapseToIdentity();
	try
	{
		ROIfilter->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}
	typename T::Pointer out = ROIfilter->GetOutput();
	out->DisconnectPipeline();
	std::cout<<" done";
	return out;
}

template <typename TIN,typename TOUT>
typename TOUT::Pointer readerINcasterOUT( const char* name )
{
	typename TIN::Pointer inputImage = readImage<TIN>(name);
	typedef typename itk::CastImageFilter< TIN,TOUT > CasterFilterType;
	typename CasterFilterType::Pointer caster = CasterFilterType::New();
	caster->SetInput(inputImage);
	caster->Update();
	return caster->GetOutput();
}




typedef unsigned char inputPixelType;
typedef itk::Image<inputPixelType,  3> rawImageType;
typedef itk::Image<unsigned char,  3> rawImageType_8bit;
typedef itk::Image<unsigned short,  3> rawImageType_16bit;
typedef itk::Image<unsigned short, 3> LabelType;
typedef itk::Image<unsigned int, 3> montageLabelType;

typedef MultipleNeuronTracer::ImageType3D gfpImageType; // Float Image
typedef itk::ImageDuplicator< rawImageType >  rawDuplicatorType;
typedef itk::ImageDuplicator< rawImageType_16bit >  rawDuplicatorType_16;
typedef itk::ImageDuplicator< montageLabelType >  LabelDuplicatorType;
typedef itk::ImageDuplicator< gfpImageType >  gfpDuplicatorType;
typedef itk::RegionOfInterestImageFilter< rawImageType, rawImageType > rawROIFilterType;
typedef itk::RegionOfInterestImageFilter< rawImageType_16bit, rawImageType_16bit > rawROIFilterType_16;
typedef itk::RegionOfInterestImageFilter< gfpImageType, gfpImageType > gfpROIFilterType;
typedef itk::RegionOfInterestImageFilter< montageLabelType, montageLabelType > LabelROIFilterType;
typedef std::pair< LabelType::Pointer, vtkSmartPointer< vtkTable > > segResultsType;

typedef itk::ImageFileReader<rawImageType_16bit> rawReaderType;
typedef itk::ImageFileWriter<montageLabelType> labelWriterType;
typedef itk::ImageFileReader<montageLabelType> labelReaderType;
typedef itk::RescaleIntensityImageFilter< rawImageType_16bit, rawImageType_8bit > RescaleFilterType;
typedef itk::RescaleIntensityImageFilter< rawImageType_16bit, rawImageType_16bit > RescaleFilterType_16to16;
typedef itk::RescaleIntensityImageFilter< gfpImageType, gfpImageType > RescaleFilterType_floToflo;
typedef itk::CastImageFilter< rawImageType_16bit, gfpImageType > CasterFilterType_16Tofloat;
typedef itk::MedianImageFilter<gfpImageType, gfpImageType> MedianFilterType;
typedef itk::ImageFileWriter<rawImageType_8bit> rawWriterType_8;


inline LabelType::Pointer RunNuclearSegmentation(rawImageType::Pointer, const char*);
vtkSmartPointer<vtkTable> ComputeFeaturesAndAssociations(rawImageType::Pointer, rawImageType::Pointer, rawImageType::Pointer, LabelType::Pointer, const char* );
void RunEverything(rawImageType::RegionType, rawImageType_8bit::Pointer, rawImageType_8bit::Pointer, rawImageType_8bit::Pointer, std::vector< LabelType::Pointer> &, std::vector < vtkSmartPointer< vtkTable > >&, std::map< unsigned int, itk::Index<3> >*, const char*, const char*, int col);
std::map< unsigned int, itk::Index<3> > GetLabelToCentroidMap( vtkSmartPointer< vtkTable > );
void WriteCenterTrace(vtkSmartPointer< vtkTable >, int, int, int, std::string);
gfpImageType::Pointer preprocessingMNT( gfpImageType::Pointer );
// itk::SizeValueType readSizeFromFile(const char*);
// typedef itk::Size<3> SizeType2;
itk::Size<3>  readSizeFromFile(const char*);

rawImageType_8bit::Pointer readAndRescale_16to8(const char*);
rawImageType_8bit::Pointer readAndRescale_16to8(const char*, int debugCopy);
rawImageType_16bit::Pointer readAndRescale_16to16(const char* nameInput);

int main(int argc, char* argv[])
{
// 	std::cout<<std::endl<<"Num_threads: "<<itk::MultiThreader::GetGlobalDefaultNumberOfThreads();
// 	std::cout<<std::endl<<"Num_threads: "<<itk::MultiThreader::GetGlobalMaximumNumberOfThreads();
// 	
// 	itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);
// 	itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);
// 	
// 	std::cout<<std::endl<<"Num_threads: "<<itk::MultiThreader::GetGlobalDefaultNumberOfThreads();
// 	std::cout<<std::endl<<"Num_threads: "<<itk::MultiThreader::GetGlobalMaximumNumberOfThreads();
// 	
// 	std::cout<<std::endl<<"TEST_v2"<<std::flush;
// 	itk::MultiThreader::SetGlobalDefaultNumberOfThreads(80); // This one can not be changed
// 	itk::MultiThreader::SetGlobalMaximumNumberOfThreads(80); // This one can chenga
// 	
// 	std::cout<<std::endl<<"Num_threads: "<<itk::MultiThreader::GetGlobalDefaultNumberOfThreads();
// 	std::cout<<std::endl<<"Num_threads: "<<itk::MultiThreader::GetGlobalMaximumNumberOfThreads();
	
// 	bool flagSmall;
	int flagSmall = atoi(argv[7]);
	
// 	if( atoi()
// 	= false;
	
	if( flagSmall == 0 )
	{
		omp_set_nested(1);
		omp_set_max_active_levels(2);
	}
	else
	{
		omp_set_nested(1);
		omp_set_max_active_levels(5);
	}
	
	int num_threads = 1;
	omp_set_num_threads(num_threads);

	if(argc < 5)
	{
		std::cout<<"usage1: darpa_tracer <dapi_montage_file> <gfp_montage_file> <cy5_montage_file> <seg_params_file> <project_definition_file> tasktodo flagforsmallorbig\n";
		return 0;
	}

// 	std::string MyName = argv[0];

// 	std::string nucFileName(argv[1]);
// 	std::string filePath = ftk::GetFilePath(nucFileName);

// 	std::string temp = nucFileName;
// 	string::iterator it;
	
	// IMPROVE THIS TIPE OF GETTING THE EXTENTION 5 FOR NRRD
// 	it = temp.end() - 5;
// 	temp.erase(it, it+5);
	
	
	std::string gfpFileName(argv[2]);
	std::string filePath = ftk::GetFilePath(gfpFileName);

	std::string gfpFileNameNoExt = gfpFileName;
	string::iterator it2;
	
	// IMPROVE THIS TIPE OF GETTING THE EXTENTION 5 FOR MHD
	it2 = gfpFileNameNoExt.end() - 4;
	gfpFileNameNoExt.erase(it2, it2+4);
	
	std::string temp = gfpFileNameNoExt;
	

	vtkSmartPointer<vtkTable> somaCentroidsTable = NULL;
// 	montageLabelType::Pointer somaMontage;
// 	rawImageType_8bit::Pointer montage_gfp;

	int onlyTrace = atoi(argv[6]);
	
	if( onlyTrace == 0 )
	{
// 		#pragma omp 
		rawImageType_8bit::Pointer montage_nuc;
		rawImageType_8bit::Pointer montage_gfp;
		rawImageType_8bit::Pointer montage_cy5;
		
		itk::SizeValueType size_nuc_montage[3];
		itk::SizeValueType size_gfp_montage[3];
		itk::SizeValueType size_cy5_montage[3];
		
		#pragma omp parallel for num_threads(2)
		for (unsigned int i = 0; i < 2; ++i)
		{
			if( i == 0 )
			{
				montage_nuc = readAndRescale_16to8(argv[1]);
				size_nuc_montage[0] = montage_nuc->GetLargestPossibleRegion().GetSize()[0];
				size_nuc_montage[1] = montage_nuc->GetLargestPossibleRegion().GetSize()[1];
				size_nuc_montage[2] = montage_nuc->GetLargestPossibleRegion().GetSize()[2];
			}
			if( i == 1 )
			{
				montage_gfp = readAndRescale_16to8(argv[2]);
				size_gfp_montage[0] = montage_gfp->GetLargestPossibleRegion().GetSize()[0];
				size_gfp_montage[1] = montage_gfp->GetLargestPossibleRegion().GetSize()[1];
				size_gfp_montage[2] = montage_gfp->GetLargestPossibleRegion().GetSize()[2];
			}
// 			if( i == 2 )
// 			{
// 				montage_cy5 = readAndRescale_16to8(argv[3]);
// 				size_cy5_montage[0] = montage_cy5->GetLargestPossibleRegion().GetSize()[0];
// 				size_cy5_montage[1] = montage_cy5->GetLargestPossibleRegion().GetSize()[1];
// 				size_cy5_montage[2] = montage_cy5->GetLargestPossibleRegion().GetSize()[2];
// 			}
		}

		//#####################################################################################################################
		//	NUCLEAR SEGMENT THE MONTAGE TILE BY TILE AND STITCH THE RESULT TILES TOGETHER
		//#####################################################################################################################
		int numThreadsRows;
		int numThreadsRows_split;;
		int numThreadsCols;;
		int numThreadsCols_split;
		if( flagSmall == 0 )
		{
			numThreadsRows = 20;
			numThreadsRows_split = 10; //10
			numThreadsCols = 8;
			numThreadsCols_split = 8; //4
		}
		else
		{
			numThreadsRows = 1;
			numThreadsRows_split = 1;
			numThreadsCols = 1;
			numThreadsCols_split = 1;
		}
		unsigned long long rowDivisor = ceil((double)size_nuc_montage[1]/numThreadsRows);//400;//861;
		unsigned long long colDivisor = ceil((double)size_nuc_montage[0]/numThreadsCols);//400;//640;
		unsigned long long num_rows = (unsigned long long)ceil((double)size_nuc_montage[1]/(double)rowDivisor);
		unsigned long long num_cols = (unsigned long long)ceil((double)size_nuc_montage[0]/(double)colDivisor);
		std::cout << "Row: " << size_nuc_montage[1] << ", Col: " << size_nuc_montage[0]<<", Stack: " <<size_nuc_montage[2];
		std::cout << "Image divided into " << num_rows << " rows and " << num_cols << " columns\n"<<std::flush;

		//##################	INITIALIZING VARIABLES FOR ROWS AND ROW_BORDERS	  ###################

		std::vector< LabelType::Pointer > Label_Rows;
		Label_Rows.resize(num_rows);
		std::vector< LabelType::Pointer > Label_RowBorders;
		Label_RowBorders.resize(num_rows-1);
		std::vector< vtkSmartPointer< vtkTable > > Table_Rows;
		Table_Rows.resize(num_rows);
		std::vector< vtkSmartPointer< vtkTable > > Table_RowBorders;
		Table_RowBorders.resize(num_rows-1);
		std::vector< std::map< unsigned int, itk::Index<3> > > Centroids_Rows;
		Centroids_Rows.resize(num_rows);
		std::vector< std::map< unsigned int, itk::Index<3> > > Centroids_RowBorders;
		Centroids_RowBorders.resize(num_rows-1);
		
		
		
		if( flagSmall == 0 )
		{
			itk::MultiThreader::SetGlobalDefaultNumberOfThreads(2); // This one can not be changed
			itk::MultiThreader::SetGlobalMaximumNumberOfThreads(2); // This one can chenga
		}
		else
		{
			itk::MultiThreader::SetGlobalDefaultNumberOfThreads(80); // This one can not be changed
			itk::MultiThreader::SetGlobalMaximumNumberOfThreads(80); // This one can chenga
		}
		//##################	SEGMENTING EACH ROW IN THE MONTAGE	  ###################
#pragma omp parallel for num_threads(numThreadsRows_split) schedule(dynamic, 1)
		for(int row=0; row<num_rows; ++row)
		{
#pragma omp critical
			{	
				std::cout<<std::endl<<"\t Row " << row<<std::flush;
			}

			//##################	INITIALIZING VARIABLES FOR TILES AND TILE_BORDERS	  ###################

			std::vector< LabelType::Pointer > Label_Tiles;
			Label_Tiles.resize(num_cols);
			std::vector< LabelType::Pointer > Label_TileBorders;
			Label_TileBorders.resize(num_cols-1);
			std::vector< vtkSmartPointer< vtkTable > > Table_Tiles;
			Table_Tiles.resize(num_cols);
			std::vector< vtkSmartPointer< vtkTable > > Table_TileBorders;
			Table_TileBorders.resize(num_cols-1);
			std::vector< std::map< unsigned int, itk::Index<3> > > Centroids_Tiles;
			Centroids_Tiles.resize(num_cols);
			std::vector< std::map< unsigned int, itk::Index<3> > > Centroids_TileBorders;
			Centroids_TileBorders.resize(num_cols-1);

			//##################	SEGMENTING EACH TILE IN A ROW    ###################
#pragma omp parallel for num_threads(numThreadsCols_split) schedule(dynamic, 1)
			for(unsigned int col=0; col<num_cols; ++col)
			{
				rawImageType::IndexType start_tile;
				start_tile[0] = col*colDivisor;
				start_tile[1] = row*rowDivisor;
				start_tile[2] = 0;

				rawImageType::SizeType size_tile;
				size_tile[0] = (((col+1)*colDivisor) < size_nuc_montage[0]) ? colDivisor : (size_nuc_montage[0] - (col*colDivisor));
				size_tile[1] = (((row+1)*rowDivisor) < size_nuc_montage[1]) ? rowDivisor : (size_nuc_montage[1] - (row*rowDivisor));
				size_tile[2] = size_nuc_montage[2];

				rawImageType::RegionType region_tile;
				region_tile.SetSize(size_tile);
				region_tile.SetIndex(start_tile);
				
				RunEverything(region_tile, montage_nuc, montage_gfp, montage_cy5, Label_Tiles, Table_Tiles, &(Centroids_Tiles[col]), argv[4], argv[5], col);
				
				//##################	EXTRACT A TILE_BORDER AND START SEGMENTING	  ###################

				if(col != 0)
				{
					// 			std::cout<<"Extracting tile border X = " << col*rowDivisor << ", Row = " << row << "\n";
					rawImageType::IndexType start_tileBorder;
					start_tileBorder[0] = (col*colDivisor) - (2*25);
					start_tileBorder[1] = (row*rowDivisor);
					start_tileBorder[2] = 0;

					rawImageType::SizeType size_tileBorder;
					size_tileBorder[0] = (4*25);
					size_tileBorder[1] = (((row+1)*rowDivisor) < size_nuc_montage[1]) ? rowDivisor : (size_nuc_montage[1] - (row*rowDivisor));
					size_tileBorder[2] = size_nuc_montage[2];

					rawImageType::RegionType region_tileBorder;
					region_tileBorder.SetSize(size_tileBorder);
					region_tileBorder.SetIndex(start_tileBorder);

					RunEverything(region_tileBorder, montage_nuc, montage_gfp, montage_cy5, Label_TileBorders, Table_TileBorders, &(Centroids_TileBorders[col-1]), argv[4], argv[5],col-1);
				}
				std::cout<<std::endl<<"\t\t\t\t ---->>>> Done with nucleus segmentation for Tile " << col << "_" << row << "\n";
			}

			std::cout<<"Stitching all tiles in Row " << row << "...";

			//##################	ALLOCATING MEMORY AND REGISTERING FOR A SEGMENTED ROW_BORDER	  ###################

			LabelType::Pointer rowSeg = LabelType::New();
			LabelType::PointType origin_row;
			origin_row[0] = 0; 
			origin_row[1] = 0;
			origin_row[2] = 0;
			rowSeg->SetOrigin( origin_row );
			LabelType::IndexType start_row;
			start_row[0] = 0;
			start_row[1] = 0;
			start_row[2] = 0;
			LabelType::SizeType size_row;
			size_row[0] = size_nuc_montage[0];
			size_row[1] = (((row+1)*rowDivisor) < size_nuc_montage[1]) ? rowDivisor : (size_nuc_montage[1] - (row*rowDivisor));
			size_row[2] = size_nuc_montage[2];
			LabelType::RegionType region_row;
			region_row.SetSize ( size_row  );
			region_row.SetIndex( start_row );
			rowSeg->SetRegions( region_row );
			rowSeg->Allocate();
			rowSeg->FillBuffer(0);
			rowSeg->Update();

			vtkSmartPointer< vtkTable > rowTable = vtkSmartPointer<vtkTable>::New();
			rowTable->Initialize();
			for(int c=0; c<(int)Table_Tiles[0]->GetNumberOfColumns(); ++c)
			{
				vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
				column->SetName( Table_Tiles[0]->GetColumnName(c) );
				rowTable->AddColumn(column);
			}
			//#############################################################

			LabelType::PixelType * rowSegArray = rowSeg->GetBufferPointer();
			unsigned long long row_slice_size = size_row[1] * size_row[0];
			unsigned long long row_row_size = size_row[0];		

			unsigned short max_value = 0, current_max = 0;
			for(int m=0; m<(int)Label_Tiles.size(); ++m)
			{
				std::cout<<std::endl<<m<<" "<<(int)Label_Tiles.size()<<" asdfasdf2";
				if(m != 0)
				{
					//##################	STITCHING THE TILE_BORDER INTO THE ROW	  ###################

					LabelType::Pointer myTileBorder = Label_TileBorders[m-1];
					LabelType::PixelType * myTileBorderArray = myTileBorder->GetBufferPointer();
					itk::Size<3> tileBorder_size = myTileBorder->GetLargestPossibleRegion().GetSize();
					unsigned long long tileBorder_slice_size = tileBorder_size[1] * tileBorder_size[0];
					unsigned long long tileBorder_row_size = tileBorder_size[0];	
					unsigned long long x_offset_tb = (m*colDivisor) - (2*25);
					for(unsigned long long z=0; z<tileBorder_size[2]; ++z)
					{
						for(unsigned long long y=0; y<tileBorder_size[1]; ++y)
						{
							for(unsigned long long x=0; x<tileBorder_size[0]; ++x)
							{
								unsigned short value = myTileBorderArray[(tileBorder_slice_size*z) + (tileBorder_row_size*y) + (x)];
								if(value == 0) 
									continue;
								unsigned long long lab_cen_x = Centroids_TileBorders[m-1][value][0];
								if((lab_cen_x < 25) || (lab_cen_x >= (3*25))) 
									continue;
								rowSegArray[(row_slice_size*z) + (row_row_size*y) + (x_offset_tb + x)] = max_value + value;
								//if((max_value + value) > current_max)
								//	current_max = max_value + value;
							}
						}
					}

					//##################	STITCHING THE TILE_BORDER_TABLE INTO THE ROW_TABLE	  ###################
					if((int)Table_TileBorders[m-1]->GetNumberOfRows() != 0)
					{
						for(int r=0; r<(int)Table_TileBorders[m-1]->GetNumberOfRows(); ++r)
						{
							if((Table_TileBorders[m-1]->GetValue(r,1).ToInt() < 25) || (Table_TileBorders[m-1]->GetValue(r,1).ToInt() >= (3*25)))
								continue;
							vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
							for(int c=0; c<(int)Table_TileBorders[m-1]->GetNumberOfColumns(); ++c)
							{
								if(c == 0)
									model_data1->InsertNextValue(vtkVariant(Table_TileBorders[m-1]->GetValue(r,c).ToUnsignedShort() + max_value));
								else if(c == 1)
									model_data1->InsertNextValue(vtkVariant(Table_TileBorders[m-1]->GetValue(r,c).ToInt() + x_offset_tb));
								else
									model_data1->InsertNextValue(Table_TileBorders[m-1]->GetValue(r,c));
							}
							rowTable->InsertNextRow(model_data1);
						}
						//max_value = current_max;
						max_value = rowTable->GetValue((int)rowTable->GetNumberOfRows()-1, 0).ToUnsignedShort();
					}
				}

				//##################	STITCHING THE TILE INTO THE ROW	  ###################

				LabelType::Pointer myTile = Label_Tiles[m];
				LabelType::PixelType * myTileArray = myTile->GetBufferPointer();
				itk::Size<3> tile_size = myTile->GetLargestPossibleRegion().GetSize();
				unsigned long long tile_slice_size = tile_size[1] * tile_size[0];
				unsigned long long tile_row_size = tile_size[0];	
				unsigned long long x_offset_t = m*colDivisor;
				for(unsigned long long z=0; z<tile_size[2]; ++z)
				{
					for(unsigned long long y=0; y<tile_size[1]; ++y)
					{
						for(unsigned long long x=0; x<tile_size[0]; ++x)
						{
							unsigned short value = myTileArray[(tile_slice_size*z) + (tile_row_size*y) + x];
							if(value == 0) continue;
							unsigned long long lab_cen_x = Centroids_Tiles[m][value][0];
							if((m != 0) && (lab_cen_x < 25)) continue;
							if((m != (Label_Tiles.size()-1)) && (lab_cen_x >= (tile_size[0]-25))) continue;
							rowSegArray[(row_slice_size*z) + (row_row_size*y) + (x_offset_t + x)] = max_value + value;
							//if((max_value + value) > current_max)
							//	current_max = max_value + value;
						}
					}
				}

				//##################	STITCHING THE TILE_TABLE INTO THE ROW_TABLE	  ###################
				if((unsigned long long)Table_Tiles[m]->GetNumberOfRows() != 0)
				{
					for(unsigned long long r=0; r<(unsigned long long)Table_Tiles[m]->GetNumberOfRows(); ++r)
					{
						if((m != 0) && (Table_Tiles[m]->GetValue(r,1).ToInt() < 25)) continue;
						if((m != (Label_Tiles.size()-1)) && (Table_Tiles[m]->GetValue(r,1).ToInt() >= (tile_size[0]-25))) continue;
						vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
						for(unsigned long long c=0; c<(unsigned long long)Table_Tiles[m]->GetNumberOfColumns(); ++c)
						{
							if(c == 0)
								model_data1->InsertNextValue(vtkVariant(Table_Tiles[m]->GetValue(r,c).ToUnsignedShort() + max_value));
							else if(c == 1)
								model_data1->InsertNextValue(vtkVariant(Table_Tiles[m]->GetValue(r,c).ToInt() + x_offset_t));
							else
								model_data1->InsertNextValue(Table_Tiles[m]->GetValue(r,c));
						}
						rowTable->InsertNextRow(model_data1);
					}
					//max_value = current_max;
					std::cout<<std::endl<<"MAX VAL ROW: "<<(int)rowTable->GetNumberOfRows()-1;
					max_value = rowTable->GetValue((int)rowTable->GetNumberOfRows()-1, 0).ToUnsignedShort();
				}
			}

			std::cout<<"Max Value in Row " << row << "..." << current_max;
			std::cout<<std::endl<<"asdfasdf1 ";

// 			Table_Tiles.clear();
// 			Label_TileBorders.clear();

			//#############################################################
			Label_Rows[row] = rowSeg;
			Table_Rows[row] = rowTable;
			Centroids_Rows[row] = GetLabelToCentroidMap(rowTable);

			//##################	EXTRACT A ROW_BORDER AND START SEGMENTING	  ###################

			if(row != 0)
			{
				rawImageType::IndexType start_rowBorder;
				start_rowBorder[0] = 0;
				start_rowBorder[1] = (row*rowDivisor) - (2*25);
				start_rowBorder[2] = 0;

				rawImageType::SizeType size_rowBorder;
				size_rowBorder[0] = size_nuc_montage[0];
				size_rowBorder[1] = (4*25);
				size_rowBorder[2] = size_nuc_montage[2];

				rawImageType::RegionType region_rowBorder;
				region_rowBorder.SetSize(size_rowBorder);
				region_rowBorder.SetIndex(start_rowBorder);

				RunEverything(region_rowBorder, montage_nuc, montage_gfp, montage_cy5, Label_RowBorders, (Table_RowBorders), &(Centroids_RowBorders[row-1]), argv[4], argv[5],row-1);
			}
			std::cout<<std::endl<<"\t\t Row "<<row<<" done with all"<<std::flush;
		}

		num_threads = 80;
		omp_set_num_threads(num_threads);
		
		std::cout<<"Stitching all rows ..."<<std::flush;

		//##################	ALLOCATING MEMORY FOR THE SEGMENTED MONTAGE	  ###################

		montageLabelType::Pointer montageSeg = montageLabelType::New();
		montageLabelType::PointType origin_montage;
		origin_montage[0] = 0; 
		origin_montage[1] = 0;
		origin_montage[2] = 0;
		montageSeg->SetOrigin( origin_montage );
		montageLabelType::IndexType start_montage;
		start_montage[0] = 0;
		start_montage[1] = 0;
		start_montage[2] = 0;
		montageLabelType::SizeType size_montage;
		size_montage[0] = size_nuc_montage[0];
		size_montage[1] = size_nuc_montage[1];
		size_montage[2] = size_nuc_montage[2];
		montageLabelType::RegionType region_montage;
		region_montage.SetSize ( size_montage  );
		region_montage.SetIndex( start_montage );
		montageSeg->SetRegions( region_montage );
		montageSeg->Allocate();
		montageSeg->FillBuffer(0);
		montageSeg->Update();

		vtkSmartPointer< vtkTable > montageTable = vtkSmartPointer<vtkTable>::New();
		montageTable->Initialize();
		for(unsigned long long c=0; c<(unsigned long long)Table_Rows[0]->GetNumberOfColumns(); ++c)
		{
			vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
			column->SetName( Table_Rows[0]->GetColumnName(c) );
			montageTable->AddColumn(column);
		}

		//#############################################################

		montageLabelType::PixelType * montageSegArray = montageSeg->GetBufferPointer();
		unsigned long long montage_slice_size = size_montage[1] * size_montage[0];
		unsigned long long montage_row_size = size_montage[0];	

		unsigned int max_value_1 = 0, current_max_1 = 0;
		for(unsigned long long n=0; n<(int)Label_Rows.size(); ++n)
		{
			if(n != 0)
			{

				//##################	STITCHING THE ROW_BORDER INTO THE MONTAGE	  ###################

				// 			std::cout << "stitchin the row border\n" ;
				LabelType::Pointer myRowBorder = Label_RowBorders[n-1];
				LabelType::PixelType * myRowBorderArray = myRowBorder->GetBufferPointer();
				itk::Size<3> rowBorder_size = myRowBorder->GetLargestPossibleRegion().GetSize();
				unsigned long long rowBorder_slice_size = rowBorder_size[1] * rowBorder_size[0];
				unsigned long long rowBorder_row_size = rowBorder_size[0];	
				unsigned long long y_offset_rb = (n*rowDivisor) - (2*25);
				for(unsigned long long z=0; z<rowBorder_size[2]; ++z)
				{
					for(unsigned long long y=0; y<rowBorder_size[1]; ++y)
					{
						for(unsigned long long x=0; x<rowBorder_size[0]; ++x)
						{
							unsigned int value_1 = myRowBorderArray[(rowBorder_slice_size*z) + (rowBorder_row_size*y) + (x)];
							if(value_1 == 0) continue;
							unsigned long long lab_cen_y = Centroids_RowBorders[n-1][value_1][1];
							if((lab_cen_y < 25) || (lab_cen_y >= (3*25))) continue;
							montageSegArray[(montage_slice_size*z) + (montage_row_size*(y_offset_rb + y)) + x] = max_value_1 + value_1;
							//if((max_value_1 + value_1) > current_max_1)
							//	current_max_1 = max_value_1 + value_1;
						}
					}
				}
				//Label_RowBorders[n-1]->UnRegister();

				//##################	STITCHING THE ROW_BORDER_TABLE INTO THE MONTAGE_TABLE	  ###################

				// 			std::cout << "stitchin the row border table\n" ;
				if((unsigned long long)Table_RowBorders[n-1]->GetNumberOfRows() != 0)
				{
					for(unsigned long long r=0; r<(unsigned long long)Table_RowBorders[n-1]->GetNumberOfRows(); ++r)
					{
						if((Table_RowBorders[n-1]->GetValue(r,2).ToInt() < 25) || (Table_RowBorders[n-1]->GetValue(r,2).ToInt() >= (3*25)))
							continue;
						vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
						for(unsigned long long c=0; c<(unsigned long long)Table_RowBorders[n-1]->GetNumberOfColumns(); ++c)
						{
							if(c == 0)
								model_data1->InsertNextValue(vtkVariant(Table_RowBorders[n-1]->GetValue(r,c).ToUnsignedInt() + max_value_1));
							else if(c == 2)
								model_data1->InsertNextValue(vtkVariant(Table_RowBorders[n-1]->GetValue(r,c).ToInt() + y_offset_rb));
							else
								model_data1->InsertNextValue(Table_RowBorders[n-1]->GetValue(r,c));
						}
						montageTable->InsertNextRow(model_data1);
					}
					max_value_1 = montageTable->GetValue((int)montageTable->GetNumberOfRows()-1, 0).ToUnsignedInt();
				}
			}

			//##################	STITCHING THE ROW INTO THE MONTAGE	  ###################

			LabelType::Pointer myRow = Label_Rows[n];
			LabelType::PixelType * myRowArray = myRow->GetBufferPointer();
			itk::Size<3> row_size = myRow->GetLargestPossibleRegion().GetSize();
			unsigned long long row_slice_size_1 = row_size[1] * row_size[0];
			unsigned long long row_row_size_1 = row_size[0];	
			unsigned long long y_offset_r = n*rowDivisor;
			for(unsigned long long z=0; z<row_size[2]; ++z)
			{
				for(unsigned long long y=0; y<row_size[1]; ++y)
				{
					for(unsigned long long x=0; x<row_size[0]; ++x)
					{
						unsigned int value_1 = myRowArray[(row_slice_size_1*z) + (row_row_size_1*y) + (x)];					
						if(value_1 == 0) continue;
						unsigned long long lab_cen_y = Centroids_Rows[n][value_1][1];
						if((n != 0) && (lab_cen_y < 25)) continue;
						if((n != (Label_Rows.size()-1)) && (lab_cen_y >= (row_size[0]-25))) continue;
						montageSegArray[(montage_slice_size*z) + (montage_row_size*(y_offset_r + y)) + x] = max_value_1 + value_1;
						//if((max_value_1 + value_1) > current_max_1)
						//	current_max_1 = max_value_1 + value_1;
					}
				}
			}

			//##################	STITCHING THE ROW_TABLE INTO THE MONTAGE_TABLE	  ###################

			if((unsigned long long)Table_Rows[n]->GetNumberOfRows() != 0)
			{
				for(unsigned long long r=0; r<(unsigned long long)Table_Rows[n]->GetNumberOfRows(); ++r)
				{
					if((n != 0) && (Table_Rows[n]->GetValue(r,2).ToInt() < 25)) continue;
					if((n != (Label_Rows.size()-1)) && (Table_Rows[n]->GetValue(r,2).ToInt() >= (row_size[0]-25))) continue;
					vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
					for(unsigned long long c=0; c<(unsigned long long)Table_Rows[n]->GetNumberOfColumns(); ++c)
					{
						if(c == 0)
							model_data1->InsertNextValue(vtkVariant(Table_Rows[n]->GetValue(r,c).ToUnsignedInt() + max_value_1));
						else if(c == 2)
							model_data1->InsertNextValue(vtkVariant(Table_Rows[n]->GetValue(r,c).ToInt() + y_offset_r));
						else
							model_data1->InsertNextValue(Table_Rows[n]->GetValue(r,c));
					}
					montageTable->InsertNextRow(model_data1);
				}
				max_value_1 = montageTable->GetValue((int)montageTable->GetNumberOfRows()-1, 0).ToUnsignedInt();
			}
		}

		std::cout << "everything done\n";
		// Not sure this is working correctly
		Label_Rows.clear();
		Label_RowBorders.clear();

		//#############################################################
		std::cout<<"done !!\n\n";

		labelWriterType::Pointer writer = labelWriterType::New();
		writer->SetFileName(temp + "_label.mhd");
		writer->SetInput(montageSeg);
		writer->Update();

		ftk::SaveTable(temp + "_table.txt", montageTable);

		//#####################################################################################################################
		//	EXTRACTING SOMA and SOMA CENTROIDS
		//#####################################################################################################################

		montageLabelType::Pointer somaMontage = montageLabelType::New();
		itk::Size<3> im_size = montageSeg->GetBufferedRegion().GetSize();
		montageLabelType::IndexType start;
		start[0] =   0;  // first index on X
		start[1] =   0;  // first index on Y    
		start[2] =   0;  // first index on Z  
		montageLabelType::PointType origin;
		origin[0] = 0; 
		origin[1] = 0;    
		origin[2] = 0;    
		somaMontage->SetOrigin( origin );
		montageLabelType::RegionType region;
		region.SetSize( im_size );
		region.SetIndex( start );
		somaMontage->SetRegions( region );
		somaMontage->Allocate();
		somaMontage->FillBuffer(0);
		somaMontage->Update();
		montageLabelType::PixelType * somaArray = somaMontage->GetBufferPointer();
		montageLabelType::PixelType * labelArray = montageSeg->GetBufferPointer();

		unsigned long long slice_size = im_size[1] * im_size[0];
		unsigned long long row_size = im_size[0];
		std::map<unsigned int, int> classMap;
		for(int row=0; row<(int)montageTable->GetNumberOfRows(); ++row)
		{
			classMap[montageTable->GetValue(row,0).ToUnsignedInt()] = montageTable->GetValueByName(row, "prediction_active_mg").ToInt();
		}

		for(int i=0; i<im_size[2]; ++i)
		{
			for(int j=0; j<im_size[1]; ++j)
			{
				for(int k=0; k<im_size[0]; ++k)
				{
					unsigned long long offset = (i*slice_size)+(j*row_size)+k;
					if(classMap[labelArray[offset]] == 1)
						somaArray[offset] = labelArray[offset];
				}
			}
		}	

		writer = labelWriterType::New();
		writer->SetFileName(temp + "_soma_montage.mhd");
		writer->SetInput(somaMontage);
		writer->Update();

		somaCentroidsTable = vtkSmartPointer<vtkTable>::New();
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

		for(int row = 0; row < (int)montageTable->GetNumberOfRows(); ++row)
		{
			if(montageTable->GetValueByName(row, "prediction_active_mg").ToInt() == 1)
			{
				vtkSmartPointer<vtkVariantArray> centroid_data = vtkSmartPointer<vtkVariantArray>::New();
				centroid_data->InsertNextValue(montageTable->GetValue(row,1));
				centroid_data->InsertNextValue(montageTable->GetValue(row,2));
				centroid_data->InsertNextValue(montageTable->GetValue(row,3));
				somaCentroidsTable->InsertNextRow(centroid_data);			
			}
		}

		for(int row = 0; row < (int)montageTable->GetNumberOfRows(); ++row)
		{
			if(montageTable->GetValueByName(row, "prediction_active_mg").ToInt() != 1)
			{
				montageTable->RemoveRow(row);
				--row;
			}
		}
		ftk::SaveTable(temp + "_soma_table.txt", montageTable);
		ftk::SaveTable(temp + "_soma_centroids_table.txt", somaCentroidsTable);
	}

	if( onlyTrace == 3 )
	{
		std::cout<<std::endl<<"OnlyTrace 3";
		somaCentroidsTable = ftk::LoadTable(temp + "_soma_centroids_table.txt");
		std::cout<<std::endl<<temp + "_soma_centroids_table.txt";
		std::string object_file = argv[8];
		
		rawImageType_8bit::Pointer binaryImage = readImage< rawImageType_8bit >(argv[8]);
		
		typedef itk::SignedMaurerDistanceMapImageFilter< rawImageType_8bit, gfpImageType > DanielssonFilterType;
		DanielssonFilterType::Pointer danielssonFilter = DanielssonFilterType::New();	
		danielssonFilter->SetInput( binaryImage );
// 		danielssonFilter->InputIsBinaryOn();
		std::cout<<std::endl<<"Start Danielson"<<std::flush;
		danielssonFilter->Update();
		gfpImageType::Pointer distMap = danielssonFilter->GetOutput();
		
		
// 		std::string nic = temp + "_Distance.mhd";
// 		gfpImageType::Pointer binaryImage2 = readImage< gfpImageType >(nic.c_str());
// 		gfpImageType::Pointer distMap = binaryImage2;
		
		

		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( "Dist_to_Object" );
		column->SetNumberOfValues(somaCentroidsTable->GetNumberOfRows());
		somaCentroidsTable->AddColumn(column);

		for(int row=0; row<(int)somaCentroidsTable->GetNumberOfRows(); ++row)
		{
			std::cout<<std::endl<<row;
			gfpImageType::IndexType indx;
			indx[0] = somaCentroidsTable->GetValue(row, 0).ToInt();
			indx[1] = somaCentroidsTable->GetValue(row, 1).ToInt();
			indx[2] = somaCentroidsTable->GetValue(row, 2).ToInt();
			std::cout<<" "<<indx[0]<<" "<<indx[1]<<" "<<indx[2]<<std::flush;
			std::cout<<" "<<indx[0]<<" "<<indx[1]<<" "<<indx[2]<<std::flush;
			float dist = distMap->GetPixel(indx);
// 			if( dist < std::numeric_limits<float>::max() )
// 			{
// 			float dist = 0;
				somaCentroidsTable->SetValueByName(row, "Dist_to_Object", vtkVariant(dist));
// 			}
			std::cout<<" "<<dist<<std::flush;
		}

		ftk::SaveTable(temp + "_soma_centroids_table_out.txt", somaCentroidsTable);
	}

	if( onlyTrace == 1 )
	{	
		std::cout<<std::endl<<"TRACE_PREPROCES";
// 		itk::MultiThreader::SetGlobalDefaultNumberOfThreads(80); // This one can not be changed
// 		itk::MultiThreader::SetGlobalMaximumNumberOfThreads(80); // This one can chenga
		
// 		rawImageType_16bit::Pointer montage_gfp = readAndRescale_16to16(argv[2]);
		
		gfpImageType::Pointer montage_gfp = readImage< gfpImageType >(argv[2]);
		
		itk::SizeValueType size_gfp_montage[3];
		size_gfp_montage[0] = montage_gfp->GetLargestPossibleRegion().GetSize()[0];
		size_gfp_montage[1] = montage_gfp->GetLargestPossibleRegion().GetSize()[1];
		size_gfp_montage[2] = montage_gfp->GetLargestPossibleRegion().GetSize()[2];

		//#####################################################################################################################
		//	NEW MULTIPLE NEURON TACER
		//#####################################################################################################################
		
		// Do the same preprocessing as in MNT
		gfpImageType::Pointer img_trace = preprocessingMNT(montage_gfp);

		typedef itk::LaplacianRecursiveGaussianImageFilter< gfpImageType , gfpImageType> GFilterType;
		GFilterType::Pointer gauss = GFilterType::New();
		
		float sigmas[] =  { 2.0f, 2.8284f, 4.0f, 5.6569f, 8.0f, 11.31f };	//LoG scales
		for (unsigned int i = 0; i < 6; ++i)
		{
			stringstream out;
			out<<i;
			string s = out.str();
			std::string tempFileName_3 = gfpFileNameNoExt + "_LOG_" + s + ".mhd";
			
			gauss->SetInput( img_trace );
			gauss->SetSigma( sigmas[i] );
			gauss->SetNormalizeAcrossScale(false);
			gauss->Update();
			
			writeImage<gfpImageType>(gauss->GetOutput(),tempFileName_3.c_str());
			std::cout<<std::endl<<"Done Writing "<<tempFileName_3.c_str()<<std::flush;

			std::cout<<std::endl<<"Done scale: "<<i<<std::flush;
		}
		std::string tempFileName_4 = gfpFileNameNoExt + "_PREP_MNT.mhd";
		writeImage<gfpImageType>(img_trace,tempFileName_4.c_str());
	}
	
	
	if( onlyTrace == 2 )
		
	{
		std::cout<<std::endl<<"TRACE_TOTAL";
		
		// to change 
		itk::Size<3> size_gfp_montage;
		size_gfp_montage = readSizeFromFile(argv[2]);
		
		//#####################################################################################################################
		//	MULTIPLE NEURON TACER
		//#####################################################################################################################

		somaCentroidsTable = ftk::LoadTable(temp + "_soma_centroids_table.txt");
		std::cout<<std::endl<<temp + "_soma_montage.mhd";

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

		std::cout << "Number of cells to be traced : " << centroid_list.size() << "\n";

		std::string SWCFilename = filePath + "/TracesAndSomas/OnlySWC.xml";
		std::ofstream outfile;
		outfile.open(SWCFilename.c_str());
		outfile << "<?xml\tversion=\"1.0\"\t?>\n";
		outfile << "<Source>\n\n";

		omp_set_nested(0);
		omp_set_max_active_levels(1);
		itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1); // This one can not be changed
		itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1); // This one can chenga
		
		std::cout<<std::endl<<"Num_threads: "<<itk::MultiThreader::GetGlobalDefaultNumberOfThreads();
		std::cout<<std::endl<<"Num_threads: "<<itk::MultiThreader::GetGlobalMaximumNumberOfThreads();
		
		num_threads = 80;
		omp_set_num_threads(num_threads);

		int counterCentro = 0;
		int tileSizeX = 500;
		int tileSizeY = 500;
		int tileSizeZ = 600;
		
		/////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////
		int numDivisionsInRowLOG = (int)floor((double)size_gfp_montage[1]/(2*tileSizeY+10)); // Minimum value 3, otherwise does not make sense
		if( numDivisionsInRowLOG < 3 )
			numDivisionsInRowLOG = 2;
		int numDivisionsInRowCEN = numDivisionsInRowLOG - 1;
		
		std::vector< itk::Size<3> > sizeOfBigTilesLOG;
		sizeOfBigTilesLOG.resize(numDivisionsInRowLOG);
		std::vector< itk::Size<3> > sizeOfBigTilesCEN;
		sizeOfBigTilesCEN.resize(numDivisionsInRowLOG-1);
		for( unsigned int jj=0; jj<numDivisionsInRowLOG; ++jj )
		{
			sizeOfBigTilesLOG[jj][0] = size_gfp_montage[0];
			sizeOfBigTilesLOG[jj][1] = (unsigned long long)floor((double)size_gfp_montage[1]/(double)numDivisionsInRowLOG);
			sizeOfBigTilesLOG[jj][2] = size_gfp_montage[2];
		}
		sizeOfBigTilesLOG[numDivisionsInRowLOG-1][1] = size_gfp_montage[1] - (numDivisionsInRowLOG-1)*sizeOfBigTilesLOG[numDivisionsInRowLOG-1][1];
		if( sizeOfBigTilesLOG[0][1] <2*tileSizeY)
		{
			std::cout<<std::endl<<"MISTAKE_1";
		}
		sizeOfBigTilesCEN[0][0] = size_gfp_montage[0];
		sizeOfBigTilesCEN[0][1] = sizeOfBigTilesLOG[0][1]+(unsigned long long)floor((double)sizeOfBigTilesLOG[0][1]/2);
		sizeOfBigTilesCEN[0][2] = size_gfp_montage[2];
		for( unsigned int jj=1; jj<numDivisionsInRowCEN; ++jj )
		{
			sizeOfBigTilesCEN[jj][0] = size_gfp_montage[0];
			sizeOfBigTilesCEN[jj][1] = sizeOfBigTilesLOG[jj][1];
			sizeOfBigTilesCEN[jj][2] = size_gfp_montage[2];
		}
		if( numDivisionsInRowCEN == 1 )
		{
			sizeOfBigTilesCEN[numDivisionsInRowCEN-1][1] = size_gfp_montage[1];
		}
		else
		{
			sizeOfBigTilesCEN[numDivisionsInRowCEN-1][1] = size_gfp_montage[1] -sizeOfBigTilesCEN[0][1] -  (numDivisionsInRowCEN-2)*sizeOfBigTilesLOG[numDivisionsInRowLOG-2][1];
		}
		
		
		unsigned long long total = 0;
		for( unsigned int jj=0; jj<numDivisionsInRowLOG; ++jj )
		{
			std::cout<<std::endl<<"Row Size Big"<<sizeOfBigTilesLOG[jj][1];
			total = total + sizeOfBigTilesLOG[jj][1];
		}
		if( total != size_gfp_montage[1] )
		{
			std::cout<<std::endl<<"MISTAKE_2";
		}
		std::cout<<std::endl<<"\t\t"<<total<<" "<<size_gfp_montage[1]<<" "<<size_gfp_montage[0];
		
		unsigned long long total2 = 0;
		for( unsigned int jj=0; jj<numDivisionsInRowCEN; ++jj )
		{
			std::cout<<std::endl<<"Col Size Big"<<sizeOfBigTilesCEN[jj][1];
			total2 = total2 + sizeOfBigTilesCEN[jj][1];
		}
		if( total2 != size_gfp_montage[1] )
		{
			std::cout<<std::endl<<"MISTAKE_3";
		}
		std::cout<<std::endl<<"\t\t"<<total2<<" "<<size_gfp_montage[1]<<" "<<size_gfp_montage[0];
		
		std::vector< itk::Index<3> > initialBigTileLOG;
		initialBigTileLOG.resize(numDivisionsInRowLOG);
		initialBigTileLOG[0][0] = 0;
		initialBigTileLOG[0][1] = 0;
		initialBigTileLOG[0][2] = 0;
		for( unsigned int jj=1; jj<numDivisionsInRowLOG; ++jj )
		{
			initialBigTileLOG[jj][0] = 0;
			initialBigTileLOG[jj][1] = sizeOfBigTilesLOG[jj-1][1] + initialBigTileLOG[jj-1][1];
			initialBigTileLOG[jj][2] = 0;
		}
		
		std::vector< itk::Index<3> > initialBigTileCEN;
		initialBigTileCEN.resize(numDivisionsInRowCEN);
		initialBigTileCEN[0][0] = 0;
		initialBigTileCEN[0][1] = 0;
		initialBigTileCEN[0][2] = 0;
		for( unsigned int jj=1; jj<numDivisionsInRowCEN; ++jj )
		{
			initialBigTileCEN[jj][0] = 0;
			initialBigTileCEN[jj][1] = sizeOfBigTilesCEN[jj-1][1] + initialBigTileCEN[jj-1][1];
			initialBigTileCEN[jj][2] = 0;
		}
		
		
		for( unsigned int bigTile = 0; bigTile<numDivisionsInRowCEN ; ++bigTile )
		{
			std::cout<<std::endl<<"ITER";
			
			stringstream out34;
			out34<<bigTile;
			string srr = out34.str();
		std::string SWCFilenameDivided = filePath + "/TracesAndSomasDivided/OnlySWC_" + srr +".xml";
		std::ofstream outfileDivided;
		outfileDivided.open(SWCFilenameDivided.c_str());
		outfileDivided << "<?xml\tversion=\"1.0\"\t?>\n";
		outfileDivided << "<Source>\n\n";

			itk::Index<3> initialBigIndexLOG;
			itk::Size<3> sizeOfTheRegionLOG;
			
			initialBigIndexLOG[0] = initialBigTileLOG[bigTile][0];
			initialBigIndexLOG[1] = initialBigTileLOG[bigTile][1];
			initialBigIndexLOG[2] = initialBigTileLOG[bigTile][2];
			
			sizeOfTheRegionLOG[0] = sizeOfBigTilesLOG[bigTile][0];
			sizeOfTheRegionLOG[1] = sizeOfBigTilesLOG[bigTile][1]+sizeOfBigTilesLOG[bigTile+1][1];
			sizeOfTheRegionLOG[2] = sizeOfBigTilesLOG[bigTile][2];
			
			gfpImageType::RegionType desiredRegionBigTileLOG;
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
			
			
// 			gfpImageType::RegionType desiredRegionBigTileCEN;
// 			desiredRegionBigTileCEN.SetSize(sizeOfTheRegionCEN);
// 			desiredRegionBigTileCEN.SetIndex(initialBigIndexCEN);
// 			
// 
// 			montageLabelType::RegionType desiredRegionBigTileMON;
// 			desiredRegionBigTileMON.SetSize(sizeOfTheRegionCEN);
// 			desiredRegionBigTileMON.SetIndex(initialBigIndexCEN);
			
			
			std::vector< gfpImageType::Pointer > LoGDesiredRegion;
			LoGDesiredRegion.resize(6);
			
			
			montageLabelType::Pointer somaMontageDesiredRegion;
			gfpImageType::Pointer img_traceDesiredRegion;
			
			#pragma omp parallel for num_threads(8)
			for (unsigned int i = 0; i < 8; ++i)
			{
			//	if( i<6 )
			//	{
			//		stringstream out;
			//		out<<i;
			//		string s = out.str();
			//		std::string tempFileName_3 = gfpFileNameNoExt + "_LOG_" + s + ".mhd";
			//		
			//		LoGDesiredRegion[i] = readImageRegion< gfpImageType >( tempFileName_3.c_str(), desiredRegionBigTileLOG );
			//	}
				
				if( i==6)
				{
					std::string tempFileName_42 = gfpFileNameNoExt + "_PREP_MNT.mhd";
					img_traceDesiredRegion = readImageRegion< gfpImageType >( tempFileName_42.c_str(), desiredRegionBigTileLOG );
					
					std::string tempFileName_46 = gfpFileNameNoExt + ".mhd";
					rawImageType_16bit::Pointer img_traceDesiredRegion_tif = readImageRegion< rawImageType_16bit >( tempFileName_46.c_str(), desiredRegionBigTileLOG );
					
					RescaleFilterType::Pointer rescaleFilter_bigt = RescaleFilterType::New();
					rescaleFilter_bigt->SetOutputMinimum(0);
					rescaleFilter_bigt->SetOutputMaximum(std::numeric_limits<unsigned char>::max());
					rescaleFilter_bigt->SetInput(img_traceDesiredRegion_tif);
					rescaleFilter_bigt->Update();
					
					std::string tempFileName_47 = filePath + "/BigTile_" + srr + "_GFP.tif";
					writeImage< rawImageType_8bit >(rescaleFilter_bigt->GetOutput(),tempFileName_47.c_str());
					
					

					
				}
// 			
				if( i ==7 )
				{
					std::string somaName = temp + "_soma_montage.mhd";
					somaMontageDesiredRegion = readImageRegion< montageLabelType >( somaName.c_str(), desiredRegionBigTileLOG );
// 					std::string tempFileName_45 = filePath + "/BigTile_" + srr + "_SOMA_PRE.mhd";
// 					writeImage< montageLabelType >(somaMontageDesiredRegion,tempFileName_45.c_str());
				}

				
// 				stringstream out2;
// 				out2<<i<<"_"<<bigTile;
// 				string s2 = out2.str();
// 				std::string tempFileName_32 = gfpFileNameNoExt + "_LOG_" + s2 + ".mhd";
// 				
// 				writeImage< gfpImageType >(LoGDesiredRegion[i],tempFileName_32.c_str());
			}

		
		#pragma omp parallel for num_threads(80) schedule(dynamic, 1)
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
				if(x >= tileSizeX/2)
					ssx_off << x - tileSizeX/2;
				else 
					ssx_off << 0;
				if(y >= tileSizeY/2)
					ssy_off << y - tileSizeY/2;
				else 
					ssy_off << 0;
				if(size_gfp_montage[2] <= tileSizeZ)
					ssz_off << 0;
				else
				{	
					if(z >= tileSizeZ/2)
						ssz_off << z - tileSizeZ/2;
					else 
						ssz_off << 0;
				}
				
				int x_local = x;
				int y_local = y - initialBigTileLOG[bigTile][1];
				int z_local = z;
				
				if(x_local >= tileSizeX/2)
					ssx_offBig << x_local - tileSizeX/2;
				else 
					ssx_off << 0;
				if(y_local >= tileSizeY/2)
					ssy_offBig << y_local - tileSizeY/2;
				else 
					ssy_off << 0;
				if(size_gfp_montage[2] <= tileSizeZ)
					ssz_offBig << 0;
				else
				{	
					if(z_local >= tileSizeZ/2)
						ssz_offBig << z_local - tileSizeZ/2;
					else 
						ssz_offBig << 0;
				}

				//########    CROP THE DESIRED DICE FROM THE GFP AND SOMA MONTAGES   ########
				montageLabelType::IndexType start;
				start[0] = ((x - tileSizeX/2)>0) ? (x - tileSizeX/2):0;
				start[1] = ((y - tileSizeY/2)>0) ? (y - tileSizeY/2):0;
				if(size_gfp_montage[2] <= tileSizeZ)
					start[2] = 0;
				else
					start[2] = ((z - tileSizeZ/2)>0) ? (z - tileSizeZ/2):0;

				rawImageType::IndexType start2;
				start2[0] = ((x - tileSizeX/2)>0) ? (x - tileSizeX/2):0;
				start2[1] = ((y - tileSizeY/2)>0) ? (y - tileSizeY/2):0;
				if(size_gfp_montage[2] <= tileSizeZ)
					start2[2] = 0;
				else
					start2[2] = ((z - tileSizeZ/2)>0) ? (z - tileSizeZ/2):0;

				montageLabelType::SizeType size;
				size[0] = ((x+tileSizeX/2)<size_gfp_montage[0]) ? tileSizeX : (tileSizeX/2+size_gfp_montage[0]-x-1); 
				size[1] = ((y+tileSizeY/2)<size_gfp_montage[1]) ? tileSizeY : (tileSizeY/2+size_gfp_montage[1]-y-1);
				if(size_gfp_montage[2] <= tileSizeZ)
					size[2] = size_gfp_montage[2];
				else
					size[2] = ((z+tileSizeZ/2)<size_gfp_montage[2]) ? tileSizeZ : (tileSizeZ/2+size_gfp_montage[2]-z-1);

				rawImageType::SizeType size2;
				size2[0] = ((x+tileSizeX/2)<size_gfp_montage[0]) ? tileSizeX : (tileSizeX/2+size_gfp_montage[0]-x-1);
				size2[1] = ((y+tileSizeY/2)<size_gfp_montage[1]) ? tileSizeY : (tileSizeY/2+size_gfp_montage[1]-y-1);
				if(size_gfp_montage[2] <= tileSizeZ)
					size2[2] = size_gfp_montage[2];
				else
					size2[2] = ((z+tileSizeZ/2)<size_gfp_montage[2]) ? tileSizeZ : (tileSizeZ/2+size_gfp_montage[2]-z-1);

				montageLabelType::RegionType desiredRegion;
				desiredRegion.SetSize(size);
				desiredRegion.SetIndex(start);

				rawImageType::RegionType desiredRegion2;
				desiredRegion2.SetSize(size2);
				desiredRegion2.SetIndex(start2);
				
				rawImageType_16bit::RegionType desiredRegion2_16bits;
				desiredRegion2_16bits.SetSize(size2);
				desiredRegion2_16bits.SetIndex(start2);
				
				gfpImageType::RegionType desiredRegion2_float;
				desiredRegion2_float.SetSize(size2);
				desiredRegion2_float.SetIndex(start2);
				

				LabelROIFilterType::Pointer ROIfilter3 = LabelROIFilterType::New();
				ROIfilter3->SetRegionOfInterest(desiredRegion);
				ROIfilter3->SetInput(somaMontageDesiredRegion);
		#pragma omp critical
				ROIfilter3->Update();
				LabelDuplicatorType::Pointer LabelDuplicator = LabelDuplicatorType::New();
				LabelDuplicator->SetInputImage(ROIfilter3->GetOutput());
		#pragma omp critical
				LabelDuplicator->Update();
				montageLabelType::Pointer img_soma = LabelDuplicator->GetOutput();
// 
				gfpROIFilterType::Pointer ROIfilter2 = gfpROIFilterType::New();
				ROIfilter2->SetRegionOfInterest(desiredRegion2_float);
				ROIfilter2->SetInput(img_traceDesiredRegion);
		#pragma omp critical
				ROIfilter2->Update();
				gfpDuplicatorType::Pointer floatDuplicator = gfpDuplicatorType::New();
				floatDuplicator->SetInputImage(ROIfilter2->GetOutput());
		#pragma omp critical
				floatDuplicator->Update();

				gfpImageType::Pointer img_trace = floatDuplicator->GetOutput();
				
				std::cout<<std::cout<<"end of the criticals";
				
				//########    FETCH ALL CENTROIDS THAT FALL WITHIN THE DICE    ########
				itk::Index<3> centroid;
				centroid[0] = ((x - tileSizeX/2)>0) ? tileSizeX/2:x; 
				centroid[1] = ((y - tileSizeY/2)>0) ? tileSizeY/2:y;
				if(size_gfp_montage[2] <= tileSizeZ)
					centroid[2] = z;
				else
					centroid[2] = ((z - tileSizeZ/2)>0) ? tileSizeZ/2:z;
// 
				std::vector< itk::Index<3> > soma_Table;      
				for(int ctr =0; ctr<centroid_list.size() ; ++ctr)
				{
					itk::Index<3> cen = centroid_list[ctr];
// 					//if(abs((double)(cen[0]-x))<=200 && abs((double)(cen[1]-y))<=200 && abs((double)(cen[2]-z))<=rowDivisor )
					if( (cen[0]>=start2[0]) && (cen[0]<(start2[0]+size2[0])) && (cen[1]>=start2[1]) && (cen[1]<(start2[1]+size2[1])) && (cen[2]>=start2[2]) && (cen[2]<(start2[2]+size2[2])) )
					{
						itk::Index<3> centroid2;
						centroid2[0] = centroid[0] + cen[0] - x;
						centroid2[1] = centroid[1] + cen[1] - y;
						centroid2[2] = centroid[2] + cen[2] - z;
						soma_Table.push_back(centroid2);
					}
				}
				
// 				#pragma omp critical
// 				{
// 				stringstream out2;
// 				out2<<"gfp2_"<<a<<"_";
// 				string s2 = out2.str();
// 				std::string tempFileName_32 = gfpFileNameNoExt + "_LOGP_" + s2 + ".mhd";
// // 				
// 				writeImage< gfpImageType >(img_trace,tempFileName_32.c_str());
// // 				writeImage< gfpImageType >(LoGDesiredRegion[0],tempFileName_32.c_str());
// 				
// 				}
// 				#pragma omp critical
// 				{
// 				stringstream out2;
// 				out2<<"soma";
// 				string s2 = out2.str();
// 				std::string tempFileName_32 = gfpFileNameNoExt + "_LOGP_" + s2 + ".mhd";
// // 				
// 				writeImage< montageLabelType >(img_soma,tempFileName_32.c_str());
// 				}
// 				
// 				int ytt;
// 				std::cin>>ytt;
				
// 				// DEBUG
// 	// 			stringstream out;
// 	// 			out<<a;
// 	// 			string s = out.str();
// 	// 			std::string temp2 = "/data/nicolas/data/0113/Debug/imagetrace" + s + ".mhd";
// 	// 			writeImage<gfpImageType>(img_trace,temp2.c_str());
// 	// 			
// 	// 			std::string temp3 = "/data/nicolas/data/0113/Debug/imagesoma" + s + ".mhd";
// 	// 			writeImage<montageLabelType>(img_soma,temp3.c_str());
// 	// 			
// 	// 			std::string temp5 = "/data/nicolas/data/0113/Debug/imagesoma_centroids" + s + ".txt";
// 	// 			std::ofstream of;
// 	// 			of.open(temp5.c_str());
// 	// 			for(int i=0; i<(int)soma_Table.size(); ++i)
// 	// 			{
// 	// 				of << soma_Table[0] << "\t" << soma_Table[1] << "\t" << soma_Table[2] << "\n";
// 	// 			}
// 	// 			of.close();
// 
// 	// 			double stopingTime = 10;
// 	// 			double curScaling = 0.5;
// 	// 			double rmsThrehold = 0.02;
// 	// 			SomaExtractor *Somas = new SomaExtractor();
// 	// 			montageLabelType::Pointer  img_soma_yan = Somas->SegmentSoma(img_trace, soma_Table,stopingTime,curScaling,rmsThrehold);
// 
// 				//########    RUN TRACING    ########
				MultipleNeuronTracer * MNT = new MultipleNeuronTracer();
				
			//	MNT->setLogScale(LoGDesiredRegion[0],0);
			//	MNT->setLogScale(LoGDesiredRegion[1],1);
			//	MNT->setLogScale(LoGDesiredRegion[2],2);
			//	MNT->setLogScale(LoGDesiredRegion[3],3);
			//	MNT->setLogScale(LoGDesiredRegion[4],4);
			//	MNT->setLogScale(LoGDesiredRegion[5],5);
				
// 				MNT->setNDXImage(NDXImage);
				
				
// 				
				MNT->setDiceSize(size);
				MNT->setDiceIndex(start);
// 				
// 				MNT->LoadCurvImage_1(img_trace, 0);
				MNT->LoadCurvImage_2(img_trace);
				MNT->ReadStartPoints_1(soma_Table, 0);
				MNT->SetCostThreshold(1000);
				MNT->LoadSomaImage_1(img_soma);
// 	// 			MNT->LoadSomaImage_1(img_soma_yan);
				bool flagLog = false;
				MNT->setFlagOutLog(flagLog);
				MNT->RunTracing();
// 
				x = min(tileSizeX/2, x);
				y = min(tileSizeY/2, y);
				if(size_gfp_montage[2] > tileSizeZ)
					z = min(tileSizeZ/2, z);
// 
//
				vtkSmartPointer< vtkTable > swcTable = MNT->GetSWCTable(0);
				std::string swcFilename = filePath + "/TracesAndSomas/Trace_" + ssx.str() + "_" + ssy.str() + "_" + ssz.str() + "_ANT.swc";
				WriteCenterTrace(swcTable, x, y, z, swcFilename);
// 				
				#pragma omp critical
				{
					outfile << "\t<File\tFileName=\"Trace_" << ssx.str() << "_" << ssy.str() << "_" << ssz.str() << "_ANT.swc\"\tType=\"Trace\"\ttX=\"" << ssx_off.str() << "\"\ttY=\"" << ssy_off.str() << "\"\ttZ=\"" << ssz_off.str() << "\"/>\n";
				}
				
// 				vtkSmartPointer< vtkTable > swcTable = MNT->GetSWCTable(0);
				std::string swcFilenameDivided = filePath + "/TracesAndSomasDivided/Trace_BigTile_" + srr + "_" + ssx.str() + "_" + ssy.str() + "_" + ssz.str() + "_ANT.swc";
				WriteCenterTrace(swcTable, x, y, z, swcFilenameDivided);
// 				
				#pragma omp critical
				{
					outfileDivided << "\t<File\tFileName=\"Trace_BigTile_" << srr << "_" << ssx.str() << "_" << ssy.str() << "_" << ssz.str() << "_ANT.swc\"\tType=\"Trace\"\ttX=\"" << ssx_offBig.str() << "\"\ttY=\"" << ssy_offBig.str() << "\"\ttZ=\"" << ssz_offBig.str() << "\"/>\n";
				}

				
				delete MNT;
			}
			
			outfileDivided << "\n</Source>";
			outfileDivided.close();
		}
		outfile << "\n</Source>";
		outfile.close();
	}
	return 0;
}


//#############################################################################################################
//#############################################################################################################
//  ADDITIONAL FINCTIONS
//#############################################################################################################
//#############################################################################################################



//#############################################################################################################
//  NUCLEAR SEGMENTATION
//#############################################################################################################
LabelType::Pointer RunNuclearSegmentation(rawImageType::Pointer inpImg, const char* fileName)
{
	rawImageType::PointType origin_row;
	origin_row[0] = 0; 
	origin_row[1] = 0;
	origin_row[2] = 0;
	inpImg->SetOrigin( origin_row );

	itk::SizeValueType size[3];
	size[0] = inpImg->GetLargestPossibleRegion().GetSize()[0];
	size[1] = inpImg->GetLargestPossibleRegion().GetSize()[1];
	size[2] = inpImg->GetLargestPossibleRegion().GetSize()[2];

	unsigned char *in_Image;
	in_Image = (unsigned char *) malloc (size[0]*size[1]*size[2]);
	if( in_Image == NULL )
	{
		std::cerr<<"Nucleus Seg failed because malloc failed\n";
		return NULL;
	}

	//memset(in_Image/*destination*/,0/*value*/,size[0]*size[1]*size[2]*sizeof(unsigned char)/*num bytes to move*/);
	typedef itk::ImageRegionIterator< rawImageType > IteratorType22;
	IteratorType22 pix_buf( inpImg, inpImg->GetRequestedRegion() );
	itk::SizeValueType ind=0;
	for ( pix_buf.GoToBegin(); !pix_buf.IsAtEnd(); ++pix_buf, ++ind )
	{
		in_Image[ind]=(pix_buf.Get());
	}
	yousef_nucleus_seg *NucleusSeg = new yousef_nucleus_seg();
	NucleusSeg->readParametersFromFile(fileName);
	NucleusSeg->setDataImage(in_Image,size[0],size[1],size[2],"test.mhd");
	NucleusSeg->runBinarization();

	std::cout<<std::endl << "\t\t Done binarization:" <<std::flush;

	try 
	{
		NucleusSeg->runSeedDetection();

		std::cout<<std::endl << "\t\t\t Done Seed Detect: " <<std::flush;

	}
	catch( bad_alloc & excp )
	{
		std::cout<<"You have requested more memory than "
			<<"what is currently available in this "
			<<"system, please try again with a smaller "
			<<"input image\n";
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cout<<"Error: " << excp <<std::endl;
	}

	NucleusSeg->runClustering();

	std::cout<<std::endl << "\t\t\t\t Done Clus: ";
	// 	}
	unsigned short *output_img;

	output_img=NucleusSeg->getClustImage();

	LabelType::Pointer image = LabelType::New();
	LabelType::PointType origin;
	origin[0] = 0; //Is this OK?
	origin[1] = 0;
	origin[2] = 0;
	image->SetOrigin( origin );
	LabelType::IndexType start1;
	start1[0] = 0;
	start1[1] = 0;
	start1[2] = 0;
	LabelType::SizeType size1;
	size1[0] = size[0];
	size1[1] = size[1];
	size1[2] = size[2];
	LabelType::RegionType region;
	region.SetSize ( size1  );
	region.SetIndex( start1 );
	image->SetRegions( region );
	image->Allocate();
	image->FillBuffer(0);
	image->Update();

	typedef itk::ImageRegionIteratorWithIndex< LabelType > IteratorType;
	IteratorType iterator1(image,image->GetRequestedRegion());
	for(itk::SizeValueType i=0; i<(size[0]*size[1]*size[2]); ++i)
	{
		//unsigned short val = (unsigned short)output_img[i];
		iterator1.Set(output_img[i]);
		++iterator1;
	}

	delete NucleusSeg;
	free(in_Image);
// 	delete output_img;
	std::cout<<std::endl << "\t\t\t\t HERE2: Done SEG: ";

	/*
	// DELETING ALL LABELS WITH LESS THAN 8610 PIXELS OF VOLUME
	typedef itk::LabelGeometryImageFilter< LabelType, rawImageType > LabelGeometryType;
	LabelGeometryType::Pointer labelGeometryFilter = LabelGeometryType::New();	

	labelGeometryFilter->SetInput( image );
	labelGeometryFilter->SetIntensityInput( inpImg );

	labelGeometryFilter->CalculatePixelIndicesOff();
	labelGeometryFilter->CalculateOrientedBoundingBoxOff();
	labelGeometryFilter->CalculateOrientedLabelRegionsOff();
	labelGeometryFilter->CalculateOrientedIntensityRegionsOff();

	try
	{
	labelGeometryFilter->Update();
	}
	catch (itk::ExceptionObject & e) 
	{
	std::cerr << "Exception in Dirk's Label Geometry Filter: " << e << std::endl;
	return false;
	}

	std::vector< unsigned short > labels, del_labels;
	labels.clear();
	std::vector< LabelGeometryType::LabelPixelType > ls = labelGeometryFilter->GetLabels();
	for (int i = 0; i < (int)ls.size(); ++i)
	{
	unsigned short l = ls.at(i);
	labels.push_back( (unsigned short)l );
	}

	for (int i = 0; i < (int)labels.size(); ++i)
	{
	if( labelGeometryFilter->GetVolume( labels[i] ) < 8610 )
	del_labels.push_back( labels[i] );
	}

	std::cout<<std::endl<<"We have: " << (int)labels.size()-(int)del_labels.size()<<" seeds";

	LabelType::PixelType * imageArray = image->GetBufferPointer();
	itk::Size<3> size_image = image->GetLargestPossibleRegion().GetSize();
	unsigned long long image_slice_size = size_image[1] * size_image[0];
	unsigned long long image_row_size = size_image[0];	

	for(int i=0; i<size_image[2]; ++i)
	{
	for(int j=0; j<size_image[1]; ++j)
	{
	for(int k=0; k<size_image[0]; ++k)
	{
	unsigned short value = imageArray[(i*image_slice_size)+(j*image_row_size)+k];
	std::vector<unsigned short>::iterator posn1 = find(del_labels.begin(), del_labels.end(), value);
	if(posn1 != del_labels.end())
	{
	imageArray[(i*image_slice_size)+(j*image_row_size)+k] = 0;
	}
	}
	}
	}*/

	return image;
}


//#############################################################################################################
//  COMPUTE FEATURES
//#############################################################################################################

vtkSmartPointer<vtkTable> ComputeFeaturesAndAssociations(rawImageType::Pointer nucImg, rawImageType::Pointer gfpImg, rawImageType::Pointer cy5Img, LabelType::Pointer labImg, const char* fileName)
{
	ftk::ProjectDefinition projectDef;
	projectDef.Load( fileName );

	vtkSmartPointer<vtkTable> table = NULL;
	ftk::Image::Pointer nucSegRes = NULL;
	ftk::ProjectProcessor * myProc = new ftk::ProjectProcessor();

	itk::SizeValueType tileSize[3];
	tileSize[0] = nucImg->GetLargestPossibleRegion().GetSize()[0];
	tileSize[1] = nucImg->GetLargestPossibleRegion().GetSize()[1];
	tileSize[2] = nucImg->GetLargestPossibleRegion().GetSize()[2];
	
	std::cout<<std::endl<<"\t\t TRTRTR: "<<fileName<<" "<<tileSize[0]<<" "<<tileSize[1]<<" "<<tileSize[2];

	ftk::Image::Pointer sourceImages = ftk::Image::New();
	std::vector< unsigned char > color;
	color.push_back(255); color.push_back(0); color.push_back(0);
	sourceImages->AppendChannelFromData3D( nucImg->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), tileSize[0], tileSize[1], tileSize[2], "dapi", color, true);
	color[0] = 0; color[1] = 255; color[2] = 0;
	sourceImages->AppendChannelFromData3D( gfpImg->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), tileSize[0], tileSize[1], tileSize[2], "gfp", color, true);
// 	color[0] = 255; color[1] = 255; color[2] = 0;
// 	sourceImages->AppendChannelFromData3D( cy5Img->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), tileSize[0], tileSize[1], tileSize[2], "cy5", color, true);

	ftk::Image::Pointer labelImage = ftk::Image::New();
	color[0] = 255; color[1] = 255; color[2] = 255;
	labelImage->AppendChannelFromData3D( labImg->GetBufferPointer(), itk::ImageIOBase::USHORT, sizeof(unsigned short), tileSize[0], tileSize[1], tileSize[2], "dapi", color, true);

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


//#############################################################################################################
//  LABEL TO CENTROID MAP
//#############################################################################################################
std::map< unsigned int, itk::Index<3> > GetLabelToCentroidMap( vtkSmartPointer< vtkTable > table)
{
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

//#############################################################################################################
//  WRITE ONLY THE CENTER TRACE OF THE DICE
//#############################################################################################################
void WriteCenterTrace(vtkSmartPointer< vtkTable > swcNodes, int x, int y, int z, std::string filename)
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

void RunEverything(rawImageType::RegionType region, rawImageType::Pointer m_nuc, rawImageType::Pointer m_gfp, rawImageType::Pointer m_cy5, std::vector<LabelType::Pointer> &myLabel, std::vector < vtkSmartPointer< vtkTable > > &myTable, std::map< unsigned int, itk::Index<3> > *myCentroids, const char* segParams, const char* projectDef, int col)
{
	rawImageType::Pointer tile_nuc;
	rawImageType::Pointer tile_gfp;
	rawImageType::Pointer tile_cy5;

	rawROIFilterType::Pointer ROIfilter_tile_nuc = rawROIFilterType::New();
	ROIfilter_tile_nuc->SetRegionOfInterest(region);
	ROIfilter_tile_nuc->SetInput(m_nuc);
#pragma omp critical
	ROIfilter_tile_nuc->Update();
	rawDuplicatorType::Pointer rawDuplicator_tile_nuc = rawDuplicatorType::New();
	rawDuplicator_tile_nuc->SetInputImage(ROIfilter_tile_nuc->GetOutput());
#pragma omp critical
	rawDuplicator_tile_nuc->Update();
	tile_nuc = rawDuplicator_tile_nuc->GetOutput();

	rawROIFilterType::Pointer ROIfilter_tile_gfp = rawROIFilterType::New();
	ROIfilter_tile_gfp->SetRegionOfInterest(region);
	ROIfilter_tile_gfp->SetInput(m_gfp);
#pragma omp critical
	ROIfilter_tile_gfp->Update();
	rawDuplicatorType::Pointer rawDuplicator_tile_gfp = rawDuplicatorType::New();
	rawDuplicator_tile_gfp->SetInputImage(ROIfilter_tile_gfp->GetOutput());
#pragma omp critical
	rawDuplicator_tile_gfp->Update();
	tile_gfp = rawDuplicator_tile_gfp->GetOutput();

// 	rawROIFilterType::Pointer ROIfilter_tile_cy5 = rawROIFilterType::New();
// 	ROIfilter_tile_cy5->SetRegionOfInterest(region);
// 	ROIfilter_tile_cy5->SetInput(m_cy5);
// #pragma omp critical
// 	ROIfilter_tile_cy5->Update();
// 	rawDuplicatorType::Pointer rawDuplicator_tile_cy5 = rawDuplicatorType::New();
// 	rawDuplicator_tile_cy5->SetInputImage(ROIfilter_tile_cy5->GetOutput());
// #pragma omp critical
// 	rawDuplicator_tile_cy5->Update();
// 	tile_cy5 = rawDuplicator_tile_cy5->GetOutput();
	//	}

	myLabel[col] = RunNuclearSegmentation(tile_nuc, segParams);
	myTable[col] = ComputeFeaturesAndAssociations(tile_nuc, tile_gfp, tile_cy5, myLabel[col], projectDef );
	*myCentroids = GetLabelToCentroidMap(myTable[col]);
}

rawImageType_8bit::Pointer readAndRescale_16to8(const char* nameInput){

	rawReaderType::Pointer reader_nuc = rawReaderType::New();
	reader_nuc->SetFileName( nameInput );
	reader_nuc->Update();

	RescaleFilterType::Pointer rescaleFilter_nuc = RescaleFilterType::New();
	rescaleFilter_nuc->SetOutputMinimum(0);
	rescaleFilter_nuc->SetOutputMaximum(std::numeric_limits<unsigned char>::max());
	rescaleFilter_nuc->SetInput(reader_nuc->GetOutput());
	rescaleFilter_nuc->Update();
	rawImageType_8bit::Pointer montage_nuc = rescaleFilter_nuc->GetOutput();
	montage_nuc->DisconnectPipeline();

	return montage_nuc;
}


rawImageType_16bit::Pointer readAndRescale_16to16(const char* nameInput){

	rawReaderType::Pointer reader_nuc = rawReaderType::New();
	reader_nuc->SetFileName( nameInput );
	reader_nuc->Update();

	RescaleFilterType_16to16::Pointer rescaleFilter_nuc = RescaleFilterType_16to16::New();
	rescaleFilter_nuc->SetOutputMinimum(0);
	rescaleFilter_nuc->SetOutputMaximum(std::numeric_limits<unsigned short>::max());
	rescaleFilter_nuc->SetInput(reader_nuc->GetOutput());
	rescaleFilter_nuc->Update();
	rawImageType_16bit::Pointer montage_nuc = rescaleFilter_nuc->GetOutput();
	montage_nuc->DisconnectPipeline();

	return montage_nuc;
}


rawImageType_8bit::Pointer readAndRescale_16to8(const char* nameInput, int debugCopy){

	rawReaderType::Pointer reader_nuc = rawReaderType::New();
	reader_nuc->SetFileName( nameInput );
	reader_nuc->Update();

	RescaleFilterType::Pointer rescaleFilter_nuc = RescaleFilterType::New();
	rescaleFilter_nuc->SetOutputMinimum(0);
	rescaleFilter_nuc->SetOutputMaximum(std::numeric_limits<unsigned char>::max());
	rescaleFilter_nuc->SetInput(reader_nuc->GetOutput());
	rescaleFilter_nuc->Update();
	rawImageType_8bit::Pointer montage_nuc = rescaleFilter_nuc->GetOutput();
	montage_nuc->DisconnectPipeline();
	
	std::string nucFileName(nameInput);
	std::string temp = nucFileName;
	string::iterator it;
	it = temp.end() - 5;
	temp.erase(it, it+5);
	
	stringstream out;
	out<<debugCopy;
	string s = out.str();
	
	std::string temp2 = temp + s + ".mhd";
	
	rawWriterType_8::Pointer writer34 = rawWriterType_8::New();
	writer34->SetFileName(temp2.c_str());
	writer34->SetInput(montage_nuc);
	writer34->Update();
	
	return montage_nuc;
}


gfpImageType::Pointer preprocessingMNT( gfpImageType::Pointer montage_gfp )
{
// 	CasterFilterType_16Tofloat::Pointer caster = CasterFilterType_16Tofloat::New();
// 	caster->SetInput(montage_gfp);
// 	caster->Update();
// 	std::cout<<std::endl<<"Caster"<<std::flush;
// 	gfpImageType::Pointer montage_gfp_float = caster->GetOutput();

	RescaleFilterType_floToflo::Pointer rescaler = RescaleFilterType_floToflo::New();
	rescaler->SetOutputMinimum(0.0);
	rescaler->SetOutputMaximum(1.0);
	rescaler->SetInput(montage_gfp);
// 	rescaler->Update();
	std::cout<<std::endl<<"Rescaler"<<std::flush;

	MedianFilterType::Pointer medfilt = MedianFilterType::New();
	medfilt->SetInput(rescaler->GetOutput());
	gfpImageType::SizeType rad = { {1, 1, 1} };
	medfilt->SetRadius(rad);
	medfilt->Update();
	std::cout<<std::endl<<"MedFilter"<<std::flush;
	gfpImageType::Pointer CurvImage = medfilt->GetOutput();
	
	return CurvImage;
}



itk::Size<3> readSizeFromFile(const char* name)
{
	itk::Size<3> size_gfp_montage;
	rawImageType_16bit::Pointer montage_gfp;
	montage_gfp = readImage< rawImageType_16bit >(name);
	size_gfp_montage[0] = montage_gfp->GetLargestPossibleRegion().GetSize()[0];
	size_gfp_montage[1] = montage_gfp->GetLargestPossibleRegion().GetSize()[1];
	size_gfp_montage[2] = montage_gfp->GetLargestPossibleRegion().GetSize()[2];
	return size_gfp_montage;
}
