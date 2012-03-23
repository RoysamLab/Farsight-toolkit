// make -j80;./exe/darpa_tracer_w_seg ../../data/0103/DAPI_montage.mhd ../../data/0103/GFP_montage.mhd ../../data/0103/Cy5_montage_BS.mhd ../../data/0103/Seg_Params.ini ../../data/0103/ProjectDefinition.xml
// make -j80;((time ./exe/darpa_tracer_w_seg ../../data/0103/DAPI_montage.mhd ../../data/0103/GFP_montage.mhd ../../data/0103/Cy5_montage_BS.mhd ../../data/0103/Seg_Params.ini ../../data/0103/ProjectDefinition.xml) 2>&1) >> salida.log

#include "itkImageFileReader.h"
#include "itkImage.h"
#include "itkIntTypes.h"
#include "boost/tokenizer.hpp"
#include <fstream>
#include "vul/vul_file.h"
#include "iostream"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "../Tracing/MultipleNeuronTracer/MultipleNeuronTracer.h"
#include "../NuclearSegmentation/exe/SomaExtraction.h"


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

typedef unsigned char inputPixelType;
typedef itk::Image<inputPixelType,  3> rawImageType;
typedef itk::Image<unsigned char,  3> rawImageType_8bit;
typedef itk::Image<unsigned short,  3> rawImageType_16bit;
typedef itk::Image<unsigned short, 3> LabelType;
typedef itk::Image<unsigned int, 3> montageLabelType;

typedef MultipleNeuronTracer::ImageType3D gfpImageType; // Float Image
typedef itk::ImageDuplicator< rawImageType >  rawDuplicatorType;
typedef itk::ImageDuplicator< montageLabelType >  LabelDuplicatorType;
typedef itk::ImageDuplicator< gfpImageType >  gfpDuplicatorType;
typedef itk::RegionOfInterestImageFilter< rawImageType, rawImageType > rawROIFilterType;
typedef itk::RegionOfInterestImageFilter< gfpImageType, gfpImageType > gfpROIFilterType;
typedef itk::RegionOfInterestImageFilter< montageLabelType, montageLabelType > LabelROIFilterType;
typedef std::pair< LabelType::Pointer, vtkSmartPointer< vtkTable > > segResultsType;

typedef itk::ImageFileReader<rawImageType_16bit> rawReaderType;
typedef itk::ImageFileWriter<montageLabelType> labelWriterType;
typedef itk::ImageFileReader<montageLabelType> labelReaderType;
typedef itk::RescaleIntensityImageFilter< rawImageType_16bit, rawImageType_8bit > RescaleFilterType;
typedef itk::ImageFileWriter<rawImageType_8bit> rawWriterType_8;


//segResultsType RunNuclearSegmentation( rawImageType::Pointer, rawImageType::Pointer, const char*);
inline LabelType::Pointer RunNuclearSegmentation(rawImageType::Pointer, const char*);
vtkSmartPointer<vtkTable> ComputeFeaturesAndAssociations(rawImageType::Pointer, rawImageType::Pointer, rawImageType::Pointer, LabelType::Pointer, const char* );
void RunEverything(rawImageType::RegionType, rawImageType_8bit::Pointer, rawImageType_8bit::Pointer, rawImageType_8bit::Pointer, std::vector< LabelType::Pointer> &, std::vector < vtkSmartPointer< vtkTable > >&, std::map< unsigned int, itk::Index<3> >*, const char*, const char*, int col);
std::map< unsigned int, itk::Index<3> > GetLabelToCentroidMap( vtkSmartPointer< vtkTable > );
void WriteCenterTrace(vtkSmartPointer< vtkTable >, int, int, int, std::string);

rawImageType_8bit::Pointer readAndRescale_16to8(const char*);
rawImageType_8bit::Pointer readAndRescale_16to8(const char*, int debugCopy);

//typedef struct{ double x_d; double y_d; double z_d; } Table;
bool myfunction (itk::SizeValueType i,itk::SizeValueType j) { return (i<j); }

int main(int argc, char* argv[])
{

	itk::MultiThreader::SetGlobalDefaultNumberOfThreads(80);
	std::cout<<std::endl<<"TEST"<<std::flush;
	omp_set_nested(1);

	omp_set_max_active_levels(2);

	int num_threads = 1;
	omp_set_num_threads(num_threads);


	int counterTiles = 0;
	int counterTiles2 = 0;
	int counterTiles3 = 0;
	if(argc < 4)
	{
		std::cout<<"usage1: darpa_tracer <dapi_montage_file> <gfp_montage_file> <cy5_montage_file> <seg_params_file> <project_definition_file> \n";
		return 0;
	}

	//argv[1] = "C:/Data/Darpa/TEST_FOR_PIPELINE/small_dapi.tif";
	//argv[2] = "C:/Data/Darpa/TEST_FOR_PIPELINE/small_gfp.tif";	
	//argv[3] = "C:/Data/Darpa/TEST_FOR_PIPELINE/ProjectDefinition.xml";

	//Input: Somas_file, nucleus_montage, trace_channel_montage, nuc_seg_parameters
	//std::fstream soma_centroid_file;
	//soma_centroid_file.open(argv[1], std::fstream::in);
	//std::vector<Table> centroid_list;
	//while (soma_centroid_file.good()){
	//	Table new_centroid;
	//	soma_centroid_file >> new_centroid.x_d >> new_centroid.y_d >> new_centroid.z_d;
	//	centroid_list.push_back( new_centroid );
	//}
	//soma_centroid_file.close();

	std::string MyName = argv[0];

	std::string nucFileName(argv[1]);
	std::string filePath = ftk::GetFilePath(nucFileName);

	std::string temp = nucFileName;
	string::iterator it;
	
	// IMPROVE THIS TIPE OF GETTING THE EXTENTION 5 FOR NRRD
	it = temp.end() - 5;
	temp.erase(it, it+5);

	
	// 	typedef itk::ImageFileReader<gfpImageType> gfpReaderType;




	vtkSmartPointer<vtkTable> somaCentroidsTable = NULL;
	montageLabelType::Pointer somaMontage;
	rawImageType_8bit::Pointer montage_gfp;
	itk::SizeValueType size_gfp_montage[3];

	int onlyTrace = atoi(argv[6]);
	
// 	if( onlyTrace == 2 )
// 	{
// 		rawImageType_8bit::Pointer montage_nuc = readAndRescale_16to8(argv[1],1);
// 		rawImageType_8bit::Pointer montage_gfp = readAndRescale_16to8(argv[2],2);
// 		rawImageType_8bit::Pointer montage_cy5 = readAndRescale_16to8(argv[3],3);
// 	}
	
	if( onlyTrace == 0 )
	{
		rawImageType_8bit::Pointer montage_nuc = readAndRescale_16to8(argv[1],1);
		rawImageType_8bit::Pointer montage_gfp = readAndRescale_16to8(argv[2],2);
		rawImageType_8bit::Pointer montage_cy5 = readAndRescale_16to8(argv[3],3);
		

		itk::SizeValueType size_nuc_montage[3];
		size_nuc_montage[0] = montage_nuc->GetLargestPossibleRegion().GetSize()[0];
		size_nuc_montage[1] = montage_nuc->GetLargestPossibleRegion().GetSize()[1];
		size_nuc_montage[2] = montage_nuc->GetLargestPossibleRegion().GetSize()[2];

		itk::SizeValueType size_gfp_montage[3];
		size_gfp_montage[0] = montage_gfp->GetLargestPossibleRegion().GetSize()[0];
		size_gfp_montage[1] = montage_gfp->GetLargestPossibleRegion().GetSize()[1];
		size_gfp_montage[2] = montage_gfp->GetLargestPossibleRegion().GetSize()[2];


		itk::SizeValueType size_cy5_montage[3];
		size_cy5_montage[0] = montage_cy5->GetLargestPossibleRegion().GetSize()[0];
		size_cy5_montage[1] = montage_cy5->GetLargestPossibleRegion().GetSize()[1];
		size_cy5_montage[2] = montage_cy5->GetLargestPossibleRegion().GetSize()[2];


		//#####################################################################################################################
		//	NUCLEAR SEGMENT THE MONTAGE TILE BY TILE AND STITCH THE RESULT TILES TOGETHER
		//#####################################################################################################################
		unsigned long long rowDivisor = 861;
		unsigned long long colDivisor = 640;
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


		itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);
		//itk::MultiThreader::SetGlobalDefaultNumberOfThreads(80); // JUST TO TEST
		//##################	SEGMENTING EACH ROW IN THE MONTAGE	  ###################
#pragma omp parallel for num_threads(7) schedule(dynamic, 1)
		for(int row=0; row<num_rows; ++row)
		{
			omp_set_nested(1);

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
#pragma omp parallel for num_threads(11) schedule(dynamic, 1)
			for(unsigned int col=0; col<num_cols; ++col)
			{

				int tid = omp_get_thread_num();
				//stringstream out;
				//out<<tid;
				//string s = out.str();

				//#pragma omp critical
				//{
				//	counterTiles++;
				//}

				//##################	EXTRACT A TILE AND START SEGMENTING	  ###################

				//#pragma omp critical
				//{
				//	++counterTiles;
				//	myfile<<"Extracting Tile " << col*rowDivisor << "_" << row*rowDivisor << " ThreadID " << counterTiles << "\n"<<std::flush;
				//}

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
				std::cout<<std::endl<<"\t\t\t\t ---->>>> Done with nucleus segmentation for Tile " << col << "_" << row << "\n";



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

			}
			omp_set_nested(0);
			// 		myfile.close();

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
			//rowSeg->Register();

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
					unsigned long long x_offset = (m*colDivisor) - (2*25);
					for(unsigned long long z=0; z<tileBorder_size[2]; ++z)
					{
						for(unsigned long long y=0; y<tileBorder_size[1]; ++y)
						{
							for(unsigned long long x=0; x<tileBorder_size[0]; ++x)
							{
								unsigned short value = myTileBorderArray[(tileBorder_slice_size*z) + (tileBorder_row_size*y) + (x)];
								if(value == 0) continue;
								unsigned long long lab_cen_x = Centroids_TileBorders[m-1][value][0];
								if((lab_cen_x < 25) || (lab_cen_x >= (3*25))) continue;
								rowSegArray[(row_slice_size*z) + (row_row_size*y) + (x_offset + x)] = max_value + value;
								if((max_value + value) > current_max)
									current_max = max_value + value;
							}
						}
					}
					//Label_TileBorders[m-1]->UnRegister();

					//##################	STITCHING THE TILE_BORDER_TABLE INTO THE ROW_TABLE	  ###################

					for(int r=0; r<(int)Table_TileBorders[m-1]->GetNumberOfRows(); ++r)
					{
						vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
						for(int c=0; c<(int)Table_TileBorders[m-1]->GetNumberOfColumns(); ++c)
						{
							if(c == 0)
								model_data1->InsertNextValue(vtkVariant(Table_TileBorders[m-1]->GetValue(r,c).ToInt() + max_value));
							else if(c == 1)
								model_data1->InsertNextValue(vtkVariant(Table_TileBorders[m-1]->GetValue(r,c).ToInt() + x_offset));
							else
								model_data1->InsertNextValue(Table_TileBorders[m-1]->GetValue(r,c));
						}
						rowTable->InsertNextRow(model_data1);
					}

					max_value = current_max;
				}

				//##################	STITCHING THE TILE INTO THE ROW	  ###################

				LabelType::Pointer myTile = Label_Tiles[m];
				LabelType::PixelType * myTileArray = myTile->GetBufferPointer();
				itk::Size<3> tile_size = myTile->GetLargestPossibleRegion().GetSize();
				unsigned long long tile_slice_size = tile_size[1] * tile_size[0];
				unsigned long long tile_row_size = tile_size[0];	
				unsigned long long x_offset = m*colDivisor;
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
							rowSegArray[(row_slice_size*z) + (row_row_size*y) + (x_offset + x)] = max_value + value;
							if((max_value + value) > current_max)
								current_max = max_value + value;
						}
					}
				}
				//Label_Tiles[m]->UnRegister();

				//##################	STITCHING THE TILE_TABLE INTO THE ROW_TABLE	  ###################

				for(unsigned long long r=0; r<(unsigned long long)Table_Tiles[m]->GetNumberOfRows(); ++r)
				{
					vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
					for(unsigned long long c=0; c<(unsigned long long)Table_Tiles[m]->GetNumberOfColumns(); ++c)
					{
						if(c == 0)
							model_data1->InsertNextValue(vtkVariant(Table_Tiles[m]->GetValue(r,c).ToInt() + max_value));
						else if(c == 1)
							model_data1->InsertNextValue(vtkVariant(Table_Tiles[m]->GetValue(r,c).ToInt() + x_offset));
						else
							model_data1->InsertNextValue(Table_Tiles[m]->GetValue(r,c));
					}
					rowTable->InsertNextRow(model_data1);
				}

				max_value = current_max;
			}

			std::cout<<std::endl<<"asdfasdf1 ";

			Table_Tiles.clear();
			Label_TileBorders.clear();

			//#############################################################
			// 		std::cout<<"done !!\n\n";
			Label_Rows[row] = rowSeg;
			Table_Rows[row] = rowTable;
			Centroids_Rows[row] = GetLabelToCentroidMap(rowTable);

			//##################	EXTRACT A ROW_BORDER AND START SEGMENTING	  ###################

			if(row != 0)
			{

				// 		std::cout<<"Extracting row border Y = " << row*rowDivisor << "\n";
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

		}

		std::cout<<"Stitching all rows ...";

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

		unsigned short max_value_1 = 0, current_max_1 = 0;
		for(unsigned long long n=0; n<(int)Label_Rows.size(); ++n)
		{
			// 		std::cout << "stitching row " << n << "\n";
			if(n != 0)
			{

				//##################	STITCHING THE ROW_BORDER INTO THE MONTAGE	  ###################

				// 			std::cout << "stitchin the row border\n" ;
				LabelType::Pointer myRowBorder = Label_RowBorders[n-1];
				LabelType::PixelType * myRowBorderArray = myRowBorder->GetBufferPointer();
				itk::Size<3> rowBorder_size = myRowBorder->GetLargestPossibleRegion().GetSize();
				unsigned long long rowBorder_slice_size = rowBorder_size[1] * rowBorder_size[0];
				unsigned long long rowBorder_row_size = rowBorder_size[0];	
				unsigned long long y_offset = (n*rowDivisor) - (2*25);
				for(unsigned long long z=0; z<rowBorder_size[2]; ++z)
				{
					for(unsigned long long y=0; y<rowBorder_size[1]; ++y)
					{
						for(unsigned long long x=0; x<rowBorder_size[0]; ++x)
						{
							unsigned short value_1 = myRowBorderArray[(rowBorder_slice_size*z) + (rowBorder_row_size*y) + (x)];
							if(value_1 == 0) continue;
							unsigned long long lab_cen_y = Centroids_RowBorders[n-1][value_1][1];
							if((lab_cen_y < 25) || (lab_cen_y >= (3*25))) continue;
							montageSegArray[(montage_slice_size*z) + (montage_row_size*(y_offset + y)) + x] = max_value_1 + value_1;
							if((max_value_1 + value_1) > current_max_1)
								current_max_1 = max_value_1 + value_1;
						}
					}
				}
				//Label_RowBorders[n-1]->UnRegister();

				//##################	STITCHING THE ROW_BORDER_TABLE INTO THE MONTAGE_TABLE	  ###################

				// 			std::cout << "stitchin the row border table\n" ;
				for(unsigned long long r=0; r<(unsigned long long)Table_RowBorders[n-1]->GetNumberOfRows(); ++r)
				{
					vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
					for(unsigned long long c=0; c<(unsigned long long)Table_RowBorders[n-1]->GetNumberOfColumns(); ++c)
					{
						if(c == 0)
							model_data1->InsertNextValue(vtkVariant(Table_RowBorders[n-1]->GetValue(r,c).ToInt() + max_value_1));
						else if(c == 2)
							model_data1->InsertNextValue(vtkVariant(Table_RowBorders[n-1]->GetValue(r,c).ToInt() + y_offset));
						else
							model_data1->InsertNextValue(Table_RowBorders[n-1]->GetValue(r,c));
					}
					montageTable->InsertNextRow(model_data1);
				}

				max_value_1 = current_max_1;
			}

			//##################	STITCHING THE ROW INTO THE MONTAGE	  ###################

			// 		std::cout << "stitchin the row\n" ;
			LabelType::Pointer myRow = Label_Rows[n];
			LabelType::PixelType * myRowArray = myRow->GetBufferPointer();
			itk::Size<3> row_size = myRow->GetLargestPossibleRegion().GetSize();
			unsigned long long row_slice_size_1 = row_size[1] * row_size[0];
			unsigned long long row_row_size_1 = row_size[0];	
			unsigned long long y_offset = n*rowDivisor;
			for(unsigned long long z=0; z<row_size[2]; ++z)
			{
				for(unsigned long long y=0; y<row_size[1]; ++y)
				{
					for(unsigned long long x=0; x<row_size[0]; ++x)
					{
						unsigned short value_1 = myRowArray[(row_slice_size_1*z) + (row_row_size_1*y) + (x)];					
						if(value_1 == 0) continue;
						unsigned long long lab_cen_y = Centroids_Rows[n][value_1][1];
						if((n != 0) && (lab_cen_y < 25)) continue;
						if((n != (Label_Rows.size()-1)) && (lab_cen_y >= (row_size[0]-25))) continue;
						montageSegArray[(montage_slice_size*z) + (montage_row_size*(y_offset + y)) + x] = max_value_1 + value_1;
						if((max_value_1 + value_1) > current_max_1)
							current_max_1 = max_value_1 + value_1;
					}
				}
			}
			//Label_Rows[n]->UnRegister();

			//##################	STITCHING THE ROW_TABLE INTO THE MONTAGE_TABLE	  ###################

			// 		std::cout << "stitchin the row table\n" ;
			for(unsigned long long r=0; r<(unsigned long long)Table_Rows[n]->GetNumberOfRows(); ++r)
			{
				vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
				for(unsigned long long c=0; c<(unsigned long long)Table_Rows[n]->GetNumberOfColumns(); ++c)
				{
					if(c == 0)
						model_data1->InsertNextValue(vtkVariant(Table_Rows[n]->GetValue(r,c).ToInt() + max_value_1));
					else if(c == 2)
						model_data1->InsertNextValue(vtkVariant(Table_Rows[n]->GetValue(r,c).ToInt() + y_offset));
					else
						model_data1->InsertNextValue(Table_Rows[n]->GetValue(r,c));
				}
				montageTable->InsertNextRow(model_data1);
			}

			max_value_1 = current_max_1;
		}

		std::cout << "everything done\n";
		Label_Rows.clear();
		Label_RowBorders.clear();

		//#############################################################
		std::cout<<"done !!\n\n";

		labelWriterType::Pointer writer = labelWriterType::New();
		writer->SetFileName(temp + "_label.mhd");
		writer->SetInput(montageSeg);
		writer->Update();

		ftk::SaveTable(temp + "_table.txt", montageTable);


		//ftk::Image::Pointer nucSegImage = NULL;
		//nucSegImage->AppendImageFromData3D( montageSeg->GetBufferPointer(), itk::ImageIOBase::USHORT, 16, size_montage[0], size_montage[1], size_montage[2], "", false);


		//#####################################################################################################################
		//	EXTRACTING SOMA and SOMA CENTROIDS
		//#####################################################################################################################

// 		somaMontage = montageLabelType::New();
// 		itk::Size<3> im_size = montageSeg->GetBufferedRegion().GetSize();
// 		montageLabelType::IndexType start;
// 		start[0] =   0;  // first index on X
// 		start[1] =   0;  // first index on Y    
// 		start[2] =   0;  // first index on Z  
// 		montageLabelType::PointType origin;
// 		origin[0] = 0; 
// 		origin[1] = 0;    
// 		origin[2] = 0;    
// 		somaMontage->SetOrigin( origin );
// 		montageLabelType::RegionType region;
// 		region.SetSize( im_size );
// 		region.SetIndex( start );
// 		somaMontage->SetRegions( region );
// 		somaMontage->Allocate();
// 		somaMontage->FillBuffer(0);
// 		somaMontage->Update();
// 		montageLabelType::PixelType * somaArray = somaMontage->GetBufferPointer();
// 		montageLabelType::PixelType * labelArray = montageSeg->GetBufferPointer();
// 
// 		unsigned long long slice_size = im_size[1] * im_size[0];
// 		unsigned long long row_size = im_size[0];
// 		std::map<unsigned short, int> classMap;
// 		for(int row=0; row<(int)montageTable->GetNumberOfRows(); ++row)
// 		{
// 			classMap[montageTable->GetValue(row,0).ToUnsignedShort()] = montageTable->GetValueByName(row, "prediction_active_mg").ToInt();
// 		}
// 
// 		for(int i=0; i<im_size[2]; ++i)
// 		{
// 			for(int j=0; j<im_size[1]; ++j)
// 			{
// 				for(int k=0; k<im_size[0]; ++k)
// 				{
// 					unsigned long long offset = (i*slice_size)+(j*row_size)+k;
// 					if(classMap[labelArray[offset]] == 1)
// 						somaArray[offset] = labelArray[offset];
// 				}
// 			}
// 		}	
// 
// 		writer = labelWriterType::New();
// 		writer->SetFileName(temp + "_soma_montage.mhd");
// 		writer->SetInput(somaMontage);
// 		writer->Update();

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
		
// 		std::cout<<std::endl<<"Segmentatino Finished";
// 		std::cout<<std::endl<<"Segmentatino Finished";

	}

	if( onlyTrace == 1 )
	{
		
		rawImageType_8bit::Pointer montage_gfp = readAndRescale_16to8(argv[2]);
		

		itk::SizeValueType size_gfp_montage[3];
		size_gfp_montage[0] = montage_gfp->GetLargestPossibleRegion().GetSize()[0];
		size_gfp_montage[1] = montage_gfp->GetLargestPossibleRegion().GetSize()[1];
		size_gfp_montage[2] = montage_gfp->GetLargestPossibleRegion().GetSize()[2];


// 		labelReaderType::Pointer reader1 = labelReaderType::New();
// 		reader1->SetFileName(temp + "_soma_montage.mhd");
// // 		std::cout<<std::endl<<temp + "_soma_montage.mhd";
// 		reader1->Update();

		// 	rescaleFilter->SetInput(reader1->GetOutput());
		// 	rescaleFilter->Update();

// 		somaMontage = reader1->GetOutput();

		somaCentroidsTable = ftk::LoadTable(temp + "_soma_centroids_table.txt");

// 		std::cout<<std::endl<<temp + "_soma_montage.mhd";


		



		//#####################################################################################################################
		//	MULTIPLE NEURON TACER
		//#####################################################################################################################

		std::vector< itk::Index<3> > centroid_list;
		for(int r=0; r<(int)somaCentroidsTable->GetNumberOfRows(); ++r)
		{
			int cx = somaCentroidsTable->GetValue(r, 0).ToInt();
			int cy = somaCentroidsTable->GetValue(r, 1).ToInt();
			int cz = somaCentroidsTable->GetValue(r, 2).ToInt();

			itk::Index<3> cen;
			cen[0] = cx; cen[1] = cy; cen[2] = cz; 
			centroid_list.push_back(cen);
			//std::cout<<cx<<" "<<cy<<" "<<cz;
		}

		std::cout << "Number of cells to be traced : " << centroid_list.size() << "\n";

		//itk::CastImageFilter< rawImageType, gfpImageType >::Pointer caster = itk::CastImageFilter< rawImageType, gfpImageType>::New();
		//caster->SetInput(montage_gfp);
		//caster->Update();
		//gfpImageType::Pointer montage_for_tracing = caster->GetOutput();

		std::string SWCFilename = filePath + "/TracesAndSomas/OnlySWC.xml";
		std::ofstream outfile;
		outfile.open(SWCFilename.c_str());
		outfile << "<?xml\tversion=\"1.0\"\t?>\n";
		outfile << "<Source>\n\n";


		omp_set_nested(0);

		omp_set_max_active_levels(1);


		itk::MultiThreader::SetGlobalDefaultNumberOfThreads(2);


		num_threads = 40;
		omp_set_num_threads(num_threads);


		diceSizeX = 400;
		diceSizeY = 400;
		
		int counterCentro = 0;
	#pragma omp parallel for
		for( long int a=0; a<centroid_list.size(); ++a )
		{


	#pragma omp critical
			{
				counterCentro++;
				std::cout<<std::endl<<"\t\t\t\t asdfasdf ----->>>>> " << counterCentro << ", of " << centroid_list.size();
				std::cout<<centroid_list[a][0]<<" "<<centroid_list[a][1]<<" "<<centroid_list[a][2];
			}


			int x, y, z;
			std::stringstream ssx, ssy, ssz;

			x = centroid_list[a][0];
			y = centroid_list[a][1];
			z = centroid_list[a][2];


			ssx << x; ssy << y; ssz << z;


			std::cout << "Tracing Dice " << ssx.str() << "_" << ssy.str() << "_" << ssz.str() << "\n";

			std::stringstream ssx_off, ssy_off, ssz_off;		
			if(x >= 200)
				ssx_off << x - 200;
			else 
				ssx_off << 0;
			if(y >= 200)
				ssy_off << y - 200;
			else 
				ssy_off << 0;
			if(z >= 100)
				ssz_off << z - 100;
			else 
				ssz_off << 0;

			//#pragma omp barrier

			//########    CROP THE DESIRED DICE FROM THE GFP AND SOMA MONTAGES   ########
			montageLabelType::IndexType start;
			start[0] = ((x - 200)>0) ? (x - 200):0;
			start[1] = ((y - 200)>0) ? (y - 200):0;
			if(size_gfp_montage[2] <= 200)
				start[2] = 0;
			else
				start[2] = ((z - 100)>0) ? (z - 100):0;

			rawImageType::IndexType start2;
			start2[0] = ((x - 200)>0) ? (x - 200):0;
			start2[1] = ((y - 200)>0) ? (y - 200):0;
			if(size_gfp_montage[2] <= 200)
				start2[2] = 0;
			else
				start2[2] = ((z - 100)>0) ? (z - 100):0;

			montageLabelType::SizeType size;
			size[0] = ((x+200)<size_gfp_montage[0]) ? 400 : (200+size_gfp_montage[0]-x-1); 
			size[1] = ((y+200)<size_gfp_montage[1]) ? 400 : (200+size_gfp_montage[1]-y-1);
			if(size_gfp_montage[2] <= 200)
				size[2] = size_gfp_montage[2];
			else
				size[2] = ((z+100)<size_gfp_montage[2]) ? 200 : (100+size_gfp_montage[2]-z-1);

			rawImageType::SizeType size2;
			size2[0] = ((x+200)<size_gfp_montage[0]) ? 400 : (200+size_gfp_montage[0]-x-1);
			size2[1] = ((y+200)<size_gfp_montage[1]) ? 400 : (200+size_gfp_montage[1]-y-1);
			if(size_gfp_montage[2] <= 200)
				size2[2] = size_gfp_montage[2];
			else
				size2[2] = ((z+100)<size_gfp_montage[2]) ? 200 : (100+size_gfp_montage[2]-z-1);




			montageLabelType::RegionType desiredRegion;
			desiredRegion.SetSize(size);
			desiredRegion.SetIndex(start);

			rawImageType::RegionType desiredRegion2;
			desiredRegion2.SetSize(size2);
			desiredRegion2.SetIndex(start2);


			rawROIFilterType::Pointer ROIfilter2 = rawROIFilterType::New();
			ROIfilter2->SetRegionOfInterest(desiredRegion2);
			ROIfilter2->SetInput(montage_gfp);
	#pragma omp critical
			ROIfilter2->Update();
			rawImageType::Pointer img_gfp = ROIfilter2->GetOutput();
			rawDuplicatorType::Pointer rawDuplicator = rawDuplicatorType::New();
			rawDuplicator->SetInputImage(ROIfilter2->GetOutput());
	#pragma omp critical
			rawDuplicator->Update();

			itk::CastImageFilter< rawImageType, gfpImageType >::Pointer caster = itk::CastImageFilter< rawImageType, gfpImageType>::New();
			caster->SetInput(rawDuplicator->GetOutput());
			caster->Update();
			gfpImageType::Pointer img_trace = caster->GetOutput();


			//########    FETCH ALL CENTROIDS THAT FALL WITHIN THE DICE    ########
			itk::Index<3> centroid;
			centroid[0] = ((x - 200)>0) ? 200:x; 
			centroid[1] = ((y - 200)>0) ? 200:y;
			if(size_gfp_montage[2] <= 200)
				centroid[2] = z;
			else
				centroid[2] = ((z - 100)>0) ? 100:z;

			std::vector< itk::Index<3> > soma_Table;      
			for(int ctr =0; ctr<centroid_list.size() ; ++ctr)
			{
				itk::Index<3> cen = centroid_list[ctr];
				//if(abs((double)(cen[0]-x))<=200 && abs((double)(cen[1]-y))<=200 && abs((double)(cen[2]-z))<=rowDivisor )
				if( (cen[0]>=start2[0]) && (cen[0]<(start2[0]+size2[0])) && (cen[1]>=start2[1]) && (cen[1]<(start2[1]+size2[1])) && (cen[2]>=start2[2]) && (cen[2]<(start2[2]+size2[2])) )
				{
					itk::Index<3> centroid2;
					centroid2[0] = centroid[0] + cen[0] - x; //Is there a reason why x and y are flipped?
					centroid2[1] = centroid[1] + cen[1] - y;
					centroid2[2] = centroid[2] + cen[2] - z;
					soma_Table.push_back(centroid2);
				}
			}


			double stopingTime = 10;
			double curScaling = 0.5;
			double rmsThrehold = 0.02;
			SomaExtractor *Somas = new SomaExtractor();
			montageLabelType::Pointer  img_soma = Somas->SegmentSoma(img_trace, soma_Table,stopingTime,curScaling,rmsThrehold);
	
// 			montageLabelType::Pointer img_soma;// = Somas->SegmentSoma(img_trace, soma_Table,curScaling,stopingTime,rmsThrehold);

			//########    RUN TRACING    ########
			MultipleNeuronTracer * MNT = new MultipleNeuronTracer();
			MNT->LoadCurvImage_1(img_trace, 0);
			MNT->ReadStartPoints_1(soma_Table, 0);
			MNT->SetCostThreshold(700);
			MNT->LoadSomaImage_1(img_soma);
			MNT->RunTracing();

			x = min(200, x);
			y = min(200, y);
			if(size_gfp_montage[2] > 200)
				z = min(100, z);

			//MNT->WriteSWCFile(filePath + "/____Trace_" + ssx.str() + "_" + ssy.str() + "_" + ssz.str() + "_ANT.swc", 0);

			vtkSmartPointer< vtkTable > swcTable = MNT->GetSWCTable(0);
			std::string swcFilename = filePath + "/TracesAndSomas/Trace_" + ssx.str() + "_" + ssy.str() + "_" + ssz.str() + "_ANT.swc";
			WriteCenterTrace(swcTable, x, y, z, swcFilename);
			//bool ok = MNT->WriteSWCforDice(filePath + "/Trace_" + ssx.str() + "_" + ssy.str() + "_" + ssz.str() + "_ANT.swc", x, y, z, 0);
			outfile << "\t<File\tFileName=\"Trace_" << ssx.str() << "_" << ssy.str() << "_" << ssz.str() << "_ANT.swc\"\tType=\"Trace\"\ttX=\"" << ssx_off.str() << "\"\ttY=\"" << ssy_off.str() << "\"\ttZ=\"" << ssz_off.str() << "\"/>\n";

			delete MNT;
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



//######	#######################################################################################################
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
	// 	#pragma omp critical
	NucleusSeg->runBinarization();

	std::cout<<std::endl << "\t\t Done binarization:" <<std::flush;



	try 
	{
		// 		std::cout<<"Starting seed detection\n";
		NucleusSeg->runSeedDetection();
		// 	#pragma omp critical
		// 	#pragma omp barrier
		// 	{
		std::cout<<std::endl << "\t\t\t Done Seed Detect: " <<std::flush;
		// 	}
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

	// 	#pragma omp critical
	// 	{
	std::cout<<std::endl << "\t\t\t\t Done Clus: ";
	// 	}
	unsigned short *output_img;

	//if(NucleusSeg->isSegmentationFinEnabled())
	//{
	//NucleusSeg->runAlphaExpansion();
	//output_img=NucleusSeg->getSegImage();
	//} 
	//else

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

	//std::cout<<std::endl<<nucImg->GetLargestPossibleRegion().GetSize()[0];

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
	color[0] = 255; color[1] = 255; color[2] = 0;
	sourceImages->AppendChannelFromData3D( cy5Img->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), tileSize[0], tileSize[1], tileSize[2], "cy5", color, true);

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

	projectDef.Write(fileName);	
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

	std::cout << " done! " << std::endl;

}

void RunEverything(rawImageType::RegionType region, rawImageType::Pointer m_nuc, rawImageType::Pointer m_gfp, rawImageType::Pointer m_cy5, std::vector<LabelType::Pointer> &myLabel, std::vector < vtkSmartPointer< vtkTable > > &myTable, std::map< unsigned int, itk::Index<3> > *myCentroids, const char* segParams, const char* projectDef, int col)
{
	rawImageType::Pointer tile_nuc;
	rawImageType::Pointer tile_gfp;
	rawImageType::Pointer tile_cy5;
	//#pragma omp critical
	//	{
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

	rawROIFilterType::Pointer ROIfilter_tile_cy5 = rawROIFilterType::New();
	ROIfilter_tile_cy5->SetRegionOfInterest(region);
	ROIfilter_tile_cy5->SetInput(m_cy5);
#pragma omp critical
	ROIfilter_tile_cy5->Update();
	rawDuplicatorType::Pointer rawDuplicator_tile_cy5 = rawDuplicatorType::New();
	rawDuplicator_tile_cy5->SetInputImage(ROIfilter_tile_cy5->GetOutput());
#pragma omp critical
	rawDuplicator_tile_cy5->Update();
	tile_cy5 = rawDuplicator_tile_cy5->GetOutput();
	//	}

	myLabel[col] = RunNuclearSegmentation(tile_nuc, segParams);
	//temp_Tile = tileSegResults.first;
	//Table_Tiles[col] = tileSegResults.second;
	myTable[col] = ComputeFeaturesAndAssociations(tile_nuc, tile_gfp, tile_cy5, myLabel[col], projectDef );
	//Centroids_Tiles[col] = GetLabelToCentroidMap(tileSegResults.second);
	*myCentroids = GetLabelToCentroidMap(myTable[col]);

}

rawImageType_8bit::Pointer readAndRescale_16to8(const char* nameInput){

// 		std::cout<<std::endl<<"Reading nuc"<<std::flush;
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


rawImageType_8bit::Pointer readAndRescale_16to8(const char* nameInput, int debugCopy){

// 		std::cout<<std::endl<<"Reading nuc"<<std::flush;
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

	

