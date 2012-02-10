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
#include "vtkTable.h"
#include "ftkUtils.h"
#include "ftkImage.h"

//#include "../NuclearSegmentation/yousef_core/yousef_seg.h"
//#include "../ftkFeatures/ftkLabelImageToFeatures.h"
#include "../NuclearSegmentation/NucleusEditor/ftkProjectProcessor.h"
#include <iostream>
#include <sstream>
#include <time.h>
#include <limits.h>
#include <math.h>
#include <algorithm>

#include "omp.h"

typedef itk::Image<unsigned char,  3> rawImageType;
typedef itk::Image<unsigned int, 3> LabelType;
typedef MultipleNeuronTracer::ImageType3D gfpImageType;
typedef itk::ImageDuplicator< rawImageType >  rawDuplicatorType;
typedef itk::ImageDuplicator< LabelType >  LabelDuplicatorType;
typedef itk::ImageDuplicator< gfpImageType >  gfpDuplicatorType;
typedef itk::RegionOfInterestImageFilter< rawImageType, rawImageType > rawROIFilterType;
typedef itk::RegionOfInterestImageFilter< gfpImageType, gfpImageType > gfpROIFilterType;
typedef itk::RegionOfInterestImageFilter< LabelType, LabelType > LabelROIFilterType;
typedef std::pair< LabelType::Pointer, vtkSmartPointer< vtkTable > > segResultsType;

//segResultsType RunNuclearSegmentation( rawImageType::Pointer, rawImageType::Pointer, const char*);
LabelType::Pointer RunNuclearSegmentation(rawImageType::Pointer, const char* );
vtkSmartPointer<vtkTable> ComputeFeatures(rawImageType::Pointer, rawImageType::Pointer, LabelType::Pointer, const char* );
std::map< unsigned int, itk::Index<3> > GetLabelToCentroidMap( vtkSmartPointer< vtkTable > );
void WriteCenterTrace(vtkSmartPointer< vtkTable >, int, int, int, std::string);

//typedef struct{ double x_d; double y_d; double z_d; } Table;
bool myfunction (itk::SizeValueType i,itk::SizeValueType j) { return (i<j); }

int main(int argc, char* argv[])
{
	omp_set_nested(1);
	int counterTiles = 0;
	if(argc < 4)
	{
		std::cout<<"usage1: darpa_tracer <dapi_montage_file> <gfp_montage_file> <project_definition_file> \n";
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

	typedef itk::ImageFileReader<rawImageType> rawReaderType;
	typedef itk::ImageFileReader<gfpImageType> gfpReaderType;
	typedef itk::ImageFileWriter<LabelType> labelWriterType;

	rawReaderType::Pointer reader_nuc = rawReaderType::New();
	reader_nuc->SetFileName(argv[1]);
	reader_nuc->Update();
	rawImageType::Pointer montage_nuc = reader_nuc->GetOutput();
	itk::SizeValueType size_nuc_montage[3];
	size_nuc_montage[0] = montage_nuc->GetLargestPossibleRegion().GetSize()[0];
	size_nuc_montage[1] = montage_nuc->GetLargestPossibleRegion().GetSize()[1];
	size_nuc_montage[2] = montage_nuc->GetLargestPossibleRegion().GetSize()[2];

	rawReaderType::Pointer reader_gfp = rawReaderType::New();
	reader_gfp->SetFileName(argv[2]);
	reader_gfp->Update();
	rawImageType::Pointer montage_gfp = reader_gfp->GetOutput();
	itk::SizeValueType size_gfp_montage[3];
	size_gfp_montage[0] = montage_gfp->GetLargestPossibleRegion().GetSize()[0];
	size_gfp_montage[1] = montage_gfp->GetLargestPossibleRegion().GetSize()[1];
	size_gfp_montage[2] = montage_gfp->GetLargestPossibleRegion().GetSize()[2];

	//ftk::Image::Pointer sourceImages = NULL;
	//std::vector< unsigned char > color;
	//color.push_back(255); color.push_back(0); color.push_back(0);
	//sourceImages->AppendChannelFromData3D( montage_nuc->GetBufferPointer(), itk::ImageIOBase::UCHAR, 8, size_nuc_montage[0], size_nuc_montage[1], size_nuc_montage[2], std::string(argv[1]), color, false);
	//color[0] = 0; color[1] = 255; color[2] = 0;
	//sourceImages->AppendChannelFromData3D( montage_gfp->GetBufferPointer(), itk::ImageIOBase::UCHAR, 8, size_nuc_montage[0], size_nuc_montage[1], size_nuc_montage[2], std::string(argv[2]), color, false);


	//#####################################################################################################################
	//	NUCLEAR SEGMENT THE MONTAGE TILE BY TILE AND STITCH THE RESULT TILES TOGETHER
	//#####################################################################################################################
	
	unsigned long long num_rows = (unsigned long long)ceil((double)size_nuc_montage[1]/1000);
	unsigned long long num_cols = (unsigned long long)ceil((double)size_nuc_montage[0]/1000);
	std::cout << "Image divided into " << num_rows << " rows and " << num_cols << " columns\n";

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

	//##################	SEGMENTING EACH ROW IN THE MONTAGE	  ###################
	//#pragma omp parallel for
	for(int row=0; row<num_rows; ++row)
	{

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
		#pragma omp parallel for num_threads(16)
		for(int col=0; col<num_cols; ++col)
		{

			//##################	EXTRACT A TILE AND START SEGMENTING	  ###################

			std::cout<<"Extracting Tile " << col*1000 << "_" << row*1000 << "\n";			
			rawImageType::IndexType start_tile;
			start_tile[0] = col*1000;
			start_tile[1] = row*1000;
			start_tile[2] = 0;

			rawImageType::SizeType size_tile;
			size_tile[0] = (((col+1)*1000) < size_nuc_montage[0]) ? 1000 : (size_nuc_montage[0] - (col*1000));
			size_tile[1] = (((row+1)*1000) < size_nuc_montage[1]) ? 1000 : (size_nuc_montage[1] - (row*1000));
			size_tile[2] = size_nuc_montage[2];

			rawImageType::RegionType region_tile;
			region_tile.SetSize(size_tile);
			region_tile.SetIndex(start_tile);

			rawROIFilterType::Pointer ROIfilter_tile_nuc = rawROIFilterType::New();
			ROIfilter_tile_nuc->SetRegionOfInterest(region_tile);
			ROIfilter_tile_nuc->SetInput(montage_nuc);
			#pragma omp critical
				ROIfilter_tile_nuc->Update();
			rawDuplicatorType::Pointer rawDuplicator_tile_nuc = rawDuplicatorType::New();
			rawDuplicator_tile_nuc->SetInputImage(ROIfilter_tile_nuc->GetOutput());
			#pragma omp critical
				rawDuplicator_tile_nuc->Update();
			rawImageType::Pointer tile_nuc = rawDuplicator_tile_nuc->GetOutput();

			rawROIFilterType::Pointer ROIfilter_tile_gfp = rawROIFilterType::New();
			ROIfilter_tile_gfp->SetRegionOfInterest(region_tile);
			ROIfilter_tile_gfp->SetInput(montage_gfp);
			#pragma omp critical
				ROIfilter_tile_gfp->Update();
			rawDuplicatorType::Pointer rawDuplicator_tile_gfp = rawDuplicatorType::New();
			rawDuplicator_tile_gfp->SetInputImage(ROIfilter_tile_gfp->GetOutput());
			#pragma omp critical
				rawDuplicator_tile_gfp->Update();
			rawImageType::Pointer tile_gfp = rawDuplicator_tile_gfp->GetOutput();

			//##################	ALLOCATING MEMORY AND REGISTERING FOR A SEGMENTED TILE	  ###################
			
			LabelType::Pointer temp_Tile = LabelType::New();
			LabelType::PointType origin_Tile;
			origin_Tile[0] = 0; 
			origin_Tile[1] = 0;
			origin_Tile[2] = 0;
			temp_Tile->SetOrigin( origin_Tile );
			LabelType::IndexType start_Tile;
			start_Tile[0] = 0;
			start_Tile[1] = 0;
			start_Tile[2] = 0;
			LabelType::SizeType size_Tile;
			size_Tile[0] = size_tile[0];
			size_Tile[1] = size_tile[1];
			size_Tile[2] = size_tile[2];
			LabelType::RegionType region_Tile;
			region_Tile.SetSize ( size_Tile  );
			region_Tile.SetIndex( start_Tile );
			temp_Tile->SetRegions( region_Tile );
			temp_Tile->Allocate();
			temp_Tile->FillBuffer(0);
			temp_Tile->Update();
			//temp_Tile->Register();
			
			std::cout<<"Starting segmentation for Tile " << col*1000 << "_" << row*1000 << "\n";
			//segResultsType tileSegResults;
			//tileSegResults = RunNuclearSegmentation(tile_nuc, tile_gfp, argv[3]);
			temp_Tile = RunNuclearSegmentation(tile_nuc, argv[3]);
			std::cout<<"Done with nucleus segmentation for Tile " << col*1000 << "_" << row*1000 << "\n\n";
			//temp_Tile = tileSegResults.first;
			Label_Tiles[col] = temp_Tile;
			//Table_Tiles[col] = tileSegResults.second;
			Table_Tiles[col] = ComputeFeatures(tile_nuc, tile_gfp, temp_Tile, argv[4] );
			//Centroids_Tiles[col] = GetLabelToCentroidMap(tileSegResults.second);
			Centroids_Tiles[col] = GetLabelToCentroidMap(Table_Tiles[col]);

			
			//##################	EXTRACT A TILE_BORDER AND START SEGMENTING	  ###################
			
			if(col == 0) continue;

			std::cout<<"Extracting tile border X = " << col*1000 << ", Row = " << row << "\n";
			rawImageType::IndexType start_tileBorder;
			start_tileBorder[0] = (col*1000) - (2*25);
			start_tileBorder[1] = (row*1000);
			start_tileBorder[2] = 0;

			rawImageType::SizeType size_tileBorder;
			size_tileBorder[0] = (4*25);
			size_tileBorder[1] = (((row+1)*1000) < size_nuc_montage[1]) ? 1000 : (size_nuc_montage[1] - (row*1000));
			size_tileBorder[2] = size_nuc_montage[2];

			rawImageType::RegionType region_tileBorder;
			region_tileBorder.SetSize(size_tileBorder);
			region_tileBorder.SetIndex(start_tileBorder);

			rawROIFilterType::Pointer ROIfilter_tileBorder_nuc = rawROIFilterType::New();
			ROIfilter_tileBorder_nuc->SetRegionOfInterest(region_tileBorder);
			ROIfilter_tileBorder_nuc->SetInput(montage_nuc);
			#pragma omp critical
				ROIfilter_tileBorder_nuc->Update();
			rawDuplicatorType::Pointer rawDuplicator_tileBorder_nuc = rawDuplicatorType::New();
			rawDuplicator_tileBorder_nuc->SetInputImage(ROIfilter_tileBorder_nuc->GetOutput());
			#pragma omp critical
				rawDuplicator_tileBorder_nuc->Update();
			rawImageType::Pointer tileBorder_nuc = rawDuplicator_tileBorder_nuc->GetOutput();

			rawROIFilterType::Pointer ROIfilter_tileBorder_gfp = rawROIFilterType::New();
			ROIfilter_tileBorder_gfp->SetRegionOfInterest(region_tileBorder);
			ROIfilter_tileBorder_gfp->SetInput(montage_gfp);
			#pragma omp critical
				ROIfilter_tileBorder_gfp->Update();
			rawDuplicatorType::Pointer rawDuplicator_tileBorder_gfp = rawDuplicatorType::New();
			rawDuplicator_tileBorder_gfp->SetInputImage(ROIfilter_tileBorder_gfp->GetOutput());
			#pragma omp critical
				rawDuplicator_tileBorder_gfp->Update();
			rawImageType::Pointer tileBorder_gfp = rawDuplicator_tileBorder_gfp->GetOutput();

			//##################	ALLOCATING MEMORY AND REGISTERING FOR A SEGMENTED TILE_BORDER	  ###################
			
			LabelType::Pointer temp_TileBorder = LabelType::New();
			LabelType::PointType origin_TileBorder;
			origin_TileBorder[0] = 0; 
			origin_TileBorder[1] = 0;
			origin_TileBorder[2] = 0;
			temp_TileBorder->SetOrigin( origin_TileBorder );
			LabelType::IndexType start_TileBorder;
			start_TileBorder[0] = 0;
			start_TileBorder[1] = 0;
			start_TileBorder[2] = 0;
			LabelType::SizeType size_TileBorder;
			size_TileBorder[0] = size_tileBorder[0];
			size_TileBorder[1] = size_tileBorder[1];
			size_TileBorder[2] = size_tileBorder[2];
			LabelType::RegionType region_TileBorder;
			region_TileBorder.SetSize ( size_TileBorder  );
			region_TileBorder.SetIndex( start_TileBorder );
			temp_TileBorder->SetRegions( region_TileBorder );
			temp_TileBorder->Allocate();
			temp_TileBorder->FillBuffer(0);
			temp_TileBorder->Update();
			//temp_TileBorder->Register();

			std::cout<<"Starting segmentation for tile border X = " << col*1000 << ", Row = " << row << "\n";
			//segResultsType tileBorderSegResults;
			//tileBorderSegResults = RunNuclearSegmentation(tileBorder_nuc, tileBorder_gfp, argv[3]);
			temp_TileBorder = RunNuclearSegmentation(tileBorder_nuc, argv[3]);
			std::cout<<"Done with nucleus segmentation for tile border X = " << col*1000 << ", Row = " << row << "\n\n";
			//temp_TileBorder = tileBorderSegResults.first;
			Label_TileBorders[col-1] = temp_TileBorder;
			//Table_TileBorders[col-1] = tileBorderSegResults.second;
			Table_TileBorders[col-1] = ComputeFeatures(tileBorder_nuc, tileBorder_gfp, temp_TileBorder, argv[4] );
			//Centroids_TileBorders[col-1] = GetLabelToCentroidMap(tileBorderSegResults.second);
			Centroids_TileBorders[col-1] = GetLabelToCentroidMap(Table_TileBorders[col-1]);

			#pragma omp critical
			{
				counterTiles++;
				std::cout << std::endl << "\t\t\t ->-> Tile " << counterTiles << " of " << num_cols*num_rows;
				std::cout << std::endl;
			}
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
		size_row[1] = (((row+1)*1000) < size_nuc_montage[1]) ? 1000 : (size_nuc_montage[1] - (row*1000));
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
			if(m != 0)
			{

				//##################	STITCHING THE TILE_BORDER INTO THE ROW	  ###################

				LabelType::Pointer myTileBorder = Label_TileBorders[m-1];
				LabelType::PixelType * myTileBorderArray = myTileBorder->GetBufferPointer();
				itk::Size<3> tileBorder_size = myTileBorder->GetLargestPossibleRegion().GetSize();
				unsigned long long tileBorder_slice_size = tileBorder_size[1] * tileBorder_size[0];
				unsigned long long tileBorder_row_size = tileBorder_size[0];	
				unsigned long long x_offset = (m*1000) - (2*25);
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
			unsigned long long x_offset = m*1000;
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
		Table_Tiles.clear();
		Label_TileBorders.clear();

		//#############################################################
		std::cout<<"done !!\n\n";
		Label_Rows[row] = rowSeg;
		Table_Rows[row] = rowTable;
		Centroids_Rows[row] = GetLabelToCentroidMap(rowTable);

		//##################	EXTRACT A ROW_BORDER AND START SEGMENTING	  ###################

		if(row == 0) continue;

		std::cout<<"Extracting row border Y = " << row*1000 << "\n";
		rawImageType::IndexType start_rowBorder;
		start_rowBorder[0] = 0;
		start_rowBorder[1] = (row*1000) - (2*25);
		start_rowBorder[2] = 0;

		rawImageType::SizeType size_rowBorder;
		size_rowBorder[0] = size_nuc_montage[0];
		size_rowBorder[1] = (4*25);
		size_rowBorder[2] = size_nuc_montage[2];

		rawImageType::RegionType region_rowBorder;
		region_rowBorder.SetSize(size_rowBorder);
		region_rowBorder.SetIndex(start_rowBorder);

		rawROIFilterType::Pointer ROIfilter_rowBorder_nuc = rawROIFilterType::New();
		ROIfilter_rowBorder_nuc->SetRegionOfInterest(region_rowBorder);
		ROIfilter_rowBorder_nuc->SetInput(montage_nuc);
		#pragma omp critical
			ROIfilter_rowBorder_nuc->Update();
		rawDuplicatorType::Pointer rawDuplicator_rowBorder_nuc = rawDuplicatorType::New();
		rawDuplicator_rowBorder_nuc->SetInputImage(ROIfilter_rowBorder_nuc->GetOutput());
		#pragma omp critical
			rawDuplicator_rowBorder_nuc->Update();
		rawImageType::Pointer rowBorder_nuc = rawDuplicator_rowBorder_nuc->GetOutput();

		rawROIFilterType::Pointer ROIfilter_rowBorder_gfp = rawROIFilterType::New();
		ROIfilter_rowBorder_gfp->SetRegionOfInterest(region_rowBorder);
		ROIfilter_rowBorder_gfp->SetInput(montage_gfp);
		#pragma omp critical
			ROIfilter_rowBorder_gfp->Update();
		rawDuplicatorType::Pointer rawDuplicator_rowBorder_gfp = rawDuplicatorType::New();
		rawDuplicator_rowBorder_gfp->SetInputImage(ROIfilter_rowBorder_gfp->GetOutput());
		#pragma omp critical
			rawDuplicator_rowBorder_gfp->Update();
		rawImageType::Pointer rowBorder_gfp = rawDuplicator_rowBorder_gfp->GetOutput();

		//##################	ALLOCATING MEMORY AND REGISTERING FOR A SEGMENTED ROW_BORDER	  ###################

		LabelType::Pointer temp_RowBorder = LabelType::New();
		LabelType::PointType origin_RowBorder;
		origin_RowBorder[0] = 0; 
		origin_RowBorder[1] = 0;
		origin_RowBorder[2] = 0;
		temp_RowBorder->SetOrigin( origin_RowBorder );
		LabelType::IndexType start_RowBorder;
		start_RowBorder[0] = 0;
		start_RowBorder[1] = 0;
		start_RowBorder[2] = 0;
		LabelType::SizeType size_RowBorder;
		size_RowBorder[0] = size_rowBorder[0];
		size_RowBorder[1] = size_rowBorder[1];
		size_RowBorder[2] = size_rowBorder[2];
		LabelType::RegionType region_RowBorder;
		region_RowBorder.SetSize ( size_RowBorder  );
		region_RowBorder.SetIndex( start_RowBorder );
		temp_RowBorder->SetRegions( region_RowBorder );
		temp_RowBorder->Allocate();
		temp_RowBorder->FillBuffer(0);
		temp_RowBorder->Update();
		//temp_RowBorder->Register();

		std::cout<<"Starting segmentation for row border Y = " << row*1000 << "\n";
		//segResultsType rowBorderSegResults;
		//rowBorderSegResults = RunNuclearSegmentation(rowBorder_nuc, rowBorder_gfp, argv[3]);
		temp_RowBorder = RunNuclearSegmentation(rowBorder_nuc, argv[3]);
		std::cout<<"Done with nucleus segmentation for border Y = " << row*1000 << "\n\n";
		//temp_RowBorder = rowBorderSegResults.first;
		Label_RowBorders[row-1] = temp_RowBorder;
		//Table_RowBorders[row-1] = rowBorderSegResults.second;
		Table_RowBorders[row-1] = ComputeFeatures(rowBorder_nuc, rowBorder_gfp, temp_RowBorder, argv[4] );
		//Centroids_RowBorders[row-1] = GetLabelToCentroidMap(rowBorderSegResults.second);
		Centroids_RowBorders[row-1] = GetLabelToCentroidMap(Table_RowBorders[row-1]);

	}

	std::cout<<"Stitching all rows ...";

	//##################	ALLOCATING MEMORY FOR THE SEGMENTED MONTAGE	  ###################

	LabelType::Pointer montageSeg = LabelType::New();
	LabelType::PointType origin_montage;
	origin_montage[0] = 0; 
	origin_montage[1] = 0;
	origin_montage[2] = 0;
	montageSeg->SetOrigin( origin_montage );
	LabelType::IndexType start_montage;
	start_montage[0] = 0;
	start_montage[1] = 0;
	start_montage[2] = 0;
	LabelType::SizeType size_montage;
	size_montage[0] = size_nuc_montage[0];
	size_montage[1] = size_nuc_montage[1];
	size_montage[2] = size_nuc_montage[2];
	LabelType::RegionType region_montage;
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

	LabelType::PixelType * montageSegArray = montageSeg->GetBufferPointer();
	unsigned long long montage_slice_size = size_montage[1] * size_montage[0];
	unsigned long long montage_row_size = size_montage[0];	

	unsigned short max_value_1 = 0, current_max_1 = 0;
	for(unsigned long long n=0; n<(int)Label_Rows.size(); ++n)
	{
		std::cout << "stitching row " << n << "\n";
		if(n != 0)
		{

			//##################	STITCHING THE ROW_BORDER INTO THE MONTAGE	  ###################

			std::cout << "stitchin the row border\n" ;
			LabelType::Pointer myRowBorder = Label_RowBorders[n-1];
			LabelType::PixelType * myRowBorderArray = myRowBorder->GetBufferPointer();
			itk::Size<3> rowBorder_size = myRowBorder->GetLargestPossibleRegion().GetSize();
			unsigned long long rowBorder_slice_size = rowBorder_size[1] * rowBorder_size[0];
			unsigned long long rowBorder_row_size = rowBorder_size[0];	
			unsigned long long y_offset = (n*1000) - (2*25);
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

			std::cout << "stitchin the row border table\n" ;
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

		std::cout << "stitchin the row\n" ;
		LabelType::Pointer myRow = Label_Rows[n];
		LabelType::PixelType * myRowArray = myRow->GetBufferPointer();
		itk::Size<3> row_size = myRow->GetLargestPossibleRegion().GetSize();
		unsigned long long row_slice_size_1 = row_size[1] * row_size[0];
		unsigned long long row_row_size_1 = row_size[0];	
		unsigned long long y_offset = n*1000;
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

		std::cout << "stitchin the row table\n" ;
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
			
	std::string temp = nucFileName;
	string::iterator it;
	it = temp.end() - 4;
	temp.erase(it, it+4);

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

	LabelType::Pointer somaMontage = LabelType::New();
	itk::Size<3> im_size = montageSeg->GetBufferedRegion().GetSize();
	LabelType::IndexType start;
	start[0] =   0;  // first index on X
	start[1] =   0;  // first index on Y    
	start[2] =   0;  // first index on Z  
	LabelType::PointType origin;
	origin[0] = 0; 
	origin[1] = 0;    
	origin[2] = 0;    
	somaMontage->SetOrigin( origin );
	LabelType::RegionType region;
	region.SetSize( im_size );
	region.SetIndex( start );
	somaMontage->SetRegions( region );
	somaMontage->Allocate();
	somaMontage->FillBuffer(0);
	somaMontage->Update();
	LabelType::PixelType * somaArray = somaMontage->GetBufferPointer();
	LabelType::PixelType * labelArray = montageSeg->GetBufferPointer();

	unsigned long long slice_size = im_size[1] * im_size[0];
	unsigned long long row_size = im_size[0];
	std::map<unsigned short, int> classMap;
	for(int row=0; row<(int)montageTable->GetNumberOfRows(); ++row)
	{
		classMap[montageTable->GetValue(row,0).ToUnsignedShort()] = montageTable->GetValueByName(row, "prediction_active_mg").ToInt();
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
	}

	std::cout << "Number of cells to be traced : " << centroid_list.size() << "\n";

	//itk::CastImageFilter< rawImageType, gfpImageType >::Pointer caster = itk::CastImageFilter< rawImageType, gfpImageType>::New();
	//caster->SetInput(montage_gfp);
	//caster->Update();
	//gfpImageType::Pointer montage_for_tracing = caster->GetOutput();

	std::string SWCFilename = filePath + "/OnlySWC.xml";
	std::ofstream outfile;
	outfile.open(SWCFilename.c_str());
	outfile << "<?xml\tversion=\"1.0\"\t?>\n";
	outfile << "<Source>\n\n";

	#pragma omp parallel for
	for( long int a=0; a<centroid_list.size(); ++a )
	{
		int x, y, z;
		x = centroid_list[a][0];
		y = centroid_list[a][1];
		z = centroid_list[a][2];

		std::stringstream ssx, ssy, ssz;
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

		//########    CROP THE DESIRED DICE FROM THE GFP AND SOMA MONTAGES   ########
		LabelType::IndexType start;
		start[0] = ((x - 200)>0) ? (x - 200):0;
		start[1] = ((y - 200)>0) ? (y - 200):0;
		if(size_nuc_montage[2] <= 200)
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

		LabelType::SizeType size;
		size[0] = ((x+200)<size_nuc_montage[0]) ? 400 : (200+size_nuc_montage[0]-x-1); 
		size[1] = ((y+200)<size_nuc_montage[1]) ? 400 : (200+size_nuc_montage[1]-y-1);
		if(size_nuc_montage[2] <= 200)
			size[2] = size_nuc_montage[2];
		else
			size[2] = ((z+100)<size_nuc_montage[2]) ? 200 : (100+size_nuc_montage[2]-z-1);

		rawImageType::SizeType size2;
		size2[0] = ((x+200)<size_gfp_montage[0]) ? 400 : (200+size_gfp_montage[0]-x-1);
		size2[1] = ((y+200)<size_gfp_montage[1]) ? 400 : (200+size_gfp_montage[1]-y-1);
		if(size_gfp_montage[2] <= 200)
			size2[2] = size_gfp_montage[2];
		else
			size2[2] = ((z+100)<size_gfp_montage[2]) ? 200 : (100+size_gfp_montage[2]-z-1);

		LabelType::RegionType desiredRegion;
		desiredRegion.SetSize(size);
		desiredRegion.SetIndex(start);

		rawImageType::RegionType desiredRegion2;
		desiredRegion2.SetSize(size2);
		desiredRegion2.SetIndex(start2);

		LabelROIFilterType::Pointer ROIfilter3 = LabelROIFilterType::New();
		ROIfilter3->SetRegionOfInterest(desiredRegion);
		ROIfilter3->SetInput(somaMontage);
		#pragma omp critical
			ROIfilter3->Update();
		LabelDuplicatorType::Pointer LabelDuplicator = LabelDuplicatorType::New();
		LabelDuplicator->SetInputImage(ROIfilter3->GetOutput());
		#pragma omp critical
			LabelDuplicator->Update();
		LabelType::Pointer img_soma = LabelDuplicator->GetOutput();

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
			//if(abs((double)(cen[0]-x))<=200 && abs((double)(cen[1]-y))<=200 && abs((double)(cen[2]-z))<=100 )
			if( (cen[0]>=start2[0]) && (cen[0]<(start2[0]+size2[0])) && (cen[1]>=start2[1]) && (cen[1]<(start2[1]+size2[1])) && (cen[2]>=start2[2]) && (cen[2]<(start2[2]+size2[2])) )
			{
				itk::Index<3> centroid2;
				centroid2[0] = centroid[0] + cen[0] - x; //Is there a reason why x and y are flipped?
				centroid2[1] = centroid[1] + cen[1] - y;
				centroid2[2] = centroid[2] + cen[2] - z;
				soma_Table.push_back(centroid2);
			}
		}


		//########    RUN TRACING    ########
		MultipleNeuronTracer * MNT = new MultipleNeuronTracer();
		MNT->LoadCurvImage_1(img_trace, 0);
		MNT->ReadStartPoints_1(soma_Table, 0);
		MNT->SetCostThreshold(1000);
		MNT->LoadSomaImage_1(img_soma);
		MNT->RunTracing();
		
		x = min(200, x);
		y = min(200, y);
		if(size_gfp_montage[2] > 200)
			z = min(100, z);

		//MNT->WriteSWCFile(filePath + "/____Trace_" + ssx.str() + "_" + ssy.str() + "_" + ssz.str() + "_ANT.swc", 0);
		
		vtkSmartPointer< vtkTable > swcTable = MNT->GetSWCTable(0);
		std::string swcFilename = filePath + "/Trace_" + ssx.str() + "_" + ssy.str() + "_" + ssz.str() + "_ANT.swc";
		WriteCenterTrace(swcTable, x, y, z, swcFilename);
		//bool ok = MNT->WriteSWCforDice(filePath + "/Trace_" + ssx.str() + "_" + ssy.str() + "_" + ssz.str() + "_ANT.swc", x, y, z, 0);
		outfile << "\t<File\tFileName=\"Trace_" << ssx.str() << "_" << ssy.str() << "_" << ssz.str() << "_ANT.swc\"\tType=\"Trace\"\ttX=\"" << ssx_off.str() << "\"\ttY=\"" << ssy_off.str() << "\"\ttZ=\"" << ssz_off.str() << "\"/>\n";

		delete MNT;
	}

	outfile << "\n</Source>";
	outfile.close();


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

	memset(in_Image/*destination*/,0/*value*/,size[0]*size[1]*size[2]*sizeof(unsigned char)/*num bytes to move*/);
	typedef itk::ImageRegionConstIterator< rawImageType > ConstIteratorType;
	ConstIteratorType pix_buf( inpImg, inpImg->GetRequestedRegion() );
	itk::SizeValueType ind=0;
	for ( pix_buf.GoToBegin(); !pix_buf.IsAtEnd(); ++pix_buf, ++ind )
	{
		in_Image[ind]=(pix_buf.Get());
	}
	yousef_nucleus_seg *NucleusSeg = new yousef_nucleus_seg();
	NucleusSeg->readParametersFromFile(fileName);
	NucleusSeg->setDataImage(in_Image,size[0],size[1],size[2],"/home/gramak/Desktop/11_10_tracing/montage_kt11306_w410DAPIdsu.mhd");
	NucleusSeg->runBinarization();

	try 
	{
		std::cout<<"Starting seed detection\n";
		NucleusSeg->runSeedDetection();
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
		unsigned short val = (unsigned short)output_img[i];
		iterator1.Set(val);
		++iterator1;
	}

	delete NucleusSeg;

	return image;

}


//#############################################################################################################
//  COMPUTE FEATURES
//#############################################################################################################

vtkSmartPointer<vtkTable> ComputeFeatures(rawImageType::Pointer nucImg, rawImageType::Pointer gfpImg, LabelType::Pointer labImg, const char* fileName)
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

	ftk::Image::Pointer sourceImages = ftk::Image::New();
	std::vector< unsigned char > color;
	color.push_back(255); color.push_back(0); color.push_back(0);
	sourceImages->AppendChannelFromData3D( nucImg->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), tileSize[0], tileSize[1], tileSize[2], "dapi", color, true);
	color[0] = 0; color[1] = 255; color[2] = 0;
	sourceImages->AppendChannelFromData3D( gfpImg->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), tileSize[0], tileSize[1], tileSize[2], "gfp", color, true);

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
		if( (swcNodes->GetValue(i,2).ToInt() != x) || (swcNodes->GetValue(i,3).ToInt() != y) || (swcNodes->GetValue(i,4).ToInt() != z) )//if( ((int)swcNodes[i][2] + x > 1800) && ((int)swcNodes[i][2] + x < 3600) && ((int)swcNodes[i][3] + y < 8000) )
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

