#include "itkImageFileReader.h"
#include "itkImage.h"
#include "boost/tokenizer.hpp"
#include <fstream>
#include "vul/vul_file.h"
#include "iostream"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "../Tracing/MultipleNeuronTracer/MultipleNeuronTracer.h"
#include "vtkTable.h"
#include "ftkUtils.h"

#include "../NuclearSegmentation/yousef_core/yousef_seg.h"
#include <iostream>
#include <sstream>
#include <time.h>
#include <limits.h>
#include <math.h>
#include <algorithm>

#include "omp.h"

typedef itk::Image<unsigned char,  3> nucImageType;
typedef itk::Image<unsigned int, 3> LabelType;
typedef MultipleNeuronTracer::ImageType3D gfpImageType;

//typedef struct{ double x_d; double y_d; double z_d; } centroids;
bool myfunction (uint64_t i,uint64_t j) { return (i<j); }

int main(int argc, char* argv[])
{
	//if(argc < 5){
	//	std::cout<<"Usage1: darpa_tracer <Global_Centroids_List> <DAPI_Montage_File> <GFP_Montage_File> <Seg_Params_File> \n";
	//	return 0;
	//}
	//Input: Somas_file, nucleus_montage, trace_channel_montage, nuc_seg_parameters
	//std::fstream soma_centroid_file;
	//soma_centroid_file.open(argv[1], std::fstream::in);
	//std::vector<centroids> centroid_list;
	//while (soma_centroid_file.good()){
	//	centroids new_centroid;
	//	soma_centroid_file >> new_centroid.x_d >> new_centroid.y_d >> new_centroid.z_d;
	//	centroid_list.push_back( new_centroid );
	//}
	//soma_centroid_file.close();

	std::vector< itk::Index<3> > centroid_list;
	std::vector< itk::Index<3> > all_centroid_list;
	//vtkSmartPointer<vtkTable> global_centroids = ftk::LoadTable(argv[1]);
	vtkSmartPointer<vtkTable> global_centroids = ftk::LoadTable("D:\\11_10_selective_tracing\\Global_centroids.txt");
	for(int r=0; r<(int)global_centroids->GetNumberOfRows(); ++r)
	{
		int cx = global_centroids->GetValue(r, 0).ToInt();
		int cy = global_centroids->GetValue(r, 1).ToInt();
		int cz = global_centroids->GetValue(r, 2).ToInt();

		itk::Index<3> all_cen;
		all_cen[0] = cx; all_cen[1] = cy; all_cen[2] = cz; 
		all_centroid_list.push_back(all_cen);
		if( (fmod((double)cx,1050)>1000) || (fmod((double)cx,1050)<50) || (fmod((double)cy,1050)>1000) || (fmod((double)cy,1050)<50))
		{			
			itk::Index<3> cen;
			cen[0] = cx; cen[1] = cy; cen[2] = cz; 
			centroid_list.push_back(cen);			

		}

	}

	std::cout << "Number of cells to be retraced : " << centroid_list.size() << "\n";

	typedef itk::ImageFileReader<nucImageType> nucReaderType;
	typedef itk::ImageFileReader<gfpImageType> gfpReaderType;
	//nucReaderType::Pointer reader_nuc   = nucReaderType::New();
	//reader_nuc->SetFileName(argv[2]);
	//reader_nuc->Update();
	//nucImageType::Pointer img_nuc   = reader_nuc->GetOutput();
	//gfpReaderType::Pointer reader_trace = gfpReaderType::New();
	//reader_trace->SetFileName(argv[3]);
	//reader_trace->Update();
	//gfpImageType::Pointer img_trace = reader_trace->GetOutput();
	//uint64_t size_trace[3],size_nuc[3];
	//size_nuc[0] = img_nuc->GetLargestPossibleRegion().GetSize()[0];
	//size_nuc[1] = img_nuc->GetLargestPossibleRegion().GetSize()[1];
	//size_nuc[2] = img_nuc->GetLargestPossibleRegion().GetSize()[2];
	//size_trace[0] = img_trace->GetLargestPossibleRegion().GetSize()[0];
	//size_trace[1] = img_trace->GetLargestPossibleRegion().GetSize()[1];
	//size_trace[2] = img_trace->GetLargestPossibleRegion().GetSize()[2];
	//if( size_trace[0]!=size_nuc[0] || size_trace[1]!=size_nuc[1] || size_trace[2]!=size_nuc[2] ) return EXIT_FAILURE;

//	//long int a=0;
//	//while( a<centroid_list.size() ){
//	#pragma omp parallel for num_threads(12)
//	for( long int a=0; a<centroid_list.size(); ++a ){
//	//for( long int a=0; a<5; ++a ){
//		int x, y, z;
//		x = centroid_list[a][0];
//		y = centroid_list[a][1];
//		z = centroid_list[a][2];
//
//		nucImageType::IndexType start;
//		start[0] = ((x - 200)>0) ? (x - 200):0;
//		start[1] = ((y - 200)>0) ? (y - 200):0;
//		start[2] = ((z - 100) >0) ? (z - 100) :0;
//
//		gfpImageType::IndexType start2;
//		start2[0] = ((x - 200)>0) ? (x - 200):0;
//		start2[1] = ((y - 200)>0) ? (y - 200):0;
//		start2[2] = ((z - 100) >0) ? (z - 100) :0;
//
//		nucImageType::SizeType size;
//		size[0] = ((x+200)<size_nuc[0]) ? 400 : (200+size_trace[0]-x-1); 
//		size[1] = ((y+200)<size_nuc[1]) ? 400 : (200+size_trace[1]-y-1);
//		size[2] = ((z+100)<size_nuc[2]) ? 200 : (100+size_trace[2]-z-1);
//
//		gfpImageType::SizeType size2;
//		size2[0] = ((x+200)<size_trace[0]) ? 400 : (200+size_trace[0]-x-1);
//		size2[1] = ((y+200)<size_trace[1]) ? 400 : (200+size_trace[1]-y-1);
//		size2[2] = ((z+100)<size_trace[2]) ? 200 : (100+size_trace[2]-z-1);
//
//		std::ostringstream output_filename_stream;
//
//		nucImageType::RegionType desiredRegion;
//		desiredRegion.SetSize(size);
//		desiredRegion.SetIndex(start);
//
//		gfpImageType::RegionType desiredRegion2;
//		desiredRegion2.SetSize(size2);
//		desiredRegion2.SetIndex(start2);
//
//		typedef itk::RegionOfInterestImageFilter< nucImageType, nucImageType > ROIFilterType;
//		typedef itk::RegionOfInterestImageFilter< gfpImageType, gfpImageType > ROIFilterType2;
//		ROIFilterType::Pointer ROIfilter = ROIFilterType::New();
//		ROIfilter->SetRegionOfInterest(desiredRegion);
//		ROIfilter->SetInput(img_nuc);
//		ROIfilter->Update();
//		nucImageType::Pointer img = ROIfilter->GetOutput();
//
//		ROIFilterType2::Pointer ROIfilter2 = ROIFilterType2::New();
//		ROIfilter2->SetRegionOfInterest(desiredRegion2);
//		ROIfilter2->SetInput(img_trace);
//		ROIfilter2->Update();
//		gfpImageType::Pointer img_tr = ROIfilter2->GetOutput();
//
//		//itk::CastImageFilter<gfpImageType, nucImageType>::Pointer caster = itk::CastImageFilter<gfpImageType, nucImageType>::New();
//		//caster->SetInput(img_tr);
//		//itk::ImageFileWriter<nucImageType>::Pointer writer2 = itk::ImageFileWriter<nucImageType>::New();
//		itk::ImageFileWriter<gfpImageType>::Pointer writer2 = itk::ImageFileWriter<gfpImageType>::New();
//		std::stringstream ss;
//		ss << "gfp_image" << x << "_" << y << "_" << z << ".mhd";
//		std::string file_output;
//		file_output = ss.str();
//		writer2->SetFileName(file_output.c_str());
//		writer2->SetInput(img_tr);
//		writer2->Update();
//
//		//Run Nucleus Segmentation
//		clock_t startTimer = clock();
//		std::cout<<"Starting segmentation\n";
//		unsigned char *in_Image;
//		in_Image = (unsigned char *) malloc (size[0]*size[1]*size[2]);
//		if( in_Image == NULL ){
//			std::cerr<<"Nucleus Seg failed because malloc failed\n";
//			continue;
//		}
//		memset(in_Image/*destination*/,0/*value*/,size[0]*size[1]*size[2]*sizeof(unsigned char)/*num bytes to move*/);
//		typedef itk::ImageRegionConstIterator< nucImageType > ConstIteratorType;
//		ConstIteratorType pix_buf( img, img->GetRequestedRegion() );
//		uint64_t ind=0;
//		for ( pix_buf.GoToBegin(); !pix_buf.IsAtEnd(); ++pix_buf, ++ind )
//			in_Image[ind]=(pix_buf.Get());
//		yousef_nucleus_seg *NucleusSeg = new yousef_nucleus_seg();
//		NucleusSeg->readParametersFromFile(argv[4]);
//		NucleusSeg->setDataImage(in_Image,size[0],size[1],size[2],"/home/gramak/Desktop/11_10_tracing/montage_kt11306_w410DAPIdsu.mhd");
//		NucleusSeg->runBinarization();
//		try 
//		{
//			std::cout<<"Starting seed detection\n";
//			NucleusSeg->runSeedDetection();
//		}
//		catch( bad_alloc & excp )
//		{
//			std::cout<<"You have requested more memory than "
//			<<"what is currently available in this "
//			<<"system, please try again with a smaller "
//			<<"input image\n";
//		}
//		catch( itk::ExceptionObject & excp )
//		{
//			std::cout<<"Error: " << excp <<std::endl;
//		}
//		NucleusSeg->runClustering();
//		unsigned short *output_img;
///*		if(NucleusSeg->isSegmentationFinEnabled())
//		{
//			NucleusSeg->runAlphaExpansion();
//			output_img=NucleusSeg->getSegImage();
//		} 
//		else
//*/
//		output_img=NucleusSeg->getClustImage();
//		std::cout<<"Done with nucleus segmentation for image "<<a<<std::endl;
//
//		LabelType::Pointer image = LabelType::New();
//		LabelType::PointType origin;
//		origin[0] = start[0]; //Is this OK?
//		origin[1] = start[1];
//		origin[2] = start[2];
//		image->SetOrigin( origin );
//		LabelType::IndexType start1;
//		start1[0] = 0;
//		start1[1] = 0;
//		start1[2] = 0;
//		LabelType::SizeType size1;
//		size1[0] = size[0];
//		size1[1] = size[1];
//		size1[2] = size[2];
//		LabelType::RegionType region;
//		region.SetSize ( size1  );
//		region.SetIndex( start1 );
//		image->SetRegions( region );
//		image->Allocate();
//		image->FillBuffer(0);
//		image->Update();
//
//		typedef itk::ImageRegionIteratorWithIndex< LabelType > IteratorType;
//		IteratorType iterator1(image,image->GetRequestedRegion());
//		for(uint64_t i=0; i<(size[0]*size[1]*size[2]); ++i)
//		{
//			unsigned short val = (unsigned short)output_img[i];
//			iterator1.Set(val);
//			++iterator1;
//		}
//
//		std::stringstream ss1;
//		ss1 << "soma_label"<< x << "_" << y << "_" << z << ".tif";
//		std::string file_output1;
//		file_output1 = ss1.str();
//		itk::ImageFileWriter<LabelType>::Pointer writer = itk::ImageFileWriter<LabelType>::New();
//		writer->SetFileName(file_output1);
//		writer->SetInput(image);
//		writer->Update();
//		delete NucleusSeg;
//	}

	long int a=0;
	int num_to_be_traced = centroid_list.size();
	std::cout << "Deleting duplicated seeds....  \n";
	while( a<centroid_list.size() )
	{
		//Run Multiple Neuron Tracer
		int x, y, z;
		x = centroid_list[a][0];
		y = centroid_list[a][1];
		z = centroid_list[a][2];
		if(x == -1)
		{
			++a;
			continue;
		}
		itk::Index<3> centroid;
		centroid[0] = ((x - 200)>0) ? 200:x; 
		centroid[1] = ((y - 200)>0) ? 200:y;
		centroid[2] = ((z - 100) >0) ? 100:z;
		std::stringstream ss1;
		ss1 << "D:\\11_10_selective_tracing\\soma_label"<< centroid_list.at(a)[0] << "_" << centroid_list.at(a)[1] << "_" << centroid_list.at(a)[2] << ".tif";
		std::string file_output1;
		file_output1 = ss1.str();
		itk::ImageFileReader<LabelType>::Pointer reader = itk::ImageFileReader<LabelType>::New();
		reader->SetFileName(file_output1);
		reader->Update();
		LabelType::Pointer image = reader->GetOutput();
		std::cout << image->GetBufferedRegion().GetSize() << "\n";
		typedef itk::ImageRegionIteratorWithIndex< LabelType > IteratorType;
		IteratorType iterator1(image,image->GetRequestedRegion());

		std::vector< itk::Index<3> > soma_centroids;
		std::vector< itk::Index<3> > soma_centroids_to_be_del;
		std::vector<unsigned short> label_vec;
		std::vector<uint64_t> delete_index;
		for(uint64_t ctr =0; ctr<all_centroid_list.size() ; ++ctr)
		{
			itk::Index<3> cen =  all_centroid_list[ctr];
			if(cen[0] == -1) continue;
			if(abs((double)(cen[0]-x))<200 && abs((double)(cen[1]-y))<200 && abs((double)(cen[2]-z))<100 )
			{
				bool centroid_not_found = true;
				itk::Index<3> centroid2;
				centroid2[0] = centroid[0] + cen[0] - x; 
				centroid2[1] = centroid[1] + cen[1] - y;
				centroid2[2] = centroid[2] + cen[2] - z;
				LabelType::IndexType centroid_index;
				centroid_index[0] = centroid2[0]; centroid_index[1] = centroid2[1]; centroid_index[2] = centroid2[2];
				iterator1.SetIndex( centroid_index );
				unsigned short centroid_label = iterator1.Get();
				for(int b=0; b<label_vec.size(); ++b)
				{
					if( label_vec.at(b)==centroid_label )
					{
						centroid_not_found = false;
						break;
					}
				}
				if( centroid_not_found )
				{
					label_vec.push_back( centroid_label );
					soma_centroids.push_back( centroid2 );
				}
				if( !centroid_not_found )
				{
					delete_index.push_back( ctr );
					soma_centroids_to_be_del.push_back( cen );
				}
			}
		}

		if( delete_index.size() )
		{
			std::sort(delete_index.begin(), delete_index.end(), myfunction);
			for( uint64_t b=0; b<delete_index.size(); ++b )
				all_centroid_list[delete_index[b]].Fill(-1);
			delete_index.clear();
			for( uint64_t b=0; b<soma_centroids_to_be_del.size(); ++b )
			{
				itk::Index<3> centroid2 = soma_centroids_to_be_del.at(b);
				for( uint64_t c=0; c<centroid_list.size(); ++c )
				{
					itk::Index<3> centroid3 = centroid_list.at(c);
					if( centroid2[0]==centroid3[0] && centroid2[1]==centroid3[1] && centroid2[2]==centroid3[2] )
					{
						delete_index.push_back( c );
						break;
					}
				}
			}
			std::sort(delete_index.begin(), delete_index.end(), myfunction);
			for( uint64_t b=0; b<delete_index.size(); ++b )
			{
				centroid_list[delete_index[b]].Fill(-1);
				--num_to_be_traced;
			}
			std::cout<<"Number of cells to be retraced : "<<num_to_be_traced << " after "<<a<<" processed\n";
		}
		++a;
	}
	std::cout << "Deleting duplicated seeds done !! \n";

		//std::ostringstream swc_filename_stream;
		//swc_filename_stream << vul_file::strip_extension(argv[3]) << "_" << x << "_" << y << "_" << z << "_ANT.swc";


	//ofstream outFile; 
	//outFile.open("D:\\11_10_selective_tracing\\append_to_swc.txt", ios::out | ios::trunc );
	//if ( !outFile.is_open() )
	//{
	//	std::cerr << "Failed to Load Document: " << outFile << std::endl;
	//	return false;
	//}
	////Write out the features:
	//for( a=0; a<centroid_list.size(); ++a )
	//{
	//	int x, y, z;
	//	x = centroid_list[a][0];
	//	y = centroid_list[a][1];
	//	z = centroid_list[a][2];
	//	if(x == -1) continue;
	//	std::stringstream ssx, ssy, ssz;
	//	ssx << x; ssy << y; ssz << z;
	//	outFile << "_Trace_" + ssx.str() + "_" + ssy.str() + "_" + ssz.str() + "_ANT.swc\n";
	//}
	//outFile.close();


	#pragma omp parallel for num_threads(12)
	for( a=0; a<centroid_list.size(); ++a )
	{
		int x, y, z;
		x = centroid_list[a][0];
		y = centroid_list[a][1];
		z = centroid_list[a][2];
		if(x == -1) continue;
		std::cout<< "Tracing from seed point " << x << "_" << y << "_" << z << "\n";
		itk::Index<3> centroid;
		centroid[0] = ((x - 200)>0) ? 200:x; 
		centroid[1] = ((y - 200)>0) ? 200:y;
		centroid[2] = ((z - 100) >0) ? 100:z;
		std::stringstream ss1;
		ss1 << "D:\\11_10_selective_tracing\\soma_label"<< centroid_list.at(a)[0] << "_" << centroid_list.at(a)[1] << "_" << centroid_list.at(a)[2] << ".tif";
		std::string file_output1;
		file_output1 = ss1.str();
		itk::ImageFileReader<LabelType>::Pointer reader = itk::ImageFileReader<LabelType>::New();
		reader->SetFileName(file_output1);
		reader->Update();
		LabelType::Pointer image = reader->GetOutput();

		itk::ImageFileReader<gfpImageType>::Pointer reader2 = itk::ImageFileReader<gfpImageType>::New();
		std::stringstream ss;
		ss << "D:\\11_10_selective_tracing\\gfp_image" << x << "_" << y << "_" << z << ".mhd";
		std::string file_output;
		file_output = ss.str();
		reader2->SetFileName(file_output.c_str());
		reader2->Update();
		gfpImageType::Pointer img_tr = reader2->GetOutput();

		std::vector< itk::Index<3> > soma_centroids;      
		for(int ctr =0; ctr<all_centroid_list.size() ; ++ctr)
		{
			itk::Index<3> cen =  all_centroid_list[ctr];
			if(cen[0] == -1) continue;
			if(abs((double)(cen[0]-x))<=200 && abs((double)(cen[1]-y))<=200 && abs((double)(cen[2]-z))<=100 )
			{
				itk::Index<3> centroid2;
				centroid2[0] = centroid[0] + cen[0] - x; //Is there a reason why x and y are flipped?
				centroid2[1] = centroid[1] + cen[1] - y;
				centroid2[2] = centroid[2] + cen[2] - z;
				soma_centroids.push_back(centroid2);
			}
		}

		MultipleNeuronTracer * MNT = new MultipleNeuronTracer();
		MNT->LoadCurvImage_1(img_tr, 1);
		MNT->ReadStartPoints_1(soma_centroids, 1);
		MNT->SetCostThreshold(600);
		MNT->LoadSomaImage_1(image);
		MNT->RunTracing();
		/*MNT->WriteSWCFile(std::string(swc_filename_stream.str()), 1);*/
		
		std::stringstream ssx, ssy, ssz;
		ssx << x; ssy << y; ssz << z;
		MNT->WriteSWCFile("D:\\11_10_selective_tracing\\_Trace_" + ssx.str() + "_" + ssy.str() + "_" + ssz.str() + "_ANT.swc", 1);

		delete MNT;
	}

}


