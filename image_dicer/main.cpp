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

#include "../NuclearSegmentation/yousef_core/yousef_seg.h"
#include <iostream>
#include <time.h>
#include <limits.h>

typedef itk::Image<unsigned char,  3> ImageType;
typedef itk::Image<unsigned short, 3> LabelType;

typedef struct{ double x_d; double y_d; double z_d; } centroids;

int main(int argc, char* argv[]){
	//Input: Somas_file, nucleus_montage, trace_channel_montage, nuc_seg_parameters
	std::fstream soma_centroid_file;
	soma_centroid_file.open(argv[1], std::fstream::in);
	std::vector<centroids> centroid_list;
	while (soma_centroid_file.good()){
		centroids new_centroid;
		soma_centroid_file >> new_centroid.x_d >> new_centroid.y_d >> new_centroid.z_d;
		centroid_list.push_back( new_centroid );
	}
	soma_centroid_file.close();

	typedef itk::ImageFileReader<ImageType> ReaderType;
	ReaderType::Pointer reader_nuc   = ReaderType::New();
	reader_nuc->  SetFileName(argv[2]);
	reader_nuc->  Update();
	ImageType::Pointer img_nuc   = reader_nuc->GetOutput();
	ReaderType::Pointer reader_trace = ReaderType::New();
	reader_trace->SetFileName(argv[3]);
	reader_trace->Update();
	ImageType::Pointer img_trace = reader_trace->GetOutput();
	uint64_t size_trace[3],size_nuc[3];
	size_trace[0]=img_nuc  ->GetLargestPossibleRegion().GetSize()[0];
	size_trace[1]=img_nuc  ->GetLargestPossibleRegion().GetSize()[1];
	size_trace[2]=img_nuc  ->GetLargestPossibleRegion().GetSize()[2];
	size_nuc[0]  =img_trace->GetLargestPossibleRegion().GetSize()[0];
	size_nuc[1]  =img_trace->GetLargestPossibleRegion().GetSize()[1];
	size_nuc[2]  =img_trace->GetLargestPossibleRegion().GetSize()[2];
	if( size_trace[0]!=size_nuc[0] || size_trace[1]!=size_nuc[1] || size_trace[2]!=size_nuc[2] ) return EXIT_FAILURE;

	#pragma omp parallel for  
	for( uint64_t a=0; a<centroid_list.size(); ++a ){
		uint64_t x, y, z;
		x = centroid_list.at(a).x_d;
		y = centroid_list.at(a).y_d;
		z = centroid_list.at(a).z_d;

		ImageType::IndexType start;
		start[0] = ((y - 150)>0) ? (y - 150):0; //Is there a reason why x and y are flipped?
		start[1] = ((x - 150)>0) ? (x - 150):0;
		start[2] = ((z - 50) >0) ? (z - 50) :0;

		ImageType::SizeType size;
		size[0] = ((y+300)<size_trace[0]) ? 300 : (size_trace[0]-y-1); //Is there a reason why x and y are flipped?
		size[1] = ((x+300)<size_trace[1]) ? 300 : (size_trace[1]-x-1);
		size[2] = ((z+100)<size_trace[2]) ? 100 : (size_trace[2]-z-1);

		std::ostringstream output_filename_stream;

		ImageType::RegionType desiredRegion;
		desiredRegion.SetSize(size);
		desiredRegion.SetIndex(start);

		typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > ROIFilterType;
		ROIFilterType::Pointer ROIfilter = ROIFilterType::New();
		ROIfilter->SetRegionOfInterest(desiredRegion);
		ROIfilter->SetInput(img_nuc);
		ROIfilter->Update();
		ImageType::Pointer img;

		//Run Nucleus Segmentation
		clock_t startTimer = clock();
		std::cout<<"Starting segmentation\n";
		unsigned char *in_Image;
		in_Image = (unsigned char *) malloc (size[0]*size[1]*size[2]);
		if( in_Image == NULL ) return EXIT_FAILURE;
		memset(in_Image/*destination*/,0/*value*/,size[0]*size[1]*size[2]*sizeof(unsigned char)/*num bytes to move*/);
		typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;
		ConstIteratorType pix_buf( img, img->GetRequestedRegion() );
		uint64_t ind=0;
		for ( pix_buf.GoToBegin(); !pix_buf.IsAtEnd(); ++pix_buf, ++ind )
			in_Image[ind]=(pix_buf.Get());
		yousef_nucleus_seg *NucleusSeg = new yousef_nucleus_seg();
		NucleusSeg->readParametersFromFile(argv[4]);
		NucleusSeg->setDataImage(in_Image,size[0],size[1],size[2],argv[2]);
		NucleusSeg->runBinarization();
		try {
			std::cout<<"Starting seed detection\n";
			NucleusSeg->runSeedDetection();
		}
		catch( bad_alloc & excp ){
			std::cout<<"You have requested more memory than "
			<<"what is currently available in this "
			<<"system, please try again with a smaller "
			<<"input image\n";
			return EXIT_FAILURE;
		}
		catch( itk::ExceptionObject & excp ){
			std::cout<<"Error: " << excp <<std::endl;
			return EXIT_FAILURE;
		}
		NucleusSeg->runClustering();
		unsigned short *output_img;
		if(NucleusSeg->isSegmentationFinEnabled()){
			NucleusSeg->runAlphaExpansion();
			output_img=NucleusSeg->getSegImage();
		} else
			output_img=NucleusSeg->getClustImage();

		LabelType::Pointer image = LabelType::New();
		LabelType::PointType origin;
		origin[0] = start[0]; //Is this OK?
		origin[1] = start[1];
		origin[2] = start[2];
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
		for(uint64_t i=0; i<(size[0]*size[1]*size[2]); ++i){
			unsigned short val = (unsigned short)output_img[i];
			iterator1.Set(val);
			++iterator1;
		}
		delete NucleusSeg;

		//Run Multiple Neuron Tracer

	}

}


