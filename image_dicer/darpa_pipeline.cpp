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

#ifdef _OPENMP
#include "omp.h"
#endif

typedef itk::Image<unsigned char,  3> nucImageType;
typedef itk::Image<unsigned short, 3> LabelType;
typedef MultipleNeuronTracer::ImageType3D gfpImageType;

//typedef struct{ double x_d; double y_d; double z_d; } centroids;

int main(int argc, char* argv[]){
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
	vtkSmartPointer<vtkTable> global_centroids = ftk::LoadTable("C:\\Data\\Darpa\\Global_centroids.txt");
	for(int r=0; r<(int)global_centroids->GetNumberOfRows(); ++r)
	{
		int cx = global_centroids->GetValue(r, 0).ToInt();
		int cy = global_centroids->GetValue(r, 1).ToInt();
		int cz = global_centroids->GetValue(r, 2).ToInt();
		if( (fmod((double)cx,1050)>1000) || (fmod((double)cx,1050)<50) || (fmod((double)cy,1050)>1000) || (fmod((double)cy,1050)<50))
		{
			itk::Index<3> cen;
			cen[0] = cx; cen[1] = cy; cen[2] = cz; 
			centroid_list.push_back(cen);

		}
	}

	typedef itk::ImageFileReader<nucImageType> nucReaderType;
	typedef itk::ImageFileReader<gfpImageType> gfpReaderType;
	nucReaderType::Pointer reader_nuc   = nucReaderType::New();
	reader_nuc->  SetFileName("C:\\Data\\Darpa\\montage_DAPI.tif");
	reader_nuc->  Update();
	nucImageType::Pointer img_nuc   = reader_nuc->GetOutput();
	gfpReaderType::Pointer reader_trace = gfpReaderType::New();
	reader_trace->SetFileName("C:\\Data\\Darpa\\montage_GFP.tif");
	reader_trace->Update();
	gfpImageType::Pointer img_trace = reader_trace->GetOutput();
	uint64_t size_trace[3],size_nuc[3];
	size_nuc[0] = img_nuc->GetLargestPossibleRegion().GetSize()[0];
	size_nuc[1] = img_nuc->GetLargestPossibleRegion().GetSize()[1];
	size_nuc[2] = img_nuc->GetLargestPossibleRegion().GetSize()[2];
	size_trace[0] = img_trace->GetLargestPossibleRegion().GetSize()[0];
	size_trace[1] = img_trace->GetLargestPossibleRegion().GetSize()[1];
	size_trace[2] = img_trace->GetLargestPossibleRegion().GetSize()[2];
	if( size_trace[0]!=size_nuc[0] || size_trace[1]!=size_nuc[1] || size_trace[2]!=size_nuc[2] ) return EXIT_FAILURE;

	#pragma omp parallel for  
	for( uint64_t a=0; a<centroid_list.size(); ++a )
	{
		int x, y, z;
		x = centroid_list[a][0];
		y = centroid_list[a][1];
		z = centroid_list[a][2];

		nucImageType::IndexType start;
		start[0] = ((x - 200)>0) ? (x - 200):0; //Is there a reason why x and y are flipped?
		start[1] = ((y - 200)>0) ? (y - 200):0;
		start[2] = ((z - 75) >0) ? (z - 75) :0;

		gfpImageType::IndexType start2;
		start2[0] = ((x - 200)>0) ? (x - 200):0; //Is there a reason why x and y are flipped?
		start2[1] = ((y - 200)>0) ? (y - 200):0;
		start2[2] = ((z - 75) >0) ? (z - 75) :0;

		nucImageType::SizeType size;
		size[0] = ((x+200)<size_nuc[0]) ? 400 : (200+size_trace[0]-x-1); //Is there a reason why x and y are flipped?
		size[1] = ((y+200)<size_nuc[1]) ? 400 : (200+size_trace[1]-y-1);
		size[2] = ((z+75)<size_nuc[2]) ? 150 : (75+size_trace[2]-z-1);

		gfpImageType::SizeType size2;
		size2[0] = ((x+200)<size_trace[0]) ? 400 : (200+size_trace[0]-x-1); //Is there a reason why x and y are flipped?
		size2[1] = ((y+200)<size_trace[1]) ? 400 : (200+size_trace[1]-y-1);
		size2[2] = ((z+75)<size_trace[2]) ? 150 : (75+size_trace[2]-z-1);

		itk::Index<3> centroid;
		centroid[0] = ((x - 200)>0) ? 200:x; //Is there a reason why x and y are flipped?
		centroid[1] = ((y - 200)>0) ? 200:y;
		centroid[2] = ((z - 75) >0) ? 75:z;

		std::ostringstream output_filename_stream;

		nucImageType::RegionType desiredRegion;
		desiredRegion.SetSize(size);
		desiredRegion.SetIndex(start);

		gfpImageType::RegionType desiredRegion2;
		desiredRegion2.SetSize(size2);
		desiredRegion2.SetIndex(start2);

		typedef itk::RegionOfInterestImageFilter< nucImageType, nucImageType > ROIFilterType;
		typedef itk::RegionOfInterestImageFilter< gfpImageType, gfpImageType > ROIFilterType2;
		ROIFilterType::Pointer ROIfilter = ROIFilterType::New();
		ROIfilter->SetRegionOfInterest(desiredRegion);
		ROIfilter->SetInput(img_nuc);
		ROIfilter->Update();
		nucImageType::Pointer img = ROIfilter->GetOutput();

		ROIFilterType2::Pointer ROIfilter2 = ROIFilterType2::New();
		ROIfilter2->SetRegionOfInterest(desiredRegion2);
		ROIfilter2->SetInput(img_trace);
		ROIfilter2->Update();
		gfpImageType::Pointer img_tr = ROIfilter2->GetOutput();

		itk::CastImageFilter<gfpImageType, nucImageType>::Pointer caster = itk::CastImageFilter<gfpImageType, nucImageType>::New();
		caster->SetInput(img_tr);
		//itk::ImageFileWriter<nucImageType>::Pointer writer2 = itk::ImageFileWriter<nucImageType>::New();
		//writer2->SetFileName("C:\\ROYSAMLAB\\FARSIGHT\\Farsight_DASH_bin\\exe\\Release\\gfp_File.tif");
		//writer2->SetInput(caster->GetOutput());
		//writer2->Update();

		//Run Nucleus Segmentation
		clock_t startTimer = clock();
		std::cout<<"Starting segmentation\n";
		unsigned char *in_Image;
		in_Image = (unsigned char *) malloc (size[0]*size[1]*size[2]);
		if( in_Image == NULL ) return EXIT_FAILURE;
		memset(in_Image/*destination*/,0/*value*/,size[0]*size[1]*size[2]*sizeof(unsigned char)/*num bytes to move*/);
		typedef itk::ImageRegionConstIterator< nucImageType > ConstIteratorType;
		ConstIteratorType pix_buf( img, img->GetRequestedRegion() );
		uint64_t ind=0;
		for ( pix_buf.GoToBegin(); !pix_buf.IsAtEnd(); ++pix_buf, ++ind )
			in_Image[ind]=(pix_buf.Get());
		yousef_nucleus_seg *NucleusSeg = new yousef_nucleus_seg();
		NucleusSeg->readParametersFromFile("C:\\ROYSAMLAB\\FARSIGHT\\Farsight_DASH_bin\\exe\\Release\\Seg_Params.ini");
		NucleusSeg->setDataImage(in_Image,size[0],size[1],size[2],"");
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
		for(uint64_t i=0; i<(size[0]*size[1]*size[2]); ++i)
		{
			unsigned short val = (unsigned short)output_img[i];
			iterator1.Set(val);
			++iterator1;
		}
		delete NucleusSeg;
			
		//itk::ImageFileWriter<LabelType>::Pointer writer = itk::ImageFileWriter<LabelType>::New();
		//writer->SetFileName("C:\\ROYSAMLAB\\FARSIGHT\\Farsight_DASH_bin\\exe\\Release\\Soma_File.tif");
		//writer->SetInput(image);
		//writer->Update();

		//Run Multiple Neuron Tracer
		std::vector< itk::Index<3> > soma_centroids;
		soma_centroids.push_back(centroid);
		
		std::ostringstream swc_filename_stream;
		//swc_filename_stream << vul_file::strip_extension(argv[3]) << "_" << x << "_" << y << "_" << z << "_ANT.swc";
		
		MultipleNeuronTracer * MNT = new MultipleNeuronTracer();
		MNT->LoadCurvImage_1(img_tr, 1);
		MNT->ReadStartPoints_1(soma_centroids, 1);
		MNT->SetCostThreshold(300);
		MNT->LoadSomaImage_1(image);
		MNT->RunTracing();
		/*MNT->WriteSWCFile(std::string(swc_filename_stream.str()), 1);*/
		
		std::stringstream ssx, ssy, ssz;
		ssx << x; ssy << y; ssz << z;
		MNT->WriteSWCFile("Trace_" + ssx.str() + "_" + ssy.str() + "_" + ssz.str() + "_ANT.swc", 1);

		


	}

}


