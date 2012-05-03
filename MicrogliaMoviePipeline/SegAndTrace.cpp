/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include "SegAndTrace.h"
#include "../NuclearSegmentation/yousef_core/yousef_seg.h"
#include "../NuclearSegmentation/exe/SomaExtraction.h"
#include "../ftkCommon/ftkProjectManager.h"
#include <QFileInfo>

MicrogliaMovieSegTracer::MicrogliaMovieSegTracer()
{

}

MicrogliaMovieSegTracer::~MicrogliaMovieSegTracer()
{

}

void MicrogliaMovieSegTracer::StartSegTracing(const char* filename, const char* paramFileName)
{
	ReadXMLProjectFile( filename);
	LoadOptions( paramFileName);
	std::cout<< "Image Stacks Number "<<files.size()<<std::endl;

	for( int i = 0; i < files.size(); i++)
	{
		CharImageType::Pointer inputImage = ReadImage(files[i].c_str());
		std::vector< itk::Index<3> > seedVec;
		GenerateSeedPointsForSoma(inputImage, seedVec);
		std::cout<<seedVec.size()<< " seeds generated!"<<std::endl;

		CasterTypeCharToFloat::Pointer cast2float = CasterTypeCharToFloat::New();
		cast2float->SetInput( inputImage);
		cast2float->Update();

		std::cout<<seedVec.size()<< "Segment Somas..."<<std::endl;
		LabelImageType::Pointer segImage = SegmentSoma( cast2float->GetOutput(), seedVec);
		WriteImage( "SomaImage.tif", segImage);
	}
}

void MicrogliaMovieSegTracer::WriteImage(const char* filename, LabelImageType::Pointer labelImage)
{
	CasterTypeLabelToChar::Pointer caster = CasterTypeLabelToChar::New();
	caster->SetInput( labelImage);
	CharWriterType::Pointer writer = CharWriterType::New();
	writer->SetInput( caster->GetOutput());
	writer->SetFileName( filename);
	writer->Update();
}

void MicrogliaMovieSegTracer::ReadXMLProjectFile(const char* filename)
{
	files.clear();
	ftk::ProjectManager * project = new ftk::ProjectManager(filename);
	for ( unsigned int i = 0; i < project->size(); i++)
	{ 
		bool found = false;
		std::string FileName = project->GetFileName(i);

		QFileInfo  NewFileInfo(QString(FileName.c_str()));
		if (!NewFileInfo.exists())
		{
			std::cout << "file not found " << FileName << std::endl;
		}
		else
		{
			found = true;
		}  
		if (found)
		{
			QString type = QString(project->GetFileType(i).c_str());
			if ((type == "Image"))
			{
				files.push_back(FileName);
			}//end type image
		} 
	}//end of for project size
}

template <class T>
bool MicrogliaMovieSegTracer::SetParamValue(std::map<std::string,std::string> &opts, std::string str, T &value, T defVal)
{
	std::map<std::string,std::string>::iterator mi;
	mi = opts.find(str); 
	if( mi != opts.end())
	{ 
		std::istringstream ss((*mi).second); 
		ss>>value;
		return true;
	}
	else
	{ 
		value = defVal; 
		return false;
	}
}

void MicrogliaMovieSegTracer::LoadOptions(const char* paramFileName)
{
	std::map<std::string,std::string> opts;
	ifstream fin(paramFileName); 
	std::string name; 
	fin>>name;
	while(fin.good()) 
	{
		char cont[100];	 
		fin.getline(cont, 99);
		opts[name] = std::string(cont);
		fin>>name;
	}
	fin.close();

	SetParamValue<int>(opts, "-high_sensitivity", shift, 0);
	SetParamValue<double>(opts, "-min_scale", scaleMin, 30);
	SetParamValue<double>(opts, "-max_scale", scaleMax, 35);
	SetParamValue<double>(opts, "-xy_clustering_res", regionXY, 30);
	SetParamValue<double>(opts, "-z_clustering_res", regionZ, 25);
	SetParamValue<int>(opts, "-sampling_ratio_XY_to_Z", useDistMap, 2);
	SetParamValue<int>(opts, "-Use_Distance_Map", sampling_ratio_XY_to_Z, 1);
	SetParamValue<double>(opts, "-sigmoid_alfa", alfa, 30);
	SetParamValue<double>(opts, "-sigmoid_beta", beta, 3);
	SetParamValue<int>(opts, "-time_threhold", timethreshold, 5);
	SetParamValue<double>(opts, "-curvature_scaling", curvatureScaling, 0.5);
	SetParamValue<double>(opts, "-rmsThreshold", rmsThres, 0.02);
	SetParamValue<int>(opts, "-hole_size", holeSize, 10);
	SetParamValue<int>(opts, "-min_object_size", minObjSize, 1000);
	SetParamValue<float>(opts, "-intensity_threshold", intensity_threshold, 0.005);
	SetParamValue<float>(opts, "-contrast_threshold", contrast_threshold, 0.0003);
	SetParamValue<int>(opts, "-cost_threshold", cost_threshold, 700);
	SetParamValue<float>(opts, "-debris_threshold", debris_threshold, 0.8);
	SetParamValue<int>(opts, "-offshoot", offshoot, 6);
	SetParamValue<int>(opts, "-device", device, 0);
}

MicrogliaMovieSegTracer::CharImageType::Pointer MicrogliaMovieSegTracer::ReadImage(const char * fileName)
{
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName (fileName);	
	reader->Update();
	CharImageType::Pointer outputImage = reader->GetOutput();
	return outputImage;
}

void MicrogliaMovieSegTracer::GenerateSeedPointsForSoma(CharImageType::Pointer inputImage, std::vector< itk::Index<3> > &seedVec)
{
	unsigned int size1 = inputImage->GetLargestPossibleRegion().GetSize()[0];
	unsigned int size2 = inputImage->GetLargestPossibleRegion().GetSize()[1];
	unsigned int size3 = inputImage->GetLargestPossibleRegion().GetSize()[2];
	std::cout<< "size"<<"\t"<<size1<<"\t"<<size2<<"\t"<<size3<<std::endl;
	unsigned char *in_Image = (unsigned char *) malloc( size1 * size2 * size3);

	if( in_Image == NULL )
	{
		std::cout<<"Memory allocation for image failed\n";
	}
	memset( in_Image, 0,size1 * size2 * size3 * sizeof(unsigned char));

	typedef itk::ImageRegionConstIterator< CharImageType > ConstIteratorType;
	ConstIteratorType pix_buf( inputImage, inputImage->GetRequestedRegion());

	unsigned int ind = 0;
	for ( pix_buf.GoToBegin(); !pix_buf.IsAtEnd(); ++pix_buf, ++ind)
	{
		in_Image[ind]=(pix_buf.Get());
	}

	vector< Seed> seeds;
	yousef_nucleus_seg *NucleusSeg = new yousef_nucleus_seg(); 
	NucleusSeg->setParamsForSeedDetection( shift, scaleMin, scaleMax, regionXY, regionZ, useDistMap, sampling_ratio_XY_to_Z);
	NucleusSeg->setDataImage(in_Image, size1, size2, size3, "");
	std::cout<< "Binarization..."<<std::endl;
	NucleusSeg->runBinarization(128);
	std::cout<< "Seed Detection..."<<std::endl;
	NucleusSeg->runSeedDetection();
	NucleusSeg->outputSeeds();
	seeds = NucleusSeg->getSeedsList();

	seedVec.clear();
	for( int i = 0; i < seeds.size(); i++)
	{
		itk::Index<3> index;
		index[0] = seeds[i].x();
		index[1] = seeds[i].y();
		index[2] = seeds[i].z();
		seedVec.push_back(index);
	}
	delete NucleusSeg;
}

MicrogliaMovieSegTracer::LabelImageType::Pointer MicrogliaMovieSegTracer::SegmentSoma(FloatImageType::Pointer inputImage, std::vector< itk::Index<3> > &seedVec)
{
	SomaExtractor *Somas = new SomaExtractor();
	SomaExtractor::SegmentedImageType::Pointer segImage = Somas->SegmentSoma(inputImage, seedVec, alfa, beta, timethreshold, curvatureScaling, rmsThres, holeSize, minObjSize);
	delete Somas;											
	return segImage;
}

void MicrogliaMovieSegTracer::StartTracing(FloatImageType::Pointer inputImage, std::vector< itk::Index<3> > &seedVec, LabelImageType::Pointer somaImage)
{

}

void MicrogliaMovieSegTracer::CalculateLMeasures()
{

}