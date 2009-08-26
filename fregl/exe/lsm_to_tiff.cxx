/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/

//: Save an lsm image to a list of tiff images
//
//  Only deals with single time-frame for now. The input is the lsm image,
//  and the output is a set of grayscale tiff images, each having the
//  filename of output_base_name_channel_name.tiff
//
//  lsm_file [output_tiff_base_name] [output_tiff_path]
//

#include <ftkImage/vtkLSMReader.h>
#include <string>
#include <vector>
#include <vtkImageData.h>
#include "itkImage.h"
#include "itkImageFileWriter.h"

typedef itk::Image<unsigned char,3> ImageType;

int main( int argc, char* argv[] )
{
  if (argc<2) {
	  std::cerr << "Usage: " << argv[0] << " lsm_file [output_tiff_base_name] [output_tiff_path]";
    return EXIT_FAILURE;
  }
  
  vtkLSMReader * lsmR = vtkLSMReader::New();
  std::string lsm_name = argv[1];
  lsmR->SetFileName(lsm_name.c_str());
  if (!lsmR->OpenFile()){
    std::cerr<<"\nFailed to open"<<lsm_name<<std::endl;
    std::cerr <<"\nUsage: " << argv[0] << " lsm_file [output_tiff_base_name] [output_tiff_path] ";
    return EXIT_FAILURE;
  }
  lsmR->Update();

  int numChannels = lsmR->GetNumberOfChannels();
  int numTimes = lsmR->GetNumberOfTimePoints();

  std::cout<<"numChannels = "<<numChannels<<std::endl;
  std::cout<<"numTimes = "<<numTimes<<std::endl;

  // Choose base name for output images
  std::string base_name;
  if (argc >= 3)
  {
	  base_name = argv[2];
  }
  else
  {
	  // Get the base name of the lsm image
	  const std::string slash = "\\/";
	  const std::string dot = ".";
	  std::string::size_type po = lsm_name.find_last_of(slash);
	  if (po != std::string::npos) lsm_name = lsm_name.substr(po+1,lsm_name.size()-po-1);
	  po = lsm_name.find_last_of(dot);
	  base_name = lsm_name.substr (0, po);
  }

  // Assign tiff file names
  std::vector<std::string> tiff_names;
  for (int i = 0; i<numChannels; i++) {
	//We want to remove anything after a space in the channel name:
	std::string chName = std::string(lsmR->GetChannelName(i));
	const std::string delin = " -/\\.";
	std::string::size_type po = chName.find_first_of(delin);
	chName = chName.substr(0,po);

    std::string tiff_name = base_name + std::string("_") + chName + std::string(".tiff");
    std::cout<<"tiff name = "<<tiff_name<<std::endl;
    tiff_names.push_back(tiff_name);
  }

  vtkImageData * vimdata;
  int imsize, size_x, size_y, size_z;
  int extent[6]; //recording dimensions of the images. extent[0,2,4]
                 //is the min location, and extent[1,3,5] is the max
                 //location in 3D
 
  for (int counter=0; counter < numChannels; counter++)
  {
    int time_stamp = 0;
    vimdata = lsmR->GetTimePointOutput(time_stamp,counter);
    lsmR->Update();
    vimdata->GetExtent(extent);
    size_x = extent[1]-extent[0]+1;
    size_y = extent[3]-extent[2]+1;
    size_z = extent[5]-extent[4]+1;
    imsize = (size_z)*(size_y)*(size_x);

    unsigned char *channel_data = (unsigned char*) vimdata->GetScalarPointer();

    int red = lsmR->GetChannelColorComponent(counter,0);
    int green = lsmR->GetChannelColorComponent(counter,1);
    int blue = lsmR->GetChannelColorComponent(counter,2);
    
    std::cout<<"rgb for channel "<<counter<<" = ["<<red<<","<<green<<","<<blue<<"]\n";
    
    ImageType::Pointer image = ImageType::New();
    ImageType::IndexType start;
    start[0] = 0; // first index on X
    start[1] = 0; // first index on Y
    start[2] = 0; // first index on Z
    
    ImageType::SizeType size;
    size[0] = size_x; // size along X
    size[1] = size_y; // size along Y
    size[2] = size_z; // size along Z
    
    ImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );

    image->SetRegions( region );
    image->Allocate();

    image->GetPixelContainer()->SetImportPointer( channel_data, sizeof(unsigned char), false );
    
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
    if (argc == 4) {
      writer->SetFileName(std::string(argv[3]) + tiff_names[counter]);
      std::cout<<std::string(argv[3])+tiff_names[counter]<<std::endl;
    }
    else writer->SetFileName(tiff_names[counter]);
    writer->SetInput( image );
    writer->Update();
  }
  return 0;
}




