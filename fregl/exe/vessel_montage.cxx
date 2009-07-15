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

//: Executable program to generate montaged distance map
//
//  The input is a file containing the xml file of the joint
//  transformations and the xml file containing the result sets. The
//  changes will be written back to the original xml segmentation
//  files.
//
//  vessel_montage xml_joint_transforms xml_result_sets
// 
//  where
//   xml_joint_transforms   Name of the xml file containing the joint transformations
//   xml_result_sets        Name of the xml file containing the result sets
//
//  Optional argument:
//   -anchor                Name of the anchor image.
//

#include <vector>

#include <vul/vul_arg.h>

#include <fregl/fregl_result_record.h>
#include <fregl/fregl_joint_register.h>
#include <fregl/fregl_space_transformer.h>
#include <Common/fsc_channel_accessor.h>

#include <maciej_seg/maciejSegmentation.h>
#include <maciej_seg/xml_util.h>

#include "itkTIFFImageIO.h"
#include "itkRGBAPixel.h"
#include "itkRGBPixel.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include <vul/vul_file.h>
#include <vnl/vnl_math.h>

typedef unsigned char                      InputPixelType;
typedef itk::Image< InputPixelType, 3 >    ImageType;
typedef unsigned short                     DMPixelType; //for distance map
typedef itk::Image< DMPixelType, 3 >       DMImageType;
typedef itk::RGBPixel< unsigned char >     OutputPixelType;
typedef itk::Image< OutputPixelType, 3 >   ColorImageType;
//typedef itk::ImageFileWriter< ColorImageType >  WriterType3D;
typedef itk::ImageFileWriter< ImageType >  WriterType3D;
typedef itk::DanielssonDistanceMapImageFilter< ImageType, DMImageType > DMFilterType;
typedef std::vector<vnl_vector_fixed<unsigned short, 3> > point_list;
typedef ImageType::SizeType                 SizeType;
typedef itk::RescaleIntensityImageFilter< DMImageType, ImageType > RescalerType;
typedef itk::ImageRegionIterator< ImageType > RegionIterator;
typedef itk::ImageRegionConstIterator< ImageType > RegionConstIterator;
typedef itk::ImageRegionConstIterator< DMImageType > RegionConstIterator2;

void
read_surface_points(std::string const& vessel_xml, SizeType const& size,
                    ImageType::Pointer& image, point_list & vessel_seg_result)
{
  std::ifstream in_file_str( vessel_xml.c_str() );
  if ( !in_file_str ){
    std::cerr<<"Couldn't open "<<vessel_xml<<std::endl;
    exit( 0 );
  }

  std::string line_str;
  vnl_vector_fixed<unsigned short, 3> point;
  vnl_vector_fixed<unsigned short, 3> one(1,1,1);
  while ( in_file_str ) {
    std::getline( in_file_str, line_str );
    if (line_str.length() == 0) continue;

    std::istringstream line_stream(line_str);
    line_stream>> point[2] >> point[1] >> point[0];
    //We have to subtract [1,1,1] from the position, since Arun starts
    //the index from [1,1,1], instead of [0,0,0]
    point = point-one;
    vessel_seg_result.push_back( point );
  }
  in_file_str.close();

  // Now generate the binary image
  ImageType::IndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  ImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );
  
  image->SetRegions( region );
  image->Allocate();
  image->FillBuffer(itk::NumericTraits<unsigned char>::Zero);
  
  point_list::iterator pt_iter;
  ImageType::IndexType pixelIndex;
  for (pt_iter = vessel_seg_result.begin(); pt_iter != vessel_seg_result.end();
       pt_iter++) {
    pixelIndex[0] = (*pt_iter)[0]; // x position
    pixelIndex[1] = (*pt_iter)[1]; // y position
    pixelIndex[2] = (*pt_iter)[2]; // z position
    image->SetPixel( pixelIndex, 255 );
  }
}

void
write_surface_points(std::string const & vessel_xml, ImageType::Pointer& image)
{
  std::ofstream file_str( vessel_xml.c_str() );
  if ( !file_str ){
    std::cerr<<"Couldn't open "<<vessel_xml<<" for writing"<<std::endl;
    exit( 0 );
  }
  //vnl_vector_fixed<unsigned short, 3> point;
  //vnl_vector_fixed<unsigned short, 3> one(1,1,1);
  ImageType::IndexType pixelIndex;
  char value;
  RegionConstIterator outputIt( image, image->GetRequestedRegion() );
  for ( outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt) {
    value = outputIt.Get();
    if (value > 0) {
      pixelIndex = outputIt.GetIndex();
      file_str<<pixelIndex[2]+1<<"\t"<<pixelIndex[1]+1<<"\t"<<pixelIndex[0]+1
              <<"\t0\t0\n";
    }
  }
  
}

int 
main( int argc, char* argv[] )
{
  vul_arg< vcl_string > arg_xml_xform  ( 0, "Name of the xml file containing the joint transformations" );
  vul_arg< vcl_string > arg_xml_result ( 0, "Name of the xml file containing the result sets" );
  vul_arg< vcl_string > arg_anchor ("-anchor", "Name of the anchor image");
  vul_arg< vcl_string > arg_outfile   ("-output","The name of the otuput image");
  vul_arg< bool > arg_in_anchor("-in_anchor", "Generate only the updated image in the anchor space", false);
  vul_arg< bool > arg_vessel ("-vessel","Output vessel image", false);
  vul_arg< bool > arg_no_map ("-no_map","Do not generate montaged distance map", false);

  vul_arg_parse( argc, argv );

  // Read in the xml_xform file
  fregl_joint_register::Pointer joint_register = new fregl_joint_register( arg_xml_xform() );
  std::vector<fregl_result_record::Pointer> result_records;
  result_record_read_xml(arg_xml_result(), result_records);
  
  // Set the space transformer, which allows only the anchor image and
  // the adjacent images to be operated on
  //
  fregl_space_transformer space_transformer(joint_register);
  bool overlap_only = true;
  bool in_anchor = arg_in_anchor();
  std::vector<std::string> const & image_names = joint_register->image_names();
  
  // Take each image as the anchor in turn. Load the vessel results of
  // the anchor and the overlapping images. Generate the distance
  // transform and output the montage.
  //
  for (unsigned int i = 0; i<image_names.size(); i++) {
    std::string anchor_name = image_names[i];
    if (arg_anchor.set() && anchor_name!=arg_anchor()) continue;
    std::cout<<"anchor image = "<<anchor_name<<std::endl;

    space_transformer.set_anchor( anchor_name, in_anchor, overlap_only );
    
    // load the segmentation files of the anchor and the overlapping
    // images. Store the result in a binary image, with 255 indicating
    // foreground and 0 for background.
    
    std::vector<std::string> overlapping_image_names = space_transformer.image_names();
    std::vector<SizeType> overlapping_image_sizes = space_transformer.image_sizes();
    std::vector<ImageType::Pointer> vessel_seg_images;
    std::vector<point_list> vessel_seg_results;
    std::vector<ImageType::Pointer> xformed_distance_maps;

    int anchor_index = -1;
    std::string anchor_vessel_xml;

    for (unsigned int j = 0; j<overlapping_image_names.size(); j++) {     
      point_list empty_point_list;

      for (unsigned int k = 0; k<result_records.size(); k++) { 
        if (overlapping_image_names[j] == anchor_name) {
          anchor_index = j;
          anchor_vessel_xml = result_records[k]->vessel_xml();
        }

        // Read in the segmentation result from the file
        if (result_records[k]->registration_image() == overlapping_image_names[j]) {
          std::cout<<"read surface point from "<<result_records[k]->vessel_xml()<<std::endl;
          // assuming everything to be in the current directory
          ImageType::Pointer image = ImageType::New();
          vessel_seg_results.push_back( empty_point_list );
          read_surface_points(result_records[k]->vessel_xml(), 
                              overlapping_image_sizes[j], image, 
                              vessel_seg_results[j]);
          vessel_seg_images.push_back( image );
        }
      }
    }

    // Transform the points of the adjacent images in the overlapping
    // regions to the anchor's result image, and vice versa
    for (unsigned int j = 0; j<vessel_seg_results.size(); j++) {
      if (j == anchor_index) continue;

      point_list::const_iterator pt_iter;
      ImageType::IndexType pixelIndex;
      vnl_vector_fixed< float, 3 > loc, xformed_loc;

      //adjacent to the anchor
      point_list const& list = vessel_seg_results[j];
      for (pt_iter = list.begin(); pt_iter != list.end(); pt_iter++) {
        loc[0] = (*pt_iter)[0];
        loc[1] = (*pt_iter)[1];
        loc[2] = (*pt_iter)[2];
        if (space_transformer.in_anchor(loc, j, xformed_loc)) {
          pixelIndex[0] = xformed_loc[0]; // x position
          pixelIndex[1] = xformed_loc[1]; // y position
          pixelIndex[2] = xformed_loc[2]; // z position
          vessel_seg_images[anchor_index]->SetPixel( pixelIndex, 255 );
        }
      }
      
      //anchor to the adjacent
      point_list const& a_list = vessel_seg_results[anchor_index];
      for (pt_iter = a_list.begin(); pt_iter != a_list.end(); pt_iter++) {
        loc[0] = (*pt_iter)[0];
        loc[1] = (*pt_iter)[1];
        loc[2] = (*pt_iter)[2];
        if (space_transformer.in_image(loc, j, xformed_loc)) {
          pixelIndex[0] = xformed_loc[0]; // x position
          pixelIndex[1] = xformed_loc[1]; // y position
          pixelIndex[2] = xformed_loc[2]; // z position
          vessel_seg_images[j]->SetPixel( pixelIndex, 255 );
        }
      }
    }

    //update the vessel xml file of the anchor
    std::cout<<"Updating the anchor vessel file ..."<<std::endl;
    write_surface_points(anchor_vessel_xml, vessel_seg_images[anchor_index]);

    if (arg_no_map()) return 0;

    // Construct the Danielsson distance map for each image. No
    // scaling is done on the distances. Anything higher than 255 is
    // set to 255. The map is transformed.
    for (unsigned int j=0; j<vessel_seg_images.size(); j++) {
      std::cout<<"constructing the distance map for image "
               <<overlapping_image_names[j]<<std::endl;
      DMFilterType::Pointer f_daniel =  DMFilterType::New(); 
      //RescalerType::Pointer scaler = RescalerType::New();
 
      //scaler->SetOutputMaximum( 255 );
      //scaler->SetOutputMinimum( 0 );

      f_daniel->InputIsBinaryOn();
      f_daniel->SetInput(vessel_seg_images[j]);
      //scaler->SetInput(f_daniel->GetOutput());
      try
        {
          //scaler->Update();
	  f_daniel->Update();
        }
      catch (itk::ExceptionObject & e)
        {
          std::cerr << "Exception in DanielFiltere: " << e << std::endl;
          exit(0);
        }

      // For display purpose, truncate the number at 255 and save the
      // image as unsigned char
      DMImageType::Pointer DM = f_daniel->GetOutput();
      ImageType::Pointer image = ImageType::New();
      image->SetRegions( DM->GetLargestPossibleRegion() );
      image->Allocate();

      RegionConstIterator2 inputIt( DM, DM->GetRequestedRegion() );
      RegionIterator outputIt( image, image->GetRequestedRegion() );
      for ( inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd();  
            ++inputIt, ++outputIt) {
        outputIt.Set( vnl_math_min( inputIt.Get(), 255 ) );
      }

      ImageType::Pointer xformed_image = space_transformer.transform_image(image, j, 255);
      xformed_distance_maps.push_back( xformed_image );
    }

    // Generate the montage of the transformed distance maps by taking
    // the minimum
    ImageType::Pointer final_image = xformed_distance_maps[0];
    for (unsigned int  j = 1; j<xformed_distance_maps.size(); j++) {
      ImageType::Pointer xformed_image = xformed_distance_maps[j];
      RegionConstIterator inputIt1( xformed_image, xformed_image->GetRequestedRegion() );
      RegionIterator outputIt( final_image, final_image->GetRequestedRegion() );
      for ( inputIt1.GoToBegin(), outputIt.GoToBegin(); !inputIt1.IsAtEnd();  
            ++inputIt1, ++outputIt) {
        outputIt.Set( vnl_math_min( outputIt.Get(), inputIt1.Get() ) );
      }
    }

    std::string outfile_name;
    if (!arg_outfile.set()) {
      std::string id = vul_file::strip_extension( anchor_name );
      outfile_name = id+"_DM_montage.tif";
    }
    else outfile_name = arg_outfile();
    WriterType3D::Pointer writer3D = WriterType3D::New();
    writer3D->SetFileName( outfile_name );
    writer3D->SetInput( final_image );
    writer3D->Update();

    if (arg_vessel.set()) {
      // If the surface points are needed, transform the binary image
      // and display using a different channel.
      for (unsigned int j=0; j<vessel_seg_images.size(); j++) {
        ImageType::Pointer xformed_image = 
          xformed_distance_maps[j] = space_transformer.transform_image(vessel_seg_images[j], j, 0, true);
      }
      // Generate the montage of the transformed result images by taking
      // the maximum
      final_image = xformed_distance_maps[0];
      for (unsigned int  j = 1; j<xformed_distance_maps.size(); j++) {
        ImageType::Pointer xformed_image = xformed_distance_maps[j];
        RegionConstIterator inputIt1( xformed_image, xformed_image->GetRequestedRegion() );
        RegionIterator outputIt( final_image, final_image->GetRequestedRegion() );
        for ( inputIt1.GoToBegin(), outputIt.GoToBegin(); !inputIt1.IsAtEnd();  
              ++inputIt1, ++outputIt) {
          outputIt.Set( vnl_math_max( outputIt.Get(), inputIt1.Get() ) );
        }
      }
      std::string id = vul_file::strip_extension( anchor_name );
      outfile_name = id+"_vessel_montage.tif";
      writer3D->SetFileName( outfile_name );
      writer3D->SetInput( final_image );
      writer3D->Update();
    }
  }
  return 0;
}
