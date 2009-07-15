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

//: Executable program to remove duplicated nuclei between adjacent images.
//
//  The input includes the xml file of the joint transformations and
//  the xml file containing the result sets. The changes will be
//  written back to the original xml segmentation files. Class
//  information is ignored in this process, since it is not very
//  reliable.
//
//  nuclear_montage xml_joint_transforms xml_result_sets
// 
//  where
//   xml_joint_transforms   Name of the xml file containing the joint transformations
//   xml_result_sets        Name of the xml file containing the result sets
//
//  Optional argument:
//   -image_before          The montage image without duplication removed
//   -image_after           The montage image with duplication removed
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

#include <vul/vul_file.h>

typedef unsigned char                      InputPixelType;
typedef itk::Image< InputPixelType, 3 >    ImageType;
typedef itk::RGBPixel< unsigned char >     OutputPixelType;
typedef itk::Image< OutputPixelType, 3 >   ColorImageType;
typedef itk::ImageFileWriter< ColorImageType >  WriterType3D;
typedef itk::ImageRegionConstIterator< ImageType > RegionConstIterator;
typedef itk::ImageRegionIterator< ColorImageType > RegionIterator;

ImageType::Pointer
read_image( std::string const & file_name, int channel )
{
  std::cout<<"Reading the image "<<file_name<<std::endl;

  ImageType::Pointer image;

  // Get pixel information
  itk::TIFFImageIO::Pointer io = itk::TIFFImageIO::New();
  io->SetFileName(file_name);
  io->ReadImageInformation();
  int pixel_type = (int)io->GetPixelType();
  std::cout<<"Pixel Type = "<<pixel_type<<std::endl; //1 - grayscale, 2-RGB, 3-RGBA, etc.,

  if (pixel_type == 3) { //RGBA pixel type
    typedef fsc_channel_accessor<itk::RGBAPixel<unsigned char>,3 > ChannelAccType;
    ChannelAccType channel_accessor(file_name);
    image = channel_accessor.get_channel(ChannelAccType::channel_type(channel));
  }
  else if (pixel_type == 2) { //RGA pixel type
    typedef fsc_channel_accessor<itk::RGBPixel<unsigned char>,3 > ChannelAccType;
    ChannelAccType channel_accessor(file_name);
    image = channel_accessor.get_channel(ChannelAccType::channel_type(channel));
  }
  else {// Gray image
    typedef itk::ImageFileReader< ImageType > ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( file_name );
    try {
      reader->Update();
    }
    catch(itk::ExceptionObject& e) {
      vcl_cout << e << vcl_endl;
    }
    image =  reader->GetOutput();
  }
  return image;
}

// The input image is superimposed with the segmentation result
void
superimpose_nuclei(ImageType::Pointer image, 
                   maciejSegmentation::Pointer result,
                   bool show_dup)
{
  std::vector<rich_cell::Pointer> const & cells = result->all_cells();
  typedef std::vector<rich_cell::Pointer>::const_iterator c_iter;
  typedef std::vector< vnl_vector_fixed< float, 3 > >::iterator b_iter;

  for (c_iter ci = cells.begin(); ci != cells.end(); ci++) {
    if ( !(*ci)->valid_ ) continue;
    if (!show_dup && (*ci)->dup_) continue; //duplicated

    for (b_iter bi = (*ci)->boundary_points_.begin(); bi != (*ci)->boundary_points_.end(); bi++) {
      vnl_vector_fixed< float, 3 > const & pt = *bi;
      ImageType::IndexType pos;
      pos[0] = pt[0];
      pos[1] = pt[1];
      pos[2] = pt[2];
      image->SetPixel(pos, 255);
    }
  }

}

ColorImageType::Pointer 
generate_montage(fregl_space_transformer const& space_transformer, 
                 std::vector<maciejSegmentation::Pointer> const& seg_results,
                 int anchor_index,
                 std::string const& image_path, 
                 std::vector<std::string> const& image_names, 
                 int channel, bool show_dup)
{
  // Read in the images
  std::vector<ImageType::Pointer> images;
  for (unsigned int i = 0; i<image_names.size(); i++) {
    std::string full_name = image_path+std::string("/")+image_names[i];
    ImageType::Pointer image = read_image( full_name, channel );
    images.push_back( image );
  }

  // Generate images with the segmented nuclei superimposed
  for (unsigned int i = 0; i<images.size(); i++) {
    superimpose_nuclei(images[i], seg_results[i], show_dup);
  }

  // Transform the images to the global space defined by the anchor
  // image.
  std::vector<ImageType::Pointer> xformed_images;
  for (unsigned int i = 0; i<images.size(); i++) {
    ImageType::Pointer xformed_image = space_transformer.transform_image(images[i], i);
    xformed_images.push_back( xformed_image );
    images[i] = 0; //early release of the memory
  }

  std::cout<<"Number of transformed images "<<xformed_images.size()<<std::endl;

  // Combine the images together. The anchor image is the red channel
  // and all overlapping images are green.
  ColorImageType::Pointer out_image = ColorImageType::New();
  out_image->SetRegions( xformed_images[0]->GetRequestedRegion() );
  out_image->Allocate();
  out_image->FillBuffer(itk::RGBPixel<unsigned char>(itk::NumericTraits<unsigned char>::Zero));

  for (unsigned int i = 0; i<xformed_images.size(); i++) {
    RegionConstIterator inputIt1( xformed_images[i], xformed_images[i]->GetRequestedRegion() );
    RegionIterator outputIt( out_image, out_image->GetRequestedRegion() );

    for ( inputIt1.GoToBegin(), outputIt.GoToBegin(); !inputIt1.IsAtEnd();  ++inputIt1, ++outputIt) {
      OutputPixelType pix = outputIt.Get();
      if (i == anchor_index) pix.SetRed( inputIt1.Get() );
      else pix.SetGreen( vnl_math_max(inputIt1.Get(), pix.GetGreen()) );
      outputIt.Set( pix );
    }
  }

  return out_image;
} 

int 
main( int argc, char* argv[] )
{
  vul_arg< vcl_string > arg_xml_xform  ( 0, "Name of the xml file containing the joint transformations" );
  vul_arg< vcl_string > arg_xml_result ( 0, "Name of the xml file containing the result sets" );
  vul_arg< float > arg_alpha ("-alpha", "Maximum separation between overlapping cells", 10);
  vul_arg< vcl_string > arg_anchor ("-anchor", "Name of the anchor image");
  vul_arg< vcl_string > arg_image_before ("-image_before", "Suffix for the name of the montage image without duplication removed");
  vul_arg< vcl_string > arg_image_after ("-image_after", "Suffix for the name of the montage image with duplication removed");
  vul_arg< vcl_string > arg_image_path  ("-image_path", "Path of images", "./");
  vul_arg< int > arg_channel ("-channel", "The color channel to extract", 0);
  vul_arg< bool > arg_in_anchor ("-in_anchor", "The global space is the anchor image", false);

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
  std::vector<std::string> const & image_names = joint_register->image_names();
  
  // Take each image as the anchor in turn. Load the xml files of the
  //anchor and the overlapping images. Mark duplicated nuclei and
  //write the xml files back to disk. The class information is
  //ignored.
  //
  for (unsigned int i = 0; i<image_names.size(); i++) {
    std::string anchor_name = image_names[i];
    if (arg_anchor.set() && anchor_name!=arg_anchor()) continue;

    std::cout<<"anchor image = "<<anchor_name<<std::endl;
    space_transformer.set_anchor( anchor_name, arg_in_anchor(), overlap_only );
    
    // load the xml files of the anchor and the overlapping images
    std::vector<std::string> overlapping_image_names = space_transformer.image_names();
    std::vector<maciejSegmentation::Pointer> nuclear_seg_results(overlapping_image_names.size());
    int anchor_index;

    std::vector<std::string> xml_nuclear_names;
    std::vector<std::string> xml_nuclear_image_paths;
    std::vector<std::string> xml_nuclear_image_names;
    std::string  image_path;
    std::string  image_name;
    for (unsigned int j = 0; j<overlapping_image_names.size(); j++) {

      for (unsigned int k = 0; k<result_records.size(); k++) { 
        if (overlapping_image_names[j] == anchor_name) anchor_index = j;

        if (result_records[k]->registration_image() == overlapping_image_names[j]) {
          std::cout<<"read xml for "<<result_records[k]->registration_image()<<std::endl;
          // assuming everything to be in the current directory
          nuclear_seg_results[j] = new maciejSegmentation();
          xml_util_read( "./", result_records[k]->nuclear_xml(),
                         image_path ,image_name ,*nuclear_seg_results[j] );
          xml_nuclear_names.push_back(result_records[k]->nuclear_xml());
          xml_nuclear_image_paths.push_back( image_path );
          xml_nuclear_image_names.push_back( image_name );
        }
      }
    }

    // Repeat the same process for each class. First, scan the cell
    // list of the anchor image to find the max class number.
    //
    maciejSegmentation::Pointer anchor_nuclear_seg = nuclear_seg_results[anchor_index];
    std::vector<rich_cell::Pointer> const & cells_in_anchor = anchor_nuclear_seg->all_cells();
        
    int max_class_num = 0; 
    for (unsigned int i = 0; i<cells_in_anchor.size(); i++) {
      if (cells_in_anchor[i]->class_type_ > max_class_num) 
        max_class_num = cells_in_anchor[i]->class_type_;
    }

    //for (int ci=1; ci<=max_class_num; ci++) {

    // Build the binning structure for adjacent images
    double bin_size = 20;
    //int class_label = ci;

    for (unsigned int j = 0; j<overlapping_image_names.size(); j++) {
        if (j != anchor_index)
          //nuclear_seg_results[j]->construct_bins(bin_size, class_label);
          nuclear_seg_results[j]->construct_bins(bin_size);
      }

    // Mark duplications
    for (std::vector<rich_cell::Pointer>::const_iterator ci = cells_in_anchor.begin(); ci != cells_in_anchor.end(); ci++ ) {
      //if (!(*ci)->valid_ || (*ci)->class_type_ != class_label) continue;
      if ( !(*ci)->valid_ ) continue;

      //if already duplicated, skip
      if ( (*ci)->dup_ ) continue;
      
      // the cell is valid and distinct (at the moment)
      for (unsigned int j = 0; j<overlapping_image_names.size(); j++) {
        if (j == anchor_index) continue; 
        
        vnl_vector_fixed< float, 3 > xformed_loc;
        if (!space_transformer.in_image((*ci)->center_, j, xformed_loc)) continue;
        // find nuclei Ni which are in circular region of radius =
        // alpha. If c is larget than all Ni, c is marked as
        // "duplicated", else all Ni are marked as "duplicated".
        //
        
        /*
        // compute the radius for the search range
        double radius = 0;
        double dist;
        for (rich_cell::PointVectorType::iterator bi = (*ci)->boundary_points_.begin(); bi != (*ci)->boundary_points_.end(); bi++) {
        dist = (*bi - (*ci)->center_).magnitude();
        if (radius < dist) radius = dist;
        }
        
        //  Get all nuclei in the search range of 3xc->radius_variance
        std::vector<rich_cell::Pointer>  nearby_cells;
        nuclear_seg_results[j]->cells_within_radius(nearby_cells, xformed_loc,
        radius);
        
        // refine the list again to get only nuclei for witch |c-Ni|
        // is less than the alpha distance
        std::vector<rich_cell::Pointer>  nearby_cells_refined;
        for (std::vector<rich_cell::Pointer>::iterator cj = nearby_cells.begin(); cj != nearby_cells.end(); cj++) {
        dist = (xformed_loc-(*cj)->center_).magnitude();
        if (dist < alpha) nearby_cells_refined.push_back(*cj);
        }
        if (!nearby_cells_refined.empty())
        std::cout<<"Nearby cells found!"<<std::endl;
        */
        
        std::vector<rich_cell::Pointer>  nearby_cells;
        std::vector<rich_cell::Pointer>  nearby_cells_refined;
        nuclear_seg_results[j]->cells_within_radius(nearby_cells, xformed_loc,
                                                    arg_alpha());
        
        // Remove cells marked as "dup" and belong to different classes
        for (std::vector<rich_cell::Pointer>::iterator cj = nearby_cells.begin(); cj != nearby_cells.end(); cj++) {
          if ((*cj)->dup_ )  continue;
          nearby_cells_refined.push_back(*cj);
        }
        
        // If c is larger than all nuclei in nearby_cells_refined, c
        // is marked as "duplicated", else all nuclei in
        // nearby_cells_refined are marked as "duplicated".
        bool duplicated = false;
        for (std::vector<rich_cell::Pointer>::iterator cj = nearby_cells_refined.begin(); cj != nearby_cells_refined.end(); cj++) {
          if ((*cj)->volume_ > (*ci)->volume_) {
            duplicated = true;
            break;
          }
        }
        
        edit_record e_record;
        e_record.status_ = edit_record::DUP;
        if (duplicated) { // (*ci) is duplicated  
          (*ci)->edit_history_.push_back( e_record );
          (*ci)->dup_ = true;
        }
        else {//others are duplicated
          for (std::vector<rich_cell::Pointer>::iterator cj = nearby_cells_refined.begin(); cj != nearby_cells_refined.end(); cj++) {
            (*cj)->edit_history_.push_back( e_record );
            (*cj)->dup_ = true;
          }
        }
        
        if (duplicated) break;
      }
    }
    //} //end of classes
    
    
    // Output the xml files again
    std::cout<<"Update the xml files ..."<<std::endl;
    for (unsigned int j = 0; j < nuclear_seg_results.size(); j++) {
      xml_util_write("./", xml_nuclear_names[j], xml_nuclear_image_paths[j],
                     xml_nuclear_image_names[j], *(nuclear_seg_results[j]) );
    }

    // Output images if necessary
    if (arg_image_before.set()) {
      std::cout<<"Output the image before removal of duplications ..."<<std::endl;
      bool show_dup = true;
      ColorImageType::Pointer out_image;
      out_image = generate_montage(space_transformer, nuclear_seg_results, 
                                   anchor_index, arg_image_path(), 
                                   xml_nuclear_image_names, arg_channel(), 
                                   show_dup); 
      std::string anchor_image_id =vul_file::strip_extension( xml_nuclear_image_names[anchor_index] );
      WriterType3D::Pointer writer3D = WriterType3D::New();
      writer3D->SetFileName( anchor_image_id+arg_image_before()+".tif" );
      writer3D->SetInput( out_image );
      writer3D->Update();
    }
    if (arg_image_after.set()) {
      std::cout<<"Output the image after removal of duplications ..."<<std::endl;
      bool show_dup = false;
      ColorImageType::Pointer out_image;
      out_image = generate_montage(space_transformer, nuclear_seg_results, 
                                   anchor_index, arg_image_path(), 
                                   xml_nuclear_image_names, arg_channel(), 
                                   show_dup); 
      std::string anchor_image_id =vul_file::strip_extension( xml_nuclear_image_names[anchor_index] );
      WriterType3D::Pointer writer3D = WriterType3D::New();
      writer3D->SetFileName( anchor_image_id+arg_image_after()+".tif" );
      writer3D->SetInput( out_image );
      writer3D->Update();
    }
  }
  return 0;
}
