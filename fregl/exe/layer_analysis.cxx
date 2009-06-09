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

//: Executable program to generate an output showing a layer using certain criteria
//
//  The input is a file containing the xml file of the joint
//  transformations and the xml file containing the result sets.
//  
//  layer_analysis xml_joint_transforms xml_result_sets
//
//  where
//   xml_joint_transforms   Name of the xml file containing the joint transformations
//   xml_result_sets        Name of the xml file containing the result sets
//

#include <vul/vul_arg.h>

#include <fregl/fregl_result_record.h>
#include <fregl/fregl_joint_register.h>
#include <fregl/fregl_space_transformer.h>

#include <maciej_seg/maciejSegmentation.h>
#include <maciej_seg/xml_util.h>

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkRGBPixel.h"

typedef itk::RGBPixel< unsigned char >   colorPixelType;
typedef itk::Image< colorPixelType, 2 > ImageType2D;
typedef itk::Image< colorPixelType, 3 > ImageType3D;
typedef itk::ImageFileWriter< ImageType2D >  WriterType2D;

int 
main( int argc, char* argv[] )
{
  vul_arg< vcl_string > arg_xml_xform  ( 0, "Name of the xml file containing the joint transformations" );
  vul_arg< vcl_string > arg_xml_result ( 0, "Name of the xml file containing the result sets" );
  vul_arg< vcl_string > arg_anchor (0, "Name of the anchor image. If not set, the image of the first result file is taken as the first"); 
  vul_arg< int > arg_class ("-class", "class number");
  vul_arg< float > arg_thrd ("-threshold"," threshold value for a layer",20); 
  vul_arg< vcl_string > arg_output ("-output", "Name of the output image","layer.tiff");

  vul_arg_parse( argc, argv );

  // Read in the xml files
  fregl_joint_register::Pointer joint_register = new fregl_joint_register( arg_xml_xform() );
  std::vector<fregl_result_record::Pointer> result_records;
  result_record_read_xml(arg_xml_result(), result_records);

  // Set the space transformer
  //
  fregl_space_transformer space_transformer(joint_register);
  space_transformer.set_anchor( arg_anchor() );
  std::vector<std::string> image_names = space_transformer.image_names();

  //Read in the xml files for nuclear segmentation. The files are
  //stored in the order corresponding to image_names
  std::string  image_path;
  std::string  image_name;
  std::vector<maciejSegmentation::Pointer> nuclear_seg_results(image_names.size());
  for (unsigned int img_ind = 0; img_ind<image_names.size(); img_ind++) {
    bool found = false;
    //if (image_names[img_ind] == arg_anchor.set()) anchor_index = img_ind;
    for (unsigned int k = 0; k<result_records.size(); k++) {  
      if (result_records[k]->registration_image() == image_names[img_ind]) {
        std::cout<<"read xml for "<<result_records[k]->registration_image()<<std::endl;
        // assuming everything to be in the current directory
        nuclear_seg_results[img_ind] = new maciejSegmentation();
        xml_util_read( "./", result_records[k]->nuclear_xml(),
                       image_path ,image_name ,*nuclear_seg_results[img_ind] );
        found = true;
        break;
      }
    }
    if (!found) {
      std::cerr<<"Image "<<image_names[img_ind]<<" has no segmentation result!"<<std::endl;
      return 1;
    }
  }

  // Prepare the output image
  ImageType2D::RegionType region;
  ImageType3D::RegionType::SizeType size3d;
  ImageType2D::RegionType::SizeType size2d;
  ImageType2D::RegionType::IndexType index;
  ImageType3D::PointType origin;
  ImageType2D::PointType origin2d;

  size3d = space_transformer.montage_size();
  origin = space_transformer.origin();
  size2d[0] = size3d[0];
  size2d[1] = size3d[1];
  index[0] = 0;
  index[1] = 0;
  region.SetIndex( index );
  region.SetSize(size2d);

  ImageType2D::Pointer image = ImageType2D::New();
  image->SetRegions( region );
  image->Allocate();
  image->FillBuffer(itk::RGBPixel<unsigned char>(255));
  origin2d[0] = origin[0];
  origin2d[1] = origin[1];
  image->SetOrigin( origin2d );
  colorPixelType black, pink;
  black.SetRed( 0 );
  black.SetGreen( 0 );
  black.SetBlue( 0 );
  pink.SetRed(255);
  pink.SetGreen( 0 );
  pink.SetBlue( 255 );

  // For each segmentation result, transform the boundary of each cell
  // from the z-slice with the maximum perimeter. If the cell belongs to
  // the given class and the nnd attribute is less than a threshold,
  // color the interior in pink.
  for (unsigned int i = 0; i<nuclear_seg_results.size(); i++) {
    std::vector<rich_cell::Pointer> const & cells = nuclear_seg_results[i]->all_cells();

    for (unsigned int ci = 0; ci<cells.size(); ci++) {
      if (!cells[ci]->valid_ || cells[ci]->dup_) continue;
      if ( arg_class.set() && cells[ci]->class_type_ != arg_class() ) continue;

      // Find the mid-z. The boundary points are first stored in an 2D
      // vector array, where the index for the first dimension is the
      // displacement of the z position from the min_z.
      int min_z = cells[ci]->bounding_box_.GetIndex()[2];
      int max_z = min_z + cells[ci]->bounding_box_.GetSize()[2];
      int z_size = cells[ci]->bounding_box_.GetSize()[2];
      rich_cell::PointVectorType const & b_points = cells[ci]->boundary_points_;
      std::vector< std::vector<rich_cell::PointType> > z_slices(z_size);
      for (unsigned int bi = 0; bi<b_points.size(); bi++) {
        z_slices[b_points[bi][2]-min_z].push_back( b_points[bi] );
      }
      
      int best_z = z_size/2;
      /*
      int best_z = 0;
      for (unsigned int zi = 1; zi<z_slices.size(); zi++) {
        if (z_slices[best_z].size() < z_slices[zi].size()) best_z = zi;
      }
      */

      // Transform all the points in z_slices[best_z] list to the
      // montage space
      for (unsigned int bi = 0; bi<z_slices[best_z].size(); bi++) {
        //check if the point is the first or the last of its x-scan line
        bool is_x_end = true;
        for (unsigned int xi = 0; xi<z_slices[best_z].size(); xi++) {
          if (bi!=xi && z_slices[best_z][xi][1] == z_slices[best_z][bi][1])
            if (z_slices[best_z][xi][0] > z_slices[best_z][bi][0]) {
              is_x_end = false;
              break;
            }
        }
        if (!is_x_end) {
          is_x_end = true; 
          for (unsigned int xi = 0; xi<z_slices[best_z].size(); xi++) {
          if (bi!=xi && z_slices[best_z][xi][1] == z_slices[best_z][bi][1])
            if (z_slices[best_z][xi][0] < z_slices[best_z][bi][0]) {
              is_x_end = false;
              break;
            }
          }
        }

        //check if the point is the first or the last of its y-scan line
        bool is_y_end = true;
        for (unsigned int yi = 0; yi<z_slices[best_z].size(); yi++) {
          if (bi!=yi && z_slices[best_z][yi][0] == z_slices[best_z][bi][0])
            if (z_slices[best_z][yi][1] > z_slices[best_z][bi][1]) {
              is_y_end = false;
              break;
            }
        }
        if (!is_y_end) {
          is_y_end = true; 
          for (unsigned int yi = 0; yi<z_slices[best_z].size(); yi++) {
          if (bi!=yi && z_slices[best_z][yi][0] == z_slices[best_z][bi][0])
            if (z_slices[best_z][yi][1] < z_slices[best_z][bi][1]) {
              is_y_end = false;
              break;
            }
          }
        }

        if (!is_x_end && !is_y_end) continue;

        ImageType3D::PointType loc;
        loc[0] = z_slices[best_z][bi][0];
        loc[1] = z_slices[best_z][bi][1];
        loc[2] = z_slices[best_z][bi][2];
        ImageType3D::PointType xformed_loc;
        space_transformer.in_anchor(loc,i,xformed_loc);

        ImageType2D::PointType xformed_loc_2d;
        ImageType2D::IndexType xformed_index;
        xformed_loc_2d[0] = xformed_loc[0];
        xformed_loc_2d[1] = xformed_loc[1];

        image->TransformPhysicalPointToIndex( xformed_loc_2d, xformed_index );
        image->SetPixel( xformed_index, black );
      }

      // If the cell meet the layer criterion, color the interior
      if ( cells[ci]->ave_nnbr_dist_ > arg_thrd() ) continue;

      rich_cell::PointVectorType const & i_points = cells[ci]->all_points_;;
      for (unsigned int ii = 0; ii<i_points.size(); ii++) {
        if (i_points[ii][2] == min_z+best_z) {
          ImageType3D::PointType loc;
          loc[0] = i_points[ii][0];
          loc[1] = i_points[ii][1];
          loc[2] = i_points[ii][2];
          ImageType3D::PointType xformed_loc;
          space_transformer.in_anchor(loc,i,xformed_loc);

          ImageType2D::PointType xformed_loc_2d;
          ImageType2D::IndexType xformed_index;
          xformed_loc_2d[0] = xformed_loc[0];
          xformed_loc_2d[1] = xformed_loc[1];

          image->TransformPhysicalPointToIndex( xformed_loc_2d, xformed_index );
          colorPixelType pix = image->GetPixel( xformed_index );
          if (pix.GetRed() == 0) continue;

          // color the pixel and its neighboring 4 pixels if not
          // already in black
          image->SetPixel( xformed_index, pink );
          ImageType2D::IndexType xformed_index_nb;
          xformed_index_nb = xformed_index;
          xformed_index_nb[0] = xformed_index[0]+1;
          pix = image->GetPixel( xformed_index_nb);
          if (pix.GetRed() == 255)
            image->SetPixel( xformed_index_nb, pink );

          xformed_index_nb = xformed_index;
          xformed_index_nb[0] = xformed_index[0]-1;
          pix = image->GetPixel( xformed_index_nb);
          if (pix.GetRed() == 255)
            image->SetPixel( xformed_index_nb, pink );

          xformed_index_nb = xformed_index;
          xformed_index_nb[1] = xformed_index[1]+1;
          pix = image->GetPixel( xformed_index_nb);
          if (pix.GetRed() == 255)
            image->SetPixel( xformed_index_nb, pink );

          xformed_index_nb = xformed_index;
          xformed_index_nb[1] = xformed_index[1]+1;
          pix = image->GetPixel( xformed_index_nb);
          if (pix.GetRed() == 255)
            image->SetPixel( xformed_index_nb, pink );
        }
      }
    }
  }

  // Output the 2D projected image
  WriterType2D::Pointer writer2D = WriterType2D::New();
  writer2D->SetFileName( arg_output() );
  writer2D->SetInput( image );
  writer2D->Update();
  
  
  return 0;
}
