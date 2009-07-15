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

//: Executable program to assign the average signal to the nuclei
//
//  The usage:
//  
//  assign_signals .xml iba_image Nissl_image  
//
//  Where:
//  .xml            The .xml file containing the segmentation result
//  Image           The image (including path) containing the signals.
//  Start           Starting distance (<=0) is the interior distance away 
//                  from the boundary.   
//  End             Ending distance (>=0) is the exterior distance away from 
//                  the boundary. (If start > end, the entire segmented 
//                  area is considered)  
//  output          Output file name

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegion.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h" 
#include "itkSubtractImageFilter.h"

#include <string>
#include <fstream>
#include <vnl/vnl_math.h>
#include <vnl/vnl_vector_fixed.h>
#include <maciej_seg/maciejSegmentation.h>
#include <maciej_seg/xml_util.h>
#include <maciej_seg/rich_cell.h>
#include <vul/vul_file.h>

typedef itk::Image< unsigned char, 3 > ImageType;
typedef itk::Image< unsigned short, 3 > LabelImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageRegion< 3 >         RegionType;
typedef itk::ImageRegionConstIterator< LabelImageType > LabelConstRegionIteratorType;
typedef itk::ImageRegionConstIterator< ImageType > ConstRegionIteratorType;
typedef itk::ImageRegionIterator< ImageType > RegionIteratorType;
typedef itk::BinaryBallStructuringElement< unsigned char, 3 > StructuringElementType;
typedef itk::GrayscaleDilateImageFilter< ImageType, ImageType, 
                StructuringElementType >  DilateFilterType;
typedef itk::GrayscaleErodeImageFilter< ImageType, ImageType, 
                StructuringElementType >  ErodeFilterType;
typedef itk::SubtractImageFilter< ImageType,ImageType,
                                  ImageType  >  SubFilterType;

int
main( int argc, char* argv[] )
{
  if (argc<6) {
    std::cerr << "Usage: " << argv[0]
              << " .xml signal_image start end output_name "<<std::endl;
  }

  /*
  if (argc<4) {
    std::cerr << "Usage: " << argv[0] 
              << " .xml info_list output_name "<<std::endl;
    std::cerr <<"Info_list contains the following information in each line:"<<std::endl;
    std::cerr <<"\t- Image: The image (including path) containing the signals."<<std::endl;
    std::cerr <<"\t- Starting distance (<=0) is the interior distance away from the boundary."<<std::endl;
    std::cerr <<"\t- Ending distance (>=0) is the exterior distance away from the boundary"<<std::endl;
    return EXIT_FAILURE;
  }
  */

  /*
  // count the number of signal channels
  std::string line_str;
  int count = 0;
  std::ifstream in_file_str( argv[2]);
  if ( !in_file_str ){
    std::cerr<<"Couldn't open "<<argv[2]<<std::endl;
    exit( 0 );
  }
  while ( in_file_str ) {
    std::getline(in_file_str, line_str);
     if (line_str.length() == 0) continue;

     count++;
  }
  in_file_str.close();
  */

  // read the segmentation
  maciejSegmentation maciejseg;
  std::string image_path, image_name;
  std::string output_name =  argv[5] ;
  xml_util_read( "./", argv[1], image_path ,image_name ,maciejseg );
  std::vector<rich_cell::Pointer> const & cells = maciejseg.all_cells();
  int r_interior, r_exterior;
  std::stringstream( argv[3] ) >> r_interior;
  std::stringstream( argv[4] ) >> r_exterior;

  /*
  out_file<<cells.size()<<"\t"<< count <<std::endl;
  for (int i = 0; i<count; i++)
    out_file<<"\t0.0";
  out_file<<"\n";
  */

  //if the output file already exist, record the contents of a vector
  //of strings, and write the string back with the new signal computed
  std::vector< std::string > old_stuff;
  std::string line_str;
  int cell_count, signal_count = 0;
  if ( vul_file::exists(output_name.c_str()) ) {
    std::ifstream in_file_str(output_name.c_str());
    old_stuff.reserve(cells.size());
    std::getline(in_file_str, line_str);
    std::istringstream line_stream(line_str);
    line_stream>>cell_count>>signal_count;
    assert( cell_count == cells.size());
    std::getline(in_file_str, line_str);// to get rid of the first row which only contains 0's.
    while ( in_file_str ) {
      std::getline(in_file_str, line_str);
      if (line_str.length() == 0) continue;
    
      old_stuff.push_back( line_str );
    }
  }
  signal_count++;

  vnl_vector_fixed<int,3>  image_size = maciejseg.image_size();
  std::ofstream out_file_str( output_name.c_str() );
  out_file_str<<cells.size()<<"\t"<<signal_count<<std::endl;
  for (int i = 0; i<signal_count; i++)
    out_file_str<<"\t0";
  out_file_str<<std::endl;

  // Read the signal images
  image_name = argv[2];
  ImageType::Pointer image;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( image_name );
  try {
    reader->Update();
  }
  catch(itk::ExceptionObject& e) {
    vcl_cout << e << vcl_endl;
  }
  image =  reader->GetOutput();

  if (r_exterior < r_interior)
    std::cout<<"Exterior range < interior range: Signal is computed from the entire segmented region"<<std::endl; 

  for (unsigned int i = 0; i<cells.size(); i++) {
    std::cout<<"Cell ID = "<<cells[i]->label_<<std::endl;
    float sum_interior = 0;
    float sum_exterior = 0;
    int count_interior = 0;
    int count_exterior = 0;

    if (r_exterior < r_interior) { // the entire segmented area is taken
      for (unsigned int b = 0; b<cells[i]->all_points_.size(); b++) {
        vnl_vector_fixed< float, 3 > const & pt =  cells[i]->all_points_[b];
        ImageType::IndexType pos;
        pos[0] = pt[0];
        pos[1] = pt[1];
        pos[2] = pt[2];
        sum_interior += image->GetPixel(pos);
      }
      if (!old_stuff.empty()) 
        out_file_str<<old_stuff[i];
      out_file_str<<"\t"<<sum_interior/cells[i]->all_points_.size()<<std::endl;
      continue;
    }
   
    RegionType region = cells[i]->bounding_box_;
    if (r_interior < 0) { //erode the mask 
      // Generate a mask image of the cell region. Erode the region by
      // r_interior
      RegionType::SizeType size = region.GetSize();
      RegionType::IndexType start={{0,0,0}};
      ImageType::Pointer cropped_mask = ImageType::New(); 
      RegionType mask_region;
      mask_region.SetIndex( start );
      mask_region.SetSize( size );
      cropped_mask->SetRegions( mask_region );
      cropped_mask->Allocate();
      cropped_mask->FillBuffer(0);
      LabelConstRegionIteratorType it1( maciejseg.label_image(), region);
      RegionIteratorType it2( cropped_mask, mask_region );
      for (it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2) {
        if (it1.Get() == cells[i]->label_)
          it2.Set( 255 );
      }
      
      ImageType::Pointer eroded_mask;
      ErodeFilterType::Pointer f_erode = ErodeFilterType::New();
      SubFilterType::Pointer f_sub = SubFilterType::New();
      StructuringElementType  structuringElement;
      structuringElement.SetRadius( -r_interior );
      structuringElement.CreateStructuringElement();
      f_erode->SetKernel( structuringElement );
      f_erode->SetInput(cropped_mask);
      f_sub->SetInput1( cropped_mask  );
      f_sub->SetInput2( f_erode->GetOutput() );
      try {
        f_sub->Update();
      }
      catch (itk::ExceptionObject & e) {
        std::cerr << "Exception in SubFilter: " << e << std::endl;
        exit(0);
      }
      eroded_mask = f_sub->GetOutput();
      
      // Sum the signal in the eroded region only
      ConstRegionIteratorType it3( eroded_mask, mask_region );
      ConstRegionIteratorType it4( image, region);
      for (it3.GoToBegin(), it4.GoToBegin(); !it3.IsAtEnd(); ++it1, ++it3, ++it4) {
        if (it3.Get() > 0) {
          sum_interior += it4.Get();
          count_interior ++;
        }
      }
    }
    if (r_exterior > 0) { //dilate the mask
      // enlarge the bounding box by r on each side.
      RegionType::SizeType size = region.GetSize();
      RegionType::IndexType start = region.GetIndex();
      RegionType::IndexType end;
      end[0] = vnl_math_min((int)(start[0]+size[0]+r_exterior), image_size[0]);
      end[1] = vnl_math_min((int)(start[1]+size[1]+r_exterior), image_size[1]);
      end[2] = vnl_math_min((int)(start[2]+size[2]+r_exterior), image_size[2]);
      start[0] = vnl_math_max((int)(start[0]-r_exterior), 0);
      start[1] = vnl_math_max((int)(start[1]-r_exterior), 0);
      start[2] = vnl_math_max((int)(start[2]-r_exterior), 0);
      
      size[0] = end[0] - start[0];
      size[1] = end[1] - start[1];
      size[2] = end[2] - start[2];
      region.SetSize( size );
      region.SetIndex( start );
      
      // Generate a mask image of the region just found. Dilate the
      // region defined by the segmentation by r. 
      ImageType::Pointer cropped_mask = ImageType::New(); 
      RegionType mask_region;
      start[0] = start[1] = start[2] = 0;
      mask_region.SetIndex( start );
      mask_region.SetSize( size );
      cropped_mask->SetRegions( mask_region );
      cropped_mask->Allocate();
      cropped_mask->FillBuffer(0);
      LabelConstRegionIteratorType it1( maciejseg.label_image(), region);
      RegionIteratorType it2( cropped_mask, mask_region );
      for (it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2) {
        if (it1.Get() == cells[i]->label_)
          it2.Set( 255 );
      }
      ImageType::Pointer dilated_mask;
      DilateFilterType::Pointer f_dilate = DilateFilterType::New();
      SubFilterType::Pointer f_sub = SubFilterType::New();
      StructuringElementType  structuringElement;
      structuringElement.SetRadius( r_exterior );
      structuringElement.CreateStructuringElement();
      f_dilate->SetKernel( structuringElement );
      f_dilate->SetInput(cropped_mask);
      f_sub->SetInput1( f_dilate->GetOutput() );
      f_sub->SetInput2( cropped_mask );
      
      try {
        f_sub->Update();
      }
      catch (itk::ExceptionObject & e) {
        std::cerr << "Exception in SubFilter: " << e << std::endl;
          exit(0);
      }
      dilated_mask = f_sub->GetOutput();
        
      // Sum the signal in the dilated region only
      ConstRegionIteratorType it3( dilated_mask, mask_region );
      ConstRegionIteratorType it4( image, region);
      for (it3.GoToBegin(), it4.GoToBegin(); !it3.IsAtEnd(); ++it1, ++it3, ++it4) {
        if (it3.Get() > 0) {
          sum_exterior += it4.Get();
          count_exterior ++;
        }
      }
    }

    // average the interior and exterior signals
    if (!old_stuff.empty()) 
      out_file_str<<old_stuff[i];
    out_file_str<<"\t"<<(sum_interior+sum_exterior)/float(count_interior+count_exterior)<<"\n";
  }
  out_file_str.close();
  
  return 0;
}

/*
    ///////////// old stuff ///////////////
    // Compute the average Iba1 of all cells, valid or not.
    for (unsigned int b = 0; b<cells[i]->all_points_.size(); b++) {
      vnl_vector_fixed< float, 3 > const & pt =  cells[i]->all_points_[b];
      ImageType::IndexType pos;
      pos[0] = pt[0];
      pos[1] = pt[1];
      pos[2] = pt[2];
      sum += Iba_image->GetPixel(pos);
    }
    out_file<<sum/cells[i]->all_points_.size()<<"\t";

    // enlarge the bounding box by r on each side.
    int r = 5;
    RegionType region = cells[i]->bounding_box_;
    RegionType::SizeType size = region.GetSize();
    RegionType::IndexType start = region.GetIndex();
    RegionType::IndexType end;
    end[0] = vnl_math_min((int)(start[0]+size[0]+r), image_size[0]);
    end[1] = vnl_math_min((int)(start[1]+size[1]+r), image_size[1]);
    end[2] = vnl_math_min((int)(start[2]+size[2]+r), image_size[2]);
    start[0] = vnl_math_max((int)(start[0]-r), 0);
    start[1] = vnl_math_max((int)(start[1]-r), 0);
    start[2] = vnl_math_max((int)(start[2]-r), 0);
    
    size[0] = end[0] - start[0];
    size[1] = end[1] - start[1];
    size[2] = end[2] - start[2];
    region.SetSize( size );
    region.SetIndex( start );

    // Generate a mask image of the region just found. Dilate the
    // region defined by the segmentation by r. 
    ImageType::Pointer cropped_mask = ImageType::New(); 
    RegionType mask_region;
    start[0] = start[1] = start[2] = 0;
    mask_region.SetIndex( start );
    mask_region.SetSize( size );
    cropped_mask->SetRegions( mask_region );
    cropped_mask->Allocate();
    cropped_mask->FillBuffer(0);
    LabelConstRegionIteratorType it1( maciejseg.label_image(), region);
    RegionIteratorType it2( cropped_mask, mask_region );
    for (it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2) {
      if (it1.Get() == cells[i]->label_)
        it2.Set( 255 );
    }
    ImageType::Pointer dilated_mask;
    DilateFilterType::Pointer f_dilate = DilateFilterType::New();
    SubFilterType::Pointer f_sub = SubFilterType::New();
    StructuringElementType  structuringElement;
    structuringElement.SetRadius( r );
    structuringElement.CreateStructuringElement();
    f_dilate->SetKernel( structuringElement );
    f_dilate->SetInput(cropped_mask);
    f_sub->SetInput1( f_dilate->GetOutput() );
    f_sub->SetInput2( cropped_mask );
    try
    {
      f_sub->Update();
    }
    catch (itk::ExceptionObject & e)
      {
        std::cerr << "Exception in SubFilter: " << e << std::endl;
        exit(0);
      }
    dilated_mask = f_sub->GetOutput();

    // Sum the signal in the dilated region only
    ConstRegionIteratorType it3( dilated_mask, mask_region );
    ConstRegionIteratorType it4( Nissl_image, region);
    sum = 0;
    for (it3.GoToBegin(), it4.GoToBegin(); !it3.IsAtEnd(); ++it1, ++it3, ++it4) {
      if (it3.Get() > 0) {
        sum += it4.Get();
        count ++;
      }
    }
    out_file<<sum/count<<"\n";
  }

  return 0;
}

*/
