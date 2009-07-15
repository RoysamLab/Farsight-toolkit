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

//: Executable program to update features in the xml file
//
//  The usage:
//
//  update_features .xml 
//
//  where :
//   .xml            The .xml file containing the segmentation result
//
//  Options: 
//   neuron          The image containing the stain for neurons
//   path            The path where images are found

#include <maciej_seg/xml_util.h>
#include <maciej_seg/maciejSegmentation.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include <vul/vul_arg.h>

typedef itk::Image< unsigned char, 3 > ImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;

void separate_dirname_from_filename( const std::string& input,
                                     std::string& dirname,
                                     std::string& filename)
{
  std::string base_name;
  const std::string slash = "\\/";
  std::string::size_type po = input.find_last_of(slash);
  if (po != std::string::npos) {
    filename = input.substr(po+1,filename.size()-po-1);
    dirname = input.substr(0,po+1);
  }
  else {
    filename = input;
    dirname = "";
  }
}

int
main(  int argc, char* argv[] )
{
   vul_arg< vcl_string > arg_xml     ( 0, "The xml file containing the segmentation result." );
   vul_arg< vcl_string > arg_neuron  ( "-neuron", "The name of the image containing the neuronal stain." );
   vul_arg< vcl_string > arg_path    ( "-path", "The path where all images are found." );
   vul_arg_parse( argc, argv );

   std::string path, xml_name;
   separate_dirname_from_filename(arg_xml(), path, xml_name);
   maciejSegmentation::Pointer seg_result = new maciejSegmentation();
   std::string image_path, image_name;
   xml_util_read( path, xml_name, image_path, image_name, *seg_result );
   std::vector<rich_cell::Pointer> const& cells = seg_result->all_cells();

   // read the intensity image
   ImageType::Pointer intensity_image, neuronal_image;
   ReaderType::Pointer reader = ReaderType::New();
   if (arg_path.set()) 
     image_name = arg_path() + image_name;

   reader->SetFileName( image_name );
   try {
     reader->Update();
   }
   catch(itk::ExceptionObject& e) {
     vcl_cout << e << vcl_endl;
   }
   intensity_image =  reader->GetOutput();
   intensity_image->DisconnectPipeline();

   // read the neuronal image if set
   if (arg_neuron.set()) {
     std::string neuronal_name = arg_path()+arg_neuron();
     reader->SetFileName( neuronal_name );
     try {
       reader->Update();
     }
     catch(itk::ExceptionObject& e) {
       vcl_cout << e << vcl_endl;
     }
     neuronal_image =  reader->GetOutput();
     neuronal_image->DisconnectPipeline();
   }

   // Build the binning structure for computation of the nearest
   // neighbor
   double bin_size = 20;
   seg_result->construct_bins(bin_size);

   // go through each cell to update the features
   for (unsigned int i = 0; i<cells.size(); i++) {
    std::cout<<"Cell ID = "<<cells[i]->label_<<std::endl;

    cells[i]->update_volume();

    cells[i]->update_eccentricity_and_radius();

    cells[i]->ave_intensity_ = cells[i]->compute_average_intensity( intensity_image, seg_result->label_image(), 0, -1);

    if (arg_neuron.set()) {
      cells[i]->neuronal_signal_ =  cells[i]->compute_average_intensity( neuronal_image, seg_result->label_image(), -5, 5);
    }

    cells[i]->update_texture( intensity_image );

    std::vector<rich_cell::Pointer>  nearby_cells;
    seg_result->k_nearest_cells( nearby_cells,  cells[i]->center_, 2);
    cells[i]->update_nearest_nbr(nearby_cells);
   }
   xml_util_write( path, arg_xml(), image_path, image_name, *seg_result);
   
   return 0;
}
