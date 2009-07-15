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

#include "maciejSegmentation.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <time.h>	//time execution
#include <vnl/algo/vnl_determinant.h>
#include <stdlib.h>	//for wriring/reading data
#include <stdlib.h>	//printf
#include <sstream>	//stringstream

#include "filter.h"
#include "cell_mgr.h"
#include "xml_util.h"
#include "rich_cell.h"
#include "nucleus_feature_util.h"
#include "fnsl3d_filter.h"

#include <string>

#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkImageFileReader.h"

#include "itkImageAdaptor.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkCastImageFilter.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkFlipImageFilter.h"

//#include <rsdl/rsdl_bins.h>
//#include <rsdl/rsdl_bins.txx> // to avoid explicit instantiation

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

/********************************************************************************/
/*Desc- main program for running nuclear seg, the program takes 19 inputs args. The inputs
/*are described below. The inputs are sets of consecitive numered images that are the slices
/*of cells. The images are loaded, filtered, using ITK, segmented, then merged and a final score is 
/*computd for the segmentation. There are options to seg the image form data already available or 
/*auto train.
/*
/*	19 input arguments
/*	1 - path - path to folder containing images
/*	2 - filename - name of file name of the image stacks, only the name without the numbering
/*	3 - type - type of image, ei extension tif
/*	4 - sstack - starting number of images ie image_nameXX.tif XX->sstack
/*	5 - estack - ending stack number
/*	6 - swidth - width of number ie image_nameXX.tif swidth -> 2 , image_nameXXX.tif swidth -> 3
/*	7 - pT - threshold parameter, gets divided by 1000, int needed for RRS
/*	8 - g - grid size parameter
/*	9 - rM - radius of median filter
/*	10 - rG - radius of morphological operators
/*	11 - min - should always be 0	(param for watershed)
/*	12 - max - should always be 255 (param for watershed)
/*	13 - step - param for step size of watershed, should be 1
/*	14 - gliadat - path and filename for glia train data
/*	15 - neurdat - path and filename for neuron train data
/*	16 - channel - channel of nuclei stain R-0 G-1 B-2
/*	17 - mode - 0 -> train from data, 1 -> autotrain
/*	18 - loadimages - 1 -> load preprocessed images, 0-> filter images, always filter images, preprocessed sometimes are incorrect
/*	19 - saveoutput - output segmented image->1 , no output->0
/********************************************************************************/

maciejSegmentation::
maciejSegmentation()
  :cells_valid_(false), cell_count_(0), bins_set_(false)
{} 

void
maciejSegmentation::
run( std::vector<std::string> const& parameters)
{
  typedef itk::RGBPixel<unsigned char>   PixelType;
  typedef itk::Image< PixelType, 3 >   ImageType;

  ImageType::Pointer proj_img, seg_result_img;

  run( parameters, proj_img, seg_result_img);
}

void
maciejSegmentation::
run( std::vector<std::string> const& parameters,
     filter::mInputImageFileType::Pointer & proj_img,
     filter::mInputImageFileType::Pointer & seg_result_img)
{
  assert(parameters.size()==14);

  int g = 1 ,		//amount of grids
    rM = 1 ,	//median radius
    rG = 1 ;	//morph radius
    short int min = 0 , 
    max = 254 , 
    step = 1;
  double  pT = 161;
  std::string temp = "_filt.tif";
  
  int sstack = 0, //starting stack
    estack = 1,	//end stack
    swidth = 2,	//width of numbers in filename
    channel = 0,	//channel of nuclei stain R-0 G-1 B-2
    mode = 0,	//0 -> train from data, 1 -> autotrain
    loadimages = 0,	//1 -> load preprocessed images, 0-> filter images
    saveoutput = 1;	//output segmented image->1 , no output->0
  
  std::string inputs , path , filename , type, gliadat, neurdat;
  
  std::ostringstream oss;
  
  path =  parameters[0];
  std::cout<<"path = "<<path<<std::endl;
  
  filename = parameters[1];
  std::cout<<"filename = "<<filename<<std::endl;

  type = parameters[2];
  std::cout<<"type = "<<type<<std::endl;
	
  inputs = parameters[3];
  std::stringstream( inputs ) >> sstack;
  std::cout<<"sstack = "<<sstack<<std::endl;

  inputs = parameters[4];
  std::stringstream( inputs ) >> estack;
  std::cout<<"estack = "<<estack<<std::endl;
	
  inputs = parameters[5];
  std::stringstream( inputs ) >> swidth;
  std::cout<<"swidth = "<<swidth<<std::endl;

  inputs = parameters[6];
  std::stringstream( inputs ) >> pT;
  std::cout<<"pT = "<< pT<<std::endl;

  inputs = parameters[7];
  std::stringstream( inputs ) >> g;
  std::cout<<"grid = "<< g<<std::endl;

  inputs = parameters[8];
  std::stringstream( inputs ) >> rM;
  std::cout<<"rM = "<< rM<<std::endl;

  inputs = parameters[9];
  std::stringstream( inputs ) >> rG;
  std::cout<<"rG = "<< rG<<std::endl;

  /*
  inputs = parameters[10];
  std::stringstream( inputs ) >> min;
  std::cout<<"min = "<< min<<std::endl;

  inputs = parameters[11];
  std::stringstream( inputs ) >> max;
  std::cout<<"max = "<< max<<std::endl;

  inputs = parameters[12];
  std::stringstream( inputs ) >> step;
  std::cout<<"step = "<< step<<std::endl;
  */

  gliadat = parameters[10];
  std::cout<<"gliadat = "<< gliadat<<std::endl;

  neurdat = parameters[11];
  std::cout<<"neurdat = "<<neurdat<<std::endl;

  inputs = parameters[12];
  std::stringstream( inputs ) >> channel;
  std::cout<<"channel = "<<channel<<std::endl;


  inputs = parameters[13];
  std::stringstream( inputs ) >> mode;	//0 -> train from data, 1 -> autotrain
  std::cout<<"mode= "<<mode<<std::endl;

  /*
  inputs = parameters[14];
  std::stringstream( inputs ) >> loadimages;
  std::cout<<"loadimages= "<<loadimages<<std::endl;

  inputs = parameters[15];
  std::stringstream( inputs ) >> saveoutput;
  std::cout<<"saveoutput= "<<saveoutput<<std::endl;
  */

  if( !(mode == 0 || mode == 1) ) {
    std::cerr << mode << " Wrong mode" << std::endl;
    exit(0);
  }
  
  parameters_ = parameters;

  //norm threshold param
  pT /= 1000.0f;

  oss << '_' << std::setw(2) << std::setfill('0') << g << '_' << std::setw(2) << std::setfill('0') << rM << '_' << std::setw(2) << std::setfill('0') << rG << '_' << std::setw(6) << std::setfill('0') << pT;
  
  std::cout << "Starting filter " << g << ' ' << rM << ' ' << pT << std::endl;
  filter f(g,g,rM,rG,pT,0.5,channel);	//setup filter to prefilter image
  
  if( loadimages == 1 )
    {
      std::cout << "Loading Images" << std::endl;
      f.load_images(std::string( "_filt" + oss.str() +".tif"),std::string( "_gradient" + oss.str() +".tif"),std::string( "_intensity" + oss.str() +".tif"),path,filename,"tif",sstack,estack,swidth,loadimages);
    }
  if( loadimages == 0 )
    {
      std::cout << "Filtering Images" << std::endl;
      f.run_filter(filename,path,type,sstack,estack,swidth);		//run filter with path to images and size of stack
      //never do this
      if( saveoutput == 10 )
        {
          
          temp = "_filt" + oss.str() + ".tif";
          f.save_image(f.getfilt(),temp);
          
          temp = "_intensity" + oss.str() + ".tif";
          f.save_image(f.getintensity(),temp);
          
          temp = "_gradient" + oss.str() + ".tif";
          f.save_image(f.getgradient(),temp);
        }
    }
  
  std::cout << "Starting cell" << std::endl;
  //start cell
  cell_mgr c(f.getsize(),f.getfilt(),f.getintensity(),f.getgradient());
  std::string a= oss.str();
  c.setup_cells(min,max,step,a);	//create cells from labeled image
  
  if( mode == 0 )
    {
      std::cout << "Getting train data" << std::endl;
      c.get_data(gliadat,neurdat);
    }
  else
    {
      std::cout << "Calc train data" << std::endl;
      c.calc_train_data( 0.95 , 1.0 , 0.9 , 0.2 , 1000 , 8000);
      
      //merge traning cells
      std::cout << "Build train rag" << std::endl;
      c.build_rag_train();	
      //get rid of cells touching border
      c.remove_border_cells();
      
      std::cout << "Classifying" << std::endl;
      c.classify_cells();	
      
      std::cout << "Calculating train data" << std::endl;
      //save data into matrix
      c.calc_ng_data();
      //only keep good cells as training cells
      //c.calc_mean_vect(2.0/3.0);
      //save good cells into matrix
      //c.calc_ng_data();
      
      //reset cell class labels
      std::cout << "Relabeling Orig Cells" << std::endl;
      c.relabel_orig_cells();		
      
    }
  
  //merge and classify cells from train data
  std::cout << "Building 2nd rag" << std::endl;
  a = oss.str();
  c.build_rag(mode,a);	//merge cells
  
  std::cout << "Removing border Cells" << std::endl;
  c.remove_border_cells();
  
  if( mode == 1)
	{
          std::cout << "Cleaning up data" << std::endl;
          //c.save_train_data();	
	}	
  
  c.calc_ng_data();
  
  //get rid of cells not touching median plane
  c.remove_non_median_cells();
  
  //correct any mis classfied cells
  std::cout << "Correct Errors" << std::endl;
  a= oss.str();
  c.build_rag(mode,a,2);
  
  //relabel cells with consecutive label numbers
  std::cout << "Relabeling" << std::endl;
  c.relabel_cells(true);
  c.calc_ng_data();
  
  //score the segmentation
  std::cout << "Scoring" << std::endl;
  c.score_foreground_mah( oss.str() );
  
  c.check_all_cells(0,100,100,1,200);

  if( saveoutput == 1 )
    {
      std::cout << "Labeling and Saving" << std::endl;
      //label border and number cells
      f.label_bound(c.get_cells());
      f.label_image(c.get_cells());
      //save image
      //f.save_image_rgb(std::string("_out_" + oss.str() + ".tif"));
      seg_result_img = f.save_image_rgb(std::string(filename + "_out_.tif"));
      //save cell features
      c.output_cells(std::string( filename + ".txt" ));
      //std::vector<rich_cell::Pointer> rich_cells, rich_cells_in;
      vnl_vector_fixed<int,3> image_size;
      image_size.set( c.mImageSize );
      convert_cells_to_rich_cells(c.mvCell, image_size);
      cells_valid_ = true;

      //xml_util_write("cells.xml", rich_cells, parameters, c.mImageSize);
      //xml_util_read("cells_demo.xml", rich_cells_in);
      //xml_util_write("cells_demo_2.xml", rich_cells_in, parameters, c.mImageSize);
      //project cells
      std::cout << "Projecting cells" << std::endl;
      c.project_cells();
      //save projected image
      std::cout << "Saving projected" << std::endl;
      a = filename + "_proj_.tif";
      proj_img = f.project_image(a,c.get_cells());
      
      //save the intensity image
      filter::mOutputImageType::Pointer img = f.getintensity();
      typedef itk::CastImageFilter<filter::mOutputImageType, 
        InternalIntensityImageType> castFilterType;
      castFilterType::Pointer caster = castFilterType::New();
      caster->SetInput(img);
      //intensity_image_ = caster->GetOutput();
      
    }
}

void
maciejSegmentation::
run( InternalIntensityImageType::Pointer image,
     std::vector<std::string> const& parameters)
{
  assert(parameters.size()==8);

  int g = 1 ,		//amount of grids
    rM = 1 ,	//median radius
    rG = 1 ;	//morph radius
    short int min = 0 , 
    max = 254 , 
    step = 1;
  double  pT = 161;
  std::string temp = "_filt.tif";
  
  int sstack = 0, //starting stack
    estack = 1,	//end stack
    swidth = 2,	//width of numbers in filename
    channel = 0,	//channel of nuclei stain R-0 G-1 B-2
    mode = 0,	//0 -> train from data, 1 -> autotrain
    loadimages = 0,	//1 -> load preprocessed images, 0-> filter images
    saveoutput = 1;	//output segmented image->1 , no output->0
  
  std::string inputs , path , filename , type, gliadat, neurdat;
  
  std::ostringstream oss;
  
  /*
  path =  parameters[0];
  std::cout<<"path = "<<path<<std::endl;
  
  filename = parameters[1];
  std::cout<<"filename = "<<filename<<std::endl;

  type = parameters[2];
  std::cout<<"type = "<<type<<std::endl;
	
  inputs = parameters[3];
  std::stringstream( inputs ) >> sstack;
  std::cout<<"sstack = "<<sstack<<std::endl;

  inputs = parameters[4];
  std::stringstream( inputs ) >> estack;
  std::cout<<"estack = "<<estack<<std::endl;
	
  inputs = parameters[5];
  std::stringstream( inputs ) >> swidth;
  std::cout<<"swidth = "<<swidth<<std::endl;
  */
  
  //image_id_ = parameters[0];
  
  inputs = parameters[0];
  std::stringstream( inputs ) >> channel;
  std::cout<<"channel = "<<channel<<std::endl;

  inputs = parameters[1];
  std::stringstream( inputs ) >> pT;
  std::cout<<"pT = "<< pT<<std::endl;

  inputs = parameters[2];
  std::stringstream( inputs ) >> g;
  std::cout<<"grid = "<< g<<std::endl;

  inputs = parameters[3];
  std::stringstream( inputs ) >> rM;
  std::cout<<"rM = "<< rM<<std::endl;

  inputs = parameters[4];
  std::stringstream( inputs ) >> rG;
  std::cout<<"rG = "<< rG<<std::endl;

  /*
  inputs = parameters[10];
  std::stringstream( inputs ) >> min;
  std::cout<<"min = "<< min<<std::endl;

  inputs = parameters[11];
  std::stringstream( inputs ) >> max;
  std::cout<<"max = "<< max<<std::endl;

  inputs = parameters[12];
  std::stringstream( inputs ) >> step;
  std::cout<<"step = "<< step<<std::endl;
  */
  
  inputs = parameters[5];
  std::stringstream( inputs ) >> mode;	//0 -> train from data, 1 -> autotrain
  std::cout<<"mode= "<<mode<<std::endl;

  gliadat = parameters[6];
  std::cout<<"gliadat = "<< gliadat<<std::endl;

  neurdat = parameters[7];
  std::cout<<"neurdat = "<<neurdat<<std::endl;

  /*
  inputs = parameters[12];
  std::stringstream( inputs ) >> channel;
  std::cout<<"channel = "<<channel<<std::endl;
  */

 

  /*
  inputs = parameters[14];
  std::stringstream( inputs ) >> loadimages;
  std::cout<<"loadimages= "<<loadimages<<std::endl;

  inputs = parameters[15];
  std::stringstream( inputs ) >> saveoutput;
  std::cout<<"saveoutput= "<<saveoutput<<std::endl;
  */

  if( !(mode == 0 || mode == 1) ) {
    std::cerr << mode << " Wrong mode" << std::endl;
    exit(0);
  }
  
  parameters_ = parameters;

  //norm threshold param
  pT /= 1000.0f;
  
  /*
  oss << '_' << std::setw(2) << std::setfill('0') << g << '_' << std::setw(2) << std::setfill('0') << rM << '_' << std::setw(2) << std::setfill('0') << rG << '_' << std::setw(6) << std::setfill('0') << pT;
  */

  std::cout << "Starting filter " << g << ' ' << rM << ' ' << pT << std::endl;
  //filter f(g,g,rM,rG,pT,0.5,channel);	//setup filter to prefilter image
  fnsl3d_filter f;

  /*
  std::cout << "Filtering Images" << std::endl;
  typedef itk::CastImageFilter<InternalIntensityImageType, 
    filter::mOutputImageType> castFilterType;
  castFilterType::Pointer caster = castFilterType::New();
  caster->SetInput(image);
  //f.run_filter(caster->GetOutput(), rM, rG, g, g, pT);
  */		
  //f.run_filter(imagename);
  f.run_filter(image, rM, rG, g, g, pT);

  /*
  f.save_image( f.get_filt(), "filt_image.tif");
  f.save_image( f.get_intensity(), "intensity_image.tif");
  f.save_image( f.get_gradient(), "gradient_image.tif");
  */
  std::cout << "Starting cell" << std::endl;
  //start cell
  cell_mgr c(f.get_size(),f.get_filt(),f.get_intensity(),f.get_gradient());
  //cell_mgr c(f.getsize(),f.getfilt(),f.getintensity(),f.getgradient());
  std::string a= oss.str();
  c.setup_cells(min,max,step,a);	//create cells from labeled image
  
  if( mode == 0 )
    {
      std::cout << "Getting train data" << std::endl;
      c.get_data(gliadat,neurdat);
    }
  else
    {
      std::cout << "Calc train data" << std::endl;
      c.calc_train_data( 0.95 , 1.0 , 0.9 , 0.2 , 1000 , 8000);
      
      //merge traning cells
      std::cout << "Build train rag" << std::endl;
      c.build_rag_train();
	
      //get rid of cells touching border
      c.remove_border_cells();
      
      std::cout << "Classifying" << std::endl;
      c.classify_cells();	
      
      std::cout << "Calculating train data" << std::endl;
      //save data into matrix
      c.calc_ng_data();
      //only keep good cells as training cells
      //c.calc_mean_vect(2.0/3.0);
      //save good cells into matrix
      //c.calc_ng_data();
      
      //reset cell class labels
      std::cout << "Relabeling Orig Cells" << std::endl;
      c.relabel_orig_cells();		
      
    }
  
  //merge and classify cells from train data
  std::cout << "Building 2nd rag" << std::endl;
  a = oss.str();
  c.build_rag(mode,a);	//merge cells
  
  std::cout << "Removing border Cells" << std::endl;
  c.remove_border_cells();
  
  if( mode == 1)
    {
      std::cout << "Cleaning up data" << std::endl;
      //c.save_train_data();	
    }	
  
  c.calc_ng_data();
  
  //get rid of cells not touching median plane
  c.remove_non_median_cells();
  
  //correct any mis classfied cells
  std::cout << "Correct Errors" << std::endl;
  a= oss.str();
  c.build_rag(mode,a,2);
  
  //relabel cells with consecutive label numbers
  std::cout << "Relabeling" << std::endl;
  c.relabel_cells(true);
  c.calc_ng_data();
  
  //score the segmentation
  std::cout << "Scoring" << std::endl;
  c.score_foreground_mah( oss.str() );
  
  c.check_all_cells(0,100,100,1,200);
  
  
  /*
    std::cout << "Labeling and Saving" << std::endl;
    //label border and number cells
    f.label_bound(c.get_cells());
    f.label_image(c.get_cells());
  */
  //save image
  //f.save_image_rgb(std::string("_out_" + oss.str() + ".tif"));
  //seg_result_img = f.save_image_rgb(std::string(filename + "_out_.tif"));
  //save cell features
  //c.output_cells(std::string( filename + ".txt" ));
  //std::vector<rich_cell::Pointer> rich_cells, rich_cells_in;
  
  vnl_vector_fixed<int,3> image_size;
  image_size.set( c.mImageSize );
  convert_cells_to_rich_cells(c.mvCell, image_size);
  cells_valid_ = true;
  
  //xml_util_write("cells.xml", rich_cells, parameters, c.mImageSize);
  //xml_util_read("cells_demo.xml", rich_cells_in);
  //xml_util_write("cells_demo_2.xml", rich_cells_in, parameters, c.mImageSize);
  /*
  //project cells
  std::cout << "Projecting cells" << std::endl;
  c.project_cells();
  //save projected image
  std::cout << "Saving projected" << std::endl;
  a = filename + "_proj_.tif";
  //proj_img = f.project_image(a,c.get_cells());
  */
  
  //save the intensity image
  //intensity_image_ = image;
}

void 
maciejSegmentation::
add_cell( rich_cell::Pointer new_cell )
{
  rich_cells_.push_back( new_cell ); 
}

rich_cell::Pointer
maciejSegmentation::
get_cell_at(vnl_vector_fixed<float,3> pos) const
{
  InputLabelImageType::IndexType index;
  index[0] = pos[0];
  index[1] = pos[1];
  index[2] = pos[2];

  InternalLabelImageType::PixelType label = label_image_->GetPixel( index );
  if ( label == 0)
    return 0;

  std::vector<rich_cell::Pointer>::const_iterator itr;
  itr = rich_cells_.begin();
  for (; itr!=rich_cells_.end(); itr++) {
    if ( (*itr)->label_ == label)
      return (*itr);
  }
  
  std::cerr<<"ERROR: No such label!!!"<<std::endl;
  return 0;
}

void
maciejSegmentation::
invalidate_cell( rich_cell::Pointer cell )
{
  cell->valid_ = false;
  update_label_image( cell );
}

rich_cell::Pointer 
maciejSegmentation::
merge_cells( std::vector<rich_cell::Pointer> cells )
{
  if (cells.size() < 2) {
    std::cout<<"Warning: No enough cells chosen. No cells are merged!"<<std::endl;
    return NULL;
  }

  typedef itk::ImageRegion< 3 >       RegionType;

  rich_cell::Pointer merged_cell = rich_cell::New();
  merged_cell->valid_ = true;
  merged_cell->label_ = ++cell_count_;

  
  //invalidate the old cells and keep track of the one with the
  //highest score and its class assignment
  merged_cell->score_ = cells[0]->score_;
  merged_cell->class_type_ = cells[0]->class_type_;
  for (unsigned int i = 0; i<cells.size(); i++) {
    invalidate_cell(cells[i]);
    if (cells[i]->score_ > merged_cell->score_) {
      merged_cell->score_ = cells[i]->score_;
      merged_cell->class_type_ = cells[i]->class_type_;
    }
  }
 
  // Compute the features for correct display of the
  // label_image_. First, merge the bounding box. Second, merge the
  // all_points_ from the invalidated cells. Third, update the
  // label_images_ by changing the old labels to the new label for the
  // merged cell.

  merged_cell->bounding_box_.SetIndex( cells[0]->bounding_box_.GetIndex() );
  merged_cell->bounding_box_.SetSize( cells[0]->bounding_box_.GetSize() );
  merged_cell->all_points_ = cells[0]->all_points_;
  assert( merged_cell->all_points_.size() );
 
  for (unsigned int i = 1; i<cells.size(); i++) {
    rich_cell::Pointer  cell2 = cells[i];
     RegionType::IndexType cell1_start = merged_cell->bounding_box_.GetIndex();
    RegionType::IndexType cell2_start = cell2->bounding_box_.GetIndex();
    RegionType::SizeType cell1_size = merged_cell->bounding_box_.GetSize();
    RegionType::SizeType cell2_size = cell2->bounding_box_.GetSize();
    RegionType::IndexType new_start, new_end; 
    new_start[0] = (cell1_start[0]<cell2_start[0])? cell1_start[0]:cell2_start[0];
    new_start[1]=(cell1_start[1]<cell2_start[1])? cell1_start[1]:cell2_start[1];
    new_start[2]=(cell1_start[2]<cell2_start[2])? cell1_start[2]:cell2_start[2];
    new_end[0]= (cell1_size[0]+cell1_start[0]>cell2_size[0]+cell2_start[0])? 
      cell1_size[0]+cell1_start[0]:cell2_size[0]+cell2_start[0];
    new_end[1]=(cell1_size[1]+cell1_start[1]>cell2_size[1]+cell2_start[1])? 
      cell1_size[1]+cell1_start[1]:cell2_size[1]+cell2_start[1];
    new_end[2]=(cell1_size[2]+cell1_start[2]>cell2_size[2]+cell2_start[2])? 
      cell1_size[2]+cell1_start[2]:cell2_size[2]+cell2_start[2];
    RegionType::SizeType new_size;
    new_size[0] = new_end[0]-new_start[0];
    new_size[1] = new_end[1]-new_start[1];
    new_size[2] = new_end[2]-new_start[2];
    
    merged_cell->bounding_box_.SetIndex( new_start );
    merged_cell->bounding_box_.SetSize( new_size );

    merged_cell->all_points_.insert(merged_cell->all_points_.end(), 
                                    cell2->all_points_.begin(), 
                                    cell2->all_points_.end());
  }
  

  rich_cells_.push_back( merged_cell ); 
  update_label_image( merged_cell );
  update_cell_pixels_for_one( merged_cell );
  //update_features( merged_cell );

  return merged_cell;
}

rich_cell::Pointer 
maciejSegmentation::
merge_cells(rich_cell::Pointer cell1, rich_cell::Pointer cell2)
{
  std::vector<rich_cell::Pointer> cells;
  cells.push_back( cell1 );
  cells.push_back( cell2 );

  return merge_cells( cells );

  /*
  typedef itk::ImageRegion< 3 >       RegionType;

  rich_cell::Pointer merged_cell = rich_cell::New();
  merged_cell->valid_ = true;
  merged_cell->label_ = ++cell_count_;
  invalidate_cell(cell1);
  invalidate_cell(cell2);

  // FEATURES OF THE MERGED CELL ARE NOT COMPUTED!!!

  // Compute the features for correct display of the
  // label_image_. First, merge the bounding box. Second, merge the
  // all_points_ from cell1 and cell2. Third, update the label_images_
  // by changing the old labels to the new label for the merged cell.

  RegionType::IndexType cell1_start = cell1->bounding_box_.GetIndex();
  RegionType::IndexType cell2_start = cell2->bounding_box_.GetIndex();
  RegionType::SizeType cell1_size = cell1->bounding_box_.GetSize();
  RegionType::SizeType cell2_size = cell2->bounding_box_.GetSize();
  RegionType::IndexType new_start, new_end; 
  (cell1_start[0]<cell2_start[0])? new_start[0]=cell1_start[0]:new_start[0]=cell2_start[0];
  (cell1_start[1]<cell2_start[1])? new_start[1]=cell1_start[1]:new_start[1]=cell2_start[1];
  (cell1_start[2]<cell2_start[2])? new_start[2]=cell1_start[2]:new_start[2]=cell2_start[2];
  (cell1_size[0]-1+cell1_start[0]>cell2_size[0]-1+cell2_start[0])? 
    new_end[0]=cell1_size[0]-1+cell1_start[0]:cell2_size[0]-1+cell2_start[0];
  (cell1_size[1]-1+cell1_start[1]>cell2_size[1]-1+cell2_start[1])? 
    new_end[1]=cell1_size[1]-1+cell1_start[1]:cell2_size[1]-1+cell2_start[1];
  (cell1_size[2]-1+cell1_start[2]>cell2_size[2]-1+cell2_start[2])? 
    new_end[2]=cell1_size[2]-1+cell1_start[2]:cell2_size[2]-1+cell2_start[2];

  RegionType::SizeType new_size;
  new_size[0] = new_end[0]-new_start[0]+1;
  new_size[1] = new_end[1]-new_start[1]+1;
  new_size[2] = new_end[2]-new_start[2]+1;

  merged_cell->bounding_box_.SetIndex( new_start );
  merged_cell->bounding_box_.SetSize( new_size );

  merged_cell->all_points_ = cell1->all_points_;
  assert( merged_cell->all_points_.size() );
  merged_cell->all_points_.insert(merged_cell->all_points_.end(), 
                                  cell2->all_points_.begin(), 
                                  cell2->all_points_.end());
 
  rich_cells_.push_back( merged_cell ); 
  update_label_image( merged_cell );

  return merged_cell;
  */
}

std::vector<rich_cell::Pointer> const &
maciejSegmentation::
all_cells() const
{
  return rich_cells_;
}

void
maciejSegmentation::
shift_cells(vnl_vector_fixed<float,3> shift)
{
  for (std::vector<rich_cell::Pointer>::iterator ci = rich_cells_.begin(); ci != rich_cells_.end(); ci++ ) {
      
      //skip both invalid and duplicated cells
      if ( !(*ci)->valid_ || (*ci)->dup_ ) continue;

      // shift bounding box
      typedef itk::Index<3> IndexType;
      IndexType index = (*ci)->bounding_box_.GetIndex();
      index[0] += shift[0];
      index[1] += shift[1];
      index[2] += shift[2];
      (*ci)->bounding_box_.SetIndex( index );
      
      // shift center point
      vnl_vector_fixed< float, 3 > center = (*ci)->center_;
      (*ci)->center_ += shift;
      
      // shift interior points
      rich_cell::PointVectorType::iterator pi; 
      for (pi=(*ci)->interior_points_.begin(); pi!=(*ci)->interior_points_.end(); pi++) {
        (*pi) += shift;
      } 
      for (pi = (*ci)->boundary_points_.begin(); pi!=(*ci)->boundary_points_.end(); pi++) {
        (*pi) += shift;
      } 
      
      // Combining interior_points_ & boundary_points_ to form all_points_
      (*ci)->all_points_.clear();
      (*ci)->all_points_ = (*ci)->boundary_points_;
    
      if ( !(*ci)->interior_points_.empty() ) 
        (*ci)->all_points_.insert((*ci)->all_points_.end(), 
                                  (*ci)->interior_points_.begin(), 
                                  (*ci)->interior_points_.end());

  }
}

vnl_vector_fixed<int,3> 
maciejSegmentation::
image_size() const
{
  vnl_vector_fixed<int,3> out_size(0);
  InternalLabelImageType::RegionType region = label_image_->GetLargestPossibleRegion();
  InternalLabelImageType::RegionType::SizeType size = region.GetSize();
  out_size[0] = size[0];
  out_size[1] = size[1];
  out_size[2] = size[2];
  return out_size;
} 

std::vector<std::string> const& 
maciejSegmentation::
parameters() const
{
  return parameters_;
}

void 
maciejSegmentation::
set_parameters(std::vector<std::string> const& parameters)
{
  parameters_ = parameters;
}

// THIS FUNCTION IS NOT TESTED!
void 
maciejSegmentation::
update_cell_pixels(InputLabelImageType::Pointer input_label_image)
{
  typedef itk::ImageRegionConstIterator< InputLabelImageType >  constRegionIteratorType;
  typedef itk::ImageRegionIterator< InternalLabelImageType >  RegionIteratorType;
  typedef itk::ImageRegion< 3 >         RegionType;

  // Convert the RGB image to the internal label image
  if (!label_image_) {
    label_image_ = InternalLabelImageType::New();
  }

  InternalLabelImageType::RegionType region = input_label_image->GetLargestPossibleRegion();
  label_image_->SetRegions( region );
  label_image_->Allocate();
  pixelType pixelValue;
  int newPixelValue;

  constRegionIteratorType  inputIt( input_label_image, input_label_image->GetLargestPossibleRegion() );
  RegionIteratorType  outputIt( label_image_, label_image_->GetLargestPossibleRegion() );

  for (inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd(); 
       ++inputIt, ++outputIt) {
        //conversion here
        pixelValue = inputIt.Get();
        newPixelValue = pixelValue[0]+pixelValue[1]*256+pixelValue[2]*65536;
        outputIt.Set(newPixelValue); 
  }

  //Now update the cell pixels of each cell. The interior pixels are
  //4-connected, i.e. if all for neighbors are part of the cell, the
  //pixel is marked as an interior pixel.

  typedef RegionType::IndexType IndexType;

 // The offsets for the neighboring pixels for 4-connectivity
  NeighborhoodIteratorType::OffsetType offset1 = {-1, 0, 0};
  NeighborhoodIteratorType::OffsetType offset2 = { 1, 0, 0};
  NeighborhoodIteratorType::OffsetType offset3 = { 0,-1, 0};
  NeighborhoodIteratorType::OffsetType offset4 = { 0, 1, 0};
  NeighborhoodIteratorType::OffsetType offset5 = { 0, 0,-1};
  NeighborhoodIteratorType::OffsetType offset6 = { 0, 0, 1};

  NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(1);

  for (unsigned int i = 0; i<rich_cells_.size(); i++) {
    rich_cell::Pointer cell = rich_cells_[i];
    cell->all_points_.reserve(100);
    cell->boundary_points_.reserve(100);
    cell->interior_points_.reserve(100);

    label_image_->SetRequestedRegion(cell->bounding_box_);
    NeighborhoodIteratorType it( radius, label_image_,
                                 label_image_->GetRequestedRegion());

    for ( it.GoToBegin(); !it.IsAtEnd(); ++it ) {
      if ( it.GetCenterPixel() == cell->label_ ) { // in the mask
      IndexType index = it.GetIndex();
      if ( it.GetPixel(offset1) == cell->label_ 
           && it.GetPixel(offset2) == cell->label_ 
           && it.GetPixel(offset3) == cell->label_
           && it.GetPixel(offset4) == cell->label_ 
           && it.GetPixel(offset5) == cell->label_ 
           && it.GetPixel(offset6) == cell->label_ ) {//interior point
        cell->interior_points_.push_back( rich_cell::PointType(index[0],index[1],index[2]) );
      }
      else cell->boundary_points_.push_back( rich_cell::PointType(index[0],index[1],index[2]) );
      }
    }

    // Combining interior_points_ & boundary_points_ to form all_points_
    cell->all_points_ = cell->boundary_points_;
    
    if ( !cell->interior_points_.empty() ) 
      cell->all_points_.insert(cell->all_points_.end(), 
                               cell->interior_points_.begin(), 
                               cell->interior_points_.end());
  }
  cells_valid_ = true;
}

void 
maciejSegmentation::
update_cell_pixels(std::string const & label_image_filename, int size_x, int size_y, int size_z)
{

  int array_size = size_x*size_y*size_z;
  unsigned short * memblock = new unsigned short[array_size];
  
  std::ifstream dat_file(label_image_filename.c_str(), std::ios::binary);
  if (dat_file.is_open()) {
    std::cerr<<"Opened .dat file "<<label_image_filename<<std::endl;
    dat_file.read(reinterpret_cast<char *>(memblock), sizeof(unsigned short)*array_size);
  }
  else {
    std::cerr<<"Cannot open the .dat file "<<label_image_filename<<std::endl;
    delete [] memblock;
    return;
  }
  dat_file.close();

  update_cell_pixels(memblock, size_x, size_y, size_z);
  delete [] memblock;
}

void 
maciejSegmentation::
update_cell_pixels(unsigned short * memblock, int size_x, int size_y, int size_z)
{  
  typedef itk::ImageRegion< 3 >         RegionType;
  int array_size = size_x*size_y*size_z;

  // Create the image
  if (!label_image_) {
    label_image_ = InternalLabelImageType::New();
  }
  
  InputLabelImageType::IndexType r_start;
  r_start[0] = 0; // first index on X
  r_start[1] = 0; // first index on Y
  r_start[2] = 0; // first index on Z
  
  InputLabelImageType::SizeType size;
  size[0] = size_x; // size along X
  size[1] = size_y; // size along Y
  size[2] = size_z; // size along Z
  
  InputLabelImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( r_start );

  label_image_->SetRegions( region );
  label_image_->Allocate();

  // Create the label image
  //
  unsigned short max_int = 0;

  typedef itk::ImageSliceIteratorWithIndex< InternalLabelImageType > IteratorType;
  IteratorType outputIt( label_image_, label_image_->GetRequestedRegion() );

  outputIt.SetFirstDirection( 2 ); //first direction is Z
  outputIt.SetSecondDirection( 0 ); //second direction is X
  int array_index = 0;

  outputIt.GoToBegin();
  while( !outputIt.IsAtEnd() ) {
    while ( !outputIt.IsAtEndOfSlice() ) {
      while ( !outputIt.IsAtEndOfLine() ) {
        assert( array_index<array_size );
        outputIt.Set(memblock[array_index]);
        if (memblock[array_index] > max_int) max_int = memblock[array_index];
        ++outputIt;
        array_index++;
      }
      outputIt.NextLine();
    }
    outputIt.NextSlice();
  }
 
  // label_image_ has to be flipped in y.
  typedef itk::FlipImageFilter< InternalLabelImageType > FilterType;
  FilterType::Pointer filter = FilterType::New();
  typedef FilterType::FlipAxesArrayType FlipAxesArrayType;
  FlipAxesArrayType flipArray;
  flipArray[0] = 0;
  flipArray[1] = 1;
  flipArray[2] = 0;
  filter->SetFlipAxes( flipArray );
  filter->SetInput( label_image_ );
  filter->Update();
  label_image_ = filter->GetOutput();
 
  // update the bounding boxes, the center and the volume information
  // of the cells. The way to do it is to first create an vcl_vector
  // of point list. The position of the list in the vextor indicates
  // the cell ID. Scan the image in raster order and store the pixel
  // position in the appropiate list. At the end, each list contains
  // all the pixels of one cell with the ID the same as the index in
  // the vector.
  typedef RegionType::IndexType IndexType;

  std::vector< std::vector<IndexType> > all_point_lists(max_int+1);
  typedef itk::ImageRegionConstIterator< InternalLabelImageType >  constRegionIteratorType;
  constRegionIteratorType  It( label_image_, label_image_->GetLargestPossibleRegion() );
  for (It.GoToBegin(); !It.IsAtEnd(); ++It) {
    if ( It.Get() == 0) continue;
    all_point_lists[ It.Get() ].push_back(It.GetIndex());
  }

  for (unsigned int i = 0; i<rich_cells_.size(); i++) {
    rich_cell::PointType center;
    rich_cell::Pointer cell = rich_cells_[i];
    unsigned short id = cell->label_;
    if (id<1) continue; //we're not sure if the background is set to 0

    std::vector<IndexType> const& list = all_point_lists[ id ];
    IndexType start = list[0]; 
    IndexType end = start;
    center[0] = start[0];
    center[1] = start[1];
    center[2] = start[2];

    for ( unsigned int j = 1; j<list.size(); j++ ) {
      for (unsigned int sj = 0; sj<3; sj++) {
        if (start[sj] > list[j][sj]) start[sj] = list[j][sj];
        center[sj] += list[j][sj];
      }
      for (unsigned int ej = 0; ej<3; ej++) {
        if (end[ej] < list[j][ej]) end[ej] = list[j][ej];
      }
    }

    RegionType::SizeType size;
    for (unsigned int sj = 0; sj<3; sj++) {
      size[sj] = end[sj] - start[sj] + 1;
      center[sj] /= (float)list.size();
    }

    cell->bounding_box_.SetSize( size );
    cell->bounding_box_.SetIndex( start );
    cell->volume_ = list.size();
    cell->center_ = center;
  }
  all_point_lists.clear();

  //Now update the cell pixels of each cell. The interior pixels are
  //4-connected, i.e. if all for neighbors are part of the cell, the
  //pixel is marked as an interior pixel.
  
  // The offsets for the neighboring pixels for 4-connectivity
  NeighborhoodIteratorType::OffsetType offset1 = {-1, 0, 0};
  NeighborhoodIteratorType::OffsetType offset2 = { 1, 0, 0};
  NeighborhoodIteratorType::OffsetType offset3 = { 0,-1, 0};
  NeighborhoodIteratorType::OffsetType offset4 = { 0, 1, 0};
  NeighborhoodIteratorType::OffsetType offset5 = { 0, 0,-1};
  NeighborhoodIteratorType::OffsetType offset6 = { 0, 0, 1};

  NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(1);

  for (unsigned int i = 0; i<rich_cells_.size(); i++) {
    rich_cell::Pointer cell = rich_cells_[i];
    //cell->all_points_.reserve(100);
    cell->boundary_points_.reserve(100);
    cell->interior_points_.reserve(100);

    label_image_->SetRequestedRegion(cell->bounding_box_);
    NeighborhoodIteratorType it( radius, label_image_,
                                 label_image_->GetRequestedRegion());

    for ( it.GoToBegin(); !it.IsAtEnd(); ++it ) {
      if ( it.GetCenterPixel() == cell->label_ ) { // in the mask
      IndexType index = it.GetIndex();
      if ( it.GetPixel(offset1) == cell->label_ 
           && it.GetPixel(offset2) == cell->label_ 
           && it.GetPixel(offset3) == cell->label_
           && it.GetPixel(offset4) == cell->label_ 
           && it.GetPixel(offset5) == cell->label_ 
           && it.GetPixel(offset6) == cell->label_ ) {//interior point
        cell->interior_points_.push_back( rich_cell::PointType(index[0],index[1],index[2]) );
      }
      else cell->boundary_points_.push_back( rich_cell::PointType(index[0],index[1],index[2]) );
      }
    }

    // Combining interior_points_ & boundary_points_ to form all_points_
    cell->all_points_ = cell->boundary_points_;
    
    if ( !cell->interior_points_.empty() ) 
      cell->all_points_.insert(cell->all_points_.end(), 
                               cell->interior_points_.begin(), 
                               cell->interior_points_.end());
  }
  cells_valid_ = true;
}

void 
maciejSegmentation::
update_cell_pixels_for_one(rich_cell::Pointer cell)
{
  typedef itk::ImageRegion< 3 > RegionType;
  typedef RegionType::IndexType IndexType;

  // The offsets for the neighboring pixels for 4-connectivity
  NeighborhoodIteratorType::OffsetType offset1 = {-1, 0, 0};
  NeighborhoodIteratorType::OffsetType offset2 = { 1, 0, 0};
  NeighborhoodIteratorType::OffsetType offset3 = { 0,-1, 0};
  NeighborhoodIteratorType::OffsetType offset4 = { 0, 1, 0};
  NeighborhoodIteratorType::OffsetType offset5 = { 0, 0,-1};
  NeighborhoodIteratorType::OffsetType offset6 = { 0, 0, 1};

  NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(1);

  cell->all_points_.reserve(100);
  cell->boundary_points_.reserve(100);
  cell->interior_points_.reserve(100);
  
  label_image_->SetRequestedRegion(cell->bounding_box_);
  NeighborhoodIteratorType it( radius, label_image_,
                               label_image_->GetRequestedRegion());

  for ( it.GoToBegin(); !it.IsAtEnd(); ++it ) {
    if ( it.GetCenterPixel() == cell->label_ ) { // in the mask
      IndexType index = it.GetIndex();
      if ( it.GetPixel(offset1) && it.GetPixel(offset2) &&
           it.GetPixel(offset3) && it.GetPixel(offset4) &&
           it.GetPixel(offset5) && it.GetPixel(offset6) ) {//interior point
        cell->interior_points_.push_back( rich_cell::PointType(index[0],index[1],index[2]) );
      }
      else cell->boundary_points_.push_back( rich_cell::PointType(index[0],index[1],index[2]) );
    }
  }
  
  // Combining interior_points_ & boundary_points_ to form all_points_
  cell->all_points_ = cell->boundary_points_;
  
  if ( !cell->interior_points_.empty() ) 
    cell->all_points_.insert(cell->all_points_.end(), 
                             cell->interior_points_.begin(), 
                             cell->interior_points_.end());
}

void 
maciejSegmentation::
output_cell_pixels(std::string const & filename) const
{
  // Save the labels of the nuclei as an RGB image. The first 1~255 is
  // represented by the red channel, 256~65535 is represented by the
  // red+green channel, and 65535~16777215 is by all 3 channels. the
  // entire label_image_ is copied to the output RGB image
  
  typedef itk::ImageFileWriter< InputLabelImageType >  writerType;
  typedef itk::ImageRegionConstIterator< InternalLabelImageType >  constRegionIteratorType;
  typedef itk::ImageRegionIterator< InputLabelImageType >  RegionIteratorType;


  InputLabelImageType::Pointer image = InputLabelImageType::New();
  image->SetRegions( label_image_->GetLargestPossibleRegion() );
  image->Allocate();

  unsigned short pixelValue;
  pixelType newPixelValue;
  int temp1; 
  constRegionIteratorType  inputIt( label_image_, label_image_->GetLargestPossibleRegion() );
  RegionIteratorType  outputIt( image, image->GetLargestPossibleRegion() );

  for (inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd(); 
       ++inputIt, ++outputIt) {
        //conversion here
        pixelValue = inputIt.Get();

        newPixelValue[2] = pixelValue/65536;
        temp1 = pixelValue%65536;
        newPixelValue[1] = temp1/256;
        newPixelValue[0] = temp1%256;
        outputIt.Set(newPixelValue);
  }

  writerType::Pointer writer = writerType::New();
  writer->SetFileName( filename.c_str() );
  
  writer->SetInput( image );
  try
    {
      writer->Update();
    }
  catch (itk::ExceptionObject & e)
    {
      std::cerr << "Exception in Writer: " << e << std::endl;
      exit(0);
    }
}


void 
maciejSegmentation::
construct_bins(double bin_size, int class_label, bool dup_removed)
{
  bins_set_ = true;

  // Determine the region for binning. Set the min to be (0,0,0) and
  // the max to be the image size
  //
  typedef bin_type::point_type point_type;
  point_type min;
  point_type max;
  min.fill( 0 );
  vnl_vector_fixed<int,3> size = image_size();
  max[0] = size[0];
  max[1] = size[1];
  max[2] = size[2];

  point_type bin_sizes;
  bin_sizes.fill( bin_size );
  bins_.reset( new bin_type( min, max, bin_sizes ) );

  // Add the data
  for (unsigned int i = 0; i<rich_cells_.size(); i++) {
    if (!rich_cells_[i]->valid_ || (dup_removed && rich_cells_[i]->dup_) ) continue;
    
    if ( class_label == 0 ) //ignore the class
      bins_->add_point( rich_cells_[i]->center_, rich_cells_[i] );
    else if (rich_cells_[i]->class_type_ == class_label)
      bins_->add_point( rich_cells_[i]->center_, rich_cells_[i] );
  }
}

void 
maciejSegmentation::
k_nearest_cells( std::vector<rich_cell::Pointer> & results, 
                 vnl_vector_fixed<float,3> const& center, 
                 unsigned int k ) const
{
  bins_->n_nearest( center, k, results );
}

void 
maciejSegmentation::
cells_within_radius( std::vector<rich_cell::Pointer> & results, 
                     vnl_vector_fixed<float,3> const& center,
                     double radius ) const
{
  bins_->points_within_radius( center, radius, results );
}

void 
maciejSegmentation::
convert_cells_to_rich_cells(std::vector<cell*> const &cells,
                            vnl_vector_fixed<int,3> const& image_size)
{
  rich_cells_.clear();
  label_image_ = 0;
  cell_count_ = 0;
  cells_valid_ = false;

  rich_cells_.reserve(cells.size());

  for (unsigned int i = 0; i<cells.size(); i++) {
    if (cells[i] == NULL) continue;

    cell* c = cells[i];
    rich_cell::Pointer rcell = rich_cell::New();
    rcell->valid_ = true;
    rcell->volume_ = c->mVolume;
    rcell->ave_intensity_ = c->mAvgInt;
 
    rcell->center_ = rich_cell::PointType(c->mCenterX, c->mCenterY, c->mCenterZ);
    rcell->label_ = c->mLabel;
    if ( cell_count_ <  rcell->label_ ) cell_count_ = rcell->label_;
    rcell->texture_ = c->mTexture;
    rcell->eccentricity_ = c->mEccentricity;
    rcell->score_ = c->mScore;
    rcell->class_type_ = c->mClass;

    /*
    rcell->ave_bound_gradient_ = c->mBoundGrad;
    rcell->bound_ints_ratio_ = c->mBoundIntsRatio;
    rcell->convexity_ = c->mConvexity;
    rcell->shape_factor_ = c->mShapeFact;    
    rcell->bending_energy_ = c->mBendEng;
    rcell->vol_grad_ = c->mVolGrad;
    rcell->radius_variance_ = c->mRadVar;
    */
    /*
    rcell->percent_nbr_ = c->mPerNbr;
    rcell->neighbors_ = c->mvNbrs;
    */

    /*
    rcell->avg_bound_intensity_ = c->mAvgBoundInt; //not stored in file
    rcell->major_axis_ = c->mMajorAxis; //not stored in file
    rcell->minor_axis_ = c->mMinorAxis; //not stored in file
    rcell->vol_grad_var_ = c->mVolGradVar; //not stored in file
    */

    //rcell->num_bound_pixel_ = c->mBoundPix;
    //rcell->ave_inside_gradient_ = c->mInsideGrad;
    //rcell->total_intensity_ = c->mTotInts;
    //rcell->var_rad_ = c->mVarRad;
    //rcell->num_shared_pixel_ = c->mSharedPix;
    //rcell->orientation_ = c->mOrientation;

    //Get all the points and set the bounding box
    bool min_set = false;
    double min_x= 0, min_y=0, min_z=0, max_x=0, max_y=0, max_z=0;
    for (unsigned int j = 0; j<c->mvPix.size(); j++) {
      rich_cell::PointType p(c->mvPix[j]->x_, c->mvPix[j]->y_, c->mvPix[j]->z_);
      rcell->all_points_.push_back(p);

      //set min
      if (!min_set) {
        min_x = c->mvPix[j]->x_;
        min_y = c->mvPix[j]->y_;
        min_z = c->mvPix[j]->z_;
        min_set = true;
      }
      else {
        if (c->mvPix[j]->x_ < min_x) min_x = c->mvPix[j]->x_;
        if (c->mvPix[j]->y_ < min_y) min_y = c->mvPix[j]->y_;
        if (c->mvPix[j]->z_ < min_z) min_z = c->mvPix[j]->z_;
      }
      // set max
      if (c->mvPix[j]->x_ > max_x) max_x = c->mvPix[j]->x_;
      if (c->mvPix[j]->y_ > max_y) max_y = c->mvPix[j]->y_;
      if (c->mvPix[j]->z_ > max_z) max_z = c->mvPix[j]->z_;
    }
    rich_cell::RegionType::IndexType start;
    start[0] = min_x;
    start[1] = min_y;
    start[2] = min_z;

    rich_cell::RegionType::SizeType size;
    size[0] = max_x - min_x +1;
    size[1] = max_y - min_y +1;
    size[2] = max_z - min_z +1;
    
    rcell->bounding_box_.SetSize( size );
    rcell->bounding_box_.SetIndex( start );
    
    /*
    //just for testing:
    if (rcell->label_ == 2) {
      edit_record rec;
      rec.status_ = edit_record::SPLIT;
      rec.replacements_.push_back(4);
      rec.replacements_.push_back(14);
      rec.date_ = "2006/09/28";
      rec.name_ = "Chia-Ling Tsai";
      rcell->edit_history_.push_back(rec);
      rcell->valid_ = false;
    }
    //end of testing
    */
    rich_cells_.push_back(rcell);
  }

  create_label_image( image_size );
}

void 
maciejSegmentation::
create_label_image(vnl_vector_fixed<int,3> const& image_size)
{
  label_image_ = InternalLabelImageType::New();
  InternalLabelImageType::IndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;

  InternalLabelImageType::SizeType size;
  size[0] = image_size[0]; // size along X
  size[1] = image_size[1]; // size along Y
  size[2] = image_size[2]; // size along Z

  InternalLabelImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );

  label_image_->SetRegions( region );
  label_image_->Allocate();
  label_image_->FillBuffer(0);

  for (unsigned int i = 0; i<rich_cells_.size(); i++) {
    if (rich_cells_[i]->valid_) {
      update_label_image( rich_cells_[i] );
      update_cell_pixels_for_one( rich_cells_[i] );
    }
  }
}

void 
maciejSegmentation::
update_label_image(rich_cell::Pointer cell)
  {

  InternalLabelImageType::IndexType index;
 
  if ( cell->valid_ ) { //mark the cell pixel 
	rich_cell::PointVectorType::const_iterator itr = cell->all_points_.begin();
    for (; itr!=cell->all_points_.end(); itr++) {
      index[0] = (*itr)[0];
      index[1] = (*itr)[1];
      index[2] = (*itr)[2];
	  
      label_image_->SetPixel( index, cell->label_);

    }
  }
  else { //unmark the cell pixels
	rich_cell::PointVectorType::const_iterator itr = cell->all_points_.begin();
    for (; itr!=cell->all_points_.end(); itr++) {
      index[0] = (*itr)[0];
      index[1] = (*itr)[1];
      index[2] = (*itr)[2];
      label_image_->SetPixel( index, 0);
    }
  }
}

maciejSegmentation::InternalLabelImageType::ConstPointer 
maciejSegmentation::
label_image() const
{
  return label_image_.GetPointer();
}

/*
void
maciejSegmentation::
update_features( rich_cell::Pointer cell )
{
  update_cell_pixels_for_one( cell );
  compute_volume_related_features( cell );
  //compute_intensity_related_features( intensity_image_, cell );
  //compute_neighbor_related_features( label_image_, cell );
}
*/
