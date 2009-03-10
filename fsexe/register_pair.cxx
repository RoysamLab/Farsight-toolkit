//: Executable program to register a pair of 3D image 
// 
//  The outcome is an xml file of a registration record. The images
//  can be gray, rgb color or rgba color. The input images are
//  assumed TIF images. The useage is
//
//   register_pair from_image_name to_image_name
//
//
// where 
//  from_image_name      name of the from_image
//  to_image_name        name of the to_image
// 
// optional arguments:
//  -channel        The channel to be extracted.
//  -background     Intensities below -background are ignored.
//  -smooth         Gaussian smoothing is required

#include "itkImageFileReader.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkRGBAPixel.h"
#include "itkRGBPixel.h"
#include "itkImage.h"
#include "itkTIFFImageIO.h"
#include "itkImageFileWriter.h"

#include <vul/vul_file.h>
#include <vul/vul_arg.h>

#include <fregl/fregl_pairwise_register.h>
#include <fregl/fregl_reg_record.h>
#include <fregl/fregl_util.h>
#include <Common/fsc_channel_accessor.h>

typedef unsigned char                   InputPixelType;
typedef itk::Image< InputPixelType, 3 > ImageType3D;
typedef float                           InternalPixelType;
typedef itk::Image<InternalPixelType, 3 > InternalImageType;

/*
ImageType3D::Pointer
read_image( std::string const & file_name, int channel )
{
  std::cout<<"Reading the image "<<file_name<<std::endl;

  ImageType3D::Pointer image;

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
    typedef itk::ImageFileReader< ImageType3D > ReaderType;
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
*/

ImageType3D::Pointer  
smooth_image(ImageType3D::Pointer image, int num_sub_images )
{
  typedef itk::DiscreteGaussianImageFilter< ImageType3D,InternalImageType > SmoothingFilterType;
  
  typedef itk::CastImageFilter< InternalImageType,ImageType3D > CastFilterType;
  typedef itk::RescaleIntensityImageFilter< ImageType3D, ImageType3D> RescalerType;
  typedef itk::StreamingImageFilter<ImageType3D, ImageType3D> StreamingFilterType;

  SmoothingFilterType::Pointer smoother = SmoothingFilterType::New();
  RescalerType::Pointer rescaler = RescalerType::New();
  CastFilterType::Pointer caster = CastFilterType::New();
  StreamingFilterType::Pointer streamer = StreamingFilterType::New();

  smoother->SetInput( image );
  smoother->SetVariance(1.0);
  smoother->SetMaximumKernelWidth(15);

  caster->SetInput(smoother->GetOutput());

  rescaler->SetInput( caster->GetOutput() );
  rescaler->SetOutputMinimum( 0 );
  rescaler->SetOutputMaximum( 250 );

  streamer->SetInput( rescaler->GetOutput() );
  streamer->SetNumberOfStreamDivisions(num_sub_images);
  try {
    streamer->Update();
  }
  catch(itk::ExceptionObject& e) {
    vcl_cout << e << vcl_endl;
  }
  
  return streamer->GetOutput();
}

int
main(  int argc, char* argv[] )
{
  vul_arg< vcl_string > arg_file_from  ( 0, "From image file" );
  vul_arg< vcl_string > arg_file_to    ( 0, "To image file" );
  //vul_arg< bool >      alpha           ("-alpha", "Input are color images with alpha channel. This is ignored if -color not set", false);
  vul_arg< int >       channel         ("-channel", "The color channel (0-red, 1-green, 2-blue).",0);
  vul_arg< float >     background      ("-bg", "threshold value for the background", 30);
  vul_arg< double >      smooth        ("-smooth", "If smoothing is performed", 0.5);
  //vul_arg< int >       streaming       ("-streaming","Apply streaming with the number of sub-images when smooth is set", 1);
  //vul_arg< bool >     exhaustive      ("-exhaustive", "Running exhaustive search", false);
  vul_arg< vcl_string > gdbicp ("-gdbicp","The path where the gdbicp exe can be found. The exe of gdbicp can be found at http://www.vision.cs.rpi.edu/gdbicp","");
  vul_arg< int >         slices         ("-slices","Number of slices to register. Needed when registering large image stacks on PC.",100);
  vul_arg<bool>          remove_2d      ("-remove_2d", "Remove the intermediate 2d stuff", false);
  vul_arg<bool>          scaling_arg    ("-scaling","Substantial scaling is expected.", false);
  
  vul_arg_parse( argc, argv );

  ImageType3D::Pointer from_image, to_image;

  // Read the image
  //
  //
  
  from_image = fregl_util_read_image( arg_file_from(), channel.set(), channel() );
  to_image = fregl_util_read_image( arg_file_to(), channel.set(), channel() );
  if (!from_image || !to_image) {
    vcl_cerr <<"Failed to read image(s)"<< vcl_endl;
  }
  /*
  typedef itk::ImageFileWriter< ImageType >  WriterType3D;
  WriterType3D::Pointer writer3D = WriterType3D::New();
  writer3D->SetFileName( "from.tiff" );
  writer3D->SetInput( from_image );
  writer3D->Update();
  writer3D->SetFileName( "to.tiff" );
  writer3D->SetInput( to_image );
  writer3D->Update();
  */
  
  // Apply Gaussian smoothing if demanded
  //
  /*
  if ( smooth() ) {
    std::cout<<"Perform smoothing ..."<<std::endl;
    from_image = smooth_image( from_image, streaming() );
    to_image = smooth_image( to_image, streaming() );
  }
  */

  // Perform registration
  //
  fregl_pairwise_register registor(from_image, to_image, background());
  //if (exhaustive()) registor.set_exhausive_search(true);
  if (slices.set()) registor.set_stack_size( slices() );
  registor.set_smoothing( smooth() );

  double obj_value;
  bool failed = false;
  if (registor.run(obj_value, gdbicp(), scaling_arg())) {
    
    // Set the registration result
    ImageType3D::SizeType from_image_size = from_image->GetLargestPossibleRegion().GetSize();
    ImageType3D::SizeType to_image_size = to_image->GetLargestPossibleRegion().GetSize();

    std::string from_image_id = vul_file::strip_directory( arg_file_from() );
    std::string to_image_id = vul_file::strip_directory( arg_file_to() );
    std::string from_image_id_wo_ext = vul_file::strip_extension( from_image_id );
    std::string to_image_id_wo_ext = vul_file::strip_extension( to_image_id );

    // Invert the transforn, since itk produces xform going from
    // to->from image
    typedef fregl_pairwise_register::TransformType TransformType;
    TransformType::Pointer xform = registor.transform();
    TransformType::Pointer inv_xform = TransformType::New();
    xform->GetInverse(inv_xform);

    //Create the reg_record and write to xml output
    std::cout<<"Writing the xml output"<<std::endl;
    
   
    fregl_reg_record::Pointer record = new fregl_reg_record();
    record->set_from_image_info(from_image_id, from_image_size);
    record->set_to_image_info(to_image_id, to_image_size);
    record->set_transform( inv_xform );
    record->set_obj_value( obj_value );
    record->set_overlap(1);
    std::string xml_file = from_image_id_wo_ext+std::string("_to_")+
      to_image_id_wo_ext+std::string("_transform.xml");
    record->write_xml( xml_file );
  }
  else {
    failed = true; 
    std::cout<<"Registration failed"<<std::endl;
  }

  if (remove_2d()) {
#if defined(VCL_WIN32) && !defined(__CYGWIN__)
    vcl_system("del xxx_from_image_proj.tif");
    vcl_system("del xxx_to_image_proj.tif");
    vcl_system("del mosaic_xxx_from_image_proj_to_xxx_to_image_proj.xform");
#else
    vcl_system("rm xxx_from_image_proj.tif");
    vcl_system("rm xxx_to_image_proj.tif");
    vcl_system("rm mosaic_xxx_from_image_proj_to_xxx_to_image_proj.xform");
#endif
  }

  if (failed) return 1;
  
  return 0;
}
