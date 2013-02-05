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

//: Executable program to mosaic a set of images with given transformations
//
//  The input is an xml file containing either one pairwise
//  registration result or a set of transformations from joint
//  registration. The images can be gray, rgb color or rgba color. The
//  input images are assumed TIF images. The output montage is a gray
//  image. The usage is
//  
//  mosaic_images xml_file 
//
//  where
//    xml_file      Name of the xml_file containing the xforms
//    anchor_image  Name of the anchor image
//
//  Optional arguments"
//    -path         The path of the image files.
//    -old_str      The old substr in the image names to be replaced
//    -new_str      The replacement of the old substr
//    -output       The output image name.

#include "mosaic_images_template.h"
#include "itkTimeProbe.h"
#ifdef _OPENMP
#include "omp.h"
#endif

template <typename TPixel>
int mosaic_images_template(
	vul_arg< vcl_string > arg_xml_file,
	vul_arg< vcl_string > arg_anchor,
	vul_arg< int > arg_channel,
	vul_arg< vcl_string > arg_img_path,
	vul_arg< vcl_string > arg_old_str,
	vul_arg< vcl_string > arg_new_str,
	vul_arg< bool > arg_3d,
	vul_arg< vcl_string > arg_outfile,
	vul_arg< bool > arg_in_anchor,
	vul_arg< bool > arg_overlap,
	vul_arg< bool > arg_nn,
	vul_arg< bool > arg_normalize, 
	vul_arg< vcl_string > arg_background,
	vul_arg< double > arg_sigma, 
	vul_arg< double > arg_median,
	vul_arg< int > arg_blending,
	vul_arg< bool > arg_denoise,
	vul_arg< bool > arg_write_interproj_image)
{
	typedef TPixel InputPixelType;
	typedef itk::Image< InputPixelType, 3 > ImageType;
	typedef itk::Image< InputPixelType, 2 > ImageType2D;
	typedef typename ImageType::PointType PointType; //physical space
	typedef typename ImageType::IndexType IndexType; //physical space
	typedef typename ImageType::SizeType SizeType;
	typedef itk::ImageRegionConstIterator< ImageType > RegionConstIterator;
	typedef itk::ImageRegionIterator< ImageType > RegionIterator;
	typedef itk::ImageRegionConstIterator< ImageType2D > RegionConstIterator2D;
	typedef itk::ImageLinearConstIteratorWithIndex< ImageType2D > LinearConstIteratorType2D;
	typedef itk::ImageSliceIteratorWithIndex< ImageType > SliceIteratorType;
	typedef itk::ImageSliceConstIteratorWithIndex< ImageType > SliceConstIteratorType;
	typedef itk::ImageFileWriter< ImageType2D > WriterType2D;

	// Cosntruct the graph of joint registration
    typename fregl_joint_register< InputPixelType >::Pointer joint_register = new fregl_joint_register< InputPixelType > (arg_xml_file());
    if (arg_old_str.set() && arg_new_str.set()) {
        std::cout << "Replace the name substr" << std::endl;
        joint_register->replace_image_name_substr(arg_old_str(), arg_new_str());
    }

    switch (arg_blending()) {
        case 0: std::cout << "Blending with maximum intensity values" << std::endl;
            break;
        case 1: std::cout << "Blending with even weighting" << std::endl;
            break;
        case 2: std::cout << "Blending with photopleaching weighting" << std::endl;
            break;
        case 3: std::cout << "Blending maximum intensity values (this function is under testing), reduces memory, and is much faster" << std::endl;
            break;
		case 4: std::cout << "Blending with XOR"<<std::endl;
			break;
        default: std::cout << "No such scheme defined for blending!" << std::endl;
            return 1;
    }

    // Transform the images
    //
    fregl_space_transformer< InputPixelType > space_transformer(joint_register);

    //bool in_anchor = false;
    space_transformer.set_anchor(arg_anchor(), arg_in_anchor(), arg_overlap());
    typename ImageType::Pointer final_image = ImageType::New();

    // Read in each image one by one. Transform each image and compose
    // the image with the final image.
    //
    std::vector<std::string> image_names = space_transformer.image_names();
    std::cout << "Total number of images = " << image_names.size() << std::endl;
    std::string name_prefix = std::string("montage_") + vul_file::strip_extension(arg_anchor());
    if (arg_outfile.set()) {
        name_prefix = arg_outfile();
    }

    if (arg_blending() == 2) {
        // This option takes into consideration the photobleahcing issue
        // as well. Photobleaching results in much lower intensity value
        // in an image when the region was previously imaged. This option
        // also attemps to remove the back noise. 

        // Individual weight maps are computed to factor in
        // photobleaching. The original 3D images are not saved to reduce
        // memory consumption. They will be read in again during blending
        for (unsigned int i = 0; i < image_names.size(); i++) {
            std::string image_name = arg_img_path() + std::string("/") + image_names[i];
            std::cout << "Image " << image_name << std::endl;
            typename ImageType::Pointer image = fregl_util< InputPixelType >::fregl_util_read_image(image_name, arg_channel.set(), arg_channel(), arg_denoise());
            float alpha = 5;
            space_transformer.set_individual_weight_map(i, image, alpha);
        }
        space_transformer.normalize_individual_weight_maps();

        // Now doing the blending
        std::string image_name = arg_img_path() + std::string("/") + image_names[0];
        typename ImageType::Pointer image, xformed_image;
        image = fregl_util< InputPixelType >::fregl_util_read_image(image_name, arg_channel.set(), arg_channel(), arg_denoise());
        std::cout << "Composing the final image ..." << std::endl;
        final_image = space_transformer.transform_image_weighted(image, 0, 0, arg_nn());
        for (unsigned int i = 1; i < image_names.size(); i++) {
            std::string image_name = arg_img_path() + std::string("/") + image_names[i];
            typename ImageType::Pointer image = fregl_util< InputPixelType >::fregl_util_read_image(image_name, arg_channel.set(), arg_channel(), arg_denoise());
            xformed_image = space_transformer.transform_image_weighted(image, i, 0, arg_nn());
            if (!xformed_image)
                continue;

            //fuse the image
	    typename ImageType::SizeType iamgeOutputSize = final_image->GetRequestedRegion().GetSize();
	    typename ImageType::PixelType * imageOutputArray = final_image->GetBufferPointer();
	    typename ImageType::PixelType * imageInputArray = xformed_image->GetBufferPointer();


		#if _OPENMP >= 200805L
			#pragma omp parallel for collapse(3)
		#else
			#pragma omp parallel for
		#endif
			for( int ii=0; ii<iamgeOutputSize[2]; ++ii )
            {
                for( int jj=0; jj<iamgeOutputSize[1]; ++jj )
                {
                    for( int kk=0; kk<iamgeOutputSize[0]; ++kk )
                    {
                        itk::Index<1> offset;
                        offset[0] = (ii*iamgeOutputSize[0]*iamgeOutputSize[1])+(jj*iamgeOutputSize[0])+kk;
                        imageOutputArray[offset[0]] = vnl_math_max(imageOutputArray[offset[0]],imageInputArray[offset[0]]);
					}
				}
			}

            ////fuse the image
            //RegionConstIterator inputIt(xformed_image, xformed_image->GetRequestedRegion());
            //RegionIterator outputIt(final_image, final_image->GetRequestedRegion());

            //for (inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++outputIt) {
            //    outputIt.Set(outputIt.Get() + inputIt.Get());
            //}
        }
    }
	else if (arg_blending() == 1) { //taking the average in the overlapping regions
        typename ImageType2D::Pointer weight_image = space_transformer.compute_weighted_image_2D();
        typename ImageType::IndexType start;
        start[0] = 0;
        start[1] = 0;
        start[2] = 0;
        typename ImageType::RegionType region;
        region.SetSize(space_transformer.montage_size());
        region.SetIndex(start);
        final_image->SetRegions(region);
        final_image->Allocate();
        RegionIterator outputIt(final_image, final_image->GetRequestedRegion());
        for (outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt)
            outputIt.Set(0);

        std::cout << "Composing the final image ..." << std::endl;
        for (unsigned int i = 0; i < image_names.size(); i++) {
            std::string image_name = arg_img_path() + std::string("/") + image_names[i];
            typename ImageType::Pointer image, xformed_image;
            image = fregl_util< InputPixelType >::fregl_util_read_image(image_name, arg_channel.set(), arg_channel(), arg_denoise());
            xformed_image = space_transformer.transform_image(image, i, 0, arg_nn());
            if (!xformed_image)
                continue;

            //fuse the image
            SliceIteratorType outputSliceIt(final_image, final_image->GetRequestedRegion());
            SliceConstIteratorType inputSliceIt(xformed_image, xformed_image->GetRequestedRegion());
            LinearConstIteratorType2D input2DIt(weight_image, weight_image->GetRequestedRegion());
            // Set the directions for traversal
            unsigned int direction[2];
            direction[0] = 0;
            direction[1] = 1;
            outputSliceIt.SetFirstDirection(direction[1]);
            outputSliceIt.SetSecondDirection(direction[0]);
            inputSliceIt.SetFirstDirection(direction[1]);
            inputSliceIt.SetSecondDirection(direction[0]);
            input2DIt.SetDirection(1 - direction[0]);

            outputSliceIt.GoToBegin();
            inputSliceIt.GoToBegin();
            while (!outputSliceIt.IsAtEnd()) {
                while (!outputSliceIt.IsAtEndOfSlice()) {
                    while (!outputSliceIt.IsAtEndOfLine()) {
                        if (input2DIt.Get() == 0)
                            outputSliceIt.Set(0);
                        else
                            outputSliceIt.Set(outputSliceIt.Get() + (int) (inputSliceIt.Get() / (float) input2DIt.Get()));
                        ++outputSliceIt;
                        ++inputSliceIt;
                        ++input2DIt;
                    }
                    outputSliceIt.NextLine();
                    inputSliceIt.NextLine();
                    input2DIt.NextLine();
                }
                input2DIt.GoToBegin();
                outputSliceIt.NextSlice();
                inputSliceIt.NextSlice();
            }
            /*
              for ( inputIt1.GoToBegin(), inputIt2.GoToBegin(), outputIt.GoToBegin(); 
              !inputIt1.IsAtEnd();  ++inputIt1, ++inputIt2, ++outputIt) {
              if (inputIt1.Get() == 0 ) outputIt.Set(0);
              else outputIt.Set( outputIt.Get() + (int)(inputIt2.Get()/(float)inputIt1.Get()) );
              }
             */
        }
    } 
	else if (arg_blending() == 0) 
	{//Taking the maximum
        std::string image_name = arg_img_path() + std::string("/") + image_names[0];
        typename ImageType::Pointer image, xformed_image;
        image = fregl_util< InputPixelType >::fregl_util_read_image(image_name, arg_channel.set(), arg_channel(), arg_denoise());
        std::cout << "Composing the final image arg_blending 0 ..." << std::endl;
        final_image = space_transformer.transform_image(image, 0, 0, arg_nn());
        for (unsigned int i = 1; i < image_names.size(); i++) 
		{
            image_name = arg_img_path() + std::string("/") + image_names[i];
            image = fregl_util< InputPixelType >::fregl_util_read_image(image_name, arg_channel.set(), arg_channel(), arg_denoise());
            itk::TimeProbe clock;
            clock.Start();
            xformed_image = space_transformer.transform_image(image, i, 0, arg_nn());
            clock.Stop();
            std::cout<< "Transform:"<<clock.GetTotal()<<std::endl;
            if (!xformed_image)
                continue;
        
		//fuse the image
			clock.Start();
			typename ImageType::SizeType imageOutputSize = final_image->GetRequestedRegion().GetSize();
			typename ImageType::PixelType * imageOutputArray = final_image->GetBufferPointer();
			typename ImageType::PixelType * imageInputArray = xformed_image->GetBufferPointer();

#if _OPENMP >= 200805L
	#pragma omp parallel for collapse(3)
#else
	#pragma omp parallel for
#endif
			for( int ii=0; ii<imageOutputSize[2]; ++ii )
			{
				for( int jj=0; jj<imageOutputSize[1]; ++jj )
				{
					for( int kk=0; kk<imageOutputSize[0]; ++kk )
					{
						itk::Index<1> offset;
						offset[0] = (ii*imageOutputSize[0]*imageOutputSize[1])+(jj*imageOutputSize[0])+kk;
						imageOutputArray[offset[0]] = vnl_math_max(imageOutputArray[offset[0]],imageInputArray[offset[0]]);
					}
				}
			}
			clock.Stop();
			std::cout<< "Fuse:"<<clock.GetTotal()<<std::endl;

			if( arg_write_interproj_image())
			{
				std::cout<< " 2d projection!!!"<<std::endl;
				std::stringstream strm;
				std::string num;
				strm << i;
				strm >> num;

				typename ImageType2D::Pointer image_2d = fregl_util< InputPixelType >::fregl_util_fast_max_projection(xformed_image);
				typename WriterType2D::Pointer writer2D = WriterType2D::New();
				std::string name_2d =  arg_img_path() + std::string("/") + num + std::string("_2d_tile.png");
				writer2D->SetFileName(name_2d);
				writer2D->SetInput(image_2d);
				writer2D->Update();
			}
		}
    }
	else if (arg_blending() == 3) 
	{ // THIS IS A TEST FOR THE DARPA DATASETS
		std::cout << std::endl << "FIXME2: this code is still under testing, this is for mosaic fast , also reduce the insane amount of memory used,  using max intensity, the spacing is sett up for the darpa project\n";
        std::string image_name = arg_img_path() + std::string("/") + image_names[0];
        typename ImageType::Pointer image, xformed_image;
		
        image = fregl_util< InputPixelType >::fregl_util_read_image(image_name, arg_channel.set(), arg_channel(), arg_denoise());
        std::cout << "Composing the final image arg_blending 3 ..." << image->GetLargestPossibleRegion().GetSize()[0]<< std::endl;
		
		typename ImageType::Pointer background_image = NULL;
		if(arg_background().length() > 0)
		{
			std::string backgroundimage_name = arg_img_path() + std::string("/") + arg_background();
			background_image = fregl_util< InputPixelType >::fregl_util_read_image(backgroundimage_name);
		}

		typename ImageType::SpacingType spacingImage;
		spacingImage[0] = 1;
		spacingImage[1] = 1;
		spacingImage[2] = 1.12359; // 0.3/0.267 FIXME this spacing is set up for the darpa project
		space_transformer.set_spacing( spacingImage );
	
		image->SetSpacing( spacingImage );
	
		final_image = space_transformer.transform_image(image, 0, 0, arg_nn(), arg_normalize(), background_image, arg_sigma(), arg_median());
	
		if( arg_write_interproj_image())
		{
			std::cout<< "2d projection 0."<<std::endl;
			typename ImageType2D::Pointer image_2d_0 = fregl_util< InputPixelType >::fregl_util_fast_max_projection(final_image);
			typename WriterType2D::Pointer writer2D = WriterType2D::New();
			std::string name_2d_0 =  arg_img_path() + std::string("/") + std::string("0_2d_proj.png");
			writer2D->SetFileName(name_2d_0);
			writer2D->SetInput(image_2d_0);
			writer2D->Update();
		}

        for (unsigned int i = 1; i < image_names.size(); i++) 
		{
            std::cout<< i<<std::endl;
			std::string image_name2 = arg_img_path() + std::string("/") + image_names[i];
			typename ImageType::Pointer image2;
            image2 = fregl_util< InputPixelType >::fregl_util_read_image(image_name2, arg_channel.set(), arg_channel(), arg_denoise());
			image2->SetSpacing( spacingImage );
	    
			typename ImageType::Pointer imageMontage;
			typename ImageType::PointType offsetNew;
			int fail = space_transformer.transform_image_fast(image2, imageMontage, offsetNew, i, 0, arg_nn(), arg_normalize(), background_image, arg_sigma(), arg_median());
            if (fail)
                continue;

			std::stringstream strm;
			std::string num;
			strm << i;
			strm >> num;
			
			if( arg_write_interproj_image() || i == 1)
			{
                std::cout<< " 2d projection before fuse."<<std::endl;
                typename ImageType2D::Pointer image_2d = fregl_util< InputPixelType >::fregl_util_fast_max_projection(imageMontage);
                typename WriterType2D::Pointer writer2D = WriterType2D::New();
                std::string name_2d =  arg_img_path() + std::string("/") + num + std::string("_2d_tile.png");
                writer2D->SetFileName(name_2d);
                writer2D->SetInput(image_2d);
                writer2D->Update();
			}

			//std::cout<< "offset!!!!!!"<<std::endl;
			//std::cout<< offsetNew[0]<<"\t"<<offsetNew[1]<<"\t"<<offsetNew[2]<<std::endl;

			typename ImageType::SizeType image_size_fast = imageMontage->GetLargestPossibleRegion().GetSize();
			typename ImageType::SizeType image_size_ = final_image->GetLargestPossibleRegion().GetSize();

			//std::cout<< "image size fast"<<std::endl;
			//std::cout<< image_size_fast[0]<<"\t"<<image_size_fast[1]<<"\t"<<image_size_fast[2]<<std::endl;

			unsigned long long sizeXY_small = image_size_fast[0]*image_size_fast[1];
			unsigned long long sizeX_small = image_size_fast[0];
			
			unsigned long long sizeXY_big = image_size_[0]*image_size_[1];
			unsigned long long sizeX_big = image_size_[0];
			//std::cout<< "image size big:"<<"\t"<<sizeXY_big<<std::endl;
	
			typename ImageType::PixelType * imageMontageArray = imageMontage->GetBufferPointer();
			typename ImageType::PixelType * imageResampleArray = final_image->GetBufferPointer();

// #pragma omp critical
// {
#if _OPENMP >= 200805L
	#pragma omp parallel for collapse(3)
#else
	#pragma omp parallel for
#endif
			for(int ii=0; ii<image_size_fast[2]; ++ii) // z
			{
					for(int j=0; j<image_size_fast[1]; ++j) // y
					{
							for(int k=0; k<image_size_fast[0]; ++k) // x
							{
								long long zz_big = ii + offsetNew[2];
								long long yy_big = j + offsetNew[1];
								long long xx_big = k + offsetNew[0];

								long long xx_small = k;
								long long yy_small = j;
								long long zz_small = ii;
					
								if( 0<=xx_big  && xx_big < image_size_[0] && 0<=yy_big  && yy_big < image_size_[1] && 0<=zz_big && zz_big < image_size_[2] ) 
								{
									long long offsetSmall = (zz_small*sizeXY_small)+(yy_small*sizeX_small)+xx_small;
									long long offsetBig = (zz_big*sizeXY_big)+(yy_big*sizeX_big)+xx_big;

									imageResampleArray[offsetBig] = vnl_math_max(imageResampleArray[offsetBig],imageMontageArray[offsetSmall]);
								}
							}
					}
			}
            
			if( arg_write_interproj_image() || i == 1)
			{
				std::cout<< "2d projection after fuse."<<std::endl;
	                      
				typename ImageType2D::Pointer image_2d2 = fregl_util< InputPixelType >::fregl_util_fast_max_projection(final_image);
				typename WriterType2D::Pointer writer2D2 = WriterType2D::New();
				std::string name_2d2 =  arg_img_path() + std::string("/") + num + std::string("_2d_proj.png");
				writer2D2->SetFileName(name_2d2);
				writer2D2->SetInput(image_2d2);
				writer2D2->Update();
			}
        }
    }
	else if (arg_blending() == 4) 
	{ // THIS IS A TEST FOR Time-lapse changes detection, blending with XOR
		std::cout<< "Start Mosaicing..."<<std::endl;
        std::string image_name = arg_img_path() + std::string("/") + image_names[0];
        typename ImageType::Pointer image, xformed_image;
        image = fregl_util< InputPixelType >::fregl_util_read_image(image_name, arg_channel.set(), arg_channel(), arg_denoise());

		typename ImageType::SpacingType spacingImage;
		spacingImage[0] = 1;
		spacingImage[1] = 1;
		spacingImage[2] = 1; 
		space_transformer.set_spacing( spacingImage );
	
		image->SetSpacing( spacingImage );
	
		final_image = space_transformer.transform_image(image, 0, 0, arg_nn());
	
        for (unsigned int i = 1; i < image_names.size(); i++) {

			std::string image_name2 = arg_img_path() + std::string("/") + image_names[i];
			typename ImageType::Pointer image2;
            image2 = fregl_util< InputPixelType >::fregl_util_read_image(image_name2, arg_channel.set(), arg_channel(), arg_denoise());
			image2->SetSpacing( spacingImage );
	    
			typename ImageType::Pointer imageMontage;
			typename ImageType::PointType offsetNew;
			int fail = space_transformer.transform_image_fast(image2, imageMontage, offsetNew, i, 0, arg_nn());
            if (fail)
                continue;

			typename ImageType::SizeType image_size_fast = imageMontage->GetLargestPossibleRegion().GetSize();
			typename ImageType::SizeType image_size_ = final_image->GetLargestPossibleRegion().GetSize();

			unsigned long long sizeXY_small = image_size_fast[0]*image_size_fast[1];
			unsigned long long sizeX_small = image_size_fast[0];
			
			unsigned long long sizeXY_big = image_size_[0]*image_size_[1];
			unsigned long long sizeX_big = image_size_[0];
	
			typename ImageType::PixelType * imageMontageArray = imageMontage->GetBufferPointer();
			typename ImageType::PixelType * imageResampleArray = final_image->GetBufferPointer();

// #pragma omp critical
// {
#if _OPENMP >= 200805L
	#pragma omp parallel for collapse(3)
#else
	#pragma omp parallel for
#endif
			for(int ii=0; ii<image_size_fast[2]; ++ii) // z
			{
					for(int j=0; j<image_size_fast[1]; ++j) // y
					{
							for(int k=0; k<image_size_fast[0]; ++k) // x
							{
								long long zz_big = ii + offsetNew[2];
								long long yy_big = j + offsetNew[1];
								long long xx_big = k + offsetNew[0];

								long long xx_small = k;
								long long yy_small = j;
								long long zz_small = ii;
					
								if( 0<=xx_big  && xx_big < image_size_[0] && 0<=yy_big  && yy_big < image_size_[1] && 0<=zz_big && zz_big < image_size_[2] ) 
								{
									long long offsetSmall = (zz_small*sizeXY_small)+(yy_small*sizeX_small)+xx_small;
									long long offsetBig = (zz_big*sizeXY_big)+(yy_big*sizeX_big)+xx_big;

									if(imageResampleArray[offsetBig] != imageMontageArray[offsetSmall])
									{
										imageResampleArray[offsetBig] = vnl_math_max(imageResampleArray[offsetBig],imageMontageArray[offsetSmall]);
									}
									else
									{
										imageResampleArray[offsetBig] = 0;
									}
								}
							}
					}
			}
        }
    }

    /*std::string command = std::string("mkdir ")+name_prefix;
      if(vcl_system(command.c_str()) != 0)
      {
      cerr << "mkdir returned nonzero" << endl;
      }
    
      typedef itk::NumericSeriesFileNames NameGeneratorType;
      typedef itk::ImageSeriesWriter< ImageType, ImageType2D > SeriesWriterType;
    
      NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
      SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
      seriesWriter->SetInput( final_image );
    
      int last=final_image->GetLargestPossibleRegion().GetSize()[2] - 1;
      nameGenerator->SetStartIndex( 0 );
      nameGenerator->SetEndIndex( last );
      nameGenerator->SetIncrementIndex( 1 );
    
      #if defined(VCL_WIN32) && !defined(__CYGWIN__)
      std::string name_pattern = name_prefix+std::string("\\slice")+std::string("_%03d.png");
      #else
      std::string name_pattern = name_prefix+std::string("/slice")+std::string("_%03d.png");
      #endif
      nameGenerator->SetSeriesFormat( name_pattern );
      seriesWriter->SetFileNames( nameGenerator->GetFileNames() );
      seriesWriter->Update();*/

    if (arg_3d()) {
        typename ImageType::PointType OutOrigin;
        // Make sure there are no offsets in the Mosaic file origin
	std::cout << std::endl << " \tNow it will create the images in 3d and NRRD";
        OutOrigin[0] = 0;
        OutOrigin[1] = 0;
        OutOrigin[2] = 0;
        final_image->SetOrigin(OutOrigin);
        std::string name_3d = name_prefix + std::string(".nrrd");
        typedef itk::ImageFileWriter< ImageType > WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetFileName(name_3d);
        writer->SetInput(final_image);
        writer->Update();
    }

    // doing the 2d maximum projection and dump it out
    typename ImageType2D::Pointer image_2d = fregl_util< InputPixelType >::fregl_util_fast_max_projection(final_image);
    typename WriterType2D::Pointer writer2D = WriterType2D::New();
    std::string name_2d = name_prefix + std::string("_2d_proj.png");
    writer2D->SetFileName(name_2d);
    writer2D->SetInput(image_2d);
    writer2D->Update();

    // dump the output to xml
    std::string xml_name = name_prefix + std::string(".xml");
    space_transformer.write_xml(xml_name, name_prefix, name_2d, arg_overlap(), arg_in_anchor(), arg_channel(), arg_blending(), arg_nn(), arg_denoise());

    return 0;
}

template int mosaic_images_template<unsigned char>(
	vul_arg< vcl_string > arg_xml_file,
	vul_arg< vcl_string > arg_anchor,
	vul_arg< int > arg_channel,
	vul_arg< vcl_string > arg_img_path,
	vul_arg< vcl_string > arg_old_str,
	vul_arg< vcl_string > arg_new_str,
	vul_arg< bool > arg_3d,
	vul_arg< vcl_string > arg_outfile,
	vul_arg< bool > arg_in_anchor,
	vul_arg< bool > arg_overlap,
	vul_arg< bool > arg_nn,
	vul_arg< bool > arg_normalize, 
	vul_arg< vcl_string > arg_background,
	vul_arg< double > arg_sigma, 
	vul_arg< double > arg_median,
	vul_arg< int > arg_blending,
	vul_arg< bool > arg_denoise,
	vul_arg< bool > arg_write_interproj_image
);

template int mosaic_images_template<unsigned short>(
	vul_arg< vcl_string > arg_xml_file,
	vul_arg< vcl_string > arg_anchor,
	vul_arg< int > arg_channel,
	vul_arg< vcl_string > arg_img_path,
	vul_arg< vcl_string > arg_old_str,
	vul_arg< vcl_string > arg_new_str,
	vul_arg< bool > arg_3d,
	vul_arg< vcl_string > arg_outfile,
	vul_arg< bool > arg_in_anchor,
	vul_arg< bool > arg_overlap,
	vul_arg< bool > arg_nn,
	vul_arg< bool > arg_normalize, 
	vul_arg< vcl_string > arg_background,
	vul_arg< double > arg_sigma, 
	vul_arg< double > arg_median,
	vul_arg< int > arg_blending,
	vul_arg< bool > arg_denoise,
	vul_arg< bool > arg_write_interproj_image
);
