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

#include <iostream>
#include <fregl/fregl_image_manager.h>
#include <fregl/fregl_joint_register.h>
#include <fregl/fregl_space_transformer.h>
#include <fregl/fregl_util.h>
#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageSeriesWriter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkTranslationTransform.h"
#include "itkResampleImageFilter.h"
#include "itkPasteImageFilter.h"


static std::string ToString(double val);

template < class TPixel >
fregl_image_manager< TPixel >::fregl_image_manager(std::string const & xml_filename, std::string const & image_path, std::string const & anchor_image,
bool use_NN) {
    //    std::cout << "Anchor Image " << anchor_image << std::endl;
    use_NN_interpolator = use_NN;
    global_anchor = anchor_image;
    global_image_path = image_path;
    global_joint_register = new fregl_joint_register< TPixel > (xml_filename);
    //    std::cout << "Creating Space Transformer" << std::endl;
    global_space_transformer = new fregl_space_transformer< TPixel > (global_joint_register);
    //    std::cout << "Created Space Transformer" << std::endl;
    //    std::cout << "Setting Anchor" << std::endl;
    global_space_transformer->set_anchor(global_anchor, false, false);
    //    std::cout << "Anchor Set" << std::endl;
    image_names = global_space_transformer->image_names();
    global_origin = global_space_transformer->origin();
    global_size = global_space_transformer->montage_size();
    roi_origin = global_origin;
    roi_index[0] = roi_origin[0];
    roi_index[1] = roi_origin[1];
    roi_index[2] = roi_origin[2];
    roi_size = global_size;
    //    std::cout << "GLOBAL Origin = " << global_origin[0] << "," << global_origin[1] << "," << global_origin[2] << std::endl;
    //    std::cout << "GLOBAL Size = " << global_size[0] << " x " << global_size[1] << " x " << global_size[2] << std::endl;
    global_space_transformer->set_roi(roi_origin, roi_size);
    global_channel = 0;
    global_use_channel = false;
    // Set up caching stuff
    set_cache_buffer_count(6);
    use_file_caching = false;
    QString cdir = QDir::tempPath();
    cache_dir = cdir.toStdString();
}

//: Set the Region of Interest
//
//  Index and size are updated accordingly.

template < class TPixel >
void fregl_image_manager< TPixel >::set_regionofinterest(PointType origin, SizeType size) {
    //    roi_origin = origin;
    // Convert the request from normal space(0,0,0) to anchor space
    //roi_origin is in global space

    roi_origin[0] = origin[0] + global_origin[0];
    roi_origin[1] = origin[1] + global_origin[1];
    roi_origin[2] = origin[2] + global_origin[2];

    roi_size = size;
    std::cout << "origin: " << origin[0] << " " << origin[1] << " " << origin[2] << std::endl;
    std::cout << "global_origin: " << global_origin[0] << " " << global_origin[1] << " " << global_origin[2] << std::endl;
    std::cout << "roi_origin: " << roi_origin[0] << " " << roi_origin[1] << " " << roi_origin[2] << std::endl;

    //if the roi_origin is less than the global_origin in anchor space, then we need to fit it by moving the roi_origin to the global_origin and then lowering the size
    if (roi_origin[0] < global_origin[0]) {
        roi_size[0] -= (global_origin[0] - roi_origin[0]);
        roi_origin[0] = global_origin[0];
    }
    if (roi_origin[1] < global_origin[1]) {
        roi_size[1] -= (global_origin[1] - roi_origin[1]);
        roi_origin[1] = global_origin[1];
    }

    if (roi_origin[2] < global_origin[2]) {
        roi_size[2] -= (global_origin[2] - roi_origin[2]);
        roi_origin[2] = global_origin[2];
    }

    if (roi_origin[0] > global_origin[0] + global_size[0]) roi_origin[0] = global_origin[0] + global_size[0];
    if (roi_origin[1] > global_origin[1] + global_size[1]) roi_origin[1] = global_origin[1] + global_size[1];
    if (roi_origin[2] > global_origin[2] + global_size[2]) roi_origin[2] = global_origin[2] + global_size[2];

    // Make sure size is in range or fix
    if (roi_origin[0] + roi_size[0] > global_origin[0] + global_size[0])
        roi_size[0] = (global_origin[0] + global_size[0]) - roi_origin[0];
    if (roi_origin[1] + roi_size[1] > global_origin[1] + global_size[1])
        roi_size[1] = (global_origin[1] + global_size[1]) - roi_origin[1];
    if (roi_origin[2] + roi_size[2] > global_origin[2] + global_size[2])
        roi_size[2] = (global_origin[2] + global_size[2]) - roi_origin[2];
    roi_index[0] = roi_origin[0];
    roi_index[1] = roi_origin[1];
    roi_index[2] = roi_origin[2];

    /*std::cout << "ROI origin: " << roi_origin[0] << "x" << roi_origin[1] << "x" << roi_origin[2] << std::endl;
    std::cout << "ROI size: " << roi_size[0] << "x" << roi_size[1] << "x" << roi_size[2] << std::endl;*/

    global_space_transformer->set_roi(roi_origin, roi_size);
}

//: Set the Region of Interest
//
//  Index and size are updated accordingly.

template < class TPixel >
void fregl_image_manager< TPixel >::set_regionofinterest(IndexType origin, SizeType size) {
    roi_origin[0] = origin[0];
    roi_origin[1] = origin[1];
    roi_origin[2] = origin[2];
    roi_size = size;
    set_regionofinterest(roi_origin, roi_size);
}
//:  Build the Montage based on the current Region
//   "s_origin" is the new roi origin point
//   "s_size" is the new roi size

template < class TPixel >
void fregl_image_manager< TPixel >::Update() {

    std::string image_name = global_image_path + std::string("/") + image_names[0];
    //    std::cout << "Image names list" << std::endl;
    //    for (unsigned int i = 0; i < image_names.size(); i++) {
    //        std::cout << image_names[i] << std::endl;
    //    }
    ImageTypePointer image, xformed_image;
    //    std::cout << "Composing the final image ..." << std::endl;
    //Create a blank ROI then merge in the images that are in the space
    IndexType source_index;
    typename ImageType::RegionType region;
    source_index[0] = 0;
    source_index[1] = 0;
    source_index[2] = 0.;
    region.SetIndex(source_index);
    region.SetSize(roi_size);
    montage_image = ImageType::New();
    montage_image->SetOrigin(roi_origin);
    montage_image->SetRegions(region);
    montage_image->Allocate();
    montage_image->FillBuffer(0);
    for (unsigned int i = 0; i < image_names.size(); i++) {
        if (global_space_transformer->image_in_roi(i)) {
            image_name = global_image_path + std::string("/") + image_names[i];
            ReadFileRegion(image_name, i, montage_image);
        } else
            continue;
    }
}

//: Return an ITK image pointer to the current Montage
//

template < class TPixel >
typename fregl_image_manager< TPixel >::ImageTypePointer fregl_image_manager< TPixel >::GetOutput() {
    typename ImageType::PointType new_origin;
    typedef itk::ImageDuplicator< ImageType > DuplicatorType;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(montage_image);
    duplicator->Update();
    new_origin[0] = roi_origin[0] - global_origin[0];
    new_origin[1] = roi_origin[1] - global_origin[1];
    new_origin[2] = roi_origin[2] - global_origin[2];
    ImageTypePointer image = duplicator->GetOutput();
    image->SetOrigin(new_origin);
    return image;



    /*    typedef itk::TranslationTransform<double, ImageType::ImageDimension> TranslationTransformType;
    TranslationTransformType::Pointer transform = TranslationTransformType::New();
    TranslationTransformType::OutputVectorType translation;
    std::cout << "global_origin: " << global_origin[0] << " " << global_origin[1] << " " << global_origin[2] << std::endl;
    translation[0] = global_origin[0];
    translation[1] = global_origin[1];
    translation[2] = global_origin[2];
    transform->Translate(translation);

    typedef itk::ResampleImageFilter<ImageType, ImageType, double> ResampleImageFilterType;
    ResampleImageFilterType::Pointer resampleFilter = ResampleImageFilterType::New();
    resampleFilter->SetTransform(transform.GetPointer());
    resampleFilter->SetInput(montage_image);

    ImageType::SizeType size = montage_image->GetLargestPossibleRegion().GetSize();
    resampleFilter->SetSize(size);

    try {
    resampleFilter->Update();
    } catch (itk::ExceptionObject &err) {
    std::cout << "Error in ResampleFilter" << std::endl;
    std::cout << err << std::endl;
    }


    return resampleFilter->GetOutput();*/
}

//: Return an ITK image pointer the the region of interest passed an argument
// This entry is mutex protected for multi threading operations.  It combines
// the set region, update and get output methods.

template < class TPixel >
typename fregl_image_manager< TPixel >::ImageTypePointer fregl_image_manager< TPixel >::MutexGetRegionOfInterest(PointType origin, SizeType size) {
    ImageTypePointer image;
    region_mutex.Lock();
    set_regionofinterest(origin, size);
    Update();
    image = GetOutput();
    region_mutex.Unlock();
    return image;
}


//: Return an ITK image pointer to the Region of Interest in the passed
//  file name.  Return from cache if file is cached otherwise read the file
//  into cache then get the region.

template < class TPixel >
void fregl_image_manager< TPixel >::ReadFileRegion(std::string file_name, int image_index, ImageTypePointer montage_image) {
    ImageTypePointer image, cached_image, xformed_image;
    typename ImageType::RegionType region, source_region, destination_region;
    typename ImageType::IndexType source_index;
    typename ImageType::SizeType source_size;
    typename ImageType::IndexType destination_index;
    typename ImageType::PointType source_origin, save_origin;
    if (!use_caching) {
        //If no caching then read in the image
        image = fregl_util< TPixel >::fregl_util_read_image(file_name, global_use_channel, global_channel, false);
        xformed_image = global_space_transformer->transform_image_roi(image, image_index, 0, use_NN_interpolator);
        RegionConstIterator inputIt(xformed_image, xformed_image->GetLargestPossibleRegion());
        if (!xformed_image) return;
        RegionIterator outputIt(montage_image, montage_image->GetLargestPossibleRegion());

        for (inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd();
                ++inputIt, ++outputIt) {
            outputIt.Set(vnl_math_max(outputIt.Get(), inputIt.Get()));
        }
        return;
    } else {
        //Using caching
        //        std::cout << "Here 1" << std::endl;
        cache_time++;
        if (!is_cached[image_index]) {
            // If we are file caching and its not already cached lets try to read it from the disk first.  If it is there
            // then use it otherwise read it in and transform it.
            if (use_file_caching) {
                image = cache_read_image(image_index);
                if (image) {
                    is_cached[image_index] = true;
                    is_cached_on_disk[image_index] = true;
                    cache_slot[image_index] = get_next_slot();
                    cached_images[cache_slot[image_index]] = image;
                }
            }
            // If we did not get the image from the disk then read and transform
            if (!is_cached[image_index]) {
                is_cached[image_index] = true;
                cache_slot[image_index] = get_next_slot();
                image = fregl_util< TPixel >::fregl_util_read_image(file_name, global_use_channel, global_channel, false);
                cached_images[cache_slot[image_index]] = global_space_transformer->transform_image_whole(image, image_index, 0, use_NN_interpolator);
                //        std::cout << "Here 4" << std::endl;
            }
        }
    }
    //If we are using file caching and we did not read the transformed file from
    //the disk then write it out.
    if (use_file_caching && !is_cached_on_disk[image_index])
        is_cached_on_disk[image_index] = cache_write_image(image_index, cached_images[cache_slot[image_index]]);
    cached_image = cached_images[cache_slot[image_index]];
    cache_last_used[cache_slot[image_index]] = cache_time;


    //Now paste in the part of the cached file that is in the Region of interest
    //Calculate the region of the cached image that is in the Region of interest
    source_size = cached_image->GetLargestPossibleRegion().GetSize();
    source_origin = cached_image->GetOrigin();
    save_origin = source_origin;
    //    std::cout << "Starting Source origin: " << source_origin[0] << " " << source_origin[1] << " " << source_origin[2] << std::endl;
    //    std::cout << "Calculating index and size" << std::endl;

    //If you have to move the origin clamp down the size as well.
    if (source_origin[0] < roi_origin[0]) {
        source_size[0] = (source_origin[0] + source_size[0]) - roi_origin[0];
        source_origin[0] = roi_origin[0];
    }
    if (source_origin[1] < roi_origin[1]) {
        source_size[1] = (source_origin[1] + source_size[1]) - roi_origin[1];
        source_origin[1] = roi_origin[1];
    }
    if (source_origin[2] < roi_origin[2]) {
        source_size[2] = (source_origin[2] + source_size[2]) - roi_origin[2];
        source_origin[2] = roi_origin[2];
    }
    //This next adjustment should not happen since it would mean that the image
    // is not in the roi but ...
    if (source_origin[0] > roi_origin[0] + roi_size[0]) {
        std::cout << "Should not be here 0 " << cache_slot[image_index] << std::endl;
        source_size[0] = 0;
        source_origin[0] = roi_origin[0] + roi_size[0];
    }
    if (source_origin[1] > roi_origin[1] + roi_size[1]) {
        std::cout << "Should not be here 1 " << cache_slot[image_index] << std::endl;
        source_size[1] = 0;
        source_origin[1] = roi_origin[1] + roi_size[1];
    }
    if (source_origin[2] > roi_origin[2] + roi_size[2]) {
        std::cout << "Should not be here 2 " << cache_slot[image_index] << std::endl;
        source_size[2] = 0;
        source_origin[2] = roi_origin[2] + roi_size[2];
    }
    // Make sure size is in range or fix
    //First make sure the size is with in the destination
    if (source_origin[0] + source_size[0] > roi_origin[0] + roi_size[0]) {
        //        std::cout << "Here0 " << cache_slot[image_index] << std::endl;
        source_size[0] = (roi_origin[0] + roi_size[0]) - source_origin[0];
    }
    if (source_origin[1] + source_size[1] > roi_origin[1] + roi_size[1]) {
        //        std::cout << "Here1 " << cache_slot[image_index] << std::endl;
        source_size[1] = (roi_origin[1] + roi_size[1]) - source_origin[1];
    }
    if (source_origin[2] + source_size[2] > roi_origin[2] + roi_size[2]) {
        //        std::cout << "Here2 " << cache_slot[image_index] << std::endl;
        source_size[2] = (roi_origin[2] + roi_size[2]) - source_origin[2];
    }
    source_index[0] = source_origin[0] - save_origin[0];
    source_index[1] = source_origin[1] - save_origin[1];
    source_index[2] = source_origin[2] - save_origin[2];
    //Now make sure that the size is not larger than the image
    destination_index[0] = (source_origin[0] - roi_origin[0]);
    destination_index[1] = (source_origin[1] - roi_origin[1]);
    destination_index[2] = (source_origin[2] - roi_origin[2]);
    //    std::cout << "Source Index: " << source_index[0] << " " << source_index[1] << " " << source_index[2] << std::endl;
    //    std::cout << "Source Origin: " << source_origin[0] << " " << source_origin[1] << " " << source_origin[2] << std::endl;
    //    std::cout << "Destination Index: " << destination_index[0] << " " << destination_index[1] << " " << destination_index[2] << std::endl;
    //    std::cout << "Source Size: " << source_size[0] << " " << source_size[1] << " " << source_size[2] << std::endl;
    source_region.SetIndex(source_index);
    source_region.SetSize(source_size);
    destination_region.SetIndex(destination_index);
    destination_region.SetSize(source_size);
    cached_image->SetRequestedRegion(source_region);
    montage_image->SetRequestedRegion(destination_region);

    RegionConstIterator inputIt(cached_image, source_region);
    RegionIterator outputIt(montage_image, destination_region);

    for (inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd();
            ++inputIt, ++outputIt) {
        outputIt.Set(vnl_math_max(outputIt.Get(), inputIt.Get()));
    }
    //    return PasteFilter->GetOutput();
    montage_image->SetRequestedRegion(montage_image->GetLargestPossibleRegion());
    return;
}

//: Set the image channel to montage
//

template < class TPixel >
void fregl_image_manager< TPixel >::set_channel(int channel) {
    global_channel = channel;
    global_use_channel = true;
}

//: Return the global space origin in anchor space
//

template < class TPixel >
typename fregl_image_manager< TPixel >::PointType fregl_image_manager< TPixel >::get_global_origin() {
    return global_origin;
}

//: Return return the global space size
//

template < class TPixel >
typename fregl_image_manager< TPixel >::SizeType fregl_image_manager< TPixel >::get_global_size() {
    return global_size;
}

//: Return the current Region origin in normalized (0,0,0) space
//

template < class TPixel >
typename fregl_image_manager< TPixel >::PointType fregl_image_manager< TPixel >::get_region_origin() {
    typename ImageType::PointType temp;
    temp[0] = roi_origin[0] - global_origin[0];
    temp[1] = roi_origin[1] - global_origin[1];
    temp[2] = roi_origin[2] - global_origin[2];
    return temp;
}

//: Return the current Region Size
//

template < class TPixel >
typename fregl_image_manager< TPixel >::SizeType fregl_image_manager< TPixel >::get_region_size() {
    return roi_size;
}

//: Return the anchor image name

template < class TPixel >
std::string const
fregl_image_manager< TPixel >::get_anchor_name() {
    return global_anchor;
}

//: Return a pointer the space transformer

template < class TPixel >
typename fregl_space_transformer< TPixel >::Pointer fregl_image_manager< TPixel >::get_space_transformer() {
    return global_space_transformer;
}

//: Set up the cache buffer count.  If 0 Turn off caching
// Will reset all cache buffers when called.
// Returns the state of caching (off or on)

template < class TPixel >
void fregl_image_manager< TPixel >::set_cache_buffer_count(int count) {
    cache_size = count;
    cache_time = 0;
    cached_images.clear();
    cache_last_used.clear();
    is_cached.clear();
    cache_slot.clear();
    if (count > 0) {
        use_caching = true;
        for (unsigned int i = 0; i < cache_size; i++) {
            cached_images.push_back(NULL);
            cache_last_used.push_back(0);
        }
        for (unsigned int i = 0; i < image_names.size(); i++) {
            is_cached.push_back(false);
            is_cached_on_disk.push_back(false);
            cache_slot.push_back(-1);
        }
    } else use_caching = false;
}

//: Set disk caching on or off.

template < class TPixel >
void fregl_image_manager< TPixel >::set_file_caching(bool use) {
    use_file_caching = use;
    //reset the cashed on disk information when we change the state
    for (unsigned int i = 0; i < is_cached_on_disk.size(); i++) {
        is_cached_on_disk[i] = false;
    }
}

//: Set the directory for the disk caching.  The default is current
// directory "."

template < class TPixel >
void fregl_image_manager< TPixel >::set_file_cache_dir(std::string use_dir) {
    cache_dir = use_dir;
    //reset the cashed on disk information when we change the directory since
    //it may no longer be accurate.
    for (unsigned int i = 0; i < is_cached_on_disk.size(); i++) {
        is_cached_on_disk[i] = false;
    }
}

template < class TPixel >
int fregl_image_manager< TPixel >::get_next_slot() {
    // First look for unused slots and return it if there.
    int oldest = 0;
    //        std::cout << "Here 6" << std::endl;
    for (unsigned int i = 0; i < cached_images.size(); i++) {
        //        std::cout << cache_last_used [i] << std::endl;
        if (cache_last_used[i] == 0) {
            //            std::cout << "Returning slot " << i << cache_time << std::endl;
            return i;
        } else {
            // Keep track of the slot with the oldest time;
            if (cache_last_used[oldest] > cache_last_used[i]) oldest = i;
        }
    }
    //Clear old pointer to the cache.
    for (unsigned int i = 0; i < is_cached.size(); i++) {
        if (cache_slot[i] == oldest) {
            is_cached[i] = false;
            cache_slot[i] = -1;
        }
    }
    //    std::cout << "Returning slot " << oldest << cache_time << std::endl;
    return oldest;
}

template < class TPixel >
bool fregl_image_manager< TPixel >::cache_write_image(int image_index, ImageTypePointer t_image) {
    // Write a transformed image to disk for caching.  The file name is made up of
    // "cache_anchorimagename_transformedimagename.mhd"
    std::string name = cache_dir + std::string("/cache_") + vul_file::strip_extension(global_anchor)
            + std::string("_") + vul_file::strip_extension(image_names[image_index]) + std::string(".mhd");
    typedef itk::ImageFileWriter< ImageType > WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    //    std::cout << "File name: " << name << std::endl;
    writer->SetFileName(name);
    writer->SetInput(t_image);
    try {
        writer->Update();
    } catch (itk::ExceptionObject& e) {
        vcl_cout << e << vcl_endl;
        return false;
    }
    return true;
}

template < class TPixel >
typename fregl_image_manager< TPixel >::ImageTypePointer fregl_image_manager< TPixel >::cache_read_image(int image_index) {
    // Read a transformed cache file from the disk
    // return NULL if the file is not there or
    // can not be read for some other reason.
    std::string name = cache_dir + std::string("/cache_") + vul_file::strip_extension(global_anchor)
            + std::string("_") + vul_file::strip_extension(image_names[image_index]) + std::string(".mhd");
    //    std::cout << "File name: " << name << std::endl;
    typedef itk::ImageFileReader< ImageType > ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(name);
    try {
        reader->Update();
    } catch (itk::ExceptionObject& e) {
        std::cerr << "Cached image does not exist on disk at " << name << ", transforming a new one" << std::endl;
        return NULL;
    }
    return reader->GetOutput();
}

//Explicit Instantiation
template class fregl_image_manager< unsigned char >;
template class fregl_image_manager< unsigned short >;
