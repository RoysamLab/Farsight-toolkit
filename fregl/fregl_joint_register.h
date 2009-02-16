//:
// \file
// \brief The engine for joint registration
// \author Charlene Tsai
// \date 12/05/2007
// 
//  This class makes sure there exists one transformation between
//  every image pair if one image can be transformed to another image
//  via a series of pairwise transforms.
//
#ifndef _fregl_joint_register_
#define _fregl_joint_register_

#include "itkAffineTransform.h"
#include "itkPoint.h"
#include <vector>
#include <string>

#include <fregl/fregl_reg_record.h>

#include <vbl/vbl_array_2d.h>

struct correspondence
{
  typedef itk::Point< double, 3> PointType;
  PointType from_;
  PointType to_;
};

class fregl_joint_register: public vbl_ref_count
{
public:
  typedef vbl_smart_ptr< fregl_joint_register >  Pointer;
  typedef itk::AffineTransform< double, 3>       TransformType;
  typedef itk::Size<3>                           SizeType;
  typedef std::vector< correspondence >  CorrespondenceList;

  //: Constructor taking only one xml_file 
  fregl_joint_register(std::string const & filename, double scale_multiplier = 0, double error_bound = 0);

  //: Constructor, taking a list of xml filenames of pairwise registratin records
  fregl_joint_register( std::vector<std::string> const & filenames, double scale_multiplier = 0, double error_bound = 0);
  
  //: Constructor, taking a list of pairwise reg records
  fregl_joint_register( std::vector<fregl_reg_record::Pointer> const& reg_records, double scale_multiplier = 0, double error_bound = 0 );

  //: Destructor
  ~fregl_joint_register(){};
  
  //: Construct the connected graph
  //
  //  Every image is treated as a graph, and a transform is the
  //  edge. The connected graph provides the capability to transform
  //  every image to any image in the graph
  bool build_graph(bool mutual_consistency = false);
  
  //: Compute the transformation from every image to the chosen anchor
  bool build_graph(int i, bool mutual_consistency = false); 
  
  //: Return the transform of an image pair
  //
  //  If the transform does not exist, a NULL is returned
  TransformType::Pointer get_transform(int from, int to) const;

  //: Return the list of image names
  std::vector<std::string> const & image_names() const;
  
  //: Return the name of image i
  std::string const & image_name(int i) const;

  //: Return the list of image sizes
  std::vector<SizeType> const & image_sizes() const;

  //: Return the size of image i
  SizeType const & image_size(int i) const;

  //: Return the number of images
  unsigned int number_of_images() const;

  //: Replace the sub-string of the image names
  //
  //  This function is needed when transformations of a different
  //  channel are required for mosaicking. The old_str in the image
  //  name is replace by the new_str.
  void replace_image_name_substr(std::string const & old_str, std::string  const & new_str);

  //: Write the result to an xml file
  void write_xml(std::string const & filename, float multiplier, bool mutual_consistency);

  //: Read the results from an xml file
  void read_xml(std::string const & filename);

  //: Return if two images overlap
  bool is_overlapped(int from, int to);

  //: Return the obj of an image pair
  //
  //  If the transform does not exist or doesn't overlap, a -1 is returned
  double get_obj(int from, int to) const;

  //: Get the registration record
  fregl_reg_record::Pointer get_reg_record(int from, int to);

private: 
  void initialize(std::vector<fregl_reg_record::Pointer> const & reg_records);
  //: Determine if two images overlap
  //
  //  Whether two images overlap is roughly determined by the
  //  translation in x and y, since we assume not much distortion from
  //  other parameters. The overlap property is only used to keep
  //  track of adjacency between two images.
  bool overlapping( int from, int to);
  
  // : Build the graph with mutual consistency
  bool estimate(int anchor);

  //: Build the graph without mutual consistency
  //
  //  Transformations are propagrated in breadth_first_search manner
  void breadth_first_connect( int anchor );

  //: Generate_correspondences
  //
  //  The correspondences are generated assuming that the pairwise
  //  transormations are correct. Points are sampled at a regular
  //  distance in the overlapping area to generate correspondences.
  void generate_correspondences();

private:
  vbl_array_2d<TransformType::Pointer> transforms_; // (from,to)
  vbl_array_2d<double> overlap_;
  vbl_array_2d<double> obj_;
  std::vector<std::string> image_ids_;
  std::vector< SizeType > image_sizes_;
  double scale_multiplier_;
  double error_bound_;
  bool corresp_generated_;
  vbl_array_2d< CorrespondenceList > pairwise_constraints_;
};

#endif
