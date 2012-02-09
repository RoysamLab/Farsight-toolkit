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

//:
// \file
// \brief The engine for joint registration
// \author Charlene Tsai
// \date 12/05/2007
// 
//  This class makes sure there exists one transformation between
//  every image pair if one image can be transformed to another image
//  via a series of pairwise transforms. At the end, only pairs in the
//  same joint graph have transformations. It is possible that
//  multiple sub-graph are found in the given dataset.
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

template < class TPixel >
class fregl_joint_register: public vbl_ref_count
{
public:
	typedef vbl_smart_ptr< fregl_joint_register >  Pointer;
	typedef itk::AffineTransform< double, 3>       TransformType;
	typedef itk::Size<3>                           SizeType;
	typedef std::vector< correspondence >  CorrespondenceList;

	//: Constructor taking only one xml_file 
	fregl_joint_register(std::string const & filename, double scale_multiplier = 0, double error_bound = 1);

	//: Constructor, taking a list of xml filenames of pairwise registratin records
	fregl_joint_register( std::vector<std::string> const & filenames, double scale_multiplier = 0, double error_bound = 1);

	//: Constructor, taking a list of pairwise reg records
	fregl_joint_register( std::vector<fregl_reg_record::Pointer> const& reg_records, double scale_multiplier = 0, double error_bound = 1 );

	//: Destructor
	~fregl_joint_register(){};

	//: Construct the connected graph
	//
	//  Every image is treated as a node, and a transform is the
	//  edge. The connected graph provides the capability to transform
	//  every image to any image in the graph. 
	void build_graph();

	//: Compute the transformation from every image to the chosen anchor
	bool build_graph(int i);

	//: Construct the graph without mutual consistency
	//
	//  exisiting pairwise links are copied over, and
	//  missing linkes are inferred from the exisiting pairwise transforms.
	void infer_graph();

	//: Return the transform of an image pair
	//
	//  If the transform does not exist, a NULL is returned
	TransformType::Pointer get_transform(int from_index, int to_index) const;
	TransformType::Pointer get_transform(std::string from_name, std::string to_name) const;

	//: Return the list of image names
	std::vector<std::string> const & image_names() const;

	//: Return the name of image i
	std::string const & image_name(int i) const;

	//: Return the index of the image name
	//
	//  If the image is not there, it returns -1
	int image_index(std::string name) const;

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
	void write_xml(std::string const & filename, int sub_graphs_built,
		bool gen_temp_stuff = false);

	//: Read the results from an xml file
	void read_xml(std::string const & filename);

	//: Return if two images overlap
	bool is_overlapped(int from, int to) const;

	//: Return the obj of an image pair
	//
	//  If the transform does not exist or doesn't overlap, a -1 is returned
	double get_obj(int from, int to) const;

	//: Get the registration record
	fregl_reg_record::Pointer get_reg_record(int from, int to) const;

	//: Get the error bound
	double get_error_bound() const;

	//: Get the list of adjacent images that overlap with the anchor image
	void get_adjacent_images(std::string anchor_image, std::vector<std::string>& adjacent_images) const;

	//: Return the number of subgraphs in the joint transformation file
	int number_of_subgraphs() const;

	//: Check if two images iare in the same subgraph
	bool in_same_subgraph(int image_index1, int image_index2) const;

private: 
	void initialize(std::vector<fregl_reg_record::Pointer> const & reg_records);

	// : Build the graph with mutual consistency
	bool estimate(int anchor);

	//: Build the graph without mutual consistency
	//
	//  Transformations are propagrated in breadth_first_search manner
	void breadth_first_connect( int anchor );

	//: Mark images which belong to the same sub-graph as the anchor
	//
	//  The computation is in breadth_first_search manner
	void generate_graph_indices( int anchor, int index );

	//: Generate_correspondences
	//
	//  The correspondences are generated assuming that the pairwise
	//  transormations are correct. Points are sampled at a regular
	//  distance in the overlapping area to generate correspondences.
	void generate_correspondences();

private:
	vbl_array_2d<TransformType::Pointer> transforms_; // (from,to)
	vbl_array_2d<double> overlap_; //initially pairwise, finally updated by joint
	vbl_array_2d<double> obj_; //values from the pairwise registration
	std::vector<std::string> image_ids_;
	std::vector< SizeType > image_sizes_;
	bool corresp_generated_;
	double scale_multiplier_;
	double error_bound_;
	vbl_array_2d< CorrespondenceList > pairwise_constraints_;
	std::vector<int> graph_indices_;
	int num_subgraphs_;
};

#endif
