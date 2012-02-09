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

#include "fregl_joint_register.h"
#include "fregl_util.h"

#include <vnl/vnl_math.h>
#include <vul/vul_timer.h>

#include <rrel/rrel_muset_obj.h>

#include <fstream>

static std::string ToString(double val);

template < class TPixel >
fregl_joint_register< TPixel >::
fregl_joint_register(std::string const & xml_filename, double scale_multiplier,  double error_bound )
: corresp_generated_(false), scale_multiplier_(scale_multiplier), error_bound_(error_bound)
{
	read_xml(xml_filename);
}

template < class TPixel >
fregl_joint_register< TPixel >::
fregl_joint_register( std::vector<std::string> const & xml_filenames, double scale_multiplier, double error_bound) : corresp_generated_(false), scale_multiplier_(scale_multiplier), error_bound_(error_bound)
{
	std::vector<fregl_reg_record::Pointer> reg_records;
	for (unsigned int i = 0; i<xml_filenames.size(); i++) {
		fregl_reg_record::Pointer reg_record = new fregl_reg_record();
		reg_record->read_xml(xml_filenames[i]);
		reg_records.push_back( reg_record );
	}
	initialize( reg_records );

	/*std::cout << "Printing out transforms" << std::endl;
	for (unsigned int j = 0; j<transforms_.cols(); j++) 
	{
	for (unsigned int i = 0; i<transforms_.rows(); i++) 
	{
	TransformType::ParametersType params = transforms_(i,j)->GetParameters(); //temporary
	std::cout << i << " " << j << " " << params[9] << " " << params[10] << " " << params[11] << std::endl;
	}
	}*/
}

template < class TPixel >
fregl_joint_register< TPixel >::
fregl_joint_register( std::vector<fregl_reg_record::Pointer> const & reg_records, double scale_multiplier, double error_bound) : corresp_generated_(false), scale_multiplier_(scale_multiplier), error_bound_(error_bound)
{
	initialize( reg_records );
}

template < class TPixel >
void 
fregl_joint_register< TPixel >::
initialize(std::vector<fregl_reg_record::Pointer> const & reg_records)
{
	std::cout<<"Number of records = "<<reg_records.size()<<std::endl;

	// Get the image_ids_
	for (unsigned int i = 0; i<reg_records.size(); i++) {
		//check if the image ids are already in the list
		bool to_image_existed = false;
		bool from_image_existed = false;
		for (unsigned int j = 0; j<image_ids_.size(); j++) {
			if (reg_records[i]->from_image() == image_ids_[j])
				from_image_existed = true;
			if (reg_records[i]->to_image() == image_ids_[j])
				to_image_existed = true;
		}
		if ( !to_image_existed ) {
			image_ids_.push_back( reg_records[i]->to_image() );
			image_sizes_.push_back( reg_records[i]->to_image_size());
		}
		if ( !from_image_existed ) {
			image_ids_.push_back( reg_records[i]->from_image() );
			image_sizes_.push_back( reg_records[i]->from_image_size());
		}
	}

	// debugging:
	std::cout<<"Total number of images = "<<image_ids_.size()<<std::endl;
	for (unsigned int i = 0; i< image_ids_.size(); i++) {
		std::cout<<"Image "<<image_ids_[i]<<std::endl;
	}

	// Get the transformations. element (i,j) is mapping from image i to
	// image j
	std::cout << "Number of images: " << image_ids_.size() << std::endl;
	transforms_.resize(image_ids_.size(), image_ids_.size());
	transforms_.fill( NULL );
	overlap_.resize(image_ids_.size(), image_ids_.size());
	overlap_.fill( 0 );
	obj_.resize(image_ids_.size(), image_ids_.size());
	obj_.fill( 1 );
	graph_indices_.resize(image_ids_.size(),0);

	// If scale_multiplier_ is set, run muse scale estimator to compute
	// the error scale. Set the error_bound_ to be
	// error_scale*scale_multiplier_
	//
	std::vector<double> errors;
	errors.reserve(reg_records.size());
	std::cout<<"Error_bound = "<<error_bound_<<std::endl;
	if (scale_multiplier_ > 0) {
		for (unsigned int i = 0; i<reg_records.size(); i++) {
			//obj is in the range of [0,1]. It is expected to be the value
			//computed using NormalizedCorssCorrelation metric shifted by 1
			//so that it is the range of [0,1]

			//errors.push_back( 1+reg_records[i]->obj());
			if (reg_records[i]->obj() < 1)
				errors.push_back( reg_records[i]->obj()); 
		}

		rrel_muset_obj muse_obj( 0, false );
		double scale = muse_obj.scale(errors.begin(), errors.end());
		//error_bound_ = scale*scale_multiplier_-1;
		error_bound_ = scale*scale_multiplier_;
		std::cout<<"scale_multiplier_ = "<<scale_multiplier_<<"\n";
		std::cout<<"Estimated muse scale = "<<scale<<", error_bound = "
			<<error_bound_<<std::endl;
		// Compute the average error from the accepted pairs
		float sum = 0;
		int count = 0;
		for (unsigned int i = 0; i<errors.size(); i++) { 
			if (errors[i]<error_bound_) {
				sum += errors[i];
				count++;
			}
		}
		std::cout<<"average error = "<<sum/count<<std::endl;
	}

	// Transfer the information from the reg_records to the proper data
	// structures
	//
	for (unsigned int i = 0; i<reg_records.size(); i++) {
		int to_image_index = 0;
		int from_image_index = 0;
		for (unsigned int j = 0; j<image_ids_.size(); j++) {
			if (reg_records[i]->from_image() == image_ids_[j])
				from_image_index = j;
			if (reg_records[i]->to_image() == image_ids_[j])
				to_image_index = j;
		}

		// Only image pairs with valid xforms are imported to transforms_
		// for joint registration. Those transforms define the valid links
		// in the joint graph. obj_contains the values from pairwise
		// registration. overlap_ is initialized from pairwise, but later
		// updated when joint results are available.
		if (reg_records[i]->obj()<= error_bound_)
			transforms_(from_image_index, to_image_index) = reg_records[i]->transform();
		obj_(from_image_index, to_image_index) = reg_records[i]->obj();
		overlap_(from_image_index, to_image_index) = reg_records[i]->overlap();

		if (reg_records[i]->obj() > error_bound_ && reg_records[i]->obj()< 1)
			std::cout<<"Eliminated pair "<<reg_records[i]->from_image()<<" to "<<reg_records[i]->to_image()<<" with obj="<<reg_records[i]->obj()<<std::endl;

	}

	// Set the diagonal elements to identity transformation
	for (unsigned int i = 0; i<transforms_.rows(); i++) {
		transforms_(i,i) = TransformType::New();
		transforms_(i,i)->SetIdentity();
		overlap_(i,i) = 1;
		obj_(i,i) = 0;
	}

	// Set the inverse transform of pair (i,j) to pair (j,i) if the
	// transform does not exist.
	for (unsigned int i = 0; i<transforms_.rows(); i++) {
		for (unsigned int j = 0; j<transforms_.cols(); j++) {
			/*if (!transforms_(i,j))
			std::cout << "transform " << i << " to " << j << " doesnt exist" << std::endl;*/
			if ( transforms_(i,j) && !transforms_(j,i)) {
				TransformType::Pointer new_xform = TransformType::New();
				if (!transforms_(i,j)->GetInverse(new_xform)) 
					std::cerr<<" Inverse does not exist! "<<std::endl;
				else {
					transforms_(j,i) = new_xform;
					overlap_(j,i) = overlap_(i,j);
					obj_(j,i) = obj_(i,j);
				}
			}
		}
	}

	// build the array which contains the sub_graph indices. Images
	// belonging to the same subgraph should have the same index
	num_subgraphs_ = 0;
	for (unsigned int i = 0; i<transforms_.rows(); i++) {
		if (!graph_indices_[i]) 
			generate_graph_indices(i, ++num_subgraphs_);
	}

	std::cout<<"End of Initialization"<<std::endl;
}

template < class TPixel >
void
fregl_joint_register< TPixel >::
build_graph()
{
	vul_timer timer;
	timer.mark();

	/* This piece of code is moved initialize(.)
	// Build the array which contains the sub_graph indices. Images
	// belonging to the same subgraph should have the same index
	graph_indices_.resize(image_ids_.size(),0);
	int index = 0;
	for (unsigned int i = 0; i<transforms_.rows(); i++) {
	if (!graph_indices_[i]) 
	generate_graph_indices(i, ++index);
	}
	*/

	/*
	std::cout << "Pre-corruption" << std::endl;
	for (unsigned int j = 0; j< transforms_.cols(); j++) 
	{
	for (unsigned int i = 0; i< transforms_.rows(); i++) 
	{
	if (transforms_(i,j))
	{
	TransformType::ParametersType params = transforms_(i,j)->GetParameters(); //temporary
	std::cout << i << " to " << j << ": " << params[9] << " " << params[10] << " " << params[11] << std::endl;
	}
	}
	}
	*/

	// Now compute the transformations by treating each image as the anchor. 
	for (unsigned int i = 0; i<transforms_.rows(); i++) {
		std::cout<<"Building the graph for image "<<image_ids_[i]<<std::endl;
		if (!build_graph(i))
			std::cerr<<"ERROR: No graph built for image "<<image_ids_[i]<<std::endl;
	}

	/*
	std::cout << "Post-corruption" << std::endl;
	for (unsigned int j = 0; j< transforms_.cols(); j++) 
	{
	for (unsigned int i = 0; i< transforms_.rows(); i++) 
	{
	if (transforms_(i,j))
	{
	TransformType::ParametersType params = transforms_(i,j)->GetParameters(); //temporary
	std::cout << i << " to " << j << ": " << params[9] << " " << params[10] << " " << params[11] << std::endl;
	}
	}
	}
	*/

	std::cout << "Timing: Joint registration in  ";
	timer.print( std::cout );
	std::cout<<std::endl;

}

template < class TPixel >
bool 
fregl_joint_register< TPixel >::
build_graph(int anchor)
{
	// If mutual consistency not required, simply use breadth first
	// search to propagate the transformation. Otherwise, perform bundle
	// adjustment to impose mutual consistency outside the anchor space.

	if ( !corresp_generated_ ) 
		generate_correspondences();
	if (!estimate( anchor )) 
		return false; 

	/*
	if ( mutual_consistency ) {
	if ( !corresp_generated_ ) generate_correspondences();
	if (!estimate( anchor )) return false; 
	}
	else breadth_first_connect( anchor );
	*/

	// update overlaps
	for (unsigned int from = 0; from<transforms_.rows(); from++) {
		overlap_(from, anchor) = fregl_util< TPixel >::fregl_util_overlap(transforms_(from, anchor), image_sizes_[from], image_sizes_[anchor]);
	}

	/*for (unsigned int from = 0; from<transforms_.rows(); from++) 
	{
	if (overlap_(anchor, from) != 0.0)
	std::cout << overlap_(anchor, from) << " overlap from anchor: " << anchor << " to: " << from << std::endl;
	}*/

	return true;
}

template < class TPixel >
void
fregl_joint_register< TPixel >::
infer_graph()
{
	//  Compute the transformations by treating each image as the anchor.
	for (unsigned int i = 0; i<transforms_.rows(); i++) {
		std::cout<<"Infer the transformations to image "<<image_ids_[i]<<std::endl;
		breadth_first_connect( i );
	}
}

template < class TPixel >
fregl_joint_register< TPixel >::TransformType::Pointer 
fregl_joint_register< TPixel >::
get_transform(int from, int to) const
{
	return transforms_(from, to);
}

template < class TPixel >
fregl_joint_register< TPixel >::TransformType::Pointer
fregl_joint_register< TPixel >::
get_transform(std::string from_name, std::string to_name) const
{
	int to_index = -1;
	int from_index = -1;
	for (unsigned i = 0; i<image_ids_.size(); i++){
		if (from_name == image_ids_[i]) 
			from_index = i;
		if (to_name == image_ids_[i]) 
			to_index = i;
	}
	if (to_index<0) {
		std::cerr<<"To_image is not found!!!"<<std::endl; 
		return NULL;
	}
	if (from_index<0) {
		std::cerr<<"From_image is not found!!!"<<std::endl; 
		return NULL;
	}
	return transforms_(from_index, to_index);
}

template < class TPixel >
double 
fregl_joint_register< TPixel >::get_obj(int from, int to) const
{
	return obj_(from, to);
}

template < class TPixel >
double 
fregl_joint_register< TPixel >::get_error_bound() const
{
	return error_bound_;
}

template < class TPixel >
void
fregl_joint_register< TPixel >::get_adjacent_images(std::string anchor_image, std::vector<std::string>& adjacent_images) const
{
	int index = image_index(anchor_image);
	assert(index >= 0);

	for (unsigned int i = 0; i<image_ids_.size(); i++) {
		if (overlap_[index][i] > 0 && index != (int)i) {
			adjacent_images.push_back(image_ids_[i]);
		}
	}
}

template < class TPixel >
std::vector<std::string> const &
fregl_joint_register< TPixel >::
image_names() const
{
	return image_ids_;
}

template < class TPixel >
int
fregl_joint_register< TPixel >::
image_index(std::string name) const
{
	for (unsigned int i = 0; i<image_ids_.size(); i++)
		if (name == image_ids_[i]) {
			return i;
		}
		return -1;
}

template < class TPixel >
std::vector< typename fregl_joint_register< TPixel >::SizeType > const &
fregl_joint_register< TPixel >::
image_sizes() const
{
	return image_sizes_;
}

template < class TPixel >
std::string const & 
fregl_joint_register< TPixel >::
image_name(int i) const
{
	return image_ids_[i];
}

template < class TPixel >
typename fregl_joint_register< TPixel >::SizeType const &
fregl_joint_register< TPixel >::
image_size(int i) const
{
	return image_sizes_[i];
}

template < class TPixel >
unsigned int
fregl_joint_register< TPixel >::
number_of_images() const
{
	return image_sizes_.size();
}

template < class TPixel >
void 
fregl_joint_register< TPixel >::
replace_image_name_substr(std::string const& old_str, std::string const& new_str)
{
	for (unsigned int i = 0; i<image_ids_.size(); i++) {
		size_t pos = image_ids_[i].find(old_str);
		if (pos == std::string::npos) {
			std::cerr<<old_str<<" cannot be found for replacement!"<<std::endl;
		}
		else {
			image_ids_[i].replace(image_ids_[i].find(old_str), old_str.length(), new_str);
		}
	}
}

template < class TPixel >
bool 
fregl_joint_register< TPixel >::
is_overlapped(int from, int to) const
{
	if (overlap_(from,to) > 0) return true;
	else return false;
}

template < class TPixel >
int
fregl_joint_register< TPixel >::
number_of_subgraphs() const
{
	return num_subgraphs_;
}

template < class TPixel >
bool
fregl_joint_register< TPixel >::
in_same_subgraph(int image_index1, int image_index2) const
{
	if (graph_indices_[image_index1] == graph_indices_[image_index2])
		return true;
	else return false;
}

template < class TPixel >
void 
fregl_joint_register< TPixel >::
generate_graph_indices( int anchor, int index ) 
{
	bool all_explored = false;

	// This is a breadth-first search 
	std::vector<bool> explored(transforms_.rows(), false);
	explored[anchor] =  true;
	graph_indices_[anchor] = index;
	std::vector<int> connected; //to keep track of new connections
	while (!all_explored) {
		//check for new connections which are not yet explored
		for (unsigned int i = 0; i<transforms_.rows(); i++) {
			if ((unsigned int)anchor != i && transforms_(i,anchor) && !explored[i]) {
				connected.push_back(i);
			}
		}

		if (connected.empty()) {
			all_explored = true;
			continue;
		}
		else {
			anchor = connected.back();
			graph_indices_[anchor] = index;
			connected.pop_back();
			explored[anchor] = true;
		}

		// Explore the new connections
		/*
		while (!connected.empty()) {
		anchor = connected.back();
		graph_indices_[anchor] = index;
		connected.pop_back();
		explored[anchor] = true;
		}
		*/
	}
}

// For this connection, exisiting pairwise links are copied over, and
// missing linkes are inferred from the exisiting pairwise transforms.
template < class TPixel >
void 
fregl_joint_register< TPixel >::
breadth_first_connect( int anchor ) 
{
	// check if all transforms are set for this reference image
	bool all_explored = true;
	for (unsigned int i = 0; i<transforms_.rows(); i++) {
		if (i != (unsigned int)anchor && !transforms_(i, anchor)) 
			all_explored  = false;
	}
	if (all_explored) return;

	all_explored = false;

	// This is b-first 
	std::vector<bool> explored(transforms_.rows(), false);
	explored[anchor] =  true;
	std::vector<int> connected; //to keep track of new connections
	while (!all_explored) {
		//check for new connections which are not yet explored
		for (unsigned int i = 0; i<transforms_.rows(); i++) {
			if ((unsigned int)anchor != i && transforms_(i,anchor) && !explored[i]) {
				connected.push_back(i);
			}
		}

		if (connected.empty()) {
			all_explored = true;
			continue;
		}

		// Explore the new connections 
		while (!connected.empty()) {
			int from_index = connected.back();
			connected.pop_back();
			for (unsigned int i = 0; i<transforms_.rows(); i++) {
				if ( i == (unsigned int)anchor || i == (unsigned int)from_index ) continue;

				if (transforms_(i,from_index) && !transforms_(i,anchor)) {
					transforms_(i,anchor) = TransformType::New();
					transforms_(i,anchor)->SetIdentity();
					transforms_(i,anchor)->Compose(transforms_(from_index,anchor));
					transforms_(i,anchor)->Compose(transforms_(i,from_index), true);
					transforms_(anchor,i) = TransformType::New();
					transforms_(i,anchor)->GetInverse( transforms_(anchor,i) );

					//if ( overlapping(i, anchor) )
					//  overlap_(i, anchor) = 1;
					//if ( overlapping( anchor, i) )
					//  overlap_(anchor, i) = 1;

				}
			}
			explored[from_index] = true;
		}
	}
}

template < class TPixel >
bool 
fregl_joint_register< TPixel >::
estimate(int anchor)
{
	std::vector<int> indices_for_other_images;
	indices_for_other_images.reserve(image_ids_.size());
	for (unsigned int i = 0; i<image_ids_.size();i++)
		if (i!=anchor && graph_indices_[anchor]==graph_indices_[i]) {
			indices_for_other_images.push_back(i);
		}

		//get the size of X_D (direct constraints) and X_I (indirect constraints)
		//
		unsigned int size = indices_for_other_images.size();
		int X_D_dim = 0, X_I_dim = 0;
		for (unsigned int i = 0; i<size; i++) {
			if (pairwise_constraints_[indices_for_other_images[i]][anchor].size() > 0) 
				X_D_dim += pairwise_constraints_[indices_for_other_images[i]][anchor].size();
			for (unsigned int j = i+1; j<size; j++) {
				if (pairwise_constraints_[indices_for_other_images[i]][indices_for_other_images[j]].size() > 0)
					X_I_dim += pairwise_constraints_[indices_for_other_images[i]][indices_for_other_images[j]].size();
			}
		}

		if ( X_D_dim == 0 ) {
			std::cerr<< "fregl_joint_register< TPixel >::estimate --- ERROR \n"
				<< "    no direct constraints\n"<<std::endl;
			return false;
		}

		// 1.  Center the data. Please note that the points are not centered
		// along the z-dimension since the z-dimension is often smaller than
		// 120.

		std::vector<vnl_vector_fixed<double, 3> > centers(size);
		vnl_vector_fixed<double, 3> anchor_center;
		double  to_xc = 0, to_yc = 0, to_zc = 0, to_sum = 0;
		for (unsigned int i = 0; i<size; i++) {
			double xc = 0, yc = 0, zc = 0, sum = 0;
			for (unsigned int j = 0; j<image_ids_.size(); j++) {
				if (pairwise_constraints_[indices_for_other_images[i]][j].size() == 0) continue;
				CorrespondenceList const& corresps = pairwise_constraints_[indices_for_other_images[i]][j];
				for (unsigned int k = 0; k<corresps.size(); k++) {
					correspondence const & corresp = corresps[k];
					xc += corresp.from_[0];
					yc += corresp.from_[1];
					zc += corresp.from_[2];
				}
				sum += corresps.size();

				if (j == anchor) {
					for (unsigned int k = 0; k<corresps.size(); k++) {
						correspondence const & corresp = corresps[k];
						to_xc += corresp.to_[0];
						to_yc += corresp.to_[1];
						to_zc += corresp.to_[2];
					}
					to_sum += corresps.size();
				}
			}
			centers[i][0] = xc/double(sum);
			centers[i][1] = yc/double(sum);
			centers[i][2] = 0; //zc/double(sum);
		}
		anchor_center[0] = to_xc/to_sum;
		anchor_center[1] = to_yc/to_sum;
		anchor_center[2] = 0;

		// 2. Convert X_D and W_D to sum_prod, but only compute the upper
		//  triangle, and the rhs matrix to sum_rhs

		int dim = size*4;
		vnl_matrix<double> sum_prod(dim, dim, 0.0), sum_rhs(dim, 3, 0.0);
		vnl_vector<double> cstr(4), cstr2(4);

		for (unsigned int i = 0; i<size;i++)
			if (pairwise_constraints_[indices_for_other_images[i]][anchor].size()>0) {
				int start_row = 4*i;
				CorrespondenceList const& corresps = pairwise_constraints_[indices_for_other_images[i]][anchor];

				for (unsigned int k = 0; k< corresps.size();k++) {
					correspondence::PointType p = corresps[k].from_;
					correspondence::PointType q = corresps[k].to_;
					//double s = scales_[indices_in_graph[i]][anchor];
					//double w = corresps[k]->weight()/(s*s);

					double  px = p[0] - centers[i][0];
					double  py = p[1] - centers[i][1];
					double  pz = p[2] - centers[i][2];
					cstr[0] = 1.0;
					cstr[1] = px;
					cstr[2] = py;
					cstr[3] = pz;
					//std::cout<<"px, py, pz = "<<px<<","<<py<<","<<pz<<std::endl;
					//cstr[3] = px * px;
					//cstr[4] = px * py;
					//cstr[5] = py * py;

					for (int g = 0; g < 4; g++)
						for (int h = g; h < 4; h++)
							sum_prod(g + start_row,h + start_row) += cstr[g] * cstr[h];

					//rhs matrix
					double  qx = q[0] - anchor_center[0];
					double  qy = q[1] - anchor_center[1];
					double  qz = q[2] - anchor_center[2];

					sum_rhs(start_row, 0) += cstr[0] * qx;
					sum_rhs(start_row+1, 0) += cstr[1] * qx;
					sum_rhs(start_row+2, 0) += cstr[2] * qx;
					sum_rhs(start_row+3, 0) += cstr[3] * qx;
					sum_rhs(start_row, 1) += cstr[0] * qy;
					sum_rhs(start_row+1, 1) += cstr[1] * qy;
					sum_rhs(start_row+2, 1) += cstr[2] * qy;
					sum_rhs(start_row+3, 1) += cstr[3] * qy;
					sum_rhs(start_row, 2) += cstr[0] * qz;
					sum_rhs(start_row+1, 2) += cstr[1] * qz;
					sum_rhs(start_row+2, 2) += cstr[2] * qz;
					sum_rhs(start_row+3, 2) += cstr[3] * qz;

				}
			}

			//compute X_I & W_I if indirect constraints exist
			//
			if ( X_I_dim > 0){
				for (unsigned int i = 0; i<size-1;i++)
					for (unsigned int j = i+1; j<size; j++) 
						if (pairwise_constraints_[indices_for_other_images[i]][indices_for_other_images[j]].size()>0) {
							CorrespondenceList const& corresps = pairwise_constraints_[indices_for_other_images[i]][indices_for_other_images[j]];
							int i_start = 4*i, j_start = 4*j;
							for ( unsigned int k = 0; k< corresps.size();k++) {
								correspondence::PointType p = corresps[k].from_;
								correspondence::PointType q = corresps[k].to_;
								//double s = scales_[indices_in_graph[i]][indices_in_graph[j]];
								//double w = corresps[k]->weight()/(s*s);

								double  px = p[0] - centers[i][0];
								double  py = p[1] - centers[i][1];
								double  pz = p[2] - centers[i][2];

								cstr[0] = 1.0;
								cstr[1] = px;
								cstr[2] = py;
								cstr[3] = pz;

								for (int g = 0; g < 4; g++)
									for (int h = g; h < 4; h++)
										sum_prod(g + i_start,h + i_start) += cstr[g] * cstr[h];

								double  qx = q[0] - centers[j][0];
								double  qy = q[1] - centers[j][1];
								double  qz = q[2] - centers[j][2];

								cstr2[0] = -1 * 1.0;
								cstr2[1] = -1 * qx;
								cstr2[2] = -1 * qy;
								cstr2[3] = -1 * qz;

								for (int g = 0; g < 4; g++)
									for (int h = g; h < 4; h++)
										sum_prod(g + j_start,h + j_start) += cstr2[g] * cstr2[h];

								for (int g = 0; g < 4; g++)
									for (int h = 0; h < 4; h++)
										sum_prod(g + i_start,h + j_start) += cstr[g] * cstr2[h];
							}
						}
			}

			//Fill in the lhs below the main diagonal
			for (int g = 1; g < dim; g++)
				for (int h = 0; h < g; h++)
					sum_prod(g,h) = sum_prod(h,g);

			// 3. Find the scales and scale the matrices appropriately to
			// normalize them and increase the numerical stability. However, it
			// might not be necessary since our transformation is only linear
			std::vector< double > scales(dim,0.0);
			double factor0 = 0, factor1 = 0, scale;
			for (unsigned int i = 0; i< size; i++) {
				int start_row = i*4;
				if ( sum_prod(start_row,start_row) > factor0 )
					factor0 = sum_prod(start_row,start_row);
				double max1 = vnl_math_max( vnl_math_abs(sum_prod(start_row+1,start_row+1)), vnl_math_abs(sum_prod(start_row+2,start_row+2)) );
				double max2 = vnl_math_max( vnl_math_abs(sum_prod(start_row+3,start_row+3)),max1 );
				if ( max2 > factor1 ) factor1 = max2;
			}
			scale = sqrt(factor1 / factor0);   // neither should be 0
			for (unsigned int i = 0; i<size; i++) {    
				int start_row = i*4;
				scales[start_row] = scale; 
				scales[start_row+1] = scales[start_row+2] = scales[start_row+3] = 1;
			}

			for ( int i=0; i<dim; i++ ) {
				sum_rhs(i,0) *= scales[i];
				sum_rhs(i,1) *= scales[i];
				sum_rhs(i,2) *= scales[i];
				for ( int j=0; j<dim; j++ )
					sum_prod(i,j) *= scales[i] * scales[j];
			}

			// 4. Estimate the transformation
			//
			vnl_svd<double> svd(sum_prod,1e-10);
			if ( svd.rank() < (unsigned int)dim ) {
				std::cerr<< "fregl_joint_register< TPixel >::estimate --- ERROR \n"
					<< "    sum_prod not full rank"<<std::endl;
				return false;
			}

			// COV(beta) = inv( A' * A )
			vnl_matrix<double> all_cov = svd.inverse();
			vnl_matrix<double> thetas = all_cov * sum_rhs;  // reuse inverse matrix

			// 5. Eliminate the scale
			//
			for ( int i=0; i<dim; i++ ) {
				thetas(i,0) *= scales[i];
				thetas(i,1) *= scales[i];
				thetas(i,2) *= scales[i];
			}

			// 6. Uncenter the xforms and covariance matrix 
			//   
			for (unsigned int i = 0; i<size; i++) {
				//vnl_matrix<double> c_params = thetas.extract(4,3,i*4,0);
				//vnl_matrix<double> theta(4,3,0.0);
				//vnl_matrix<double> A(4,4,0.0);
				//A.fill_diagonal( 1.0 );

				vnl_matrix<double> theta= thetas.extract(4,3,i*4,0);
				vnl_matrix<double> A = theta.extract(3,3,1,0).transpose();
				vnl_vector<double> t = theta.get_row(0);
				vnl_vector<double> uncentered_trans = t-A*centers[i]+anchor_center;
				transforms_[indices_for_other_images[i]][anchor] = TransformType::New();
				TransformType::ParametersType params(12);
				params[0] = A(0,0);
				params[1] = A(0,1);
				params[2] = A(0,2);
				params[3] = A(1,0);
				params[4] = A(1,1);
				params[5] = A(1,2);
				params[6] = A(2,0);
				params[7] = A(2,1);
				params[8] = A(2,2);
				params[9] = uncentered_trans[0];
				params[10] = uncentered_trans[1];
				params[11] = uncentered_trans[2];
				transforms_[indices_for_other_images[i]][anchor]->SetParameters(params);
			}

			return true;

}

template < class TPixel >
void 
fregl_joint_register< TPixel >::
generate_correspondences()
{
	corresp_generated_ = true;
	pairwise_constraints_.resize(image_ids_.size(), image_ids_.size());

	int space = 20;
	itk::Point<double, 3> point, xformed_pt;
	correspondence corresp;
	SizeType size_from, size_to; 
	for (unsigned int from = 0; from<image_ids_.size(); from++) {
		for (unsigned int to = from+1; to<image_ids_.size(); to++) {
			// generate the correspondences if there is a valid pairwise
			// xform between two images. 
			if ( transforms_(from,to) ) {
				size_from = image_sizes_[from];
				int z_space = vnl_math_min(10, int(size_from[2]/3));
				size_to = image_sizes_[to];
				for (unsigned int ix = 0; ix<size_from[0]; ix += space )
					for (unsigned int iy = 0; iy<size_from[1]; iy += space )
						for (unsigned int iz = 0; iz<size_from[2]; iz += z_space ) {
							point[0] = ix;
							point[1] = iy;
							point[2] = iz;

							xformed_pt = transforms_[from][to]->TransformPoint(point);
							//Check if the point is in range
							if ( xformed_pt[0]>0 && xformed_pt[0] < size_to[0] &&
								xformed_pt[1]>0 && xformed_pt[1] < size_to[1] && 
								xformed_pt[2]>0 && xformed_pt[2] < size_to[2]) {
									//std::cout<<"iz = "<<iz<<std::endl;
									corresp.from_ = point;
									corresp.to_ = xformed_pt;
									pairwise_constraints_[from][to].push_back(corresp);
									corresp.from_ = xformed_pt;
									corresp.to_ = point;
									pairwise_constraints_[to][from].push_back(corresp);
							}
						}
						std::cout<<image_ids_[from]<<" shares with "<<image_ids_[to]
						<<" "<<pairwise_constraints_[from][to].size()
							<<" correspondences."<<std::endl;

			} //if overlap
		}
	}
}

template < class TPixel >
fregl_reg_record::Pointer 
fregl_joint_register< TPixel >::
get_reg_record(int from, int to) const
{
	fregl_reg_record::Pointer reg_rec = new fregl_reg_record();
	reg_rec->set_from_image_info( image_ids_[from], image_sizes_[from] );
	reg_rec->set_to_image_info( image_ids_[to], image_sizes_[to]);
	reg_rec->set_transform( transforms_[from][to] );
	reg_rec->set_overlap( overlap_(from,to) );
	reg_rec->set_obj_value( obj_(from,to) );
	return reg_rec;
}

//#if defined(LIBXML_TREE_ENABLED) && defined(LIBXML_OUTPUT_ENABLED)
template < class TPixel >
void
fregl_joint_register< TPixel >::
write_xml(std::string const& filename, int sub_graphs_built, bool gen_temp_stuff)
{
	TiXmlDocument doc;       /* document pointer */
	std::string str;

	/* 
	* Creates a new document, a node and set it as a root node
	*/

	TiXmlElement* root_node = new TiXmlElement("Joint_Registration");
	str = ToString( transforms_.rows() );
	root_node->SetAttribute("number_of_images", str.c_str());

	str = ToString( sub_graphs_built ); 
	root_node->SetAttribute("sub_graphs_built", str);

	str = ToString( scale_multiplier_ );
	root_node->SetAttribute("error_scale_multiplier",str.c_str() );

	str = ToString( error_bound_ );
	root_node->SetAttribute("error_bound",str.c_str() );
	doc.LinkEndChild( root_node );

	// dump out the sub-graph information
	for (int index = 1; index <= sub_graphs_built; index++) {
		TiXmlElement* graph_node = new TiXmlElement("subgraph");
		for (unsigned int i = 0; i<image_ids_.size(); i++) {
			if (graph_indices_[i] == index) {
				TiXmlElement* node = new TiXmlElement("image_ID");
				node->LinkEndChild(new TiXmlText(image_ids_[i]));
				graph_node->LinkEndChild(node);
			}
		}
		root_node->LinkEndChild(graph_node);
	}

	// Set up the temp file for storing the overlap & errors
	std::ofstream temp_good, temp_bad, temp_overlap, temp_debug;
	if (gen_temp_stuff) {

		std::string filename_good = filename+"_debug_good.txt";
		std::string filename_bad = filename+"_debug_bad.txt";
		std::string filename_overlap = filename+"_debug_overlapped.txt";
		std::string filename_debug = filename+"_debug.txt";
		temp_good.open(filename_good.c_str());
		temp_good<<"from_image"<<"\t"<<"to_image"<<"\t"<<"overlap"<<"\t"<<"objective_value" <<"\n";
		temp_bad.open(filename_bad.c_str());
		temp_bad<<"from_image"<<"\t"<<"to_image"<<"\t"<<"overlap"<<"\t"<<"objective_value"<<"\n";
		temp_overlap.open(filename_overlap.c_str());
		temp_overlap<<"from_image"<<"\t"<<"to_image"<<"\t"<<"overlap"<<"\t"<<"objective_value"<<"\n";
		temp_debug.open(filename_debug.c_str());
		temp_debug<<"from_image"<<"\t"<<"to_image"<<"\t"<<"overlap"<<"\t"<<"objective_value"<<"\n";
	}

	for (unsigned int j = 0; j<transforms_.cols(); j++) {
		for (unsigned int i = 0; i<transforms_.rows(); i++) {
			if (i == j) 
				continue;

			fregl_reg_record::Pointer reg_rec = this->get_reg_record(i,j);

			/*TransformType::ParametersType params = this->transforms_(i,j)->GetParameters();

			std::cout << params[11] << std::endl;*/

			reg_rec->write_xml_node(root_node);

			// dump to the output files for debugging
			if (gen_temp_stuff && i>j) {
				if (reg_rec->obj()< 1) { //pairwise suceeded
					if (reg_rec->obj() < error_bound_) //pairwise accepted
						temp_good<<reg_rec->from_image()<<"\t"<<reg_rec->to_image()<<"\t"<<vnl_math_rnd(reg_rec->overlap()*1000)/1000.0<<"\t"<<vnl_math_rnd(reg_rec->obj()*1000)/1000.0 <<"\n";
					else  //pairwise rejected by joint
						temp_bad<<reg_rec->from_image()<<"\t"<<reg_rec->to_image()<<"\t"<<vnl_math_rnd(reg_rec->overlap()*1000)/1000.0<<"\t"<<vnl_math_rnd(reg_rec->obj()*1000)/1000.0 <<"\n";

					temp_debug<<reg_rec->from_image()<<"\t"<<reg_rec->to_image()<<"\t"<<vnl_math_rnd(reg_rec->overlap()*1000)/1000.0<<"\t"<<vnl_math_rnd(reg_rec->obj()*1000)/1000.0 <<"\n";
				}
				else if (reg_rec->overlap() > 0.01){ //pairwise failed, but images overlap
					temp_overlap<<reg_rec->from_image()<<"\t"<<reg_rec->to_image()<<"\t"<<vnl_math_rnd(reg_rec->overlap()*1000)/1000.0<<"\t"<<vnl_math_rnd(reg_rec->obj()*1000)/1000.0 <<"\n";
					temp_debug<<reg_rec->from_image()<<"\t"<<reg_rec->to_image()<<"\t"<<vnl_math_rnd(reg_rec->overlap()*1000)/1000.0<<"\t"<<vnl_math_rnd(reg_rec->obj()*1000)/1000.0 <<"\n";
				}
			}
		}
	}
	if (gen_temp_stuff) {
		temp_good.close();
		temp_bad.close();
		temp_overlap.close();
	}
	// Do the adjacency output if debug set
	if (gen_temp_stuff) {
		std::string filename_adjacency = filename+"_debug_adjacency.txt";
		std::ofstream temp_adjacency(filename_adjacency.c_str());
		for (unsigned int i = 0; i<image_ids_.size(); i++) {
			temp_adjacency<<image_ids_[i]<<":\t";
			std::vector<std::string> neighbors;
			get_adjacent_images(image_ids_[i],neighbors);
			for (unsigned int j = 0; j<neighbors.size(); j++) {
				temp_adjacency<<neighbors[j]<<"\t";
			}
			temp_adjacency<<"\n";
		}
		temp_adjacency.close();
	}

	/* 
	* Dumping document to a file
	*/
	doc.SaveFile(filename.c_str());

}

template < class TPixel >
void 
fregl_joint_register< TPixel >::
read_xml(std::string const & filename)
{
	TiXmlDocument doc;

	//Parse the resource 
	if ( !doc.LoadFile( filename.c_str() ) ) {
		vcl_cerr<<"Unable to load XML File"<<vcl_endl;
		return ;
	}

	/*Get the root element node */
	TiXmlElement* root_element = doc.RootElement();
	const char* docname = root_element->Value();
	if ( strcmp( docname, "Joint_Registration" ) != 0 &&
		strcmp( docname, "Pairwise_Registration" ) !=0) {
			vcl_cerr<<"Incorrect XML root Element: "<<vcl_endl;
			return;
	}

	/*
	if (strcmp(docname, "Joint_Registration") ==0) {
	// get the error bound
	scale_multiplier_ = atoi(root_element->Attribute("error_scale_multiplier"));
	// get the scale_multiplier
	std::stringstream(root_element->Attribute("error_bound")) >> error_bound_;
	}
	*/

	TiXmlElement* cur_node = root_element->FirstChildElement();

	// skip the information on subgraphs, found in joint registration only  
	if (strcmp(docname, "Joint_Registration") ==0) {  
		int num_graphs;
		num_graphs = atoi(root_element->Attribute("sub_graphs_built"));
		std::cout<<"number of subgraphs = "<<num_graphs<<std::endl;
		for (int i = 0; i<num_graphs; i++)
			cur_node = cur_node->NextSiblingElement();
	}

	// Reading the records
	std::vector<fregl_reg_record::Pointer> reg_records;
	for (; cur_node; cur_node = cur_node->NextSiblingElement() ) {
		fregl_reg_record::Pointer reg_rec = new fregl_reg_record();
		reg_rec->read_xml_node(cur_node);
		reg_records.push_back( reg_rec );
	}

	initialize( reg_records );
	//build_graph();

}

std::string 
ToString(double val)
{
	std::ostringstream strm;
	strm<< val<<std::endl;
	return strm.str();
}

//Explicit Template Instantiation
template class fregl_joint_register< unsigned char >;
template class fregl_joint_register< unsigned short >;
