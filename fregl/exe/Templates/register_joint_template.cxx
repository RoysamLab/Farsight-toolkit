#include "register_joint_template.h"

template <typename TPixel>
int								//return type
register_joint_template(
	vul_arg< vcl_string > arg_in_file,
	vul_arg< vcl_string > arg_xml_file,
	vul_arg< double > arg_multiplier,
	vul_arg< double > arg_error_bound,
	vul_arg< vcl_string > arg_roi_file,
	vul_arg< bool > arg_debug
)
{
	typedef TPixel PixelType;
	
	//Read in the filename
	std::ifstream in_file_str( arg_in_file().c_str() );
	std::vector<std::string> filenames;

	if ( !in_file_str ){
		std::cerr<<"Couldn't open "<<arg_in_file()<<std::endl;
		exit( 0 );
	}

	std::string filename;
	while (in_file_str) {
		vcl_getline(in_file_str, filename);
		if (filename.length() == 0) continue;

		filenames.push_back( filename );
		std::cout<<"Read "<<filename<<std::endl;
	}
	in_file_str.close();

	typename fregl_joint_register< PixelType >::Pointer jointer_register;
	if ( !arg_roi_file.set() ) {
		jointer_register = new fregl_joint_register < PixelType >( filenames, arg_multiplier(), arg_error_bound() );
	}
	else 
	{
		// if roi set, read in the image names
		std::vector<std::string> roi_imagenames;
		std::ifstream in_file_str2( arg_roi_file().c_str() );

		if ( !in_file_str2 ){
			std::cerr<<"Couldn't open "<<arg_roi_file()<<std::endl;
			exit( 0 );
		}

		while (in_file_str2) {
			vcl_getline(in_file_str2, filename);
			if (filename.length() == 0) continue;

			roi_imagenames.push_back( filename );
			std::cout<<"Image "<<filename<<" in ROI"<<std::endl;
		}
		in_file_str2.close(); 

		std::vector< fregl_reg_record::Pointer > reg_records;
		for (unsigned int i = 0; i<filenames.size(); i++) {
			fregl_reg_record::Pointer reg_record = new fregl_reg_record();
			reg_record->read_xml(filenames[i]);
			bool to_found = false;
			bool from_found = false;
			for (unsigned int j = 0; j<roi_imagenames.size(); j++) {
				if (roi_imagenames[j] == reg_record->to_image()) to_found = true;
				if (roi_imagenames[j] == reg_record->from_image()) from_found = true;
			}
			if (to_found && from_found) reg_records.push_back( reg_record );
		}

		jointer_register = new fregl_joint_register < PixelType >( reg_records, arg_multiplier(), arg_error_bound() );

	}

	/*
	bool mutual_consistency = !arg_no_mc();
	if ( arg_no_mc() ) 
	std::cout<<"Joint registration without mutual consistency ..."<<std::endl;
	else 
	std::cout<<"Joint registration mutual consistency ..."<<std::endl;
	*/
	
	jointer_register->build_graph();
        std::cout<<"subgraphs built = "<<jointer_register->number_of_subgraphs()<<std::endl;
        
	/*for (unsigned int j = 0; j<jointer_register->transforms_.cols(); j++) 
	{
		for (unsigned int i = 0; i<jointer_register->transforms_.rows(); i++) 
		{
			if (jointer_register->transforms_(i,j))
			{
				TransformType::ParametersType params = jointer_register->transforms_(i,j)->GetParameters(); //temporary
				std::cout << i << " to " << j << ": " << params[9] << " " << params[10] << " " << params[11] << std::endl;
			}
		}
	}*/
	
	jointer_register->write_xml(arg_xml_file(), jointer_register->number_of_subgraphs(), arg_debug());

	/*
	bool graph_build = jointer_register->build_graph(!arg_no_mc());
	if ( !graph_build && !arg_no_mc() ) {
	std::cout<<"Failed to build the graph with mutual consistency. The graph is  now built without mutual consistency."<<vcl_endl;
	jointer_register->build_graph(arg_no_mc());
	mutual_consistency = false;
	}
	jointer_register->write_xml(arg_xml_file(), mutual_consistency, arg_debug());
	*/

	return 0;
}

template int
	register_joint_template<unsigned char>(
	vul_arg< vcl_string > arg_in_file,
	vul_arg< vcl_string > arg_xml_file,
	vul_arg< double > arg_multiplier,
	vul_arg< double > arg_error_bound,
	vul_arg< vcl_string > arg_roi_file,
	vul_arg< bool > arg_debug
);

template int
	register_joint_template<unsigned short>(
	vul_arg< vcl_string > arg_in_file,
	vul_arg< vcl_string > arg_xml_file,
	vul_arg< double > arg_multiplier,
	vul_arg< double > arg_error_bound,
	vul_arg< vcl_string > arg_roi_file,
	vul_arg< bool > arg_debug
);