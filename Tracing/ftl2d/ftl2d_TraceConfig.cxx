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

//tracing_config.cpp
#include "ftl2d_TraceConfig.h"

#include <string>




#define MY_ENCODING "ISO-8859-1"

TraceConfig::TraceConfig()  {
     GridSpacing = 5;
     StepRatio = 0.8;
     AspectRatio = 2.5;
     SeedIntensityThreshold = 100;
     MinimumVesselWidth = 1.0;
     MaximumVesselWidth = 50.0;
     MinimumVesselLength = 10.0;
     StepSize = 0.5;
//     ParameterFileName = "Tracing_parameters.input";
}

//TraceConfig::~TraceConfig() { }

void TraceConfig::LoadParameters(char * ParameterFileName)  {
	
	char *contents;
	xmlChar* attribute;
	xmlDocPtr doc; /* the resulting document tree */
	xmlNodePtr root_element, main_node = NULL, parm_node = NULL;
	std::string whitespaces (" \t\f\v\n\r");
  
	LIBXML_TEST_VERSION;
	     
	//Parse the resource 
	doc = xmlReadFile(ParameterFileName, NULL, 0);
	if (doc == NULL) {
		std::cerr<<"Failed to parse "<<ParameterFileName<<std::endl;
	    return;
	}
	/*Get the root element node */
	main_node = xmlDocGetRootElement(doc);
		
	//main_node =  root_element->children;
	std::cout<<"Parsing input XML tag : " <<main_node->name <<std::endl;
	if (main_node->type == XML_ELEMENT_NODE && !xmlStrcmp(main_node->name, BAD_CAST "Tracing_Paramenter_2D") ){
		parm_node = main_node->children;
		
			
		for ( ; parm_node; parm_node = parm_node->next) {	
			
			//std::cout<<"Parameter name - " <<parm_node->name <<std::endl;
			if (parm_node->type != XML_ELEMENT_NODE) continue;
			 
			if ( !xmlStrcmp(parm_node->name, BAD_CAST "InputFileName") ) {

			      InputFileName = (char*)xmlNodeGetContent(parm_node);
			      //erasing whitespaces -from begining and end only
			      InputFileName = InputFileName.substr(InputFileName.find_first_not_of(whitespaces), 
			      InputFileName.find_last_not_of(whitespaces) - InputFileName.find_first_not_of(whitespaces)+1);
			      
			      std::cout<<"InputFileName = "<<InputFileName << " Size = "<<InputFileName.size() <<std::endl;
			      continue;
			   }
			
			if ( !xmlStrcmp(parm_node->name, BAD_CAST "OutputFileName") ) {
			
				  OutputFileName = (char*)xmlNodeGetContent(parm_node);
				  //erasing whitespaces -from begining and end only
				  OutputFileName = OutputFileName.substr(OutputFileName.find_first_not_of(whitespaces), 
			      OutputFileName.find_last_not_of(whitespaces) - OutputFileName.find_first_not_of(whitespaces)+1);
			      
				  std::cout<<"OutputFileName = "<<OutputFileName<<" Size = "<<OutputFileName.size() <<std::endl;
				  continue;
			   }

			if ( !xmlStrcmp(parm_node->name, BAD_CAST "MinIsolatedVesselLength") ) {
				  MinimumVesselLength = atoi((char*)xmlNodeGetContent(parm_node));
				  std::cout<<"MinimumIsolatedVesselLength = "<<MinimumVesselLength<<std::endl;
				  continue;
			
			   }
			   
			if ( !xmlStrcmp(parm_node->name, BAD_CAST "GridSpacing") ) {

				  GridSpacing = atoi((char*)xmlNodeGetContent(parm_node));
				  std::cout<<"GridSpacing = "<<GridSpacing<<std::endl;
				  continue;
			   }
			   
			if ( !xmlStrcmp(parm_node->name, BAD_CAST "StepRatio") ) {
			
				  StepRatio = atof((char*)xmlNodeGetContent(parm_node));
				  std::cout<<"StepRatio = "<<StepRatio<<std::endl;
				  continue;
			  }
			  
			if ( !xmlStrcmp(parm_node->name, BAD_CAST "AspectRatio") ) {

				  AspectRatio = atof((char*)xmlNodeGetContent(parm_node));
				  std::cout<<"AspectRatio = "<<AspectRatio<<std::endl;
				  continue;
			  }

			if ( !xmlStrcmp(parm_node->name, BAD_CAST "Sensitivity") ) {

				  PROP = atof((char*)xmlNodeGetContent(parm_node));
				  std::cout<<"Sensitivity = "<<PROP<<std::endl;			//check range
				  continue;
			  }			  
		}
			
	}
	xmlCleanupParser();
  	xmlMemoryDump();

}


