//#include <libxml/encoding.h>
//#include <libxml/xmlwriter.h>

#include "parse_xml.h"

#include <iostream>
#include <sstream>
#include <string>

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xmlreader.h>

#define M(a,b,c) matrix[(a)*pimsize+(b)*pwidth+(c)]
#define L(a,b,c) lmatrix[(a)*pimsize+(b)*pwidth+(c)]


// following code was adapted from Charlene's code

params_surface parseXML_surface(char *filename, params_surface default_params)
{

	xmlDocPtr doc; /* the resulting document tree */
	xmlNodePtr root_element, cur_node = NULL;

	LIBXML_TEST_VERSION;

	doc = xmlReadFile(filename, NULL, 0);
	if (doc == NULL) {
		std::cerr<<"Failed to parse "<<filename<<std::endl;
		return default_params;
	}

	root_element = xmlDocGetRootElement(doc);
	char * contents;

	params_surface p=default_params;
	int num_parsed=0;
	cur_node =  root_element->children;
	for ( ; cur_node; cur_node = cur_node->next) {
		if (cur_node->type == XML_ELEMENT_NODE ) {
			num_parsed ++;
			if ( !xmlStrcmp(cur_node->name, BAD_CAST "default_input_filename") ) {
				contents = (char*)xmlNodeGetContent(cur_node);
				sscanf(contents,"%s",p.fname);
				xmlFree( contents );
			}
			else if ( !xmlStrcmp(cur_node->name, BAD_CAST "skip") ) {
				contents = (char*)xmlNodeGetContent(cur_node);
				std::stringstream( contents ) >> p.skip;
				xmlFree( contents );
				
			}
			else if ( !xmlStrcmp(cur_node->name, BAD_CAST "window_z") ) {
				contents = (char*)xmlNodeGetContent(cur_node);
				std::stringstream( contents ) >> p.window1;
				xmlFree( contents );
				
			}
			else if ( !xmlStrcmp(cur_node->name, BAD_CAST "window_xy") ) {
				contents = (char*)xmlNodeGetContent(cur_node);
				std::stringstream( contents ) >> p.window;
				xmlFree( contents );
				
			}
			else if ( !xmlStrcmp(cur_node->name, BAD_CAST "foreground_width") ) {
				contents = (char*)xmlNodeGetContent(cur_node);
				std::stringstream( contents ) >> p.epsilonw;
				xmlFree( contents );
				
			}
			else if ( !xmlStrcmp(cur_node->name, BAD_CAST "do_debug") ) {
				contents = (char*)xmlNodeGetContent(cur_node);
				std::stringstream( contents ) >> p.debug;
				xmlFree( contents );
				
			}
			else if ( !xmlStrcmp(cur_node->name, BAD_CAST "alpha") ) {
				contents = (char*)xmlNodeGetContent(cur_node);
				std::stringstream( contents ) >> p.alpha1;
				xmlFree( contents );
				
			}
			else {
				num_parsed--;
			}
		}
	}

	std::cout<<p.fname<<std::endl<<p.skip<<std::endl<<p.window1<<std::endl<<p.window<<std::endl<<p.epsilonw<<std::endl<<p.debug<<std::endl<<p.alpha1;

	return p;
}

struct point{
	int x,y,z;
	int d;
	double l;
};

bool parseXML_output(char *filename, int* &matrix, double * &lmatrix, params_surface &par)
{

	xmlDocPtr doc; /* the resulting document tree */
	xmlNodePtr root_element, cur_node = NULL, child_node = NULL;

	LIBXML_TEST_VERSION;

	doc = xmlReadFile(filename, NULL, 0);
	if (doc == NULL) {
		std::cerr<<"Failed to parse "<<filename<<std::endl;
		return 0;
	}

	root_element = xmlDocGetRootElement(doc);
	char * contents;
	char * attribute1;

	int num_parsed=0;
	point point_;
	int pimsize, pwidth;
	cur_node =  root_element->children;
	for ( ; cur_node; cur_node = cur_node->next) {
		if (cur_node->type == XML_ELEMENT_NODE ) {
			num_parsed ++;
			if ( !xmlStrcmp(cur_node->name, BAD_CAST "vessel") ) {
				attribute1 = (char*)xmlGetProp(cur_node,BAD_CAST "program");
				printf("Read program = \"%s\"\n",attribute1);
				xmlFree(attribute1);
				attribute1 = (char*)xmlGetProp(cur_node,BAD_CAST "version");
				printf("Read version = \"%s\"\n",attribute1);
				xmlFree(attribute1);
			}
			else if ( !xmlStrcmp(cur_node->name, BAD_CAST "parameters") ) {
				attribute1 = (char*)xmlGetProp(cur_node,BAD_CAST "sizex");
				std::stringstream( attribute1 ) >> par.sizex;
				xmlFree( attribute1 );
				attribute1 = (char*)xmlGetProp(cur_node,BAD_CAST "sizey");
				std::stringstream( attribute1 ) >> par.sizey;
				xmlFree( attribute1 );
				attribute1 = (char*)xmlGetProp(cur_node,BAD_CAST "sizez");
				std::stringstream( attribute1 ) >> par.sizez;
				xmlFree( attribute1 );
				attribute1 = (char*)xmlGetProp(cur_node,BAD_CAST "window_z");
				std::stringstream( attribute1 ) >> par.window1;
				xmlFree( attribute1 );
				attribute1 = (char*)xmlGetProp(cur_node,BAD_CAST "window_xy");
				std::stringstream( attribute1 ) >> par.window;
				xmlFree( attribute1 );
				attribute1 = (char*)xmlGetProp(cur_node,BAD_CAST "foreground_width");
				std::stringstream( attribute1 ) >> par.epsilonw;
				xmlFree( attribute1 );
				attribute1 = (char*)xmlGetProp(cur_node,BAD_CAST "do_debug");
				std::stringstream( attribute1 ) >> par.debug;
				xmlFree( attribute1 );
				attribute1 = (char*)xmlGetProp(cur_node,BAD_CAST "alpha");
				std::stringstream( attribute1 ) >> par.alpha1;
				xmlFree( attribute1 );
				// allocate memory for matrix, lmatrix
				pimsize = (par.sizex+2*par.window)*(par.sizey+2*par.window);
				pwidth = (par.sizex+2*par.window);
				matrix = (int*) malloc(pimsize*pwidth*sizeof(unsigned char));
				lmatrix = (double*) malloc(pimsize*pwidth*sizeof(double));
			}
			else {
				num_parsed--;
			}
		}
		else{
			if( !xmlStrcmp(cur_node->name, BAD_CAST "segmentation")) {
				child_node = cur_node->children;
				for ( ; child_node; child_node = child_node->next) {
					if( !xmlStrcmp(child_node->name, BAD_CAST "point")) {
						attribute1 = (char*)xmlGetProp(child_node,BAD_CAST "x");
						std::stringstream(attribute1) >> point_.x;
						xmlFree(attribute1);
						attribute1 = (char*)xmlGetProp(child_node,BAD_CAST "y");
						std::stringstream(attribute1) >> point_.y;
						xmlFree(attribute1);
						attribute1 = (char*)xmlGetProp(child_node,BAD_CAST "z");
						std::stringstream(attribute1) >> point_.z;
						xmlFree(attribute1);
						attribute1 = (char*)xmlGetProp(child_node,BAD_CAST "d");
						std::stringstream(attribute1) >> point_.d;
						xmlFree(attribute1);
						attribute1 = (char*)xmlGetProp(child_node,BAD_CAST "l");
						std::stringstream(attribute1) >> point_.l;
						xmlFree(attribute1);
						M(point_.z+par.window1-1,point_.y+par.window-1,point_.x+par.window-1)= point_.d;
						L(point_.z+par.window1-1,point_.y+par.window-1,point_.x+par.window-1)= point_.l;
					}
				}
			}
		}
	}
	return true;
}
