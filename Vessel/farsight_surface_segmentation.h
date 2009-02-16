#ifndef _FARSIGHT_SURFACE_SEGMENTATION_H
#define _FARSIGHT_SURFACE_SEGMENTATION_H

#include <stdio.h>
#include <vector>
//#include "tiffio.h"
#include <algorithm>
#include <math.h>
#include <string.h>
#include <string>
#include <sstream>
#include <assert.h>
#include <ctime>
#include "parse_xml.h"
#include <Common/farsight_object.h>

int compare(const void *a,const void *b)
{
	return (int)(*(unsigned char*)a)-(int)(*(unsigned char*)b);
}
int compare_double(const void *a,const void *b)
{
	return int(*(double*)a - *(double*)b);
}
struct data{
	int x,y,z;
};

struct return_data{
	int M;
	double lvalue;
};


class farsight_surface_segmentation: public farsight_object<3>
{
public:
	farsight_surface_segmentation()
	{
		window = 2;
		window1 = 2;
		lambda_lower_threshold = 3;
		lambda_higher_threshold = 10;
		prune = 0;
		epsilonw = 1.5;
		alpha1 = 0.2;
	}
private:
	// variable declarations
	// window - max window width for X,Y directions
	int window;
	// window1 - max window width for Z direction
	int window1;
	// lower threshold
	double lambda_lower_threshold;
	// higher threshold
	double lambda_higher_threshold;
	//prune value
	double prune;
	double epsilonw,alpha1;
	int psize,pimsize,pwidth ,pheight,pdepth ;
	double planes[42][3];
	unsigned char * null,*fore,*back;
	unsigned char *p;
	ImagePtrType imptr;
	std::string image_filename;
	return_data function(int coz, int coy, int cox);
	int *matrix;
	int skip;
	double *lmatrix;
	bool debug;
public:
	virtual void set_image(std::string const& image_path,std::string const& image_name);
	virtual void set_parameters(std::string const& xml_filename);
	virtual bool run();
	virtual void read_xml(std::string const& file_path,std::string const& filename);
	virtual ImageConstPtrType display();
	virtual void write_xml(std::string const& file_path,std::string const& xml_filename);

	~farsight_surface_segmentation()
	{
		free(matrix);
		free(lmatrix);
		free(p);
	}
};

#endif
