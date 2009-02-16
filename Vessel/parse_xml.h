#ifndef _PARSE_XML_H
#define _PARSE_XML_H
struct params_surface{
	char fname[2048]; // bad programming, I know. likely to segfault when we handle terabytes of data with kilobytes of filanames
	int skip;
	int window1;
	int window;
	double alpha1;
	bool debug;
	double epsilonw;
	int sizex,sizey,sizez;
};


params_surface parseXML_surface(char *,params_surface);
bool parseXML_output(char *filename, int* &matrix, double * &lmatrix, params_surface &par);
#endif
