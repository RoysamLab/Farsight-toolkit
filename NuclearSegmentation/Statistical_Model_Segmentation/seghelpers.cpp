

#include "seghelpers.h"
#include "itkRegionOfInterestImageFilter.h"




template <typename T>
int writeImage(typename T::Pointer im, const char* filename)
{
	printf("Writing %s ... ",filename);
	typedef typename itk::ImageFileWriter<T> WriterType;
	
	typename WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(filename);
	writer->SetInput(im);
	try
	{
		writer->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cout << "ExceptionObject caught!" <<std::endl;
		std::cout << err << std::endl;
		return EXIT_FAILURE;
	}
	printf("Done.\n");
	return EXIT_SUCCESS;
}


template <typename T>
typename T::Pointer returnROI(typename T::Pointer image,typename T::RegionType region1)
{
	
	typedef itk::RegionOfInterestImageFilter< T,T> ROIFilterType;
	
	
	typename ROIFilterType::Pointer roifilter = ROIFilterType::New();
 	roifilter->SetInput(image);
	roifilter->SetRegionOfInterest(region1); 	
	
	try
	{
		roifilter->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cout << "ExceptionObject caught!" <<std::endl;
		std::cout << err << std::endl;
		//return EXIT_FAILURE;
	}
	printf("Done.\n");
	return roifilter->GetOutput();
	
}


bool sorthelp (double i,double j) 
{ 
	return (i<j); 
}


bool uniquehelp (unsigned short i, unsigned short j) {
	return (i==j);
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
//NORMALIZE A SET OF VECTORS BY THE STD DEVIATIONS OF EACH OF THE VARIABLES
///////////////////////////////////////////////////////////////////////////////////////////////////////
void norm_vec_std(double **mat1, long m, long n1, double **mat2, long n2) 
{
	double *ave;
	double *std;
	double anom;		//intermediate output
	
	//allocate space:
	ave=new double[m];
	std=new double[m];
	
	for (long i=0; i<m; i++) {
		ave[i]=0;
		for (long j=0; j<n1; j++) ave[i]+=mat1[j][i];
		ave[i]/=n1;
		std[i]=0;
		for (long j=0; j<n1; j++) {
			anom=mat1[j][i]-ave[i];
			std[i]+=anom*anom;
			mat1[j][i]=anom;
		}
		std[i]=sqrt(std[i]/(n1-1));
		
		//normalize the first matrix:
		for (long j=0; j<n1; j++) mat1[j][i]=mat1[j][i]/std[i];
		
		//normalize the second matrix:
		for (long j=0; j<n2; j++) mat2[j][i]=(mat2[j][i]-ave[i])/std[i];
	}
	//clean up:
	delete [] ave;
	delete [] std;
}






///////////////////////////////////////////////////////////////////////////////////////////////////////
// GET THE NUMBER OF ROWS AND COLUMNS IN THE TRAINING FILE
///////////////////////////////////////////////////////////////////////////////////////////////////////

int Num_Lines2(char *root, int mode) {
	int   row_no, column_no, zctr1 = 0, zctr2 = 0, zctr3 = 0;
	char  file_name[80], temp[80], ch;
	FILE  *fpinn;
	
	strcpy(file_name, root);
	if ((fpinn = fopen(file_name, "r")) == NULL) {
		printf("\n  The file %s does not exist !\n", file_name);
		printf("      Press <ENTER> to exit     \n\n");
		getchar();
		exit(1);
	}
	while (!feof(fpinn)) {
		fscanf(fpinn, "%c",&ch);
		zctr2++;
		//--- we don't want to count empty rows
		if (isspace(ch)) zctr2--;
		//
		if (ch=='\n' && zctr2 > 0){
			zctr3++;
			zctr2=0;
		}
		if (feof(fpinn) && zctr2 > 0){//--- this part required for linux
			zctr3++;
			zctr1++;
		}
	}
	rewind(fpinn);
	while (!feof(fpinn)){ //--- counts all data points
		fscanf(fpinn, "%s",&temp);
		zctr1++;
	}
	fclose(fpinn);
	zctr1     = zctr1-1;
	row_no    = zctr3;
	column_no = zctr1/row_no;
	//printf("\n\ncolumn no = %d  row_no= %d\n\n", column_no,   row_no);
	if (mode == 0) return row_no;
	if (mode == 1) return column_no;
	return -999;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// FILE READING/WRITING 
///////////////////////////////////////////////////////////////////////////////////////////////////////

FILE *FDeclare2(char *root, char *extension, char key) {
	char string1[80];
	FILE *fpot;
	
	strcpy(string1, root);
	if (strcmp(extension,"")!=0) strcat(string1, ".");
	strcat(string1, extension);
	if ((key != 'w')&&(key!='r')) {
		printf("either read or write, wrong key %c \n", key);getchar();
		exit(1);
	}
	if (key == 'r') {
		//if (( (FILE *) fpot = fopen(string1, "r")) == NULL) {
		//changed for DOS?UNIX compatibility
		if ((fpot = fopen(string1, "r")) == NULL) {
			printf("\n");
			printf("input file %s does not exist\n", string1);
			exit(1); getchar();
		}
	}
	if (key == 'w') {
		if ((fpot = fopen(string1, "w")) == NULL) {
			printf("\n");
			printf("output file %s does not exist\n", string1);
			exit(1);getchar();
		}
	}
	return fpot;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//CALCULATE THE AVERAGES AND STANDARD DEVIATIONS OF A SET OF VECTORS
///////////////////////////////////////////////////////////////////////////////////////////////////////

void calc_norm(double **mat, long D,long n, double *ave,double *stdev)
{
	double anom;		//intermediate output
	
	for (long i=0; i<D; i++) {
		stdev[i]=0;
		ave[i]=0;
		for (long j=0; j<n; j++) 
		{
			ave[i]+=mat[j][i];
			//std::cout<<mat[j][i] <<std::endl;
		}
		ave[i]/=n;		
		for (long j=0; j<n; j++) {
			anom=mat[j][i]-ave[i];	
			stdev[i]+=anom*anom;
		}
		stdev[i]=sqrt(stdev[i]/(n-1));
		
	}
	
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// CONVERT THE INPUT DATA TYPE TO STRING 
// NEED TO TEMPLATE
///////////////////////////////////////////////////////////////////////////////////////////////////////
std::string convert2string(unsigned long id)
{
	std::stringstream out;
	out << id;
	std::string s = out.str();
	return s;
}