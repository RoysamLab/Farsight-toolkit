
#include "mex.h"
#include "matrix.h"
#include <math.h>
#include "fftw3.h"


#include <stdio.h>
#include <string>
#include <iostream>
#include <time.h>
#include <winsock.h>
#include <iomanip> 


#define WISDOMNAME_SINGLE "myfftwWisdom"
#define WISDOMNAME_DOUBLE "myfftwWisdom_double"
#define _MYFFTW_VERBOSE

#define REAL float


#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif

struct timezone 
{
  int  tz_minuteswest; /* minutes W of Greenwich */
  int  tz_dsttime;     /* type of dst correction */
};

int gettimeofday(struct timeval *tv, struct timezone *tz)
{
  FILETIME ft;
  unsigned __int64 tmpres = 0;
  static int tzflag = 0;

  if (NULL != tv)
  {
    GetSystemTimeAsFileTime(&ft);

    tmpres |= ft.dwHighDateTime;
    tmpres <<= 32;
    tmpres |= ft.dwLowDateTime;

    tmpres /= 10;  /*convert into microseconds*/
    /*converting file time to unix epoch*/
    tmpres -= DELTA_EPOCH_IN_MICROSECS; 
    tv->tv_sec = (long)(tmpres / 1000000UL);
    tv->tv_usec = (long)(tmpres % 1000000UL);
  }

  if (NULL != tz)
  {
    if (!tzflag)
    {
      _tzset();
      tzflag++;
    }
    tz->tz_minuteswest = _timezone / 60;
    tz->tz_dsttime = _daylight;
  }

  return 0;
}


inline double
gettime(void)
{
  static struct timeval timeS;
  gettimeofday( &timeS, NULL);
  return (timeS.tv_sec + timeS.tv_usec/1000000.0);
}




inline void report_timing( std::string message)
{
  static double startTime = -1.0;
  static double lastTime = -1.0;
  if( startTime < 0.0) {
    startTime = gettime();
    lastTime = startTime;
  }
  double currentTime = gettime();

  std::cerr.setf(std::ios::fixed, std::ios::floatfield);


  std::cerr << std::setw(7)
          << std::setprecision(3) << currentTime - startTime << "s ("
          << std::setw(6) << std::setprecision(4)
          << currentTime - lastTime << "s)  "
          << message << std::endl;
  lastTime = currentTime;
}



void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
 	
    if(nrhs != 2 && nrhs != 3) {
	printf("\nUsage:\n");
    return;
	} else if(nlhs>1) {
	printf("Too many output arguments\n");
    return;
	}


    ////////////// fetching data 

    int pcnt = 0;
    const mxArray *Img;
    Img = prhs[pcnt++];       
    const int numdim = mxGetNumberOfDimensions(Img);
    const int *dims_sz = mxGetDimensions(Img);
    int *dims_sz_ = (int *)mxCalloc(numdim,sizeof(int));//[numdim];
    memcpy(dims_sz_,dims_sz,sizeof(int)*numdim);

    int totsz = 1;
    for (int k = 0; k < numdim;k++)
	totsz = dims_sz[k] * totsz;
    totsz /= 2;


    /////////////// reading dimensions along fft computed

    const mxArray *_along;
    _along = prhs[pcnt++];       

    double *along = mxGetPr(_along);
    int rank = mxGetN(_along);

    int *fftdim = (int *)mxCalloc(numdim, sizeof(int));//[numdim];
    for (int i = 0; i < numdim;i++)
	fftdim[i] = 0;
    for (int i = 0; i < rank; i++)
	fftdim[(int) along[i]-1] = 1;

    fftwf_iodim *dims = new fftwf_iodim[rank];
    fftwf_iodim *howmany_dims = new fftwf_iodim[numdim-rank];



    int howmany_rank = 0;
    rank = 0;
    int stride = 1;
    dims_sz_[0] /= 2;
    for (int i = 0; i < numdim; i++)
    {
	if (fftdim[i] == 1)
	{
		dims[rank].n = dims_sz_[i];
		dims[rank].is = stride;
		dims[rank].os = dims[rank].is;
		rank++;
	}
	else
	{
		howmany_dims[howmany_rank].n = dims_sz_[i];
		howmany_dims[howmany_rank].is = stride;
		howmany_dims[howmany_rank].os = howmany_dims[howmany_rank].is;
		howmany_rank++;
	}
	stride *= dims_sz_[i];
	
    }








    /////////////// reading planner flags

    unsigned flags =FFTW_ESTIMATE;

    if (nrhs == 3)
	{
		int planflag = (int) *mxGetPr(prhs[pcnt++]);
		switch(planflag)
		{
			case 0: flags = FFTW_ESTIMATE;   break;
			case 1: flags = FFTW_MEASURE;    break;
			case 2: flags = FFTW_PATIENT;    break;
			case 3: flags = FFTW_EXHAUSTIVE; break;			
		}
	}


   /////////////////// single support

   if (mxGetClassID(Img) == mxSINGLE_CLASS)
   {
	fftwf_forget_wisdom();
	mxArray *wisdom = mexGetVariable("global",WISDOMNAME_SINGLE);
	if (wisdom != NULL)
	{
        printf("I got the wisdom now\n");
		int buflen = mxGetN(wisdom)*sizeof(mxChar)+1;
		char *fftwwisdom = (char*)mxCalloc(buflen,sizeof(char));//[buflen];
		mxGetString(wisdom,fftwwisdom,buflen);
		int ret = fftwf_import_wisdom_from_string(fftwwisdom);
        if(ret)
            printf("Wisdom imported successfully\n");
	}


	fftwf_complex *ri = (fftwf_complex*) mxGetData(Img);
	plhs[0] = mxCreateNumericArray(numdim,dims_sz,mxGetClassID(Img),mxREAL);
	fftwf_complex *ro = (fftwf_complex*) mxGetData(plhs[0]);	
#ifdef _MYFFTW_VERBOSE
        report_timing("\n init, float ");	
#endif
	fftwf_plan plan = fftwf_plan_guru_dft(rank,(const fftwf_iodim*) dims, howmany_rank, (const fftwf_iodim*)howmany_dims,ro,ro, 1, flags);

	memcpy(ro,ri,sizeof(fftwf_complex)*totsz);
	

#ifdef _MYFFTW_VERBOSE
	report_timing("planning, float");
#endif	
	fftwf_execute_dft(plan,ro,ro);
#ifdef _MYFFTW_VERBOSE
	report_timing("execution, float ");
#endif
	fftwf_destroy_plan(plan);

        mxArray *newwisdom = mxCreateString(fftwf_export_wisdom_to_string());
        mexPutVariable("global",WISDOMNAME_SINGLE,newwisdom);

   }



   /////////////////// double support

   if (mxGetClassID(Img) == mxDOUBLE_CLASS)
   {
	fftw_forget_wisdom();
	mxArray *wisdom = mexGetVariable("global",WISDOMNAME_DOUBLE);
	if (wisdom != NULL)
	{
		int buflen = mxGetN(wisdom)*sizeof(mxChar)+1;
		char *fftwwisdom = (char*)mxCalloc(buflen, sizeof(char));//[buflen];
		mxGetString(wisdom,fftwwisdom,buflen);
		int ret = fftw_import_wisdom_from_string(fftwwisdom);
	}


	fftw_complex *ri = (fftw_complex*) mxGetData(Img);
	plhs[0] = mxCreateNumericArray(numdim,dims_sz,mxGetClassID(Img),mxREAL);
	fftw_complex *ro = (fftw_complex*) mxGetData(plhs[0]);	
#ifdef _MYFFTW_VERBOSE
        report_timing("\n init, double ");	
#endif
	const fftw_plan plan = fftw_plan_guru_dft(rank,(const fftwf_iodim*) dims, howmany_rank, (const fftwf_iodim*)howmany_dims,ro,ro, 1, flags| FFTW_PRESERVE_INPUT);

	memcpy(ro,ri,sizeof(fftw_complex)*totsz);

#ifdef _MYFFTW_VERBOSE
	report_timing("planning, double");
#endif
	fftw_execute_dft(plan,ro,ro);
#ifdef _MYFFTW_VERBOSE
	report_timing("execution, double ");
#endif
	fftw_destroy_plan(plan);

        mxArray *newwisdom = mxCreateString(fftw_export_wisdom_to_string());
        mexPutVariable("global",WISDOMNAME_DOUBLE,newwisdom);

   }


}
