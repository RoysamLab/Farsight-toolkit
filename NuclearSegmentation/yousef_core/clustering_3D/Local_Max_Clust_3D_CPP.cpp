/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

//
#include "local_max_clust_3D.h"
#include<math.h>
#include <time.h>


#ifdef _OPENMP
#include "omp.h"
#endif

//OpenCL support
#ifdef OPENCL
	#include "CL\cl.h"
	#include "CL\cl_ext.h"
	#ifndef MSTRINGIFY
		#define MSTRINGIFY(A) #A
		char* InitialClusteringKernel =
		#include "InitialClusteringKernel.cl"
	#endif
#endif

void clustering_pfn_notify(const char *errinfo, const void *private_info, size_t cb, void *user_data);
extern "C" void initialClustering_CUDA (float* im_vals, unsigned short* local_max_vals, unsigned short* max_response_r, unsigned short* max_response_c, unsigned short* max_response_z , int r, int c, int z, int scale_xy, int scale_z);


void get_maximum(float* A, int r1, int r2, int c1, int c2, int z1, int z2, int* rx, int* cx, int* zx, int R, int C, int Z)
{
	float mx = A[(z1*R*C)+(r1*C)+c1];//A[r1][c1][z1];
    rx[0] = r1;
    cx[0] = c1;
	zx[0] = z1;
    
	for(int i=r1; i<=r2; i++)
    {
        for(int j=c1; j<=c2; j++)
        {
			for(int k=z1; k<=z2; k++)
			{
				if(A[(k*R*C)+(i*C)+j]>=mx)
				{
					mx = A[(k*R*C)+(i*C)+j];//A[i][j][k];
					cx[0] = j;
					rx[0] = i;
					zx[0] = k;
				}
            }
        }
    }
}


void local_max_clust_3D(float* im_vals, unsigned short* local_max_vals, unsigned short* bImg, unsigned short* out1, int r, int c, int z, int scale_xy, int scale_z)
{  
	//im_vals is the Laplacian of Gaussian
	//local_max_vals is the seed points (local maximum) with foreground seeds assigned an id > 0 and background seeds id == -1
	// out1 will contain the clustering output

	int*** max_nghbr_im;
    
	//create max_nghbr_im and initialize it with its index (node) value
	max_nghbr_im = (int ***) malloc(r*sizeof(int**)); 
	
	#pragma omp parallel for
	for(int i=0; i<r; i++)
	{        
		max_nghbr_im[i] = (int **) malloc(c*sizeof(int*));
		for(int j=0; j<c; j++)
		{			
			max_nghbr_im[i][j] = (int *) malloc(z*sizeof(int));
			for(int k=0; k<z; k++)
			{				
				max_nghbr_im[i][j][k] = (k*r*c)+(i*c)+j;//LMX;
			}
		}
	}
	std::cout << "max_nghbr_im initialized" << std::endl;

	//In this loop we look in a local region around each point and find the maximum value in the LoG image
	//Set the value to the index of the local maximum, (so if I am a seed point do nothing).
	
	/*cerr << "sizeof(int ***) = " << sizeof(int ***) << endl;
	cerr << "sizeof(int **) = " << sizeof(int **) << endl;
	cerr << "sizeof(int *) = " << sizeof(int *) << endl;
	cerr << "sizeof(int) = " << sizeof(int) << endl;
	cerr << "Total size of array = " << (r * c * z * 3 * sizeof(int) + r * c * z * sizeof(int *) + r * c * sizeof(int **) + r * sizeof(int ***) + sizeof (int ****)) / (1024 * 1024) << " MB" << endl; //Probably incorrect */

	clock_t start_time_init_max_nghbr_im = clock();
	
	unsigned short *max_response_r;
	unsigned short *max_response_c;
	unsigned short *max_response_z;
	

	max_response_r = new unsigned short[r * c * z];
	max_response_c = new unsigned short[r * c * z];
	max_response_z = new unsigned short[r * c * z];
	
#ifdef OPENCL
	// START OPENCL BOILERPLATE ----------------------------------------------------------------------------------------------------------------
	cl_platform_id platforms[10];
	cl_uint num_platforms;
	cl_device_id device[10];
	cl_uint num_devices;
	cl_context context;
	cl_command_queue queue;
	cl_program program;
	cl_kernel kernel;
	cl_int errorcode;

	// Platform
	clGetPlatformIDs(10, platforms, &num_platforms); //Get OpenCL Platforms
	clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, 10, device, &num_devices); //Get the first platform and get a list of devices
	context = clCreateContext(0, 1, device, &clustering_pfn_notify, NULL, NULL); //Create context from first device
	queue = clCreateCommandQueue(context, device[0], 0, NULL); //Create a command queue for the first device
	
	//cout << endl << LocalMaximaKernel << endl << endl; //print out strinfified kernel
	
	program = clCreateProgramWithSource(context, 1, (const char **) &InitialClusteringKernel, 0, &errorcode); //Read in kernel and create a program
	
	if (errorcode != CL_SUCCESS)
		cout << "clCreateProgramWithSource Error code: " << errorcode << endl;
	
	errorcode = clBuildProgram(program, 0, 0, 0, 0, 0); //Build the program

	if (errorcode != CL_SUCCESS) //If there was a build error, print out build_info
	{
		cout << "clBuildProgram Error Code: " << errorcode << endl;
		char build_log[1024*1024];
		errorcode = clGetProgramBuildInfo(program, device[0], CL_PROGRAM_BUILD_LOG, sizeof(build_log), build_log, NULL); //Get the build log
		if (errorcode == CL_SUCCESS)
			cout << "Build Log:" << endl << build_log << endl;
		else
			cout << "clGetProgramBuildInfo Error Code: " << errorcode << endl;
	}

	kernel = clCreateKernel(program, "InitialClusteringKernel", &errorcode); //Create the kernel from the source code
	
	if (errorcode != CL_SUCCESS)
		cout << "clCreateKernel Error code: " << errorcode << endl;

	//END OPENCL BOILERPLATE ---------------------------------------------------------------------------------------------------------------------

	size_t cnDimension = r * c * z; //array size
	
	cl_ulong device_max_mem_alloc_size, device_global_mem_size;
	clGetDeviceInfo(device[0], CL_DEVICE_MAX_MEM_ALLOC_SIZE, 1024, &device_max_mem_alloc_size,	NULL);
	clGetDeviceInfo(device[0], CL_DEVICE_GLOBAL_MEM_SIZE, 1024, &device_global_mem_size,	NULL);
	cout << "Maximum memory allocation size: " << device_max_mem_alloc_size / (double)(1024*1024) << " MB" << endl;
	cout << "Maximum global memory allocation size: " << device_global_mem_size / (double)(1024*1024) << " MB" << endl;

	cout << "Allocating " << (sizeof(*im_vals) * cnDimension)/(double)(1024*1024) << " MB of memory on GPU for im_vals" << endl;
	cout << "Allocating " << (sizeof(*local_max_vals) * cnDimension)/(double)(1024*1024) << " MB of memory on GPU for local_max_vals" << endl;
	cout << "Allocating " << (sizeof(*max_response_r) * cnDimension)/(double)(1024*1024) << " MB of memory on GPU for max_response_r" << endl;
	cout << "Allocating " << (sizeof(*max_response_c) * cnDimension)/(double)(1024*1024) << " MB of memory on GPU for max_response_c" << endl;
	cout << "Allocating " << (sizeof(*max_response_z) * cnDimension)/(double)(1024*1024) << " MB of memory on GPU for max_response_z" << endl;
		
	//Allocate device memory
	cl_int errorcode1, errorcode2, errorcode3, errorcode4, errorcode5;	
	cl_mem device_mem_im_vals = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_float) * cnDimension, NULL, &errorcode1);
	cl_mem device_mem_local_max_vals = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_ushort) * cnDimension, NULL, &errorcode2);
	cl_mem device_mem_max_response_r = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(cl_ushort) * cnDimension, NULL, &errorcode3);
	cl_mem device_mem_max_response_c = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(cl_ushort) * cnDimension, NULL, &errorcode4);
	cl_mem device_mem_max_response_z = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(cl_ushort) * cnDimension, NULL, &errorcode5);
	
	if (errorcode1 || errorcode2 || errorcode3 || errorcode4 || errorcode5)
	{
		cout << "Failed to allocate buffer memory on GPU" << endl; 
	}	
	
	//Write memory from host to device
	clEnqueueWriteBuffer(queue, device_mem_im_vals, CL_TRUE, 0, sizeof(cl_float) * cnDimension, im_vals, NULL, NULL, NULL);
	clEnqueueWriteBuffer(queue, device_mem_local_max_vals, CL_TRUE, 0, sizeof(cl_ushort) * cnDimension, local_max_vals, NULL, NULL, NULL);

	//Set kernel arguments
	clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &device_mem_im_vals);
	clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &device_mem_local_max_vals);
	clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &device_mem_max_response_r);
	clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &device_mem_max_response_c);
	clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &device_mem_max_response_z);
	clSetKernelArg(kernel, 5, sizeof(int), &r);
	clSetKernelArg(kernel, 6, sizeof(int), &c);
	clSetKernelArg(kernel, 7, sizeof(int), &z);
	clSetKernelArg(kernel, 8, sizeof(int), &scale_xy);
	clSetKernelArg(kernel, 9, sizeof(int), &scale_z);	

	//Execute the kernel
	clEnqueueNDRangeKernel(queue, kernel, 1, 0, (const size_t *) &cnDimension, 0, 0, 0, 0);
	
	//Read the output from the device back into host memory
	clEnqueueReadBuffer(queue, device_mem_max_response_r, CL_TRUE, 0, sizeof(cl_ushort) * cnDimension, max_response_r, NULL, NULL, NULL);
	clEnqueueReadBuffer(queue, device_mem_max_response_c, CL_TRUE, 0, sizeof(cl_ushort) * cnDimension, max_response_c, NULL, NULL, NULL);
	clEnqueueReadBuffer(queue, device_mem_max_response_z, CL_TRUE, 0, sizeof(cl_ushort) * cnDimension, max_response_z, NULL, NULL, NULL);
	
	//Block till all commands are complete
	clFinish(queue);
	
	/*for (int i = 0; i < cnDimension * 3; i+=3)
		cout << i << " " << max_response[i] << " " << max_response[i+1] << " " << max_response[i+2] << endl;*/

	cout << endl;

	clReleaseKernel(kernel);
	clReleaseProgram(program);
	clReleaseMemObject(device_mem_im_vals);
	clReleaseMemObject(device_mem_local_max_vals);
	clReleaseMemObject(device_mem_max_response_r);
	clReleaseMemObject(device_mem_max_response_c);
	clReleaseMemObject(device_mem_max_response_z);
	clReleaseCommandQueue(queue);
	clReleaseContext(context);

	cout << endl;

#elif CUDA
    initialClustering_CUDA(im_vals, local_max_vals, max_response_r, max_response_c, max_response_z, r, c, z, scale_xy, scale_z);
#else
	int min_r, min_c, min_z, max_r, max_c, max_z;	

	#pragma omp parallel for private(min_r, min_c, min_z, max_r, max_c, max_z)
	for(int i=0; i<r; i++)
    {
		for(int j=0; j<c; j++)
        {
			for(int k=0; k<z; k++)
			{        				
				min_r = (int) std::max((double)(0.0),(double)(i-scale_xy));
				min_c = (int) std::max((double)(0.0),(double)(j-scale_xy));
				min_z = (int) std::max((double)(0.0),(double)(k-scale_z));
				max_r = (int) std::min((double)(r-1),(double)(i+scale_xy));
				max_c = (int) std::min((double)(c-1),(double)(j+scale_xy));                         
				max_z = (int) std::min((double)(z-1),(double)(k+scale_z));

				int R, C, Z;

				if(local_max_vals[(k*r*c)+(i*c)+j] !=0)//local_max_im[i][j][k]!=0)			
                    continue;					
				else
				{
					get_maximum(im_vals, min_r, max_r, min_c, max_c, min_z, max_z, &R, &C, &Z, r, c, z);                                              
					
					max_response_r[i * (c * z) + j * z + k] = R;
					max_response_c[i * (c * z) + j * z + k] = C;
					max_response_z[i * (c * z) + j * z + k] = Z;
				}				
			}
        }
    }

#endif // OPENCL

	/*for (int k = 0; k < r * c * z * 3; k+=3)
		cout << k << " " << max_response[k] << " " << max_response[k+1] << " " << max_response[k+2] << endl;*/


	std::cout << "Max_response array done" << std::endl;

	for(int i=0; i<r; i++)
    {
		for(int j=0; j<c; j++)
        {
			for(int k=0; k<z; k++)
			{
				if(local_max_vals[(k*r*c)+(i*c)+j] !=0)//local_max_im[i][j][k]!=0)			
                    continue;					
				else
					max_nghbr_im[i][j][k] =	max_nghbr_im[	max_response_r[i * (c * z) + j * z + k]	]
														[	max_response_c[i * (c * z) + j * z + k]	]
														[	max_response_z[i * (c * z) + j * z + k]	];
			}
		}
	}
	
	delete [] max_response_r;
	delete [] max_response_c;
	delete [] max_response_z;


	std::cout << "Initial max_nghbr_im took " << (clock() - start_time_init_max_nghbr_im)/(float)CLOCKS_PER_SEC << " seconds" << std::endl;
	
    int change = 1;
    double LM;
	std::cout << "Entering main Clustering Loop" << std::endl;
	//Now continue to update until no more changes occur, eventually will have clusters pointing to seeds	
	int iterr = 0;
    while(change)	
    {			
		//For now, limit it to a maximum of 10 iterations
		iterr++;
		if(iterr == 10)
			break;		
        change=0;
		
		
        for(int i=0; i<r; i++)
        {
            for(int j=0; j<c; j++)
            {
				for(int k=0; k<z; k++)
				{                   
					LM = max_nghbr_im[i][j][k];
					if(LM==0)
						continue;													    
					
					//Calculate coordinates of local maximum based on its index
					int rem = ((long)LM) % (r*c);
					int Z = ((int)LM-rem) / (r*c); 
				    int C = ((long)rem) % c;
					int R = (rem-C)/c;
					
                
					//Check for seed (already a local max value) or connected to seed (my local maximum is a seed point)					
					if(local_max_vals[(k*r*c)+(i*c)+j] !=0 || local_max_vals[(Z*r*c)+(R*c)+C]!=0 ) 
						continue;
					else
					{
						change++;
						max_nghbr_im[i][j][k]=max_nghbr_im[R][C][Z];
					}
				}
            }
        }
		std::cout<< "change=" << change << std::endl;
    }
    
	//cerr << "Preparing Output" << endl;
    
    for(int i=0; i<r; i++)
    {        
        for(int j=0; j<c; j++)
        {
			for(int k=0; k<z; k++)
            {
                LM = max_nghbr_im[i][j][k];
								    
                //if(local_max_vals[(int)LM] == -1 || bImg[(k*r*c)+(i*c)+j]==0)
				if(local_max_vals[(int)LM] == 65535 || bImg[(k*r*c)+(i*c)+j]==0)					
                    out1[(k*r*c)+(i*c)+j] = 0;
                else
				{
					//modified by Yousef on 8/21/2009
					//if the distance between me and my seed is more than a threshold.. then ignore me
					/*int rem = ((long)LM) % (r*c);
					int Z = (LM-rem) / (r*c); 
				    int C = ((long)rem) % c;
					int R = (rem-C)/c;
					double d = (i-R)*(i-R) + (j-C)*(j-C) + 3*(k-Z)*(k-Z);
					d = sqrt(d);
					if(d>10)
						out1[(k*r*c)+(i*c)+j] = 0;
					else*/
						out1[(k*r*c)+(i*c)+j] =local_max_vals[(int)LM];
				}
            }
        }
    }    

    #pragma omp parallel for
    for(int i=0; i<r; i++)
    {        
        for(int j=0; j<c; j++)
        {			
			free(max_nghbr_im[i][j]);
        }
		free(max_nghbr_im[i]);
    }
	free(max_nghbr_im);
}

void clustering_pfn_notify(const char *errinfo, const void *private_info, size_t cb, void *user_data)
{
	std::cerr << errinfo << std::endl;
}
 
