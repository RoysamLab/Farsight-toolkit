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

#include "gradientweighteddistancemap.h" 

//#define grad_im_2D(a,b) GRAD_IM_2D[(a-1)+(b-1)*size1]
//#define grad_imw(a,b) GRAD_IMW[(a)+(b)*(size1+2)]

//This function takes in an image with the background pixels at -inf,
// the pixels from which the gradient enhanced distance map is to be initialized at 0,
// and the pixels where they have  to be computed at +inf
//Edited - Kedar 06/30/2009


//By Yousef: define is now working, so let's try to use functions
 float grad_im_2D(int a,int b, float* GRAD_IM_2D, int size1) 
 {
	 return GRAD_IM_2D[(a-1)+(b-1)*size1];
 }
 float grad_imw(int a,int b,float* GRAD_IMW, int size1)
 {
	 return GRAD_IMW[(a)+(b)*(size1+2)];
 }

int gradient_enhanced_distance_map_2D( float *GRAD_IM_2D, float *GRAD_IMW, int size1, int size2 ){


//Copy data type floats' lower limit into a and get root2 into a local variable
	float flot_mini,root2;

	flot_mini = FLT_MIN;	
	
	root2 = sqrt(2.0);

	std::queue<im_ind> icy_needs_a_change;	

//Initialize queue with pixels that are neighbours of the the pixels at zero regions
	for( int j=2; j<size2-1; j++ )
	{
		for( int i=2; i<size1-1; i++ )
		{
			if(	( grad_im_2D(i-1,j-1, GRAD_IM_2D, size1) != flot_mini && grad_im_2D(i,j, GRAD_IM_2D, size1) > grad_im_2D(i-1,j-1, GRAD_IM_2D, size1) ) || \
				( grad_im_2D(i-1,j, GRAD_IM_2D, size1)   != flot_mini && grad_im_2D(i,j, GRAD_IM_2D, size1) > grad_im_2D(i-1,j, GRAD_IM_2D, size1)   ) || \
				( grad_im_2D(i-1,j+1, GRAD_IM_2D, size1) != flot_mini && grad_im_2D(i,j, GRAD_IM_2D, size1) > grad_im_2D(i-1,j+1, GRAD_IM_2D, size1) ) || \
				( grad_im_2D(i,j-1, GRAD_IM_2D, size1)   != flot_mini && grad_im_2D(i,j, GRAD_IM_2D, size1) > grad_im_2D(i,j-1, GRAD_IM_2D, size1)   ) || \
				( grad_im_2D(i,j+1, GRAD_IM_2D, size1)   != flot_mini && grad_im_2D(i,j, GRAD_IM_2D, size1) > grad_im_2D(i,j+1, GRAD_IM_2D, size1)   ) || \
				( grad_im_2D(i+1,j-1, GRAD_IM_2D, size1) != flot_mini && grad_im_2D(i,j, GRAD_IM_2D, size1) > grad_im_2D(i+1,j-1, GRAD_IM_2D, size1) ) || \
				( grad_im_2D(i+1,j, GRAD_IM_2D, size1)   != flot_mini && grad_im_2D(i,j, GRAD_IM_2D, size1) > grad_im_2D(i+1,j, GRAD_IM_2D, size1)   ) || \
				( grad_im_2D(i+1,j+1, GRAD_IM_2D, size1) != flot_mini && grad_im_2D(i,j, GRAD_IM_2D, size1) > grad_im_2D(i+1,j+1, GRAD_IM_2D, size1) ) )  {
					im_ind temp_inds;
					temp_inds.i = i;
					temp_inds.j = j;
					icy_needs_a_change.push( temp_inds );
			}
		}
	}
//Pop and push pixels till all the pixels are set, i.e, queue is empty
	//float temp,temp1;
	int best_val_not_found; //added by Yousef (was producing a compilation error)
	while( !icy_needs_a_change.empty() ){
		im_ind temp_inds;
		int i,j,k;//,best_val_found;
		int neigh_vals[8],neigh_vals_cpy[8];
		//float temp_hold;
		best_val_not_found = 1;
		icy_needs_a_change.pop(); //modified by Yousef (was producing a compilation error)

		temp_inds = icy_needs_a_change.front();
		i = temp_inds.i; j = temp_inds.j;

		//Calculate the relative weighted distances from all the neighbors
		neigh_vals[0] = grad_im_2D(i-1,j-1, GRAD_IM_2D, size1)+root2*fabs(grad_imw(i-1,j-1, GRAD_IMW, size1)-grad_imw(i,j, GRAD_IMW, size1));
		neigh_vals[1] = grad_im_2D(i-1,j, GRAD_IM_2D, size1)+fabs(grad_imw(i-1,j, GRAD_IMW, size1)-grad_imw(i,j, GRAD_IMW, size1));
		neigh_vals[2] = grad_im_2D(i-1,j+1, GRAD_IM_2D, size1)+root2*fabs(grad_imw(i-1,j+1, GRAD_IMW, size1)-grad_imw(i,j, GRAD_IMW, size1));
		neigh_vals[3] = grad_im_2D(i,j-1, GRAD_IM_2D, size1)+fabs(grad_imw(i,j-1, GRAD_IMW, size1)-grad_imw(i,j, GRAD_IMW, size1));
		neigh_vals[4] = grad_im_2D(i,j+1, GRAD_IM_2D, size1)+fabs(grad_imw(i,j+1, GRAD_IMW, size1)-grad_imw(i,j, GRAD_IMW, size1));
		neigh_vals[5] = grad_im_2D(i+1,j-1, GRAD_IM_2D, size1)+root2*fabs(grad_imw(i+1,j-1, GRAD_IMW, size1)-grad_imw(i,j, GRAD_IMW, size1));
		neigh_vals[6] = grad_im_2D(i+1,j, GRAD_IM_2D, size1)+fabs(grad_imw(i+1,j, GRAD_IMW, size1)-grad_imw(i,j, GRAD_IMW, size1));
		neigh_vals[7] = grad_im_2D(i+1,j+1, GRAD_IM_2D, size1)+root2*fabs(grad_imw(i+1,j+1, GRAD_IMW, size1)-grad_imw(i,j, GRAD_IMW, size1));

		//Duplicate and order the list
		for( k=0; k<8; k++ ) neigh_vals_cpy[k] = neigh_vals[k];
		int num_elements = sizeof(neigh_vals_cpy) / sizeof(neigh_vals_cpy[0]); 
		std::sort(neigh_vals_cpy, neigh_vals_cpy + num_elements);

		for( k=0; k<8; k++ ){
			//Skip value if the computed gradient enhanced value is in either the background of it is a repeated value
			if( neigh_vals_cpy[k] < 0 ) continue;
			if( k>0 ) if( neigh_vals_cpy[k-1] == neigh_vals_cpy[k] ) continue;
			//At this point the 
			if( best_val_not_found ){
				GRAD_IM_2D[(i-1)+(j-1)*size1]= neigh_vals_cpy[k];
				best_val_not_found = 0;
				continue;
			}
			if( neigh_vals_cpy[k] ==  neigh_vals[0] ){
				temp_inds.i = i-1; temp_inds.j = j-1;
				icy_needs_a_change.push( temp_inds );
			}
			if( neigh_vals_cpy[k] ==  neigh_vals[1] ){
				temp_inds.i = i-1; temp_inds.j = j;
				icy_needs_a_change.push( temp_inds );
			}
			if( neigh_vals_cpy[k] ==  neigh_vals[2] ){
				temp_inds.i = i-1; temp_inds.j = j+1;
				icy_needs_a_change.push( temp_inds );
			}
			if( neigh_vals_cpy[k] ==  neigh_vals[3] ){
				temp_inds.i = i; temp_inds.j = j-1;
				icy_needs_a_change.push( temp_inds );
			}
			if( neigh_vals_cpy[k] ==  neigh_vals[4] ){
				temp_inds.i = i; temp_inds.j = j+1;
				icy_needs_a_change.push( temp_inds );
			}
			if( neigh_vals_cpy[k] ==  neigh_vals[5] ){
				temp_inds.i = i+1; temp_inds.j = j-1;
				icy_needs_a_change.push( temp_inds );
			}
			if( neigh_vals_cpy[k] ==  neigh_vals[6] ){
				temp_inds.i = i+1; temp_inds.j = j;
				icy_needs_a_change.push( temp_inds );
			}
			if( neigh_vals_cpy[k] ==  neigh_vals[7] ){
				temp_inds.i = i+1; temp_inds.j = j+1;
				icy_needs_a_change.push( temp_inds );
			}
		}
	}

	return 1;
}

/* //Older and naieve implementation. Last three varibles pass to the function have been discarded in the above implementation
int gradient_weighted_distance_map( float *GRAD_IM_2D, float *GRAD_IMW, int size1, int size2, int count,float large_val, int scaling_pass ){
 	float temp;
	int hope_change=1;
	int counter=0;
	while((hope_change) && (counter < count)) //when there is hope and change there is still work to be done
	{
		hope_change=0; counter++;

		for( int j=1; j<size2+1; j++ ){
			for( int i=1; i<size1+1; i++ ){
				if( grad_imw(i,j) > 0 ){
					if( grad_imw(i-1,j) > 0 ){
					if( grad_im_2D(i,j) > grad_im_2D(i-1,j) ) temp = grad_im_2D(i,j)-grad_im_2D(i-1,j)+1; else temp = grad_im_2D(i-1,j)-grad_im_2D(i,j)+1;
					if( ( grad_imw(i,j)-grad_imw(i-1,j)) > temp ) { grad_imw(i,j) = grad_imw(i-1,j)+temp; hope_change++; }}

					if( grad_imw(i+1,j) > 0 ){
					if( grad_im_2D(i,j) > grad_im_2D(i+1,j) ) temp = grad_im_2D(i,j)-grad_im_2D(i+1,j)+1; else temp = grad_im_2D(i+1,j)-grad_im_2D(i,j)+1;
					if( ( grad_imw(i,j)-grad_imw(i+1,j)) > temp ) { grad_imw(i,j) = grad_imw(i+1,j)+temp; hope_change++; }}

					if( grad_imw(i,j-1) > 0 ){
					if( grad_im_2D(i,j) > grad_im_2D(i,j-1) ) temp = grad_im_2D(i,j)-grad_im_2D(i,j-1)+1; else temp = grad_im_2D(i,j-1)-grad_im_2D(i,j)+1;
					if( ( grad_imw(i,j)-grad_imw(i,j-1)) > temp ) { grad_imw(i,j) = grad_imw(i,j-1)+temp; hope_change++; }}

					if( grad_imw(i,j+1) > 0 ){
					if( grad_im_2D(i,j) > grad_im_2D(i,j+1) ) temp = grad_im_2D(i,j)-grad_im_2D(i,j+1)+1; else temp = grad_im_2D(i,j+1)-grad_im_2D(i,j)+1;
					if( ( grad_imw(i,j)-grad_imw(i,j+1)) > temp ) { grad_imw(i,j) = grad_imw(i,j+1)+temp; hope_change++; }}

					if( grad_imw(i-1,j-1) > 0 ){
					if( grad_im_2D(i,j) > grad_im_2D(i-1,j-1) ) temp = grad_im_2D(i,j)-grad_im_2D(i-1,j-1)+root2; else temp = grad_im_2D(i-1,j-1)-grad_im_2D(i,j)+root2;
					if( ( grad_imw(i,j)-grad_imw(i-1,j-1)) > temp ) { grad_imw(i,j) = grad_imw(i-1,j-1)+temp; hope_change++; }}

					if( grad_imw(i+1,j-1) > 0 ){
					if( grad_im_2D(i,j) > grad_im_2D(i+1,j-1) ) temp = grad_im_2D(i,j)-grad_im_2D(i+1,j-1)+root2; else temp = grad_im_2D(i+1,j-1)-grad_im_2D(i,j)+root2;
					if( ( grad_imw(i,j)-grad_imw(i+1,j-1)) > temp ) { grad_imw(i,j) = grad_imw(i+1,j-1)+temp; hope_change++; }}

					if( grad_imw(i-1,j+1) > 0 ){
					if( grad_im_2D(i,j) > grad_im_2D(i-1,j+1) ) temp = grad_im_2D(i,j)-grad_im_2D(i-1,j+1)+root2; else temp = grad_im_2D(i-1,j+1)-grad_im_2D(i,j)+root2;
					if( ( grad_imw(i,j)-grad_imw(i-1,j+1)) > temp ) { grad_imw(i,j) = grad_imw(i-1,j+1)+temp; hope_change++; }}

					if( grad_imw(i+1,j+1) > 0 ){
					if( grad_im_2D(i,j) > grad_im_2D(i+1,j+1) ) temp = grad_im_2D(i,j)-grad_im_2D(i+1,j+1)+root2; else temp = grad_im_2D(i+1,j+1)-grad_im_2D(i,j)+root2;
					if( ( grad_imw(i,j)-grad_imw(i+1,j+1)) > temp ) { grad_imw(i,j) = grad_imw(i+1,j+1)+temp; hope_change++; }}
				}
			}
		}

		for( int j=size2; j>0; j-- ){
			for( int i=size1; i>0; i-- ){
				if( grad_imw(i,j) > 0 ){
					if( grad_imw(i-1,j) > 0 ){
					if( grad_im_2D(i,j) > grad_im_2D(i-1,j) ) temp = grad_im_2D(i,j)-grad_im_2D(i-1,j)+1; else temp = grad_im_2D(i-1,j)-grad_im_2D(i,j)+1;
					if( ( grad_imw(i,j)-grad_imw(i-1,j)) > temp ) { grad_imw(i,j) = grad_imw(i-1,j)+temp; hope_change++; }}

					if( grad_imw(i+1,j) > 0 ){
					if( grad_im_2D(i,j) > grad_im_2D(i+1,j) ) temp = grad_im_2D(i,j)-grad_im_2D(i+1,j)+1; else temp = grad_im_2D(i+1,j)-grad_im_2D(i,j)+1;
					if( ( grad_imw(i,j)-grad_imw(i+1,j)) > temp ) { grad_imw(i,j) = grad_imw(i+1,j)+temp; hope_change++; }}

					if( grad_imw(i,j-1) > 0 ){
					if( grad_im_2D(i,j) > grad_im_2D(i,j-1) ) temp = grad_im_2D(i,j)-grad_im_2D(i,j-1)+1; else temp = grad_im_2D(i,j-1)-grad_im_2D(i,j)+1;
					if( ( grad_imw(i,j)-grad_imw(i,j-1)) > temp ) { grad_imw(i,j) = grad_imw(i,j-1)+temp; hope_change++; }}

					if( grad_imw(i,j+1) > 0 ){
					if( grad_im_2D(i,j) > grad_im_2D(i,j+1) ) temp = grad_im_2D(i,j)-grad_im_2D(i,j+1)+1; else temp = grad_im_2D(i,j+1)-grad_im_2D(i,j)+1;
					if( ( grad_imw(i,j)-grad_imw(i,j+1)) > temp ) { grad_imw(i,j) = grad_imw(i,j+1)+temp; hope_change++; }}

					if( grad_imw(i-1,j-1) > 0 ){
					if( grad_im_2D(i,j) > grad_im_2D(i-1,j-1) ) temp = grad_im_2D(i,j)-grad_im_2D(i-1,j-1)+root2; else temp = grad_im_2D(i-1,j-1)-grad_im_2D(i,j)+root2;
					if( ( grad_imw(i,j)-grad_imw(i-1,j-1)) > temp ) { grad_imw(i,j) = grad_imw(i-1,j-1)+temp; hope_change++; }}

					if( grad_imw(i+1,j-1) > 0 ){
					if( grad_im_2D(i,j) > grad_im_2D(i+1,j-1) ) temp = grad_im_2D(i,j)-grad_im_2D(i+1,j-1)+root2; else temp = grad_im_2D(i+1,j-1)-grad_im_2D(i,j)+root2;
					if( ( grad_imw(i,j)-grad_imw(i+1,j-1)) > temp ) { grad_imw(i,j) = grad_imw(i+1,j-1)+temp; hope_change++; }}

					if( grad_imw(i-1,j+1) > 0 ){
					if( grad_im_2D(i,j) > grad_im_2D(i-1,j+1) ) temp = grad_im_2D(i,j)-grad_im_2D(i-1,j+1)+root2; else temp = grad_im_2D(i-1,j+1)-grad_im_2D(i,j)+root2;
					if( ( grad_imw(i,j)-grad_imw(i-1,j+1)) > temp ) { grad_imw(i,j) = grad_imw(i-1,j+1)+temp; hope_change++; }}

					if( grad_imw(i+1,j+1) > 0 ){
					if( grad_im_2D(i,j) > grad_im_2D(i+1,j+1) ) temp = grad_im_2D(i,j)-grad_im_2D(i+1,j+1)+root2; else temp = grad_im_2D(i+1,j+1)-grad_im_2D(i,j)+root2;
					if( ( grad_imw(i,j)-grad_imw(i+1,j+1)) > temp ) { grad_imw(i,j) = grad_imw(i+1,j+1)+temp; hope_change++; }}
				}
			}
		}

		//temporaray fix for regions not connected to nuclei
		if( counter == scaling_pass ){
			for( int i=1; i<(size1+2); i++ ){
				for( int j=1; j<(size2+2); j++ ){
					if( grad_imw(i,j) == large_val ) grad_imw(i,j) = 0;
				}
			}
		}
	}//Will hope and chage be unnecessary when it is all done??
	
	//Change background pixels to 0 from their original value of -1
	//for( int i=0; i<(size1+2); i++ )
	//	for( int j=0; j<(size2+2); j++ )
	//		if( grad_imw(i,j) == -1 ) grad_imw(i,j)=0;

	//std::cout<<"Change: "<<hope_change<<std::endl;
	//std::cout<<"Count: "<<counter<<std::endl;
	
	return 1;
}
*/

int computeGradientImage(unsigned char *IM, float* gradIM, int r, int c, int z)
{
	typedef    unsigned char     InputPixelType;
	typedef    float     OutputPixelType;
	typedef itk::Image< InputPixelType,  3 >   InputImageType;
	typedef itk::Image< OutputPixelType, 3 >   OutputImageType;
	

	InputImageType::Pointer im;
	im = InputImageType::New();
	InputImageType::PointType origin;
    origin[0] = 0; 
    origin[1] = 0;    	
	origin[2] = 0;    
    im->SetOrigin( origin );

    InputImageType::IndexType start;
    start[0] =   0;  // first index on X
    start[1] =   0;  // first index on Y    	
	start[2] =   0;  // first index on Z    
    InputImageType::SizeType  size;
    size[0]  = c;  // size along X
    size[1]  = r;  // size along Y
	size[2]  = z;  // size along Z

    InputImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );

    im->SetRegions( region );
    im->Allocate();
    im->FillBuffer(0);
	im->Update();
	
	//copy the input image into the ITK image
	typedef itk::ImageRegionIteratorWithIndex< InputImageType > IteratorType;
	IteratorType iterator1(im,im->GetRequestedRegion());
	for(int i=0; i<r*c*z; i++)
	{
		iterator1.Set(IM[i]);
		++iterator1;	
	}

	typedef itk::GradientMagnitudeImageFilter<InputImageType, OutputImageType > FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput( im );
	filter->Update();

	typedef itk::ImageRegionIteratorWithIndex< OutputImageType > IteratorType2;
	IteratorType2 iterator2(filter->GetOutput(),filter->GetOutput()->GetRequestedRegion());	

	//try to write that image for now
	/*typedef itk::ImageFileWriter< MyInputImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName("grad.tif");
	writer->SetInput( filter->GetOutput() );
	writer->Update();	*/

	for(int i=0; i<r*c*z; i++)
	{
		gradIM[i] = iterator2.Get();	
		++iterator2;
	}	

	return 1;
}

void prepareInputImage(unsigned short *binImage, unsigned short* seedsImage, float* outImage, int r, int c, int z)
{
	//Mempry for outImage must be allocated outside
	if(!outImage)
		return;
	
	for(int k=0; k<z; k++)
	{
		for(int i=0; i<r; i++)
		{
			for(int j=0; j<c; j++)
			{
				if(binImage[k*r*c+i*c+j]==0)
					outImage[k*r*c+i*c+j] = FLT_MIN;
				else
				{
					if(seedsImage[k*r*c+i*c+j]>0)
						outImage[k*r*c+i*c+j] = 0 ;
					else
						outImage[k*r*c+i*c+j] = FLT_MAX ;
				}
			}
		}
	}
}
