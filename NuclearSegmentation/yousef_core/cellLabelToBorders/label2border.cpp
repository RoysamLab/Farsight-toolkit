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

#include "label2border.h"

void label2border(int* img,int c, int r, int z)
{
	const int dim = 3;
  
	typedef int PixelType;
	typedef itk::Image< PixelType, dim > ImageType;

	//read the input image
	//Create an itk image
	ImageType::Pointer myItkImg;
	myItkImg = ImageType::New();
	ImageType::PointType origin;
    origin[0] = 0.; 
    origin[1] = 0.;    
	origin[2] = 0.;    
    myItkImg->SetOrigin( origin );

    ImageType::IndexType start;
    start[0] =   0;  // first index on X
    start[1] =   0;  // first index on Y    
	start[2] =   0;  // first index on Z    
    ImageType::SizeType  size;
    size[0]  = c;  // size along X
    size[1]  = r;  // size along Y
	size[2]  = z;  // size along Z
  
    ImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );
    
    myItkImg->SetRegions( region );
    myItkImg->Allocate();
    myItkImg->FillBuffer(0);
	myItkImg->Update();
	
	//copy the input image into the ITK image
	typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
	IteratorType iterator1(myItkImg,myItkImg->GetRequestedRegion());
	for(long int i=0; i<r*c*z; i++)
	{		
		iterator1.Set(img[i]);
		++iterator1;	
	}

	typedef itk::LabelBorderImageFilter< ImageType, ImageType > FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput( myItkImg );
	filter->SetFullyConnected( false );
	filter->SetBackgroundValue( 0 );
	filter->Update();
	//filter->Print(std::cerr);

	//   Copy the resulting image into the input array
	IteratorType iterator2(filter->GetOutput(),filter->GetOutput()->GetRequestedRegion());
	//IteratorType iterator2(myItkImg,myItkImg->GetRequestedRegion());
	for(long int i=0; i<r*c*z; i++)
	{
		img[i] = iterator2.Get();
		++iterator2;
	}

	/*
	typedef itk::ImageFileWriter< ImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput( filter->GetOutput() );
	writer->SetFileName( "labelImg.tiff" );
	writer->Update();
	*/

}

//Added By Yousef Al-Kofahi on 5-29-2008
void label2border2(int* img,int* bordImg, int c, int r, int z)
{
	int v, v1, v2, v3;
	for(int i=0; i<r-1; i++)
	{
		for(int j=0; j<c-1; j++)
		{
			for(int k=0; k<z; k++)
			{				
				v=img[(k*r*c)+(i*c)+j];
				v1 = img[(k*r*c)+((i+1)*c)+j];
				v2 = img[(k*r*c)+(i*c)+j+1];
				v3 = img[(k*r*c)+((i+1)*c)+(j+1)];
				
				if(v!=v1 || v!=v2 || v!=v3)
					bordImg[(k*r*c)+(i*c)+j] = 255;
				else
					bordImg[(k*r*c)+(i*c)+j] = 0;
			}
		}
	}
}

void label2border3(int* img,int* segImg, int c, int r, int z)
{
	//int v, v1, v2, v3;
	for(int i=0; i<r-1; i++)
	{
		for(int j=0; j<c-1; j++)
		{
			for(int k=0; k<z-1; k++)
			{		
				if(segImg[(k*r*c)+(i*c)+j]==0)
					continue;

				for(int l1=0; l1<2; l1++)
				{
					for(int l2=0; l2<2; l2++)
					{
						for(int l3=0; l3<2; l3++)
						{
							if(segImg[(k*r*c)+(i*c)+j] != segImg[((k+l3)*r*c)+((i+l2)*c)+(j+l1)] && segImg[((k+l3)*r*c)+((i+l2)*c)+(j+l1)]>0)
								img[(k*r*c)+(i*c)+j] = 0;
						}				
					}
				}
			}
		}
	}
	for(int i=0; i<r; i++)
	{
		for(int j=0; j<c; j++)
		{
			for(int k=0; k<z; k++)
			{						
			if(img[(k*r*c)+(i*c)+j]>0 && segImg[(k*r*c)+(i*c)+j]>0)
				img[(k*r*c)+(i*c)+j] = 255;	
			}
		}
	}
}


