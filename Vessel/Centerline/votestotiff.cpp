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

#include <stdio.h>
#include <itkImage.h>
#include <itkImageFileWriter.h>

typedef unsigned char PixelType;
typedef itk::Image<PixelType,3> ImageType;
typedef itk::ImageFileWriter<ImageType> FileWriterType;

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

int main()
{
	std::string votesfilename = "../../half_vessel_votes2.txt";
	std::string outputfilename = "half_vessel_votes2.tif";
	FILE * fp = fopen(votesfilename.c_str(), "r");


	int x,y,z,v;
	int maxx = -1,maxy = -1, maxz = -1;
	while(fscanf(fp,"%d %d %d %d",&x,&y,&z,&v)>0)
	{
		maxx = MAX(maxx, x);
		maxy = MAX(maxy, y);
		maxz = MAX(maxz, z);
	}

	printf("maxx = %d maxy = %d maxz = %d\n",maxx,maxy,maxz);
	fclose(fp);


	fp = fopen(votesfilename.c_str(), "r");


	ImageType::Pointer im = ImageType::New();

	ImageType::SizeType size;
	ImageType::IndexType index;
	ImageType::RegionType region;

	size[0] = maxx+1;
	size[1] = maxy+1;
	size[2] = maxz+1;

	index[0] = index[1] = index[2] = 0;

	region.SetIndex(index);
	region.SetSize(size);
	im->SetRegions(region);
	im->Allocate();

	while(fscanf(fp,"%d %d %d %d",&x,&y,&z,&v)>0)
	{
		index[0] = x;
		index[1] = y;
		index[2] = z;
		im->SetPixel(index,v);
	}

	//	scanf("%*d");

	fclose(fp);

	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(outputfilename.c_str());
	writer->SetInput(im);
	writer->Update();


	return 0;
}
