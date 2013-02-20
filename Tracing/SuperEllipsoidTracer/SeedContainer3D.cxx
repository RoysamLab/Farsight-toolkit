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

#include "SeedContainer3D.h"
#include "itkImageFileWriter.h"
#include <itkImageLinearIteratorWithIndex.h>
#include <itkImageSliceConstIteratorWithIndex.h>



SeedContainer3D::SeedContainer3D()   {
}

SeedContainer3D::~SeedContainer3D() {

	std::vector<SeedPoint3D *>::iterator it;
	for (it = SeedContainer.begin(); it != SeedContainer.end(); ++it)
               delete *it;
}


void SeedContainer3D::Configure(TraceConfig::Pointer m_Config)	{
	this->GridSpacing = m_Config->getGridSpacing();
}

void SeedContainer3D::Detect(ImageType3D::Pointer im3D, ImageType2D::Pointer im2D)	{
	//Obtian the MIP
	MIPGenerator(im3D, im2D);
	Detect3Dseeds(im3D);
	//Detect2Dseeds(im2D);
	//LocateZPosition(im3D);
}


void SeedContainer3D::Detect3Dseeds(ImageType3D::Pointer image) {
	ImageType3D::IndexType glndx = {{0 , 0 , 0}};
	ImageType3D::SizeType sz = image->GetBufferedRegion().GetSize();
	for (glndx[2] = 0; glndx[2] < sz[2]-GridSpacing; glndx[2] += GridSpacing) {
		for (glndx[1] = 0; glndx[1] < sz[1]-GridSpacing; glndx[1] += GridSpacing) {
			for (glndx[0] = 0; glndx[0] < sz[0]-GridSpacing; glndx[0] += GridSpacing) {

				ImageType3D::IndexType ndx = glndx;
				ImageType3D::IndexType endx;
				endx[0] = glndx[0] + GridSpacing;	endx[0] = (endx[0]>(unsigned int)sz[0]) ? sz[0] : endx[0];
				endx[1] = glndx[1] + GridSpacing;	endx[1] = (endx[1]>(unsigned int)sz[1]) ? sz[1] : endx[1];
				endx[2] = glndx[2] + GridSpacing;	endx[2] = (endx[2]>(unsigned int)sz[2]) ? sz[2] : endx[2];

				ImageType3D::PixelType minVal = 255;
				ImageType3D::IndexType minNdx;
				for ( ndx[2] = glndx[2] ; ndx[2] < endx[2] ;  ndx[2]++)	{
					for ( ndx[1] = glndx[1] ; ndx[1] < endx[1] ;  ndx[1]++)	{
						for ( ndx[0] = glndx[0] ; ndx[0] < endx[0] ;  ndx[0]++)	{
							ImageType3D::PixelType val = image->GetPixel(ndx);
							if (val<minVal)	{
								minVal = val;
								minNdx = ndx;
								}
							}
						}
					}

				//store the index, val, and scale  here
				AddSeedPoint(minNdx, minVal, GridSpacing);
				}
			}
		}
	}

void SeedContainer3D::Detect2Dseeds(ImageType2D::Pointer image) {

	ImageType2D::RegionType::SizeType size = image->GetRequestedRegion().GetSize();
	ImageType2D::RegionType::IndexType p, pmin;
	float val, m = 1000.0f;

	std::cout << "Detecting Seeds ...";

	for ( p[0] = 0; p[0] <(unsigned int)size[0] ; p[0]+=GridSpacing)	{
		for ( p[1] = 0; p[1] <(unsigned int)size[1] ; p[1]++)	{

			val = image->GetPixel(p);
			if(val < m)	{
				m = val;
				pmin = p;
			}

			if(!(p[1]%GridSpacing) && (p[1]>0))	{
				AddSeedPoint(pmin, m);
				m = 1000.0f;
			}
		}
	}

	for ( p[1] = 0; p[1] <(unsigned int)size[1] ; p[1]+=GridSpacing)	{
			for ( p[0] = 0; p[0] <(unsigned int)size[0] ; p[0]++)	{

			val = image->GetPixel(p);
			if(val < m)	{
				m = val;
				pmin = p;
			}

			if(!(p[0]%GridSpacing) && (p[0]>0))	{
				AddSeedPoint(pmin, m);
				m = 1000.0f;
			}
		}
	}

	//std::cout <<SeedContainer.size() <<" seeds found!!"<< std::endl;
}

void SeedContainer3D::AddSeedPoint(ImageType2D::RegionType::IndexType ndx, const float& mn)	{
	SeedPoint3D *s = new SeedPoint3D();
	s->setXYPosition(ndx);
	s->setScale(static_cast<float>(GridSpacing)/4);
	this->SeedContainer.push_back(s);
	//std::cout <<"Seed detected at " << ndx << " Value" << mn <<std::endl;
}

void SeedContainer3D::AddSeedPoint(ImageType3D::IndexType ndx, ImageType3D::PixelType val, long grid)	{
	SeedPoint3D *s = new SeedPoint3D();
	s->setSeed3D(ndx, val, grid);
	this->SeedContainer.push_back(s);
	}

void SeedContainer3D::MIPGenerator(ImageType3D::Pointer im3D, ImageType2D::Pointer& im2D)	{


	ImageType3D::SizeType sz3 = im3D->GetRequestedRegion().GetSize();
	ImageType2D::SizeType sz2 = {{sz3[0], sz3[1] }};

	//im2D = ImageType2D::New();  // new is already done in
	im2D->SetRegions(sz2 );
	im2D->Allocate();


	for (unsigned long x=0; x<sz3[0]; x++)	{
		for (unsigned long y=0; y<sz3[1]; y++)	{

			//ImageType3D::PixelType maxVal = 0.0;
			ImageType3D::PixelType minVal = 255.0;

			for (unsigned long z=0; z<sz3[2]; z++)	{
				ImageType3D::IndexType nd3 = {{x,y,z}};
				ImageType3D::PixelType val = im3D->GetPixel(nd3);
			//	maxVal = (maxVal > val) ? maxVal : val;
				minVal = (minVal < val) ? minVal : val;
			}

			ImageType2D::IndexType nd2 = {{x, y}};
			//im2D->SetPixel(nd2, maxVal);
			im2D->SetPixel(nd2, minVal);
		}
	}


/*	itk::ImageFileWriter<ImageType2D>::Pointer wt = itk::ImageFileWriter<ImageType2D>::New();
	wt->SetFileName("MIP.mhd");
	wt->SetInput(im2D);
	wt->Update();
*/
	std::cout <<"MIP Generated!!" <<std::endl;
}



void SeedContainer3D::LocateZPosition(ImageType3D::Pointer im3D)	{
	std::vector <SeedPoint3D *>::iterator it;
	ImageType3D::RegionType::SizeType sz = im3D->GetRequestedRegion().GetSize();
	ImageType3D::RegionType::IndexType ndx;
	float currentValue, miniminValue = 255;
	unsigned int minimumNdx = 0;

	for (it = SeedContainer.begin(); it!= SeedContainer.end(); ++it)	{
		ndx = (*it)->getIndex();
		miniminValue = 255;
		for (unsigned int i = 0; i < sz[2]; i++)	{
			ndx[2] = i;
			currentValue = 	im3D->GetPixel(ndx);
			if (currentValue < miniminValue)	{
				miniminValue = currentValue;
				minimumNdx = i;
			}
		}
		(*it)->setZPosition(minimumNdx);
		(*it)->setIntensity(miniminValue);
	}
}
