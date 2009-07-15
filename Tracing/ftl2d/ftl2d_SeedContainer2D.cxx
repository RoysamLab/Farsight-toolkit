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

#include "ftl2d_SeedContainer2D.h"
//#include "SeedPoint2D.h"
#include <itkImageLinearConstIteratorWithIndex.h>

bool SortIntensityFunction(SeedPoint2D*  , SeedPoint2D* );


SeedContainer2D::SeedContainer2D(TraceConfig *m_Config)   {
	this->GridSpacing = m_Config->getGridSpacing();
	this->valid = 1;
	this->IntensityThreshold = m_Config->getSeedIntensityThreshold();
}

SeedContainer2D::~SeedContainer2D() {

	std::vector<SeedPoint2D *>::iterator it;
	for (it = SeedContainer.begin(); it != SeedContainer.end(); ++it)
               delete *it;
}


void SeedContainer2D::Detect(ImageType::Pointer image) {

	ImageType::RegionType::SizeType size = image->GetRequestedRegion().GetSize();
	ImageType::RegionType::IndexType p, pmin;
	float val, m = 1000.0f;
	
	std::cout << "Detecting Seeds ...";
	
	for ( p[0] = 0; p[0] <size[0] ; p[0]+=GridSpacing)	{
		for ( p[1] = 0; p[1] <size[1] ; p[1]++)	{
			
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

	for ( p[1] = 0; p[1] <size[1] ; p[1]+=GridSpacing)	{
			for ( p[0] = 0; p[0] <size[0] ; p[0]++)	{
				
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
		
	std::cout <<SeedContainer.size() <<" seeds found!!"<< std::endl;
}

void SeedContainer2D::AddSeedPoint(ImageType::RegionType::IndexType ndx, const float& mn)	{
	SeedPoint2D *s = new SeedPoint2D();
	s->setx(ndx[0]);
	s->sety(ndx[1]);
	s->setIntensity(mn);
	//s->setScale(static_cast<float>(GridSpacing)/4);
	s->setScale(2.5);
	this->SeedContainer.push_back(s);
//	std::cout <<"Pixel indices " << ndx <<std::endl;
//	std::cout <<"Pixel " << s->getx() <<" "<< s->gety() <<" "<<s->getIntensity() <<" "<< s->getScale() <<std::endl;
//	std::cin.get();
}

/*void SeedContainer2D::Detect_old(ImageType::Pointer image) {

	
	typedef itk::ImageLinearConstIteratorWithIndex <ImageType> CIterType;
	CIterType It(image, image->GetRequestedRegion() );
	
	ImageType::RegionType::SizeType size = image->GetRequestedRegion().GetSize();
	
	std::cout << "Detecting Seeds"<< std::endl;
	//Scan for the x direction
	It.SetDirection(0);
	It.GoToBegin();

	int line_count = 0;
	
	
	
	while ( !It.IsAtEnd() )	{
		
		ImageType::IndexType ndx1;
		It.GoToBeginOfLine();
		while (!It.IsAtEndOfLine() )	{
			float p, mn = 1000;
			SeedPoint2D *s = new SeedPoint2D();
			for (int i=0;i<GridSpacing;i++, ++It)	{
				p = It.Get();
				if (mn > p)	{
					mn = p;
					ndx1 = It.GetIndex();
				}

				//std::cout <<"Pixel " << ndx1 << std::endl;
			}

			s->setx(ndx1[0]);
			s->sety(ndx1[1]);
			s->setIntensity(mn);
			s->setScale(static_cast<float>(GridSpacing)/4);

//			std::cout <<"Pixel " << s->getx() <<" "<< s->gety() <<" "<<s->getIntensity() <<" "<< s->getScale() <<std::endl;			
			AddSeedPoint(s);
		}
		
		#if DEBUG_LEVEL > 1
		SeedPoint2D *s1 =  SeedContainer.back();
		std::cout <<"Pixel " <<s1->getx() <<" "<< s1->gety() <<" "<<s1->getIntensity() <<" "<< s1->getScale() << " ^ " << ndx1[1] << std::endl;		
		#endif
		
		for (int i=0;i<GridSpacing;i++)
			It.NextLine();
		line_count++;		
		if (ndx1[1] > (size[1]-GridSpacing) )
			break;
	}
	
	std::cout << "ROWWISE OPERATION DONE - "<< line_count <<std::endl;
	
	
	//Scan for the y direction
	It.SetDirection(1);
	It.GoToBegin();

	line_count = 0;
	
		
	while ( !It.IsAtEnd() )	{
		
		ImageType::IndexType ndx1;
		It.GoToBeginOfLine();
		while (!It.IsAtEndOfLine() && !It.IsAtEnd() )	{
			float p, mn = 1000;
			SeedPoint2D *s = new SeedPoint2D();
			for (int i=0;(i<GridSpacing) && !(It.IsAtEnd());i++, ++It)	{
				
				if(It.IsAtEnd()) 
					break;
				
				p = It.Get();
				if (mn > p)	{
					mn = p;
					ndx1 = It.GetIndex();
					//std::cout <<"Pixel " << ndx1 << std::endl;
				}
				
			}
			s->setx(ndx1[0]);
			s->sety(ndx1[1]);
			s->setIntensity(mn);
			s->setScale(static_cast<float>(GridSpacing)/4);

			std::cout <<"Pixel " << s->getx() <<" "<< s->gety() <<" "<<s->getIntensity() <<" "<< s->getScale() <<std::endl;	
			AddSeedPoint(s);
		}
		#if DEBUG_LEVEL > 1
		SeedPoint2D *s1 =  SeedContainer.back();
		std::cout <<"Pixel " <<s1->getx() <<" "<< s1->gety() <<" "<<s1->getIntensity() <<" "<< s1->getScale() << " ^ " << ndx1[1] << std::endl;		
		#endif


		for (int i=0;i<GridSpacing;i++)
			It.NextLine();
		if (ndx1[0] > (size[0]-GridSpacing) )
			break;
	}
	std::cout << "COLUMNWISE OPERATION DONE"<< std::endl;

}
*/
void SeedContainer2D::Validate(ImageType::Pointer im)	{

/*	std::vector<SeedPoint2D *>::iterator it;

	for (it = SeedContainer.begin(); it!=SeedContainer.end(); ++it)	{
		if (*it->getIntensity < this->IntensityThreshold)
			this->valid = 1;
		else
			this->valid = 0;
*/
}

void SeedContainer2D::SortSeeds()	{
	std::sort(SeedContainer.begin(), SeedContainer.end(), SortIntensityFunction);
//	if (SeedContainer.size() > 2000)
//		SeedContainer.erase(SeedContainer.begin()+2000,SeedContainer.end());
}

bool SortIntensityFunction( SeedPoint2D* lhs, SeedPoint2D* rhs)
{
  return (lhs->getIntensity() < rhs->getIntensity());
}

