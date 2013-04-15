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

#include <vtksys/hash_map.hxx>
#include "tinyxml/tinyxml.h"

#define VTK_CREATE(type, var) \
  vtkSmartPointer<type> var = vtkSmartPointer<type>::New()

#define MY_ENCODING "ISO-8859-1"

#include "TraceLine.h"

///////////////////////////////////////////////////////////////////////////////
TraceLine::TraceLine()
{
	this->m_parent.clear();
	this->root= -1;
	this->level = 0;
	this->m_id = -(1<<30);
	this->m_branches.clear();
	this->EuclideanD = -1;
	this->length = 0;
	this->radii = 0; 
	this->sectionArea = 0;
	this->surfaceArea = 0;
	this->volume = 0;
	this->somaSurfaceArea = 0;
	this->somaVolume = 0;
	this->BurkTaper = 0;
	this->HillmanTaper = 0;
	this->HillmanThreshold = 0;
	this->FileName = "default";

	this->daughterRatio = 0;
	this->parentDaughterRatio = 0;
	this->daughterLengthRatio = 0;
	this->partitionAsymmetry = 0;

	this->rallPower = 0;
	this->Pk = 0;
	this->Pk_2 = 0;
	this->Pk_classic = 0;

	this->BifAmplLocal = -1;
	this->BifAmplRemote = -1;
	this->BifTiltLocalAvg = -1;
	this->BifTiltRemoteAvg = -1;
	this->BifTorqueLocalAvg = -1;
	this->BifTorqueRemoteAvg = -1;

	this->BifTiltLocal.push_back(-1);
	this->BifTiltLocal.push_back(-1);
	this->BifTiltRemote.push_back(-1);
	this->BifTiltRemote.push_back(-1);

	this->BifTorqueLocal = -1;
	this->planeAngleLocal.push_back(-1);
	this->planeAngleLocal.push_back(-1);
	this->planeAngleRemote.push_back(-1);
	this->planeAngleRemote.push_back(-1);

	this->CellFeatures.clear();
//for Dist to device calculations
	this->DistanceToROI = 0;
	this->ROICoord_X = 0;
	this->ROICoord_Y = 0;
	this->ROICoord_Z = 0;
	this->ROIAzimuth = -1;
	this->ROIElevation = -1;
	this->somaROIAngle = -1;

	this->actualBifurcation = false;
	this->TraceFeatures.clear();
	modified = true;

	//// classification results
	//this->prediction = -1;
	//this->confidence = -1;
}

///////////////////////////////////////////////////////////////////////////////
TraceLine::~TraceLine()
{
}

///////////////////////////////////////////////////////////////////////////////
std::vector<TraceLine*> TraceLine::GetParents()
{
	return this->m_parent;
}

TraceLine* TraceLine::GetParent(int i)
{
	return this->m_parent.at(i);
}

int TraceLine::ParentSize()
{
	return this->m_parent.size();
}

unsigned int TraceLine::GetParentID(int i)
{
	if (!this->isParentLess())
	{
		return this->m_parent.at(i)->GetId();
	}
	else 
	{
		return -1;
	}
}

bool TraceLine::isParentLess()
{
	return (this->m_parent.empty()||this->m_parent.at(0)==NULL);
}

bool TraceLine::isMarked()
{
	return this->marked;
}

void TraceLine::MarkLine()
{
	this->marked = true;
}

void TraceLine::UnmarkLine()
{
	this->marked = false;
}

void TraceLine::RemoveParents()
{
	this->m_parent.clear();
	this->SetParent(NULL);
}

void TraceLine::RemoveParent(int i)
{
	this->m_parent.erase(m_parent.begin()+i);
}

int TraceLine::GetParentNumber(TraceLine* p)
{
	int number = 0;
	for(int i = 0; i < this->ParentSize(); i++)
	{
		if(m_parent.at(i) == p)
		{
			number = i;
			break;
		}
	}
	return number;
}

void TraceLine::SetParent(TraceLine* p)
{
	if(this->ParentSize()>0 && this->GetParent(this->ParentSize()-1)==NULL)
	{
		this->RemoveParent(this->ParentSize()-1);
		this->m_parent.push_back(p);
	}
	else
		this->m_parent.push_back(p);
}
///////////////////////////////////////////////////////////////////////////////
int TraceLine::GetRootID()
{
	return this->root;
}
int TraceLine::GetLevel()
{
	return this->level;
}
void TraceLine::setTerminalDegree(int degree)
{
	this->terminalDegree = degree;
}
bool TraceLine::isLeaf()
{
	if (this->m_branches.size() ==0)
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool TraceLine::isRoot()
{
	if (this->isParentLess())
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool TraceLine::isFree()
{
	if (this->isRoot() && this->isLeaf())
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool TraceLine::isBranch()
{
	if (this->m_branches.size() > 1)
	{
		return true;
	}
	else
	{
		return false;
	}
}
void TraceLine::setActualBifurcation(bool bifurcate)
{
	this->actualBifurcation = bifurcate;
}
bool TraceLine::isActualBifurcation()
{
	return this->actualBifurcation;
}
void TraceLine::setRoot(int RootID, int traceLevel, double parentPath)
{
	this->root = RootID;
	this->level = traceLevel;
	this->PathLength = parentPath + this->GetLength();
}
void TraceLine::setRoot(int RootID)
{
	this->root = RootID;
}
void TraceLine::calculateVol()
{
	if (this->m_trace_bits.size() >1)
	{
		double dist = 0, r = 0, Df = 0, Dl = 0;
		TraceBit pre, cur;
		TraceBitsType::iterator it = this->m_trace_bits.begin();
		pre = *it; 
		it++;
		Df = 2*this->m_trace_bits.front().r;
		Dl = 2*this->m_trace_bits.back().r;
		r = pre.r;
		for (; it != this->m_trace_bits.end(); it++)
		{
			cur = *it;
			r += cur.r;
			dist += Euclidean(pre, cur);
			pre = cur;
		}
		if (Df !=Dl)
		{
			this->BurkTaper = (Df - Dl)/dist;
			this->HillmanTaper = (Df - Dl)/Df;
		}
		else
		{
			this->BurkTaper = 0;
			this->HillmanTaper = 0;
		}
		this->length = dist;
		this->radii = r / this->m_trace_bits.size(); //ave radii of segment
		this->sectionArea = PI*pow((this->radii),2);
		this->volume = this->sectionArea*this->length;
		this->surfaceArea = 2*this->radii*PI*this->length;
	}//end size >1
	else
	{
		this->length = 0;
		this->radii = this->m_trace_bits.front().r; 
		this->sectionArea = PI*pow((this->radii),2);
		this->surfaceArea = 0;
		this->volume = 0;
		this->somaSurfaceArea = 4*this->sectionArea;
		this->somaVolume = (4.0/3.0)*PI*pow((this->radii),3);
		this->BurkTaper = 0;
		this->HillmanTaper = 0;
	}//end else
}//end vol calculation
void TraceLine::calculateBifFeatures()
{
	this->actualBifurcation = true;
	TraceBit BranchBit= this->m_trace_bits.back();
	TraceBit previousBit = this->GetBitXFromEnd(2);
	double BranchBitRadii = BranchBit.r;
	TraceLine* Daughter1 = this->GetBranch1();
	TraceBit  D1F = Daughter1->GetTraceBitsPointer()->front();
	TraceBit  D1B = Daughter1->GetTraceBitsPointer()->back();

	int n1 = Daughter1->GetTerminalDegree();
	double D1Radii = D1F.r;
	TraceLine* Daughter2 = this->GetBranch2();
	TraceBit  D2F = Daughter2->GetTraceBitsPointer()->front();
	TraceBit  D2B = Daughter2->GetTraceBitsPointer()->back();

	int n2 = Daughter2->GetTerminalDegree();
	double D2Radii = D2F.r;
	if (D1Radii > D2Radii)
	{
		this->daughterRatio = D1Radii /D2Radii;
	}
	else
	{
		this->daughterRatio = D2Radii / D1Radii;
	}
	
	Daughter1->SetParentDaughterRatio( D1Radii / BranchBitRadii);
	Daughter2->SetParentDaughterRatio( D2Radii / BranchBitRadii);
	this->HillmanThreshold = .5*BranchBitRadii + .25*D1Radii +.25*D2Radii;

	this->rallPower = this->RallPower(BranchBitRadii, D1Radii, D2Radii);
	//this->rallPower = pow ((1 + pow (ratio, power)), -1/power);	//need to solve ratio and power
	if (this->rallPower != -1)
	{
		this->Pk = this->CalculatePk(BranchBitRadii, D1Radii, D2Radii, this->rallPower);
	}
	this->Pk_classic = this->CalculatePk(BranchBitRadii, D1Radii, D2Radii, 1.5);
	this->Pk_2 = this->CalculatePk(BranchBitRadii, D1Radii, D2Radii, 2);

	if ((n1+n2) > 2)
	{
		this->partitionAsymmetry = abs(n1 - n2) / (n1 + n2 - 2);
	}else
	{
		this->partitionAsymmetry = 1;
	}

	double daughter1pathLength = this->GetBranch1()->GetLength();
	double daughter2pathLength = this->GetBranch2()->GetLength();
	if (daughter1pathLength <= daughter2pathLength)
	{
		if (daughter2pathLength != 0)
			this->daughterLengthRatio = daughter1pathLength/daughter2pathLength;
	}
	else
	{
		if (daughter1pathLength != 0)
			this->daughterLengthRatio = daughter2pathLength/daughter1pathLength;
	}

	this->BifAmplLocal = this->Angle(D1F, BranchBit, D2F);
	this->BifAmplRemote = this->Angle(D1B, BranchBit, D2B);

	//Identify bigger vs smaller daughter angle
	double tiltAngle1 = this->Angle(previousBit, BranchBit, D1F);
	double tiltAngle2 = this->Angle(previousBit, BranchBit, D2F);
	this->setBifTiltLocalAvg((tiltAngle1+tiltAngle2)/2);

	this->BifTiltLocal[0] = tiltAngle1;
	this->BifTiltLocal[1] = tiltAngle2;

	tiltAngle1 = this->Angle(previousBit, BranchBit, D1B);
	tiltAngle2 = this->Angle(previousBit, BranchBit, D2B);
	this->setBifTiltRemoteAvg((tiltAngle1+tiltAngle2)/2);

	this->BifTiltRemote[0] = tiltAngle1;
	this->BifTiltRemote[1] = tiltAngle2;
	//Daughter1->setBifTiltLocal( this->Angle(previousBit, BranchBit, D1F));
	//Daughter2->setBifTiltLocal( this->Angle(previousBit, BranchBit, D2F));

	//Daughter1->setBifTiltRemote(this->Angle(previousBit, BranchBit, D1B));
	//Daughter2->setBifTiltRemote(this->Angle(previousBit, BranchBit, D2B));

    
    double ParentPlaneLocal[3]; this->Plane(D1F, BranchBit, D2F, ParentPlaneLocal);
    double ParentPlaneRemote[3]; this->Plane(D1B, BranchBit, D2B, ParentPlaneRemote);
    
	double Daughter1PlaneLocal[3];
	double Daughter1PlaneRemote[3];
	double Daughter2PlaneLocal[3];
	double Daughter2PlaneRemote[3];
	
	double planeAnglelocal1 = -1;
	double planeAnglelocal2 = -1;
	double planeAngleremote1 = -1;
	double planeAngleremote2 = -1;

	if (!Daughter1->isLeaf() && Daughter1->GetBranchPointer()->size() == 2) //Check if Daughter1 is not trifurcation
	{
		TraceLine* GrandDaughter1of1 = Daughter1->GetBranch1();
		TraceLine* GrandDaughter2of1 = Daughter1->GetBranch2();
		TraceBit GD1F1 = GrandDaughter1of1->GetTraceBitsPointer()->front();
		TraceBit GD1B1 = GrandDaughter1of1->GetTraceBitsPointer()->back();
		TraceBit GD2F1 = GrandDaughter2of1->GetTraceBitsPointer()->front();
		TraceBit GD2B1 = GrandDaughter2of1->GetTraceBitsPointer()->back();
		this->Plane(GD1F1,D1B,GD2F1, Daughter1PlaneLocal);
		this->Plane(GD1B1,D1B,GD2B1, Daughter1PlaneRemote);
		planeAnglelocal1 = this->PlaneAngle(ParentPlaneLocal, Daughter1PlaneLocal);
		planeAngleremote1 = this->PlaneAngle(ParentPlaneRemote,Daughter1PlaneRemote);
	}
	if (!Daughter2->isLeaf() && Daughter2->GetBranchPointer()->size() == 2)
	{
		TraceLine* GrandDaughter1of2 = Daughter2->GetBranch1();
		TraceLine* GrandDaughter2of2 = Daughter2->GetBranch2();
		TraceBit GD1F2 = GrandDaughter1of2->GetTraceBitsPointer()->front();
		TraceBit GD1B2 = GrandDaughter1of2->GetTraceBitsPointer()->back();
		TraceBit GD2F2 = GrandDaughter2of2->GetTraceBitsPointer()->front();
		TraceBit GD2B2 = GrandDaughter2of2->GetTraceBitsPointer()->back();
        this->Plane(GD1F2,D2B,GD2F2, Daughter1PlaneRemote);
		this->Plane(GD1B2,D2B,GD2B2, Daughter2PlaneRemote);
		planeAnglelocal2 = this->PlaneAngle(ParentPlaneLocal, Daughter2PlaneLocal);
		planeAngleremote2 = this->PlaneAngle(ParentPlaneRemote,Daughter2PlaneRemote);
	}

	if (planeAnglelocal1 != planeAnglelocal1) //quickfix for invalid numbers
		planeAnglelocal1 = -1;
	if (planeAnglelocal2 != planeAnglelocal2)
		planeAnglelocal2 = -1;
	if (planeAngleremote1 != planeAngleremote1)
		planeAngleremote1 = -1;
	if (planeAngleremote2 != planeAngleremote2)
		planeAngleremote2 = -1;

	planeAngleLocal[0] = planeAnglelocal1;
	planeAngleLocal[1] = planeAnglelocal2;
	planeAngleRemote[0] = planeAngleremote1;
	planeAngleRemote[1] = planeAngleremote2;
	
	if (planeAnglelocal1 != -1 || planeAnglelocal2 != -1)
	{
		if (planeAnglelocal1 != -1 && planeAnglelocal2 != -1)
		{
			this->setBifTorqueLocalAvg((planeAnglelocal1+planeAnglelocal2)/2);
		}
		else if (planeAnglelocal1 != -1)
		{
			this->setBifTorqueLocalAvg(planeAnglelocal1);
		}
		else
			this->setBifTorqueLocalAvg(planeAnglelocal2);
	}
	if (planeAngleremote1 != -1 || planeAngleremote2 != -1)
	{
		if (planeAngleremote1 != -1 && planeAngleremote2 != -1)
		{
			this->setBifTorqueRemoteAvg((planeAngleremote1+planeAngleremote2)/2);
		}
		else if (planeAngleremote1 != -1)
		{
			this->setBifTorqueRemoteAvg(planeAngleremote1);
		}
		else
			this->setBifTorqueRemoteAvg(planeAngleremote2);
	}

}
void TraceLine::setTraceBitIntensities(vtkSmartPointer<vtkImageData> imageData, std::string ImageName)
{
	TraceBit curBit;
	TraceBitsType::iterator it = this->m_trace_bits.begin();
	double totalIntensity = 0; 
	for (; it != this->m_trace_bits.end(); it++)
	{
		int lx = 0, ly = 0, lz = 0; 
		curBit = *it;
		lx = (int) floor(curBit.x + 0.5);
		ly = (int) floor(curBit.y + 0.5);
		lz = (int) floor(curBit.z + 0.5);
		(*it).I = imageData->GetScalarComponentAsDouble(lx,ly,lz,0);
		totalIntensity += (*it).I;
		//std::cout<< "\nid:\t"<< (*it).id << "\tI:\t" << (*it).I; //for checking data
	}
	vtkVariant aveIntensity = vtkVariant(totalIntensity/this->m_trace_bits.size());
	this->modified = true;
	this->SetTraceFeature(ImageName, aveIntensity);
	//should have ended
}
void TraceLine::setTraceBitWeightedIntensities(ImageType::Pointer input_image, std::string ImageName)
{
	/*!
	 * a circle kernel is used to evaluate the voxel intensities along the centerline
	 * @author Audrey Cheong
	 * @param input_image itk image
	 * @param ImageName feature header
	 */

	int totalIntensity = 0;
	if (this->m_trace_bits.size()>1)
	{
		//FILE*fp = fopen("C:/TestDelete/intensity.txt","a");
		//fprintf(fp,"\nRootID	SegmentID	BitID	Pathlength	Intensity\n");
		double dist = 0;

		TraceBit curBit, nextBit;
		TraceBitsType::iterator it = this->m_trace_bits.begin();
		curBit = *it;
		it++;
		for (; it != this->m_trace_bits.end(); it++)
		{
			int firstBit[3];
			firstBit[0] = (int) floor(curBit.x + 0.5);
			firstBit[1] = (int) floor(curBit.y + 0.5);
			firstBit[2] = (int) floor(curBit.z + 0.5);
			//int radius1 = (int) floor(curBit.r + 0.5);
			double radius1 = curBit.r;
			double new_radius = radius1;

			//std::cout << "First bit: " << curBit.x << "," << curBit.y << "," << curBit.z << std::endl;
			//std::cout << "First bit: " << firstBit[0] << "," << firstBit[1] << "," << firstBit[2] << std::endl;

			int secondBit[3];
			nextBit = *it;
			secondBit[0] = (int) floor(nextBit.x + 0.5);
			secondBit[1] = (int) floor(nextBit.y + 0.5);
			secondBit[2] = (int) floor(nextBit.z + 0.5);
			//int radius2 = (int) floor(nextBit.r + 0.5);
			double radius2 = nextBit.r;

			//std::cout << "Second bit: " << nextBit.x << "," << nextBit.y << "," << nextBit.z << std::endl;
			//std::cout << "First bit: " << secondBit[0] << "," << secondBit[1] << "," << secondBit[2] << std::endl;

			double distance[3]; //actually an integer but pow needs double
			int imageCenter[3];
			for (int i = 0; i < 3; i++)
			{
				distance[i] = secondBit[i] - firstBit[i];
				imageCenter[i] = firstBit[i] + distance[i]/2;
			}
			int radius = (int) floor(new_radius+0.5);

			dist += Euclidean(curBit, nextBit);
			//int length = (int) floor(sqrt(pow(distance[0],2)+pow(distance[1],2)+pow(distance[2],2))+0.5);
			//std::cout << "Length: " << length << std::endl;

			double azimuth = AzimuthAngle(curBit,nextBit)/180*PI;
			double elevation = ElevationAngle(curBit,nextBit)/180*PI;

			//number of voxels along the length
			int num_of_voxels = (int) distance[0];
			for (int i = 1; i < 3; i++)
			{
				if (distance[i] > num_of_voxels)
					num_of_voxels = (int) distance[i];
			}

			if (num_of_voxels == 0) //do not evaluate if the 1st and 2nd coordinates are the same
			{
				continue;
			}

			double radius_increment = (double) (radius2 - radius1)/num_of_voxels;
			double line_increment[3];
			int center_voxel[3] = { radius, radius, radius };
			for (int i = 0; i < 3; i++)
			{
				line_increment[i] = distance[i]/num_of_voxels;
			}

			//image size
			ImageType::RegionType input_volume_region = input_image->GetLargestPossibleRegion().GetSize();
			ImageType::SizeType input_size = input_volume_region.GetSize();

			//crop image
			ImageType::SizeType size;
			int boxSide = radius*2+1;
			size[0] = boxSide;
			size[1] = boxSide;
			size[2] = boxSide;

			//std::cout << "Image size: " << size[0] << " " <<size[1] << " " << size[2] <<std::endl;

			//fix for out of bounds value
			ImageType::IndexType start_index;
			start_index[0] = firstBit[0]-radius;
			start_index[1] = firstBit[1]-radius;
			start_index[2] = firstBit[2]-radius;

			double moving_index[3];
			moving_index[0] = curBit.x-curBit.r;
			moving_index[1] = curBit.y-curBit.r;
			moving_index[2] = curBit.z-curBit.r;

			int localIntensity = 0;
			for (int i = 0; i < num_of_voxels; i++)
			{
				//boundary conditions
				for (int j = 0; j < 3; j++)
				{
					//if start_index is less than 0, the center voxel needs to shift
					if (start_index[j] < 0)
					{
						size[j] += start_index[j];
						center_voxel[j] += start_index[j];
						start_index[j] = 0;
					}
					int max_size = input_size.GetSize()[j];
					int cur_size = start_index[j]+boxSide;
					if (cur_size > max_size)
					{
						size[j] = max_size - start_index[j];
					}
				}

				//std::cout << "Start index: " << start_index[0] << " " <<start_index[1] << " " << start_index[2] <<std::endl;

				ImageType::RegionType volume_region;
				volume_region.SetSize( size );
				volume_region.SetIndex( start_index );

				VolumeOfInterestFilterType::Pointer volumeOfInterestFilter = VolumeOfInterestFilterType::New();
				volumeOfInterestFilter->SetInput( input_image );
				volumeOfInterestFilter->SetRegionOfInterest( volume_region );
				volumeOfInterestFilter->Update();

				ImageType::Pointer crop_input_image = volumeOfInterestFilter->GetOutput();

				ImageType::Pointer mask = ImageType::New();

				StructuredObject * circleMask = new StructuredObject();
				circleMask->circleKernel(crop_input_image,mask,center_voxel,radius,azimuth,elevation);

				MaskFilterType::Pointer maskFilter = MaskFilterType::New();
				maskFilter->SetInput(crop_input_image);
				maskFilter->SetMaskImage(mask);

				//mask->Print(std::cout);

				ImageType::Pointer crop_output_image = maskFilter->GetOutput();
				crop_output_image->Update();


				//Validate accuracy
			/*	WriterType::Pointer writer1 = WriterType::New();
				writer1->SetFileName( "Input.tif" );
				writer1->SetInput( crop_input_image );*/

				//WriterType::Pointer writer2 = WriterType::New();
				//writer2->SetFileName( "Output.tif" );
				//writer2->SetInput( crop_output_image );
				
				//try
				//{
				////	//writer1->Update();
				////	//writer2->Update();
				//}
				//catch( itk::ExceptionObject & excp )
				//{
				//	std::cerr << excp << std::endl;
				//}
				//scanf("%*d");

				ConstIteratorType imageIterator( crop_output_image,crop_output_image->GetRequestedRegion() );
				imageIterator.GoToBegin();
				while( !imageIterator.IsAtEnd() ) //imageIterator not working? sometimes - only work for 8 bit images
				{
					// Get value of current pixel
					unsigned char pixel_value = imageIterator.Get();
					totalIntensity += (int) pixel_value;
					localIntensity += (int) pixel_value;
					//std::cout << "Pixel intensity: " << (int) pixel_value << std::endl;
					++imageIterator;

				}//end imageIterator
				//std::cout << "Total intensity: " << totalIntensity << std::endl;

				//next image settings
				new_radius += radius_increment;
				radius = (int) floor( new_radius + 0.5 );
				boxSide = radius*2+1;
				for (int i = 0; i < 3; i++)
				{
					moving_index[i] += line_increment[i];
					start_index[i] = (int) floor( moving_index[i] + 0.5 );
					size[i] = boxSide;
				}
			}//endfor num_of_voxels
				
			//fprintf(fp,"%d\t%d\t%d\t%d\t%f\t%d\n", GetRootID(), GetId(), GetLevel(), curBit.id, dist, localIntensity);
			
			curBit = nextBit;
		}//endfor m_trace_bits
		//fclose(fp);
	}//endif size
	vtkVariant Intensity = vtkVariant(totalIntensity);
	this->modified = true;
	this->SetTraceFeature(ImageName, Intensity);

}
double TraceLine::GetEuclideanLength()
{
	if (this->m_trace_bits.size() <2)
	{
		this->EuclideanD =0;
	}
	else
	{
		TraceBit front = this->m_trace_bits.front();
		TraceBit back  = this->m_trace_bits.back();
		this->EuclideanD = this->Euclidean(front, back);
	}
	return this->EuclideanD;
}
double TraceLine::GetBitDensity()
{
	if (this->GetSize() >1)
	{
		this->BitDensity = this->GetSize() / this->GetLength();
	}
	else
	{
		this->BitDensity = 1;
	}
	return this->BitDensity;
}
double TraceLine::GetDistToParent()
{
	if (!this->isParentLess())
	{
		this->DistToParent = this->Euclidean(this->m_trace_bits.front(), //
			this->m_parent.at(0)->m_trace_bits.back());
		if (this->m_trace_bits.size()>1)
		{
			double Leading = 0;//, dist =0;
			TraceBit pre, cur;
			TraceBitsType::iterator it = this->m_trace_bits.begin();
			pre = *it; 
			it++;
			cur = *it;
			Leading = Euclidean(pre, cur);
			if (Leading > 2*this->DistToParent)
			{
				this->DistToParent = Leading;
			}
		}
		return this->DistToParent;
	}
	else
	{
		this->DistToParent = 0;
		return -1;
	}
}
double TraceLine::GetFragmentationSmoothness()
{
	if (!(this->EuclideanD > -1))
	{
		this->GetEuclideanLength();
	}
	double t = -1;
	if ( this->m_trace_bits.size() > 1)
	{
		t = this->length/this->EuclideanD;
	}
	return t;
}
///////////////////////////////////////////////////////////////////////////////
void TraceLine::AddBranch(TraceLine* b)
{ 
  this->m_branches.push_back(b);
}

///////////////////////////////////////////////////////////////////////////////
TraceLine* TraceLine::GetBranch1()
{
  return this->m_branches[0];
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::SetBranch1(TraceLine* b0)
{
  this->m_branches[0] = b0;
}

///////////////////////////////////////////////////////////////////////////////
TraceLine* TraceLine::GetBranch2()
{
  return this->m_branches[1];
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::SetBranch2(TraceLine* b1)
{
  this->m_branches[1] = b1;
}

///////////////////////////////////////////////////////////////////////////////
unsigned char TraceLine::GetType()
{
  return this->m_type;
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::SetType(unsigned char t)
{
  this->m_type = t;
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::AddTraceBit(TraceBit tbit)
{
  this->m_trace_bits.push_back(tbit);
}
void TraceLine::ExtendTrace(TraceBit tbit)
{
	bool front = false;
	if (this->m_trace_bits.size()>1)
	{
		if (!this->Orient(tbit))
		{
			this->m_trace_bits.push_front(tbit);
			front = true;
		}
	}
	if (!front)
	{
		this->m_trace_bits.push_back(tbit);
	}
}
TraceBit TraceLine::removeLastBit()
{
	TraceBit lastBit = this->m_trace_bits.back();
	if (this->m_trace_bits.size() > 1)
	{
		this->m_trace_bits.pop_back();
	}
	return lastBit;
}
TraceBit TraceLine::removeFirstBit()
{
	TraceBit firstBit = this->m_trace_bits.front();
	if (this->m_trace_bits.size() > 1)
	{
		this->m_trace_bits.pop_front();
	}
	return firstBit;
}

//This is an error fix
bool TraceLine::removeLeadingBit()
{
	if (this->m_trace_bits.size() < 3)
	{
		return false;
	}else
	{		
		//std::vector<double> dist;// = 0, r = 0;
		double Leading = 0, dist =0;
		TraceBit pre, cur;
		TraceBitsType::iterator it = this->m_trace_bits.begin();
		pre = *it; 
		it++;
		cur = *it;
		Leading = Euclidean(pre, cur);
		it++;
		for (; it != this->m_trace_bits.end(); it++)
		{
			cur = *it;
			dist += this->Euclidean(pre, cur);
			pre = cur;
		}//end for
		if(Leading > (dist/ (this->m_trace_bits.size() -1)))
		{
			this->m_trace_bits.pop_front();
			modified = true;
			return true;
		}//end dist
		return false;
	}//end else size
}
///////////////////////////////////////////////////////////////////////////////
TraceLine::TraceBitsType::iterator TraceLine::GetTraceBitIteratorBegin()
{
  return this->m_trace_bits.begin();
}
TraceBit TraceLine::GetBitXFromBegin(int x)
{
	TraceBit curBit;
	TraceBitsType::iterator iter = this->m_trace_bits.begin();
	TraceBitsType::iterator iterend = this->m_trace_bits.end();
	if (3 > (int)this->m_trace_bits.size())
	{
		--iterend;
		curBit = * iterend;
		//return  curBit;
	}
	else
	{
		int i = 0;
		while(( i< (int)(this->m_trace_bits.size()-2))&&( i < x ))
		{
			++iter;
			i++;
		}
		curBit = *iter;
	}
	return curBit;
}
///////////////////////////////////////////////////////////////////////////////
TraceLine::TraceBitsType::iterator TraceLine::GetTraceBitIteratorEnd()
{
  return this->m_trace_bits.end();
}

TraceBit TraceLine::GetBitXFromEnd(int x)
{
	TraceBit curBit;
	TraceBitsType::iterator iter = this->m_trace_bits.begin();
	TraceBitsType::iterator iterend = this->m_trace_bits.end();
	if (3 >(int)this->m_trace_bits.size())
	{
		curBit = * iter;
		
	}
	else
	{
		int i = 0;
		while(( i<(int)(this->m_trace_bits.size()-2))&&( i < x ))
		{
			--iterend;
			i++;
		}
		curBit = *iterend;
	}
	return curBit;
}
///////////////////////////////////////////////////////////////////////////////
TraceLine::TraceBitsType* TraceLine::GetTraceBitsPointer()
{
  return &m_trace_bits;
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::SetId(unsigned int lid){ m_id = lid;}

///////////////////////////////////////////////////////////////////////////////
unsigned int TraceLine::GetId()
{ 
  return m_id;
};

///////////////////////////////////////////////////////////////////////////////
int TraceLine::GetSize()
{
  return (int)m_trace_bits.size();
};

///////////////////////////////////////////////////////////////////////////////
void TraceLine::Print(std::ostream &c,int indent)
{
  for(int counter=0; counter<indent; counter++)
	  c<<" ";
  c<<"TraceLine: "<<std::endl;
  for(int counter=0; counter<indent; counter++)
	  c<<" ";
  c<<"Size: "<<m_trace_bits.size()<<std::endl;
  for(unsigned int counter=0; counter< m_branches.size(); counter++)
  {
	  m_branches[counter]->Print(std::cout,indent+4);
  }
}

///////////////////////////////////////////////////////////////////////////////
std::vector<unsigned int> * TraceLine::GetMarkers()
{
  return &(this->m_markers);
}

///////////////////////////////////////////////////////////////////////////////
std::vector<TraceLine*> * TraceLine::GetBranchPointer()
{
  return &(this->m_branches);
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::setTraceColor(double newColor)
{
  this->traceColor = newColor;
};

///////////////////////////////////////////////////////////////////////////////
double TraceLine::getTraceColor()
{
  return this->traceColor;
};

///////////////////////////////////////////////////////////////////////////////
TraceLine::TraceLine(const TraceLine &t)
{
  this->traceColor = t.traceColor;
  this->m_id = t.m_id;
  this->m_markers = t.m_markers;
  this->m_type = t.m_type;
  this->m_parent.clear();
  modified = true;
  for(unsigned int counter=0; counter< t.m_branches.size(); counter++)
  {
    TraceLine *temp = new TraceLine();
    temp->SetParent(this);
    this->AddBranch(temp);
    *temp = *t.m_branches[counter];
  }
  this->m_trace_bits = t.m_trace_bits;
  this->FileName = t.FileName;
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::Getstats()
{
  double XF, XB, YF, YB, ZF, ZB;
  XF = m_trace_bits.front().x;
  XB = m_trace_bits.back().x;
  YF= m_trace_bits.front().y;
  YB= m_trace_bits.back().y;
  ZF= m_trace_bits.front().z;
  ZB= m_trace_bits.back().z;
  printf("Trace # %d \t Trace Size: \t %d ", m_id, (int)m_trace_bits.size());
  printf("First bit x: %4.2f y: %4.2f z: %4.2f \t", XF, YF, ZF); 
  printf("Endt bit x: %4.2f y: %4.2f z: %4.2f \n", XB, YB, ZB); 
}


vtkVariant TraceLine::GetTraceFeature(std::string FeatureName)
{
	/*!
	* return -pi if feature not computed
	* segment level features
	*/
	std::map<std::string ,vtkVariant>::iterator find = this->TraceFeatures.find(FeatureName);
	if (find != this->TraceFeatures.end())
	{
		return (*find).second;
	}
	else
	{
		vtkVariant value = -PI;
		return value;
	}
}

void TraceLine::SetTraceFeature(std::string FeatureName,vtkVariant FeatureValue)
{
	/*!
	* Set or update segment level features 
	* creates if not computed
	*/

	std::map<std::string ,vtkVariant>::iterator find = this->TraceFeatures.find(FeatureName);
	if (find != this->TraceFeatures.end())
	{
		(*find).second = FeatureValue;
	}
	else
	{
		this->TraceFeatures[FeatureName] = FeatureValue;
	}

}

vtkVariant TraceLine::GetCellFeature(std::string FeatureName)
{
	/*!
	* return -pi if feature not computed
	* cell level features
	*/
	std::map<std::string ,vtkVariant>::iterator find = this->CellFeatures.find(FeatureName);
	if (find != this->CellFeatures.end())
	{
		return (*find).second;
	}
	else
	{
		vtkVariant value = -PI;
		return value;
	}
}

void TraceLine::SetCellFeature(std::string FeatureName,vtkVariant FeatureValue)
{
	/*!
	* Set or update cell level features 
	* creates if not computed
	*/
	std::map<std::string ,vtkVariant>::iterator find = this->CellFeatures.find(FeatureName);
	if (find != this->CellFeatures.end())
	{
		(*find).second = FeatureValue;
	}
	else
	{
		this->CellFeatures[FeatureName] = FeatureValue;
	}

}

///////////////////////////////////////////////////////////////////////////////
double TraceLine::Euclidean(TraceBit bit1, TraceBit bit2)
{
	double distance, x, y, z;
	x = pow((bit1.x -bit2.x),2);
	y = pow((bit1.y -bit2.y),2);
	z = pow((bit1.z -bit2.z),2);
	distance = sqrt(x +y +z);
	return distance;

}
double TraceLine::Angle(TraceBit bit1f, TraceBit bit1b, TraceBit bit2f, TraceBit bit2b)
{
	double delX1, delX2, delY1, delY2, delZ1, delZ2,  norm1, norm2, angle;
	//delta x,y,z 
	angle = 0;
	delX1= fabs (bit1f.x-bit1b.x);	delX2= fabs (bit2f.x-bit2b.x);
	delY1= fabs (bit1f.y-bit1b.y);	delY2= fabs (bit2f.y-bit2b.y);
	delZ1= fabs (bit1f.z-bit1b.z);	delZ2= fabs (bit2f.z-bit2b.z);
	norm1=sqrt(pow((delX1),2)+ pow((delY1),2)+ pow((delZ1),2));
	norm2=sqrt(pow((delX2),2)+ pow((delY2),2)+ pow((delZ2),2));
	angle = acos(((delX1 * delX2) + (delY1 *delY2) + (delZ1 * delZ2)) /
          (norm2 * norm1));
	if (!(angle >= 0))
	{//prevent nan rounding errors
		angle = 0;
	}
	return angle;
}
double TraceLine::Angle(TraceBit bit1, TraceBit vertex, TraceBit bit2)
{
	double delX1, delX2, delY1, delY2, delZ1, delZ2,  norm1, norm2;
	double NewAngle = 0;
	//delta x,y,z 	
	delX1= fabs (vertex.x-bit1.x);	delX2= fabs (vertex.x-bit2.x);
	delY1= fabs (vertex.y-bit1.y);	delY2= fabs (vertex.y-bit2.y);
	delZ1= fabs (vertex.z-bit1.z);	delZ2= fabs (vertex.z-bit2.z);
	norm1=sqrt(pow((delX1),2)+ pow((delY1),2)+ pow((delZ1),2));
	norm2=sqrt(pow((delX2),2)+ pow((delY2),2)+ pow((delZ2),2));
	NewAngle = acos(((delX1 * delX2) + (delY1 *delY2) + (delZ1 * delZ2)) /
          (norm2 * norm1));
	if (!(NewAngle >= 0))
	{//prevent nan rounding errors
		NewAngle = 0;
	}
	return (NewAngle *180 )/PI;
}
double TraceLine::AzimuthAngle(TraceBit vertex, TraceBit bit1)
{
	double delX, delY;
	double NewAngle = 0;
	//delta x,y
	delX = bit1.x - vertex.x;
	delY = bit1.y - vertex.y;

	NewAngle = atan2(delY,delX);
	//gives positive angle if counterclockwise and negative angle if clockwise

	return (NewAngle *180 )/PI;
}
double TraceLine::ElevationAngle(TraceBit vertex, TraceBit bit1)
{
	double delX, delY, delZ;
	double NewAngle = 0;
	//delta x,y,z
	delX = bit1.x - vertex.x;
	delY = bit1.y - vertex.y;
	delZ = bit1.z - vertex.z;
	
	double hypotenuse = sqrt(pow(delX,2) + pow(delY,2));
	NewAngle = atan2(delZ,hypotenuse);
	//gives positive angle if counterclockwise and negative angle if clockwise

	return (NewAngle *180 )/PI;
}
double TraceLine::RallPower(double diamParent, double diamD1, double diamD2)
{
	double m = .001;
	double min = 1000; //min should be less than max
	while (m < 5)
	{
		double a1 = pow(diamParent, m);
		double a2 = pow(diamD1, m) + pow(diamD2, m);
		double a3 = fabs(a1 - a2);
		if (a3 <= .001)
		{
			return m;
		}
		m += .001;
	}
	return -1;
}
bool TraceLine::EndPtDist(TraceLine *Trace2, int &dir1, int &dir2, double &dist,
                          double &maxdist, double &angle) 
{
  //int center1=this->GetSize()/2, center2=Trace2->GetSize()/2;
  double min, distances[4];

  int Xbits = 10; //small step from end to determine 
  //compute the endpt distances
  distances[0]=Euclidean(m_trace_bits.front(), Trace2->m_trace_bits.front());//0 F-F
  distances[1]=Euclidean(m_trace_bits.front(), Trace2->m_trace_bits.back());//1 F-B
  distances[2]=Euclidean(m_trace_bits.back(), Trace2->m_trace_bits.front());//2 B-F
  distances[3]=Euclidean(m_trace_bits.back(), Trace2->m_trace_bits.back());//3 B-B

  //determine minimum spacing
  min = distances[0];
  int i, mark=0;
  for (i = 1; i<4; i++)
    {
    if (min > distances[i])
      {
      min= distances[i];
      mark=i;
      }
    } 
// std::cout << mark << "=mark\n";
  // from min determine orientation and return distance
  if (mark ==0)
    {//F-F
		//std::cout <<"FF\n";
	if (this->m_branches.size() > 0 && Trace2->m_branches.size()> 0)
		{return false;}
	if(this->GetParentID(0)!=-1 || Trace2->GetParentID(0)!=-1)
	{return false;}
    dist = distances[0];
    dir1= m_trace_bits.front().marker;
    dir2= Trace2->m_trace_bits.front().marker; 
	angle = Angle(m_trace_bits.front(),this->GetBitXFromBegin(Xbits), 
		Trace2->m_trace_bits.front(),Trace2->GetBitXFromBegin(Xbits));
    //angle=PI-angle;
    maxdist= distances[3];
    }
  else if (mark ==1)
    {//1 F-B
		//std::cout <<"FB\n";
	if ((this->GetParentID(0) != -1)||!Trace2->isLeaf())
		{return false;}
    dist = distances[1];
    dir1= m_trace_bits.front().marker;
    dir2= Trace2->m_trace_bits.back().marker; 
	angle = Angle(m_trace_bits.front(),this->GetBitXFromBegin(Xbits), 
		Trace2->GetBitXFromEnd(Xbits),Trace2->m_trace_bits.back());
    maxdist= distances[2];
    }
  else if (mark ==2)
    {//2 B-F
		//std::cout <<"BF\n";
		if (!this->isLeaf() || (Trace2->GetParentID(0)!=-1))
		{return false;}
    dist = distances[2];
    dir1= m_trace_bits.back().marker;
    dir2= Trace2->m_trace_bits.front().marker; 
	angle = Angle(this->GetBitXFromEnd(Xbits),m_trace_bits.back(), 
		Trace2->m_trace_bits.front(),Trace2->GetBitXFromBegin(Xbits));
    maxdist= distances[1];
    }
  else
    {//3 B-B
		//std::cout <<"BB\n";
	if (!this->isLeaf() || !Trace2->isLeaf()) 
		{return false;}
    dist = distances[3];
    dir1= m_trace_bits.back().marker;
    dir2= Trace2->m_trace_bits.back().marker;
	angle = Angle(this->GetBitXFromEnd(Xbits),m_trace_bits.back(), 
		Trace2->GetBitXFromEnd(Xbits),Trace2->m_trace_bits.back());
    //angle=PI-angle;
    maxdist= distances[0];
    }
  return true;
}

void TraceLine::EndPtDistVessel(TraceLine *Trace2, TraceBit &dir1, TraceBit &dir2, double &dist, double &maxdist, double &angle){

#ifdef USE_BALL_TRACER

  //int center1=this->GetSize()/2, center2=Trace2->GetSize()/2;
  double min, distances[4];

  int Xbits = 10; //small step from end to determine 
  //compute the endpt distances
  distances[0]=Euclidean(m_trace_bits.front(), Trace2->m_trace_bits.front());//0 F-F
  distances[1]=Euclidean(m_trace_bits.front(), Trace2->m_trace_bits.back());//1 F-B
  distances[2]=Euclidean(m_trace_bits.back(), Trace2->m_trace_bits.front());//2 B-F
  distances[3]=Euclidean(m_trace_bits.back(), Trace2->m_trace_bits.back());//3 B-B

  //determine minimum spacing
  min = distances[0];
  int i, mark=0;
  for (i = 1; i<4; i++)
  {
    if (min > distances[i])
    {
      min = distances[i];
      mark = i;
    }
  } 
  
  // std::cout << mark << "=mark\n";
  // from min determine orientation and return distance
  if (mark ==0)
    {//F-F
		//std::cout <<"FF\n";
	if (this->m_branches.size() > 0 && Trace2->m_branches.size()> 0)
		{//return false;
	}
	if(this->GetParentID(0)!=-1 || Trace2->GetParentID(0)!=-1)
	{//return false;
	}
    
	dist = distances[0];
    dir1 = m_trace_bits.front(); //.marker;
    dir2 = Trace2->m_trace_bits.front(); //.marker; 
	angle = Angle(m_trace_bits.front(),this->GetBitXFromBegin(Xbits), 
		Trace2->m_trace_bits.front(),Trace2->GetBitXFromBegin(Xbits));
    //angle=PI-angle;
    maxdist= distances[3];
    }
  else if (mark ==1)
    {//1 F-B
		//std::cout <<"FB\n";
	if ((this->GetParentID(0) != -1)||!Trace2->isLeaf())
	{//return false;
	}

    dist = distances[1];
    dir1 = m_trace_bits.front(); //.marker;
    dir2 = Trace2->m_trace_bits.back(); //.marker; 
	angle = Angle(m_trace_bits.front(),this->GetBitXFromBegin(Xbits), 
		Trace2->GetBitXFromEnd(Xbits),Trace2->m_trace_bits.back());
    maxdist= distances[2];
    }
  else if (mark ==2)
    {//2 B-F
		//std::cout <<"BF\n";
		if (!this->isLeaf() || (Trace2->GetParentID(0)!=-1))
		{//return false;
		}
    dist = distances[2];
    dir1 = m_trace_bits.back(); //.marker;
    dir2 = Trace2->m_trace_bits.front(); //.marker; 
	angle = Angle(this->GetBitXFromEnd(Xbits),m_trace_bits.back(), 
		Trace2->m_trace_bits.front(),Trace2->GetBitXFromBegin(Xbits));
    maxdist= distances[1];
    }
  else
    {//3 B-B
		//std::cout <<"BB\n";
	if (!this->isLeaf() || !Trace2->isLeaf()) 
	{//return false;
	}

    dist = distances[3];
    dir1 = m_trace_bits.back(); //.marker;
    dir2 = Trace2->m_trace_bits.back(); //.marker;
	angle = Angle(this->GetBitXFromEnd(Xbits),m_trace_bits.back(), 
		Trace2->GetBitXFromEnd(Xbits),Trace2->m_trace_bits.back());
    //angle=PI-angle;
    maxdist= distances[0];
    }
  //return true;

#endif 
}

double TraceLine::CalculatePk(double Dp, double Da, double Db, double n)
{
	return (pow(Da, n) + pow(Db, n))/pow(Dp, n);
}
bool TraceLine::Orient(TraceLine * Trunk)
{
	double distances[2];
	//compute the endpt distances
	distances[0]= Euclidean(m_trace_bits.front(),	Trunk->m_trace_bits.back());// F-B
	distances[1]= Euclidean(m_trace_bits.back(),	Trunk->m_trace_bits.back());// B-B
	if(distances[0] < distances[1])
	{
		return true;	//oriented correctly
	}
	return false;		//needs to be flipped
}
bool TraceLine::Orient(TraceBit bit)
{
	double distances[2];
	//compute the endpt distances
	distances[0]= Euclidean(m_trace_bits.front(),	bit);// F-Bit
	distances[1]= Euclidean(m_trace_bits.back(),	bit);// B-Bit
	if(distances[0] > distances[1])
	{
		return true;	//oriented correctly
	}
	return false;		//needs to be flipped
}
///////////////////////////////////////////////////////////////////////////////
void TraceLine::SetFileName(char *newFileName)
{
	this->FileName = newFileName;
}
const char * TraceLine::GetFileName()
{
	return this->FileName.c_str();
}
void TraceLine::getEndPtBounds(double bounds[])
{
	//x
	if ( this->m_trace_bits.front().x < this->m_trace_bits.back().x)
	{
		bounds[0] = this->m_trace_bits.front().x;
		bounds[1] = this->m_trace_bits.back().x;
	}
	else
	{
		bounds[1] = this->m_trace_bits.front().x;
		bounds[0] = this->m_trace_bits.back().x;
	}
	//y
	if ( this->m_trace_bits.front().y < this->m_trace_bits.back().y)
	{
		bounds[2] = this->m_trace_bits.front().y;
		bounds[3] = this->m_trace_bits.back().y;
	}
	else
	{
		bounds[3] = this->m_trace_bits.front().y;
		bounds[2] = this->m_trace_bits.back().y;
	}
	//z
	if ( this->m_trace_bits.front().z < this->m_trace_bits.back().z)
	{
		bounds[4] = this->m_trace_bits.front().z;
		bounds[5] = this->m_trace_bits.back().z;
	}
	else
	{
		bounds[5] = this->m_trace_bits.front().z;
		bounds[4] = this->m_trace_bits.back().z;
	}
	//debug
	for (int i = 0; i < 6; i++)
	{
		std::cout<< i <<" Value " << bounds[i];
	}
}
std::string TraceLine::stats()
{
	std::stringstream thisStats;
	thisStats << this->GetId();
	thisStats << "\t";
	thisStats << (int)this->GetType();
	thisStats << "\t" ;
	thisStats << this->GetSize();
	thisStats << "\t";
	thisStats << this->GetLength();
	thisStats << "\t" ;
	thisStats << this->GetEuclideanLength();
	thisStats << "\t" ;
	thisStats << this->GetRadii();
	thisStats << "\t";
	thisStats << this->GetFragmentationSmoothness() ;
	thisStats << "\t";
	for(int i = 0; i<this->ParentSize();i++)
	{
		thisStats << "\t";	
		thisStats << this->GetParentID(i);
		
	}
	return thisStats.str();
}
std::string TraceLine::RootCoord()
{
	std::stringstream thisStats;
	thisStats << this->m_trace_bits.front().x;
	thisStats << "_";
	thisStats << this->m_trace_bits.front().y;
	thisStats << "_";
	thisStats << this->m_trace_bits.front().z;
	return thisStats.str();

}
std::string TraceLine::statHeaders()
{
	std::stringstream thisStatsHeaders;
	thisStatsHeaders <<"ID"
		<<"\tType"
		<<"\tSize"
		<<"\tLength"
		<<"\tEuclidean L"
		<<"\tRadii"
		<<"\tContraction"
		<<"\tParent ID";
	return thisStatsHeaders.str();
}

void TraceLine::Plane(TraceBit bit1, TraceBit vertex, TraceBit bit2, double vector[])
{
	vector[0] = bit1.x - vertex.x;
	vector[1] = bit1.y - vertex.y;
	vector[2] = bit1.z - vertex.z;
}

double TraceLine::PlaneAngle(double *plane1, double *plane2)
{
	double top = plane1[0]*plane2[0] + plane1[1]*plane2[1] + plane1[2]*plane2[2];
	double bottom = sqrt(plane1[0]*plane1[0] + plane1[1]*plane1[1]+ plane1[2]*plane1[2])*sqrt(plane2[0]*plane2[0] + plane2[1]*plane2[1]+ plane2[2]*plane2[2]);
	return std::acos(top/bottom)*180/PI;
}

double TraceLine::GetCompartmentCurvature()
{
	double angle = 0;
	if (!this->m_parent.empty())
	{
		if (this->m_trace_bits.size()>2)
		{
			TraceBit pre, cur, next;
			TraceBitsType::iterator it = this->m_trace_bits.begin();
			pre = *it;
			it++;
			cur = *it;
			it++;
			next = *it;
			angle = this->Angle(pre,cur,next);
		}
		return angle;
	}
	else
		return -1;
}

double TraceLine::GetAzimuth()
{
	double newAzimuthAngle = 0;
	if (!this->m_parent.empty())
	{
		if (this->m_trace_bits.size()>1)
		{
			//TraceBit pre, cur;
			//TraceBitsType::iterator it = this->m_trace_bits.begin();
			//pre = *it;
			//it++;
			//cur = *it;
			//newAzimuthAngle = this->AzimuthAngle(pre,cur);

			TraceBit leadBit = this->m_trace_bits.front();
			TraceBit endBit = this->m_trace_bits.back();
			newAzimuthAngle = this->AzimuthAngle(leadBit,endBit);
		}
		return newAzimuthAngle;
	}
	else
	{
		return -1;
	}
}
double TraceLine::GetElevation()
{
	double newElevationAngle = 0;
	if (!this->m_parent.empty())
	{
		if (this->m_trace_bits.size()>1)
		{
			//TraceBit pre, cur;
			//TraceBitsType::iterator it = this->m_trace_bits.begin();
			//pre = *it;
			//it++;
			//cur = *it;
			//newElevationAngle = this->ElevationAngle(pre,cur);
			
			TraceBit leadBit = this->m_trace_bits.front();
			TraceBit endBit = this->m_trace_bits.back();
			newElevationAngle = this->ElevationAngle(leadBit,endBit);
		}
		return newElevationAngle;
	}
	else
	{
		return -1;
	}
}

double TraceLine::GetAngle(TraceBit bit1f, TraceBit bit1b, TraceBit bit2f, TraceBit bit2b){

	return Angle(bit1f, bit1b, bit2f, bit2b);
}

void TraceLine::CalculateDirectionToROI(TraceBit tipsPt)
{
	TraceBit somaPt = this->m_trace_bits.front();
	TraceBit ROIpt;
	ROIpt.x = this->ROICoord_X;
	ROIpt.y = this->ROICoord_Y;
	ROIpt.z = this->ROICoord_Z;

	this->ROIAzimuth = AzimuthAngle(somaPt,ROIpt);
	this->ROIElevation = ElevationAngle(somaPt,ROIpt);
	this->somaROIAngle = Angle(ROIpt,somaPt,tipsPt);
}