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

#include "TraceBit.h"
#include "TraceLine.h"
#include "CellTrace.h"
#include "SuperEllipsoid.h"

static double DefaultValue[115] ={-1,-1,-1,-1,-1,-1,-1,0,0,0,-1,-1,-1,-1,0,0,0,0,100,0,0,100,0,0,
0,0,1000,0,0,100,0,0,0,100,0,100,0,0,100,0,0,0,0,0,0,-PI,-PI,-PI,100,-PI,0,0,100,0,0,100,0,0,100,
0,0,100,0,0,100,-1,-1,100,-1,-1,100,0,0,100,100,0,-PI,180,0,-PI,180,0,0,180,0,0,180,0,0,180,0,0,180,
0,100,100,0,0,0,0,0,100,0,0,0,0,100,0,0,100,100,-PI,0,-1,-1};

CellTrace::CellTrace()
{
	this->clearAll();
	CellData = vtkSmartPointer<vtkVariantArray>::New();
}
CellTrace::CellTrace(std::vector<TraceLine*> Segments)
{
	this->clearAll();
	CellData = vtkSmartPointer<vtkVariantArray>::New();
	this->setTraces(Segments);
}
void CellTrace::setTraces(std::vector<TraceLine*> Segments)//& Segments)
{
	this->segments = Segments;
	unsigned int i = 0;
	this->NumSegments = (int) this->segments.size();
	this->stems = (int) this->segments[0]->GetBranchPointer()->size();
	/*this->prediction = this->segments[0]->getPrediction();
	this->confidence = this->segments[0]->getConfidence();*/
	if (this->stems > 0)
	{
		for ( int j = 0; j < this->stems; j++)
		{
			this->MaxMin(this->segments[0]->GetBranchPointer()->at(j)->GetDistToParent(), 
				this->TotalStemDistance, this->MinStemDistance, this->MaxStemDistance);
			if(this->segments[0]->GetBranchPointer()->at(j)->GetTerminalDegree() > 1)
				this->branchingStem++;

		}
		this->EstimatedSomaRadius = this->TotalStemDistance/(double)this->stems;
	}
	else
	{
		this->MinStemDistance = -PI;
		this->MaxStemDistance = -PI;
		this->TotalStemDistance = -PI;
		this->EstimatedSomaRadius = -PI;
	}
	TraceBit rootBit = this->segments[0]->GetTraceBitsPointer()->front();
	//set the soma point as the intital bounds
	this->somaX = rootBit.x;
	this->somaY = rootBit.y;
	this->somaZ = rootBit.z;
	this->maxX = rootBit.x;
	this->maxY = rootBit.y;
	this->maxZ = rootBit.z;
	this->minX = rootBit.x;
	this->minY = rootBit.y;
	this->minZ = rootBit.z;

	tips.push_back(rootBit);

	this->somaSurface = this->segments[0]->GetSurfaceArea();
	this->somaVolume = this->segments[0]->GetVolume();
	this->somaRadii = this->segments[0]->GetRadii();
	this->DeviceDistance = this->segments[0]->GetDistanceToROI();

	for(i = 0; i < this->segments.size(); i++)
	{
		this->IDs.insert(this->segments[i]->GetId());
		this->MaxMin(this->segments[i]->GetLength(), this->PathLengthTotal, this->PathLengthMin, this->PathLengthMax);
		this->MaxMin(this->segments[i]->GetVolume(), this->TotalVolume, this->SegmentVolumeMin, this->SegmentVolumeMax);
		this->MaxMin(this->segments[i]->GetEuclideanLength(),this->TotalEuclideanPath, this->MinEuclideanPath, this->MaxEuclideanPath);
		this->MaxMin(this->segments[i]->GetSurfaceArea(), this->surfaceAreaTotal, this->SurfaceAreaMin, this->SurfaceAreaMax);
		this->MaxMin(this->segments[i]->GetSectionArea(), this->sectionAreaTotal, this->SectionAreaMin, this->SectionAreaMax);
		this->MaxMin(this->segments[i]->GetSize(), this->FragmentationTotal, this->FragmentationMin, this->FragmentationMax);
		this->MaxMin(this->segments[i]->GetBurkTaper(), this->BurkTaperTotal, this->BurkTaperMin, this->BurkTaperMax);
		this->MaxMin(this->segments[i]->GetHillmanTaper(), this->HillmanTaperTotal, this->HillmanTaperMin, this->HillmanTaperMax);

		double diam = 2*this->segments[i]->GetRadii();
		this->MaxMin(diam, this->DiameterTotal, this->DiameterMin, this->DiameterMax);
		double diamPower = pow(diam,1.5);
		this->MaxMin(diamPower, this->DiameterPowerTotal, this->DiameterPowerMin, this->DiameterPowerMax);

		this->MaxMin(this->segments[i]->GetFragmentationSmoothness(), this->ContractionTotal, this->ContractionMin, this->ContractionMax);
		int tempLevel = this->segments[i]->GetLevel();
		if (this->segments[i]->isLeaf())
		{
			this->terminalTips++;

			TraceBit leafBit = this->segments[i]->GetTraceBitsPointer()->back();
			TraceBit leadBit = this->segments[i]->GetTraceBitsPointer()->front();
			this->MaxMin(this->segments[i]->GetSize(), this->TerminalSegmentTotal, this->TerminalSegmentMin, this->TerminalSegmentMax);
			double tipEucDistanceToSoma = this->segments[i]->Euclidean(leafBit, rootBit);
			this->MaxMin(tipEucDistanceToSoma, this->TipToSomaEucDisTotal, this->TipToSomaEucDisMin, this->TipToSomaEucDisMax);
			if(!this->segments[i]->isRoot())
			{
				TraceBit parentBit = this->segments[i]->GetParent()->GetTraceBitsPointer()->back();
				double parentDiam = 2*parentBit.r;
				this->TotalLastParentDiam += parentDiam;
				if (this->LastParentDiamMax < parentDiam)
				{
					this->LastParentDiamMax = parentDiam;
				}else if(this->LastParentDiamMin > parentDiam)
				{
					this->LastParentDiamMin = parentDiam;
				}
			}//fails if non branching tree
			double DiamThreshold = 2*leadBit.r;
			this->DiamThresholdTotal += DiamThreshold;
			if (this->DiamThresholdMin > DiamThreshold)
			{
				this->DiamThresholdMin = DiamThreshold;
			}else if(this->DiamThresholdMax < DiamThreshold)
			{
				this->DiamThresholdMax = DiamThreshold;
			}
			float lx = leafBit.x;
			float ly = leafBit.y;
			float lz = leafBit.z;
			float total; //just so things work
			this->MaxMin(lx, total, this->minX, this->maxX);
			this->MaxMin(ly, total, this->minY, this->maxY);
			this->MaxMin(lz, total, this->minZ, this->maxZ);
			this->MaxMin(tempLevel, this->SumTerminalLevel, this->MinTerminalLevel, this->MaxTerminalLevel);
			this->MaxMin(this->segments[i]->GetPathLength(), this->TerminalPathLength, this->TerminalPathLengthMin, this->TerminalPathLengthMax);
			this->totalTipX += (lx - this->somaX);
			this->totalTipY += (ly - this->somaY);
			this->totalTipZ += (lz - this->somaZ);

			tips.push_back(leafBit);
		}//end if leaf
		else if(!this->segments[i]->isRoot())
		{
			if (this->segments[i]->isActualBifurcation())
			{
				this->actualBifurcations++;
			}
			this->branchPoints++;
			TraceBit leadBit = this->segments[i]->GetTraceBitsPointer()->front();
			double BifToSomaEucDistance = this->segments[i]->Euclidean(leadBit,rootBit);
			this->MaxMin(BifToSomaEucDistance, this->BranchPtToSomaEucDisTotal, this->BranchPtToSomaEucDisMin, this->BranchPtToSomaEucDisMax);
			this->MaxMin(this->segments[i]->GetpartitionAsymmetry(), this->partitionAsymmetry, this->partitionAsymmetryMin, this->partitionAsymmetryMax);
			this->MaxMin(this->segments[i]->GetdaughterRatio(), this->daughterRatio, this->daughterRatioMin, this->daughterRatioMax);
			this->MaxMin(this->segments[i]->GetparentDaughterRatio(), this->parentDaughterRatio, this->parentDaughterRatioMin, this->parentDaughterRatioMax);
			this->MaxMin(this->segments[i]->GetRallPower(), this->rallPower, this->rallPowerMin, this->rallPowerMax);
			this->MaxMin(this->segments[i]->GetPk(), this->Pk, this->PkMin, this->PkMax);
			this->MaxMin(this->segments[i]->GetPk_2(), this->Pk_2, this->Pk_2Min, this->Pk_2Max);
			this->MaxMin(this->segments[i]->GetPk_classic(), this->Pk_classic, this->Pk_classicMin, this->Pk_classicMax);
			this->MaxMin(this->segments[i]->GetBifAmplLocal(), this->BifAmplLocal, this->BifAmplLocalMin, this->BifAmplLocalMax);
			this->MaxMin(this->segments[i]->GetBifAmpRemote(), this->BifAmpRemote, this->BifAmpRemoteMin, this->BifAmpRemoteMax);
			//this->MaxMin(this->segments[i]->GetBranch1()->GetBifTiltLocal(), this->BifTiltLocal, this->BifTiltLocalMin, this->BifTiltLocalMax);
			this->MaxMin(this->segments[i]->GetBifTiltLocalSmall(), this->BifTiltLocalSmall, this->BifTiltLocalSmallMin, this->BifTiltLocalSmallMax);
			this->MaxMin(this->segments[i]->GetBifTiltLocalBig(), this->BifTiltLocalBig, this->BifTiltLocalBigMin, this->BifTiltLocalBigMax);
			this->MaxMin(this->segments[i]->GetBifTiltLocalTwoDaughter(), this->BifTiltLocalTwoDaughter, this->BifTiltLocalTwoDaughterMin, this->BifTiltLocalTwoDaughterMax);
			this->MaxMin(this->segments[i]->GetBifTiltRemoteTwoDaughter(), this->BifTiltRemoteTwoDaughter, this->BifTiltRemoteTwoDaughterMin, this->BifTiltRemoteTwoDaughterMax);
			this->MaxMin(this->segments[i]->GetBifTorqueLocalTwoDaughter(), this->BifTorqueLocalTwoDaughter, this->BifTorqueLocalTwoDaughterMin, this->BifTorqueLocalTwoDaughterMax);
			this->MaxMin(this->segments[i]->GetBifTorqueRemoteTwoDaughter(), this->BifTorqueRemoteTwoDaughter, this->BifTorqueRemoteTwoDaughterMin, this->BifTorqueRemoteTwoDaughterMax);
			//this->MaxMin(this->segments[i]->GetBranch1()->GetBifTorqueLocal(), this->BifTorqueLocal, this->BifTorqueLocalMin, this->BifTorqueLocalMax);
			if(this->segments[i]->isBranch())
			{
			this->MaxMin(this->segments[i]->GetBifTiltRemoteSmall(), this->BifTiltRemoteSmall, this->BifTiltRemoteSmallMin, this->BifTiltRemoteSmallMax);
			this->MaxMin(this->segments[i]->GetBifTiltRemoteBig(), this->BifTiltRemoteBig, this->BifTiltRemoteBigMin, this->BifTiltRemoteBigMax);
			this->MaxMin(this->segments[i]->GetBifTorqueLocalSmall(), this->BifTorqueLocalSmall, this->BifTorqueLocalSmallMin, this->BifTorqueLocalSmallMax, this->TwoDaughterTorqueCount);
			this->MaxMin(this->segments[i]->GetBifTorqueLocalBig(), this->BifTorqueLocalBig, this->BifTorqueLocalBigMin, this->BifTorqueLocalBigMax, this->TwoDaughterTorqueCount);
			if (this->segments[i]->GetTerminalDegree() > 2)
			{
				this->MaxMin(this->segments[i]->GetBifTorqueRemoteBig(), this->BifTorqueRemoteBig, this->BifTorqueRemoteBigMin, this->BifTorqueRemoteBigMax);
				//make test if there is a small angle (maybe if getTorque = -1)
				this->MaxMin(this->segments[i]->GetBifTorqueRemoteSmall(), this->BifTorqueRemoteSmall, this->BifTorqueRemoteSmallMin, this->BifTorqueRemoteSmallMax);
			}

        //this->MaxMin(this->segments[i]->GetBranch2()->GetBifTiltRemote(), this->BifTiltRemote, this->BifTiltRemoteMin, this->BifTiltRemoteMax);
        //this->MaxMin(this->segments[i]->GetBranch2()->GetBifTorqueRemote(), this->BifTorqueRemote, this->BifTorqueRemoteMin, this->BifTorqueRemoteMax);
      }
			if (this->segments[i]->GetLevel() == 1)
			{
				this->MaxMin(this->segments[i]->GetAzimuth(), this->Azimuth, this->AzimuthMin, this->AzimuthMax);
				this->MaxMin(this->segments[i]->GetElevation(), this->Elevation, this->ElevationMin, this->ElevationMax);
			}
			if (this->segments[i]->GetTerminalDegree() ==2)
			{
				this->terminalBifCount++;
				this->MaxMin(this->segments[i]->GetHillmanThreshold(), this->HillmanThreshTotal, this->HillmanThreshMin, this->HillmanThreshMax);
			}
		}//end bifurcation features
	}//end for segment size
	this->skewnessX = (somaX - (minX + (maxX - minX)/2)); /// ((float)(maxX - minX)/2);
	this->skewnessY = (somaY - (minY + (maxY - minY)/2)); /// ((float)(maxY - minY)/2);
	this->skewnessZ = (somaZ - (minZ + (maxZ - minZ)/2)); /// ((float)(maxZ - minZ)/2);
	this->euclideanSkewness = sqrt(pow(skewnessX, 2) + pow(skewnessY, 2) + pow(skewnessZ, 2)) / sqrt(3.0);

	if (tips.size() > 0) //do superellipsoid here
	{
		this->tipMagnitude = sqrt(pow(totalTipX,2)+pow(totalTipY,2)+pow(totalTipZ,2));
		this->tipAzimuth = atan2(totalTipY,totalTipX)*180/PI;
		double hypotenuse = sqrt(pow(totalTipX,2)+pow(totalTipY,2));
		this->tipElevation = atan2(totalTipZ,hypotenuse)*180/PI;
	}
}
void CellTrace::setFileName(std::string newFileName)
{
	this->FileName = newFileName;
}
std::string CellTrace::GetFileName()
{
	std::string outputFileName;
	if (this->FileName != "file")
	{
		outputFileName = this->FileName;
	}
	else
	{
		std::stringstream newCellName; 
		newCellName<<"cell" << this->segments[0]->GetId();
		outputFileName = newCellName.str();
	}
	return outputFileName;
}
void CellTrace::getSomaCoord(double xyz[])
{
	xyz[0] = this->somaX;
	xyz[1] = this->somaY;
	xyz[2] = this->somaZ;
}
void CellTrace::getCellBounds(double bounds[])
{//min then max of x, y, z
	bounds[0] = this->minX;
	bounds[1] = this->maxX;
	bounds[2] = this->minY;
	bounds[3] = this->maxY;
	bounds[4] = this->minZ;
	bounds[5] = this->maxZ;
}
void CellTrace::setDistanceToROI(double newDistance, double Coord_X , double Coord_Y, double Coord_Z)
{
	this->segments[0]->SetDistanceToROI(newDistance);
	this->segments[0]->SetDistanceToROICoord_X(Coord_X);
	this->segments[0]->SetDistanceToROICoord_Y(Coord_Y);
	this->segments[0]->SetDistanceToROICoord_Z(Coord_Z);
}
void CellTrace::SetClassifcation(int predicCol, double prediction, int confCol,double confidence)
{
	//this->segments[0]->SetClassification(prediction, confidence);
	this->segments[0]->editCellFeature(prediction, predicCol);
	this->segments[0]->editCellFeature(confidence, confCol);
}
void CellTrace::addNewFeature(vtkVariant nextFeature)
{
	this->segments[0]->addCellFeature(nextFeature);
}
void CellTrace::clearAll()
{
	this->segments.clear();
	this->NumSegments = 0;
	this->stems = 0;
	this->branchPoints = 0;
	this->terminalTips = 0;
	this->actualBifurcations = 0;
	this->branchingStem = 0;
	//this->TriCount = 0;
	//this->TerminalTriCount = 0;
	//this->notTerminalTriCount = 0;

	this->MinTerminalLevel = 100; //something large for initial value
	this->MaxTerminalLevel = 0;
	this->SumTerminalLevel = 0;

	this->TerminalPathLengthMin = 100;
	this->TerminalPathLengthMax=0;
	this->TerminalPathLength = 0;

	this->FragmentationTotal = 0;
	this->FragmentationMin = 100;
	this->FragmentationMax = 0;

	this->HillmanTaperTotal = 0;
	this->HillmanTaperMin = 100;
	this->HillmanTaperMax = 0;

	this->BurkTaperTotal = 0;
	this->BurkTaperMin = 100;
	this->BurkTaperMax = 0;

	this->HillmanThreshTotal = 0;
	this->HillmanThreshMin = 100;
	this->HillmanThreshMax = 0;
	this->terminalBifCount = 0;

	this->ContractionTotal = 0;
	this->ContractionMin = 100;
	this->ContractionMax = 0;

	this->DiameterTotal = 0;
	this->DiameterMin = 100;
	this->DiameterMax = 0;

	this->DiameterPowerTotal = 0;
	this->DiameterPowerMin = 100;
	this->DiameterPowerMax = 0;

	this->TotalEuclideanPath = 0;
	this->MinEuclideanPath = 1000;
	this->MaxEuclideanPath = 0;

	this->PathLengthTotal = 0;
	this->PathLengthMin = 1000;
	this->PathLengthMax = 0;

	this->TerminalSegmentTotal = 0;
	this->TerminalSegmentMin = 100;
	this->TerminalSegmentMax = 0;

	this->TotalVolume = 0;
	this->SegmentVolumeMax = 0;
	this->SegmentVolumeMin = 1000;

	this->FileName = "file";

	this->somaX = 0;
	this->somaY = 0;
	this->somaZ = 0;

	this->maxX = 0;
	this->maxY = 0;
	this->maxZ = 0;

	this->minX = 100;
	this->minY = 100;
	this->minZ = 100;

	this->sectionAreaTotal = 0;
	this->SectionAreaMax = 0;
	this->SectionAreaMin = 100;

	this->surfaceAreaTotal = 0;
	this->SurfaceAreaMax = 0;
	this->SurfaceAreaMin = 100;

	this->DiamThresholdTotal = 0;
	this->DiamThresholdMax = 0;
	this->DiamThresholdMin = 100;

	this->TotalLastParentDiam = 0;
	this->LastParentDiamMax = 0;
	this->LastParentDiamMin = 100;

	this->somaVolume= 0;
	this->somaSurface = 0;
	this->somaRadii = 0;

	this->MinStemDistance = 1000;
	this->TotalStemDistance = 0;
	this->MaxStemDistance = 0; 
	this->EstimatedSomaRadius = -PI;

	this->daughterRatio = 0;
	this->daughterRatioMin = 100;
	this->daughterRatioMax = 0;

	this->parentDaughterRatio = 0;
	this->parentDaughterRatioMin = 100;
	this->parentDaughterRatioMax = 0;

	this->partitionAsymmetry = 0;
	this->partitionAsymmetryMin = 100;
	this->partitionAsymmetryMax = 0;

	this->rallPower = 0;
	this->rallPowerMin = 100;
	this->rallPowerMax = 0;

	this->Pk = 0;
	this->PkMin = 100;
	this->PkMax = 0;

	this->Pk_2 = 0;
	this->Pk_2Min = 100;
	this->Pk_2Max = 0;

	this->Pk_classic = 0;
	this->Pk_classicMin = 100;
	this->Pk_classicMax = 0;

	this->BifAmplLocal = 0;
	this->BifAmplLocalMin = 180;
	this->BifAmplLocalMax = 0;

	this->BifAmpRemote = 0;
	this->BifAmpRemoteMin = 180;
	this->BifAmpRemoteMax = 0;

	this->BifTiltLocalSmall = 0;
	this->BifTiltLocalSmallMin = 180;
	this->BifTiltLocalSmallMax = 0;
	this->BifTiltLocalBig = 0;
	this->BifTiltLocalBigMin = 180;
	this->BifTiltLocalBigMax = 0;

	this->BifTiltRemoteSmall = 0;
	this->BifTiltRemoteSmallMin = 180;
	this->BifTiltRemoteSmallMax = 0;
	this->BifTiltRemoteBig = 0;
	this->BifTiltRemoteBigMin = 180;
	this->BifTiltRemoteBigMax = 0;

	this->BifTiltLocalTwoDaughter = 0;
	this->BifTiltRemoteTwoDaughter = 0;

	this->BifTorqueLocalSmall = 0;
	this->BifTorqueLocalSmallMin = 180;
	this->BifTorqueLocalSmallMax = 0;
	this->BifTorqueLocalBig = 0;
	this->BifTorqueLocalBigMin = 180;
	this->BifTorqueLocalBigMax = 0;

	this->BifTorqueRemoteSmall = 0;
	this->BifTorqueRemoteSmallMin = 180;
	this->BifTorqueRemoteSmallMax = 0;
	this->BifTorqueRemoteBig = 0;
	this->BifTorqueRemoteBigMin = 180;
	this->BifTorqueRemoteBigMax = 0;

	this->BifTorqueLocalTwoDaughter = 0;
	this->BifTorqueRemoteTwoDaughter = 0;
	this->BifTorqueLocalTwoDaughterMax = 0;
	this->BifTorqueRemoteTwoDaughterMax = 0;
	this->TwoDaughterTorqueCount = 0;

	this->Azimuth = 0;
	this->AzimuthMin = 180;
	this->AzimuthMax = 0;

	this->Elevation = 0;
	this->ElevationMin = 180;
	this->ElevationMax = 0;

	this->TipToSomaEucDisTotal = 0;
	this->TipToSomaEucDisMin = 1000;
	this->TipToSomaEucDisMax = 0;

	this->BranchPtToSomaEucDisTotal = 0;
	this->BranchPtToSomaEucDisMin = 1000;
	this->BranchPtToSomaEucDisMax = 0;

	this->totalTipX = 0;
	this->totalTipY = 0;
	this->totalTipZ = 0;

	this->tipMagnitude = -1000;
	this->tipAzimuth = -1000;
	this->tipElevation = -1000;

	this->DeviceDistance = 0;
	this->prediction = -PI;
	this->confidence = -PI;
}
void CellTrace::MaxMin(double NewValue, double &total, double &Min, double &Max)
{
	if (NewValue != -1)
	{
		total += NewValue;
		if (NewValue > Max)
		{
			Max = NewValue;
		}
		if (NewValue < Min)
		{
			Min = NewValue;
		}
	}
}
void CellTrace::MaxMin(double NewValue, double &total, double &Min, double &Max, int &Count)
{
	if (NewValue != -1)
	{
		this->MaxMin(NewValue,total,Min,Max);
		Count++;
	}
}
void CellTrace::MaxMin(float NewValue, float &total, float &Min, float &Max)
{
	if (NewValue != -1)
	{
		total += NewValue;
		if (NewValue > Max)
		{
			Max = NewValue;
		}
		if (NewValue < Min)
		{
			Min = NewValue;
		}
	}
}
void CellTrace::MaxMin(int NewValue, int &total, int &Min, int &Max)
{
	if (NewValue != -1)
	{
		total += NewValue;
		if (NewValue > Max)
		{
			Max = NewValue;
		}
		if (NewValue < Min)
		{
			Min = NewValue;
		}
	}
}
vtkSmartPointer<vtkVariantArray> CellTrace::DataRow()
{
	CellData->Reset();
	CellData->InsertNextValue(this->segments[0]->GetId());

	CellData->InsertNextValue(this->maxX - this->minX);//Width
	CellData->InsertNextValue(this->maxY - this->minY);//Length
	CellData->InsertNextValue(this->maxZ - this->minZ);//Height
	
	CellData->InsertNextValue(this->somaX);
	CellData->InsertNextValue(this->somaY);
	CellData->InsertNextValue(this->somaZ);
	CellData->InsertNextValue(this->somaRadii);
	CellData->InsertNextValue(this->somaSurface);
	CellData->InsertNextValue(this->somaVolume);

	CellData->InsertNextValue(this->skewnessX);
	CellData->InsertNextValue(this->skewnessY);
	CellData->InsertNextValue(this->skewnessZ);
	CellData->InsertNextValue(this->euclideanSkewness);

	CellData->InsertNextValue(this->NumSegments);
	CellData->InsertNextValue(this->stems);
	//CellData->InsertNextValue(this->branchingStem);
	CellData->InsertNextValue(this->branchPoints);
	//CellData->InsertNextValue(this->actualBifurcations);
	CellData->InsertNextValue(this->terminalTips);
	if (this->branchPoints == 0) 
	{
		this->branchPoints = 1;
	}//protect from divide by zero
	if (this->actualBifurcations == 0)
	{
		this->actualBifurcations = 1;
	}

	CellData->InsertNextValue(this->DiameterMin);
	CellData->InsertNextValue(this->DiameterTotal / this->NumSegments);
	CellData->InsertNextValue(this->DiameterMax);

	CellData->InsertNextValue(this->DiameterPowerMin);
	CellData->InsertNextValue(this->DiameterPowerTotal / this->NumSegments);
	CellData->InsertNextValue(this->DiameterPowerMax);

	CellData->InsertNextValue(this->TotalVolume);
	CellData->InsertNextValue(this->SegmentVolumeMin);
	CellData->InsertNextValue(this->TotalVolume/this->NumSegments);////average segment Volume
	CellData->InsertNextValue(this->SegmentVolumeMax);
	CellData->InsertNextValue(this->surfaceAreaTotal);
	CellData->InsertNextValue(this->SurfaceAreaMin);
	CellData->InsertNextValue(this->surfaceAreaTotal/this->NumSegments);
	CellData->InsertNextValue(this->SurfaceAreaMax);
	CellData->InsertNextValue(this->sectionAreaTotal);
	CellData->InsertNextValue(this->SectionAreaMin);
	CellData->InsertNextValue(this->sectionAreaTotal/this->NumSegments);
	CellData->InsertNextValue(this->SectionAreaMax);

	//CellData->InsertNextValue(this->BurkTaperTotal);
	CellData->InsertNextValue(this->BurkTaperMin);
	CellData->InsertNextValue(this->BurkTaperTotal / this->NumSegments);
	CellData->InsertNextValue(this->BurkTaperMax);

	//CellData->InsertNextValue(this->HillmanTaperTotal);
	CellData->InsertNextValue(this->HillmanTaperMin);
	CellData->InsertNextValue(this->HillmanTaperTotal / this->NumSegments);
	CellData->InsertNextValue(this->BurkTaperMax);

	CellData->InsertNextValue(this->TotalEuclideanPath);
	CellData->InsertNextValue(this->TotalEuclideanPath/this->NumSegments);//average segment euclidean length
	CellData->InsertNextValue(this->PathLengthTotal);
	CellData->InsertNextValue(this->PathLengthMin);
	CellData->InsertNextValue(this->PathLengthTotal/this->NumSegments);//average segment length
	CellData->InsertNextValue(this->PathLengthMax);

	CellData->InsertNextValue(this->MinStemDistance);
	CellData->InsertNextValue(this->EstimatedSomaRadius);
	CellData->InsertNextValue(this->MaxStemDistance);

	CellData->InsertNextValue(this->ContractionMin);
	double aveContraction = this->ContractionTotal / this->NumSegments;
	if (aveContraction != aveContraction)
	{
		CellData->InsertNextValue(-PI);
	}
	else
	{
		CellData->InsertNextValue(aveContraction);
	}
	CellData->InsertNextValue(this->ContractionMax);

	CellData->InsertNextValue(this->FragmentationTotal);
	CellData->InsertNextValue(this->FragmentationMin);
	CellData->InsertNextValue(this->FragmentationTotal / this->actualBifurcations);
	CellData->InsertNextValue(this->FragmentationMax);

	CellData->InsertNextValue(this->daughterRatioMin);
	CellData->InsertNextValue(this->daughterRatio / this->actualBifurcations);
	CellData->InsertNextValue(this->daughterRatioMax);

	CellData->InsertNextValue(this->parentDaughterRatioMin);
	CellData->InsertNextValue(this->parentDaughterRatio/ this->actualBifurcations);
	CellData->InsertNextValue(this->parentDaughterRatioMax);

	CellData->InsertNextValue(this->partitionAsymmetryMin);
	CellData->InsertNextValue(this->partitionAsymmetry / this->actualBifurcations);
	CellData->InsertNextValue(this->partitionAsymmetryMax);

	CellData->InsertNextValue(this->rallPowerMin);
	CellData->InsertNextValue(this->rallPower / this->actualBifurcations);
	CellData->InsertNextValue(this->rallPowerMax);

	CellData->InsertNextValue(this->PkMin);
	CellData->InsertNextValue(this->Pk / this->actualBifurcations);
	CellData->InsertNextValue(this->PkMax);

	CellData->InsertNextValue(this->Pk_classicMin);
	CellData->InsertNextValue(this->Pk_classic / this->actualBifurcations);
	CellData->InsertNextValue(this->Pk_classicMax);

	CellData->InsertNextValue(this->Pk_2Min);
	CellData->InsertNextValue(this->Pk_2 / this->actualBifurcations);
	CellData->InsertNextValue(this->Pk_2Max);

	double AveAzimuth = -PI;
	double AveElevation = -PI;
	if (this->stems !=0)
	{
		AveAzimuth = this->Azimuth / this->stems;
		AveElevation = this->Elevation / this->stems;
	}
	CellData->InsertNextValue(this->AzimuthMin);
	CellData->InsertNextValue(AveAzimuth);
	CellData->InsertNextValue(this->AzimuthMax);
	CellData->InsertNextValue(this->ElevationMin);
	CellData->InsertNextValue(AveElevation);
	CellData->InsertNextValue(this->ElevationMax);

	CellData->InsertNextValue(this->BifAmplLocal / this->actualBifurcations);
	CellData->InsertNextValue(this->BifAmplLocalMin);
	CellData->InsertNextValue(this->BifAmplLocalMax);
	CellData->InsertNextValue(this->BifTiltLocalSmall/ this->actualBifurcations);
	CellData->InsertNextValue(this->BifTiltLocalSmallMin);
	CellData->InsertNextValue(this->BifTiltLocalSmallMax);
	CellData->InsertNextValue(this->BifTiltLocalBig/ this->actualBifurcations);
	CellData->InsertNextValue(this->BifTiltLocalBigMin);
	CellData->InsertNextValue(this->BifTiltLocalBigMax);

	if (this->TwoDaughterTorqueCount == 0)
	{
		this->TwoDaughterTorqueCount = 1;
	}
	CellData->InsertNextValue(this->BifTorqueLocalTwoDaughter/ this->TwoDaughterTorqueCount);
	//CellData->InsertNextValue(this->BifTorqueLocalSmall/ (this->branchPoints-2));
	CellData->InsertNextValue(this->BifTorqueLocalSmallMin);
	CellData->InsertNextValue(this->BifTorqueLocalSmallMax);
	//CellData->InsertNextValue(this->BifTorqueLocalBig/ (this->branchPoints-2));
	CellData->InsertNextValue(this->BifTorqueLocalBigMin);
	CellData->InsertNextValue(this->BifTorqueLocalBigMax);

	CellData->InsertNextValue(this->BifAmpRemote / this->actualBifurcations);
	CellData->InsertNextValue(this->BifAmpRemoteMin);
	CellData->InsertNextValue(this->BifAmpRemoteMax);
	CellData->InsertNextValue(this->BifTiltRemoteSmall / this->actualBifurcations);
	CellData->InsertNextValue(this->BifTiltRemoteSmallMin);
	CellData->InsertNextValue(this->BifTiltRemoteSmallMax);
	CellData->InsertNextValue(this->BifTiltRemoteBig / this->actualBifurcations);
	CellData->InsertNextValue(this->BifTiltRemoteBigMin);
	CellData->InsertNextValue(this->BifTiltRemoteBigMax);
	CellData->InsertNextValue(this->BifTorqueRemoteTwoDaughter/ this->TwoDaughterTorqueCount);
	//CellData->InsertNextValue(this->BifTorqueRemoteSmall / (this->branchPoints-2));
	CellData->InsertNextValue(this->BifTorqueRemoteSmallMin);
	CellData->InsertNextValue(this->BifTorqueRemoteSmallMax);
	//CellData->InsertNextValue(this->BifTorqueRemoteBig / (this->branchPoints-2));
	CellData->InsertNextValue(this->BifTorqueRemoteBigMin);
	CellData->InsertNextValue(this->BifTorqueRemoteBigMax);

	CellData->InsertNextValue(this->MinTerminalLevel);
	CellData->InsertNextValue(this->TerminalPathLengthMin);
	CellData->InsertNextValue(this->SumTerminalLevel /this->terminalTips);//average terminal level
	CellData->InsertNextValue(this->TerminalPathLength/this->terminalTips);//now average path to end
	CellData->InsertNextValue(this->MaxTerminalLevel);
	CellData->InsertNextValue(this->TerminalPathLengthMax);

	CellData->InsertNextValue(this->TerminalSegmentTotal);
	CellData->InsertNextValue(this->TerminalSegmentMin);
	CellData->InsertNextValue(this->TerminalSegmentTotal / this->terminalTips); //average
	CellData->InsertNextValue(this->TerminalSegmentMax);

	CellData->InsertNextValue(this->DiamThresholdMin);
	CellData->InsertNextValue(this->DiamThresholdTotal/this->terminalTips);
	CellData->InsertNextValue(this->DiamThresholdMax);
	CellData->InsertNextValue(this->LastParentDiamMin);
	CellData->InsertNextValue(this->TotalLastParentDiam/this->terminalTips);
	CellData->InsertNextValue(this->LastParentDiamMax);

	CellData->InsertNextValue(this->HillmanThreshMin); 
	if (this->terminalBifCount != 0)
	{
		CellData->InsertNextValue(this->HillmanThreshTotal/ this->terminalBifCount);
	}else
	{
		CellData->InsertNextValue(-PI);
	}
	CellData->InsertNextValue(this->HillmanThreshMax);

	CellData->InsertNextValue(this->BranchPtToSomaEucDisMin);
	CellData->InsertNextValue(this->BranchPtToSomaEucDisTotal/this->branchPoints);
	CellData->InsertNextValue(this->BranchPtToSomaEucDisMax);
	CellData->InsertNextValue(this->TipToSomaEucDisMin);
	CellData->InsertNextValue(this->TipToSomaEucDisTotal/this->terminalTips);
	CellData->InsertNextValue(this->TipToSomaEucDisMax);

	CellData->InsertNextValue(this->tipMagnitude);
	CellData->InsertNextValue(this->tipAzimuth);
	CellData->InsertNextValue(this->tipElevation);

	CellData->InsertNextValue(this->GetFileName().c_str());
	/*CellData->InsertNextValue(this->prediction);
	CellData->InsertNextValue(this->confidence);*/
	CellData->InsertNextValue( this->segments[0]->GetDistanceToROI());
	//std::cout << this->FileName << std::endl;
	return CellData;
}

vtkSmartPointer<vtkVariantArray> CellTrace::ConvertDefaultValueToNull(vtkSmartPointer<vtkVariantArray> row)
{
	vtkSmartPointer<vtkVariantArray> newRow = vtkSmartPointer<vtkVariantArray>::New();
	newRow->SetNumberOfValues( row->GetNumberOfValues());
	for( vtkIdType i = 0; i < row->GetNumberOfValues(); i++)
	{
		vtkVariant var = row->GetValue( i);
		if( i == row->GetNumberOfValues() - 2) // file name
		{
			newRow->SetValue( i, var);
		}
		else
		{
			double varDouble = var.ToDouble();
			if( DefaultValue[i] == -1)
			{
				newRow->SetValue( i, var);
			}
			else if( varDouble == DefaultValue[i])
			{
				//newRow->SetValue( i, NULL);
			}
			else
			{
				newRow->SetValue( i, var);
			}
		}
	}
	return newRow;
}

vtkSmartPointer<vtkVariantArray> CellTrace::GetExtendedDataRow(int CheckAddFeatures)
{
	this->DataRow();
	std::vector<vtkVariant> extendedFeatures = this->segments[0]->GetCellFeatures();

	if (extendedFeatures.size() != CheckAddFeatures)
	{
		for (int i = 0; i < CheckAddFeatures; i++)
		{
			CellData->InsertNextValue(- PI);
		}
	}
	else
	{
		for (int j = 0; j < CheckAddFeatures; j++)
		{
			CellData->InsertNextValue(extendedFeatures[j]);
		}
	}
	//std::cout << "row size " << this->CellData->GetNumberOfValues()<< std::endl;
	return this->CellData;
}
vtkSmartPointer<vtkVariantArray> CellTrace::BoundsRow()
{
	// id, somaX somaY somaZ min/max xyz, skewness xyz
	vtkSmartPointer<vtkVariantArray> CellBoundsRow = vtkSmartPointer<vtkVariantArray>::New();
	CellBoundsRow->InsertNextValue(this->segments[0]->GetId());
	CellBoundsRow->InsertNextValue(this->somaX);
	CellBoundsRow->InsertNextValue(this->somaY);
	CellBoundsRow->InsertNextValue(this->somaZ);
	CellBoundsRow->InsertNextValue(this->minX);
	CellBoundsRow->InsertNextValue(this->maxX);
	CellBoundsRow->InsertNextValue(this->minY);
	CellBoundsRow->InsertNextValue(this->maxY);
	CellBoundsRow->InsertNextValue(this->minZ);
	CellBoundsRow->InsertNextValue(this->maxZ);
	return CellBoundsRow;
}
std::string CellTrace::BasicFeatureString()
{
	std::stringstream features;
	features << this->GetFileName().c_str() << "\t";
	features << this->somaX << "\t";
	features << this->somaY << "\t";
	features << this->somaZ << "\t";
	return features.str();
}
std::set<long int> CellTrace::TraceIDsInCell()
{
	return this->IDs;
}
unsigned int CellTrace::rootID()
{
	return this->segments[0]->GetId();
}

TraceLine * CellTrace::getRootTrace()
{
	return this->segments.front();
}

std::vector<TraceLine *> CellTrace::getSegments()
{
	return this->segments;
}
