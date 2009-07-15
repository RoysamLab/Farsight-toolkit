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

void FindSeedPoints2(CImage& anImage, std::vector<CPoint> & seed_candidates)
{
	// relative axial resolution = lateral resolution / axial resolution
	float relative_axial_resolution = 1.0;
	int giMARGIN2 = 3;

	unsigned char * *ImageData = anImage.data;

	// do multi-scale edge detection
	typedef itk::Image< float, 1 > ImageType;
	ImageType::Pointer axial_intensities = ImageType::New();
	ImageType::SizeType axial_size;
	axial_size[0] = giSLICES;
	ImageType::RegionType axial_region;
	axial_region.SetSize(axial_size);
	axial_intensities->SetRegions(axial_region);
	axial_intensities->Allocate();

	ImageType::Pointer horz_intensities = ImageType::New();
	ImageType::SizeType horz_size;
	horz_size[0] = giCOLS;
	ImageType::RegionType horz_region;
	horz_region.SetSize(horz_size);
	horz_intensities->SetRegions(horz_region);
	horz_intensities->Allocate();

	ImageType::Pointer vert_intensities = ImageType::New();
	ImageType::SizeType vert_size;
	vert_size[0] = giROWS;
	ImageType::RegionType vert_region;
	vert_region.SetSize(vert_size);
	vert_intensities->SetRegions(vert_region);
	vert_intensities->Allocate();

	std::vector<CPoint> pre_seed_candidates;

	ImageType::IndexType index;
	// search axially at x,y grid intersections
	for (int x = giGRIDSIZE; x < giCOLS - giGRIDSIZE; x += giGRIDSIZE)
	{
		for (int y = giGRIDSIZE; y < giROWS - giGRIDSIZE; y += giGRIDSIZE)
		{
			for (int z = 0; z < giSLICES; z++)
			{
				index[0] = z;
				axial_intensities->SetPixel(index, static_cast<float>(The3DImage->data[z][y][x]));;
			}

			std::vector< std::pair<int, std::pair<int,int> > > candidates 
				= MatchedDetectorFor1DProfile(axial_intensities, giPadding, 1, giShiftDistance);

			for(unsigned int i = 0; i < candidates.size(); i++)
			{
				int z = candidates[i].second.first;
				int value = candidates[i].first;
				int width = candidates[i].second.second;

				CPoint candidate;
				candidate.m_iX = x;
				candidate.m_iY = y;
				candidate.m_iZ = z;
				candidate.m_iValue = value;
				candidate.m_fVWidth = static_cast<float>(width);
				pre_seed_candidates.push_back(candidate);
			}
		}
	}

	std::cout << "icandidates: " << pre_seed_candidates.size() << std::endl;

	int num_of_axial_seeds = pre_seed_candidates.size();
	std::vector<CPoint> laterally_detected_seed_candidates;
	// search for more seeds. 
	for(unsigned int i = 0; i < num_of_axial_seeds; i++)
	{
		//pre_seed_candidates[i].Print();
		std::cout << i << " of " << num_of_axial_seeds << "\r" << std::endl;
		int z = pre_seed_candidates[i].m_iZ;
		int y = pre_seed_candidates[i].m_iY;
		int x;
		// search along a horizontal line
		for (x = 0; x < giCOLS; x++)
		{
			index[0] = x;
			horz_intensities->SetPixel(index, static_cast<float>(The3DImage->data[z][y][x]));;
		}
		std::vector< std::pair<int, std::pair<int,int> > > candidates_horz 
			= MatchedDetectorFor1DProfile(horz_intensities, giMARGIN2, 1, giShiftDistance);

		int closest_x, index_x;
		float distance; 
		float closest_distance;
		index_x = -1;
		for( int a = 0; a < candidates_horz.size(); a++ )
		{
			x = candidates_horz[a].second.first;
			distance = (x - pre_seed_candidates[i].m_iX) * (x - pre_seed_candidates[i].m_iX);
			if (a == 0)
				closest_distance = distance;
			if( distance <= closest_distance )
			{
				closest_distance = distance;
				closest_x = x;
				index_x = a;
			}
			CPoint candidate(pre_seed_candidates[i]);
			candidate.m_iX = x;
			candidate.m_iValue = max(candidate.m_iValue, candidates_horz[a].first);
			candidate.m_fHWidth = static_cast<float>(candidates_horz[a].second.second);
			candidate.m_iHDir = NumOfDirections/4;
			laterally_detected_seed_candidates.push_back(candidate);

		}		

		x = pre_seed_candidates[i].m_iX;
		//search along a vertical line
		for (y = 0; y < giROWS; y++)
		{
			index[0] = y;
			vert_intensities->SetPixel(index, static_cast<float>(The3DImage->data[z][y][x]));;
		}
		std::vector< std::pair<int, std::pair<int,int> > > candidates_vert 
			= MatchedDetectorFor1DProfile(vert_intensities, giMARGIN2, 1, giShiftDistance);

		int closest_y, index_y;
		index_y = -1;
		for( int b = 0; b < candidates_vert.size(); b++ )
		{
			y = candidates_vert[b].second.first;
			distance = (y - pre_seed_candidates[i].m_iY) * (y - pre_seed_candidates[i].m_iY);
			if (b == 0)
				closest_distance = distance;
			if( distance <= closest_distance )
			{
				closest_distance = distance;
				closest_y = y;
				index_y = b;
			}
			CPoint candidate(pre_seed_candidates[i]);
			candidate.m_iY = y;
			candidate.m_iValue = max(candidate.m_iValue, candidates_vert[b].first);
			candidate.m_fHWidth = static_cast<float>(candidates_vert[b].second.second);
			candidate.m_iHDir = 0;
			laterally_detected_seed_candidates.push_back(candidate);
		}

		int response_x = (candidates_horz.size()) ? candidates_horz[index_x].first : 0;
		int response_y = (candidates_vert.size()) ? candidates_vert[index_y].first : 0;
		int response_z = pre_seed_candidates[i].m_iValue;	
		int best_x = (candidates_horz.size()) ? candidates_horz[index_x].second.first : 0;
		int best_y = (candidates_vert.size()) ? candidates_vert[index_y].second.first : 0;

		float size_ratio;

		if( response_x > response_y )
		{
			if( candidates_horz.size() )
			{
				pre_seed_candidates[i].m_iX = candidates_horz[index_x].second.first;
				pre_seed_candidates[i].m_fHWidth = static_cast<float>(candidates_horz[index_x].second.second);
				size_ratio = max(pre_seed_candidates[i].m_fHWidth, pre_seed_candidates[i].m_fVWidth) / min(pre_seed_candidates[i].m_fHWidth, pre_seed_candidates[i].m_fVWidth);
				//if (size_ratio > sqrt(2.0) * 1.2)
				//	pre_seed_candidates[i].m_iValue = 0;
				pre_seed_candidates[i].m_iValue = max(pre_seed_candidates[i].m_iValue, candidates_horz[index_x].first);
				pre_seed_candidates[i].m_iHDir = NumOfDirections/4;
				CPoint test(pre_seed_candidates[i]); 
			}
		}
		else
		{
			if( candidates_vert.size() )
			{
				pre_seed_candidates[i].m_iY = candidates_vert[index_y].second.first;
				pre_seed_candidates[i].m_fHWidth = static_cast<float>(candidates_vert[index_y].second.second);
				size_ratio = max(pre_seed_candidates[i].m_fHWidth, pre_seed_candidates[i].m_fVWidth) / min(pre_seed_candidates[i].m_fHWidth, pre_seed_candidates[i].m_fVWidth);
				//if (size_ratio > sqrt(2.0) * 1.2)
				//	pre_seed_candidates[i].m_iValue = 0;
				pre_seed_candidates[i].m_iValue = max(pre_seed_candidates[i].m_iValue, candidates_vert[index_y].first);
				pre_seed_candidates[i].m_iHDir = 0;
				CPoint test(pre_seed_candidates[i]); 
			}
		}
	}
	

	std::cout << "Found " << laterally_detected_seed_candidates.size() << " more" << std::endl;
	int where = 0;
	for(std::vector<CPoint>::iterator i = laterally_detected_seed_candidates.begin(); i != laterally_detected_seed_candidates.end(); i++)
	{
		std::cout << where++ << " of " << laterally_detected_seed_candidates.size() << "\r" << std::endl;
		int x = i->m_iX;
		int y = i->m_iY;
		for (int z = 0; z < giSLICES; z++)
		{
			index[0] = z;
			axial_intensities->SetPixel(index, static_cast<float>(The3DImage->data[z][y][x]));;
		}

		std::vector< std::pair<int, std::pair<int,int> > > candidates 
			= MatchedDetectorFor1DProfile(axial_intensities, giPadding, 1, giShiftDistance);

		if(candidates.size())
		{
			float distance; 
			float closest_distance;
			int closest_z, index_z;
			index_z = -1;
			for( int b = 0; b < candidates.size(); b++ )
			{
				int z = candidates[b].second.first;
				distance = (z - i->m_iZ) * (z - i->m_iZ);
				if (b == 0)
					closest_distance = distance;
				if( distance <= closest_distance )
				{
					closest_distance = distance;
					closest_z = z;
					index_z = b;
				}
			}
			i->m_iZ = closest_z;
			i->m_fVWidth = candidates[index_z].second.second;
			i->m_iValue = max(i->m_iValue , candidates[index_z].first);
			pre_seed_candidates.push_back(*i);
		}
	}

	for(std::vector<CPoint>::iterator i = pre_seed_candidates.begin(); i != pre_seed_candidates.end(); i++)
	{
		if(i->m_iValue)
			seed_candidates.push_back(*i);
		if(i->m_iX == 135 && i->m_iY == 541)
			i->Print();
	}	


	sort (seed_candidates.begin(), seed_candidates.end(), value_descending());
	//FilterSeedCandidates2D(seed_candidates);
}

bool VerifySeedPoint3D2(CPoint* aPoint)
{
	bool result = true;

	//int y = aPoint->m_iY;
	//int x = aPoint->m_iX;
	//
	//std::vector<unsigned char> intensities(giSLICES,0);

	//for (int z = 0; z < giSLICES; z++)
	//	intensities[z] = The3DImage->data[z][y][x];

	//std::vector< std::pair<int, std::pair<int,int> > > candidates 
	//	= MatchedDetectorFor1DProfile(intensities, giPadding, 1, giShiftDistance);
	//
	//if(candidates.size() == 0)
	//	return false;

	//aPoint->m_iZ = candidates[0].second.first;
	//aPoint->m_fVWidth = candidates[0].second.second*2.0f;

	// an array of votes for each direction
	//int HDirVotes[NumOfDirections];
	int* DirVotes = new int [NumOfDirections* NumOfDirections];
	int* Response = new int [NumOfDirections* NumOfDirections];
	//memset(HDirVotes, 0, sizeof(int)*NumOfDirections);
	memset(DirVotes, 0, sizeof(int) * NumOfDirections * NumOfDirections);
	memset(Response, 0, sizeof(int) * NumOfDirections * NumOfDirections);

	// since we call this function before we perform any tracking, 
	// we don't have to see if the candidate point lies on an already 
	// tracked vessel

	CPoint tempPoint(*aPoint);
	CPoint HLeftPoint, HRightPoint;
	CPoint VLeftPoint, VRightPoint;
	CPoint BestHLeftPoint, BestVLeftPoint;
	CPoint BestHRightPoint, BestVRightPoint;

	int VDirections[9];
	PrepareDirectionsArray(VDirections, 9, 0);
	int HDirections[9];
	PrepareDirectionsArray(HDirections, 9, aPoint->m_iHDir);

	int* HRResponse = new int [NumOfDirections* NumOfDirections];
	int* HLResponse = new int [NumOfDirections* NumOfDirections];
	int* VRResponse = new int [NumOfDirections* NumOfDirections];
	int* VLResponse = new int [NumOfDirections* NumOfDirections];
	memset(HRResponse, 0, sizeof(int) * NumOfDirections * NumOfDirections);
	memset(HLResponse, 0, sizeof(int) * NumOfDirections * NumOfDirections);
	memset(VRResponse, 0, sizeof(int) * NumOfDirections * NumOfDirections);
	memset(VLResponse, 0, sizeof(int) * NumOfDirections * NumOfDirections);

	int HRBestResponse, HLBestResponse, VRBestResponse, VLBestResponse;
	HRBestResponse = HLBestResponse = VRBestResponse = VLBestResponse = 0;
	int HRBestForegroundEstimate, HLBestForegroundEstimate, VRBestForegroundEstimate, VLBestForegroundEstimate;
	HRBestForegroundEstimate = HLBestForegroundEstimate = VRBestForegroundEstimate = VLBestForegroundEstimate = 0;
	int HRBestBackgroundEstimate, HLBestBackgroundEstimate, VRBestBackgroundEstimate, VLBestBackgroundEstimate;
	HRBestBackgroundEstimate = HLBestBackgroundEstimate = VRBestBackgroundEstimate = VLBestBackgroundEstimate = 0;
	int index = 0;

	int Hdir, Vdir;
	for (unsigned int h = 0; h < 9; h++)
	{
		Hdir = HDirections[h];
		for (unsigned int v = 0; v < 9; v++)
		{
			Vdir = VDirections[v];
			HRResponse[index] = gHRightTemplatesArray[Hdir][Vdir]->CalculateMaxResponse2(&tempPoint, & HRightPoint);
			HLResponse[index] = gHLeftTemplatesArray[Hdir][Vdir]->CalculateMaxResponse2(&tempPoint, & HLeftPoint);
			VRResponse[index] = gVRightTemplatesArray[Hdir][Vdir]->CalculateMaxResponse3(&tempPoint, & VRightPoint);
			VLResponse[index] = gVLeftTemplatesArray[Hdir][Vdir]->CalculateMaxResponse3(&tempPoint, & VLeftPoint);	

			if( HRResponse[index] > HRBestResponse)
			{
				BestHRightPoint = HRightPoint;
				HRBestForegroundEstimate = gHRightTemplatesArray[Hdir][Vdir]->m_iForegroundPixelValue;
				HRBestBackgroundEstimate = gHRightTemplatesArray[Hdir][Vdir]->m_iBackgroundPixelValue;
				HRBestResponse = HRResponse[index]; 
			}

			if( HLResponse[index] > HLBestResponse)
			{
				BestHLeftPoint = HLeftPoint;
				HLBestForegroundEstimate = gHLeftTemplatesArray[Hdir][Vdir]->m_iForegroundPixelValue;
				HLBestBackgroundEstimate = gHLeftTemplatesArray[Hdir][Vdir]->m_iBackgroundPixelValue;
				HLBestResponse = HLResponse[index];

			}

			if( VRResponse[index] > VRBestResponse)
			{
				BestVRightPoint = VRightPoint;
				VRBestForegroundEstimate = gVRightTemplatesArray[Hdir][Vdir]->m_iForegroundPixelValue;
				VRBestBackgroundEstimate = gVRightTemplatesArray[Hdir][Vdir]->m_iBackgroundPixelValue;
				VRBestResponse = VRResponse[index];
			}

			if( VLResponse[index] > VLBestResponse)
			{
				BestVLeftPoint = VLeftPoint;
				VLBestForegroundEstimate = gVLeftTemplatesArray[Hdir][Vdir]->m_iForegroundPixelValue;
				VRBestBackgroundEstimate = gVLeftTemplatesArray[Hdir][Vdir]->m_iBackgroundPixelValue;
				VLBestResponse = VLResponse[index];
			}

			index++;
		}
	}

	int BestResponse = HRBestResponse + HLBestResponse + VRBestResponse + VLBestResponse;
	int BestForegroundEstimate = HRBestForegroundEstimate + HLBestForegroundEstimate 
		+ VRBestForegroundEstimate + VLBestForegroundEstimate;
	int BestBackgroundEstimate = HRBestBackgroundEstimate + HLBestBackgroundEstimate
		+ VRBestBackgroundEstimate + VLBestBackgroundEstimate;

	BestResponse /= 4;
	BestForegroundEstimate /= 4;
	BestBackgroundEstimate /= 4;

	int threshold = 0;
	if (gf3DStdDev > 3.0)
		threshold = 3 * 3 * giUsedTemplateLength;
	else if (gf3DStdDev > 0.0)
		threshold = static_cast<int>(gf3DStdDev * 3.0 * giUsedTemplateLength);
	else
		threshold = static_cast<int>(3.0 * giUsedTemplateLength);

	threshold *= 3;

	// we assume starting points must have at least a contrast of 3
	if (BestResponse < threshold)
		result = false;

	if (BestHLeftPoint.m_iValue < threshold ||
		BestHRightPoint.m_iValue < threshold ||
		BestVRightPoint.m_iValue < threshold ||
		BestVLeftPoint.m_iValue < threshold)
		result = false;

	if (result/* && VerifySeedResponses2D(Response, DirVotes)*/)
	{
		aPoint->m_iValue = Round(static_cast<double>(BestResponse) / 4.0 );//* H_similarity * V_similarity);

		aPoint->m_iHDir = BestHRightPoint.m_iHDir;
		if (BestHLeftPoint.m_iValue > BestHRightPoint.m_iValue)
			aPoint->m_iHDir = BestHLeftPoint.m_iHDir;

		aPoint->m_iVDir = BestVRightPoint.m_iVDir;
		if (BestVLeftPoint.m_iValue > BestVRightPoint.m_iValue)
			aPoint->m_iVDir = BestVLeftPoint.m_iVDir;

		GetCenterLocation(*aPoint,
			BestHRightPoint,
			BestHLeftPoint,
			BestVRightPoint,
			BestVLeftPoint);

		if(!The3DImage->WithinImageMargin(*aPoint,giMARGIN))
			return false;
		double Xdiff = BestHLeftPoint.m_iX - BestHRightPoint.m_iX;
		double Ydiff = BestHLeftPoint.m_iY - BestHRightPoint.m_iY;
		double Zdiff = BestHLeftPoint.m_iZ - BestHRightPoint.m_iZ;

		double width = std::sqrt(Xdiff* Xdiff + Ydiff* Ydiff + Zdiff* Zdiff);

		// the following two lines were added for CANCER images
		float h_width = static_cast<float>(BestHLeftPoint.FindDistance(&BestHRightPoint));
		float v_width = static_cast<float>(BestVLeftPoint.FindDistance(&BestVRightPoint));
		gfHWidth += h_width;
		gfVWidth += v_width;
		aPoint->m_fHWidth = h_width;
		aPoint->m_fVWidth = v_width;
		gfWidthSum += static_cast<float>(width);
		giNumOfWidthSumMembers++;
	}
	else
		result = false;

	// 1. if the point is a valid point, update the foreground and background histograms
	// 2. add the point to the list of seed points
	if (result == true)
	{
		deque<CPoint> points; 
		points.push_back(*aPoint);
		points.push_back(BestHRightPoint);
		points.push_back(BestHLeftPoint);
		points.push_back(BestVRightPoint);
		points.push_back(BestVLeftPoint);
		Seeds.push_back(points);

		gaiForegroundHistogram[BestForegroundEstimate]++;
		gaiBackgroundHistogram[BestBackgroundEstimate]++;
	}

	delete [] Response;
	delete [] DirVotes;

	return result;
}

