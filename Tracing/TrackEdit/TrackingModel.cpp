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

#include "TrackingModel.h"

TrackingModel::TrackingModel()
{
	num_time_points = -1;
	//connect(this,SIGNAL(labelsChanged(int)),this,SLOT(UpdateFeaturesForTimePoint(int)));
}
void TrackingModel::saveLabels()
{
	for(unsigned int counter=0; counter<m_lfilenames.size(); counter++)
	{
		writeImage<LabelImageType>(m_limages[counter],m_lfilenames[counter].c_str());
	}
}
void TrackingModel::SetLabelImageFileNames(std::vector<std::string> &filenames)
{
	if(num_time_points==-1)
		num_time_points = filenames.size();
	else
		assert(num_time_points == filenames.size());
	this->m_lfilenames = filenames; // element-by-element copy
}

void TrackingModel::SetRawImageFileNames(std::vector<std::string> &filenames)
{
	if(num_time_points==-1)
		num_time_points = filenames.size();
	else
		assert(num_time_points == filenames.size());
	this->m_rfilenames = filenames; // element-by-element copy
}

void TrackingModel::SetTracksFileName(const char * filename)
{
	this->m_tfilename = filename;
}

void TrackingModel::ReadLabelImages()
{
	for(unsigned int counter=0; counter<m_lfilenames.size(); counter++)
	{
		m_limages.push_back(readImage<LabelImageType>(m_lfilenames[counter].c_str()));
	}
}

void TrackingModel::GenerateFeatures()
{
	features.clear();
	for(unsigned int counter=0; counter< m_limages.size(); counter++)
	{
		std::vector<FeatureType> f;
		getFeatureVectorsFarsight(m_limages[counter],m_rimages[counter],f,counter,0);
		features.push_back(f);
	}
}
void TrackingModel::UpdateFeaturesForTimePoint(int t)
{
	features[t].clear();
	getFeatureVectorsFarsight(m_limages[t],m_rimages[t],features[t],t,0);
	for(unsigned int counter=0; counter<features[t].size(); counter++)
	{
		features[t][counter].Fprintf();
	}
}
void TrackingModel::ReadRawImages()
{
	for(unsigned int counter=0; counter<m_rfilenames.size(); counter++)
	{
		m_rimages.push_back(readImage<InputImageType>(m_rfilenames[counter].c_str()));
	}
	DEBUG3("Read %d raw files\n", m_rfilenames.size());
}

void TrackingModel::ReadTracks()
{
	FILE *fp = fopen(m_tfilename.c_str(),"r");
	if(fp==NULL)
	{
		printf("Could not open %s for reading\n",m_tfilename.c_str());
		return;
	}
	TrackPoint tr;
	while(1)
	{
		fscanf(fp,"%d %d %lf %lf %lf",&tr.id,&tr.t,&tr.x,&tr.y,&tr.z);
		if(feof(fp))
			break;
		m_tracks.push_back(tr);
	}
}

std::vector<int> TrackingModel::GetLabelsAlongZ(double x_arg, double y_arg , double t_arg)
{
	int t = int(t_arg+0.5);
	int x = int(x_arg+0.5);
	int y = int(y_arg+0.5);
	
	LabelImageType::SizeType size = m_limages[t]->GetLargestPossibleRegion().GetSize();

	x = CLAMP(x,0,static_cast<int>(size[0]-1));
	y = CLAMP(y,0,static_cast<int>(size[1]-1));
	t = CLAMP(t,0,static_cast<int>(m_limages.size()-1));

//	printf("x = %d y = %d t = %d\n",x,y,t);

	std::vector<int> labels;
	
	LabelImageType::IndexType index;
	index[0]=x;
	index[1]=y;
	unsigned short last = 0;
	for(unsigned int counter=0; counter<size[2]; counter++)
	{
		index[2] = counter;
		unsigned short l = m_limages[t]->GetPixel(index);
//		printf("%d ",static_cast<int>(l));
		if(l!=0)
		{
			labels.push_back(l);
			last = l;
		}
	}
//	printf("\nAbout to sort\n");
	
	sort(labels.begin(),labels.end());
	std::vector<int>::iterator new_last = std::unique(labels.begin(),labels.end());
//	printf("About to delete duplicate elements\n");
	labels.erase(new_last,labels.end());
	
	return labels;
}

void TrackingModel::mergeTracks(TrackSelection *select, int trackid)
{
	DEBUG3("Entering mergeTracks\n");
	for(int counter=0; counter< num_time_points; counter++)
	{
		std::set<int> tracks_at_t = select->get(counter);
		std::set<int>::iterator iter = tracks_at_t.begin();
		int n1,n2,n3,n4,n5;
		bool flag = false;
		if(tracks_at_t.size()>0)
		{
			switch(tracks_at_t.size())
			{
			case 1:
				n1 = *iter;
				flag = mergeLabels1(m_limages[counter],trackid,n1);
				break;
			case 2:
				n1 = *iter;++iter;
				n2 = *iter;
				flag = mergeLabels2(m_limages[counter],trackid,n1,n2);
				break;
			case 3:
				n1 = *iter;++iter;
				n2 = *iter;++iter;
				n3 = *iter;
				flag = mergeLabels3(m_limages[counter],trackid,n1,n2,n3);
				break;
			case 4:
				n1 = *iter;++iter;
				n2 = *iter;++iter;
				n3 = *iter;++iter;
				n4 = *iter;
				flag = mergeLabels4(m_limages[counter],trackid,n1,n2,n3,n4);
				break;
			case 5:
				n1 = *iter;++iter;
				n2 = *iter;++iter;
				n3 = *iter;++iter;
				n4 = *iter;++iter;
				n5 = *iter;
				flag = mergeLabels5(m_limages[counter],trackid,n1,n2,n3,n4,n5);
				break;
			default:
				printf("Can't merge more than 5 labels at a timepoint, as of now, for efficiency(read laziness)\n");
				break;
			}
			if(flag)
			{
				UpdateFeaturesForTimePoint(counter);
				emit labelsChanged(counter);
			}
		}
	}
	DEBUG3("Leaving mergeTracks\n");
}
void TrackingModel::deleteTracks(TrackSelection *select)
{
	DEBUG3("Entered deleteTracks\n");
	for(int counter=0; counter< num_time_points; counter++)
	{
		std::set<int> tracks_at_t = select->get(counter);
		std::set<int>::iterator iter = tracks_at_t.begin();
		bool flag = false;
		int n1,n2,n3,n4,n5;
		if(tracks_at_t.size()>0) // as of now a useless check
		{
			switch(tracks_at_t.size())
			{
			case 1:
				n1 = *iter;
				flag = deleteLabels1(m_limages[counter],n1);
				break;
			case 2:
				n1 = *iter;++iter;
				n2 = *iter;
				flag = deleteLabels2(m_limages[counter],n1,n2);
				break;
			case 3:
				n1 = *iter;++iter;
				n2 = *iter;++iter;
				n3 = *iter;
				flag = deleteLabels3(m_limages[counter],n1,n2,n3);
				break;
			case 4:
				n1 = *iter;++iter;
				n2 = *iter;++iter;
				n3 = *iter;++iter;
				n4 = *iter;
				flag = deleteLabels4(m_limages[counter],n1,n2,n3,n4);
				break;
			case 5:
				n1 = *iter;++iter;
				n2 = *iter;++iter;
				n3 = *iter;++iter;
				n4 = *iter;++iter;
				n5 = *iter;
				flag = deleteLabels5(m_limages[counter],n1,n2,n3,n4,n5);
				break;
			default:
				printf("Can't delete more than 5 labels, as of now, for efficiency(read laziness)\n");
				break;
			}
			if(flag)
			{
				UpdateFeaturesForTimePoint(counter);
				emit labelsChanged(counter);
			}
		}
	}
	DEBUG3("Leaving deleteTracks\n");
}
