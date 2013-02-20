#include "TrackingDataModel.h"

void TrackingDataModel::InitializeModel()
{
	ReadLabelImages();
	ReadRawImages();
	GenerateObjects();
	printf("Finished initializing the TrackingDataModel\n");
	//FIXME: more things to be added
}
void TrackingDataModel::SetLabelImageFileNames(std::vector<std::string> &filenames)
{
	if(m_num_time_points==-1)
		m_num_time_points = filenames.size();
	else
		assert(m_num_time_points == (int)filenames.size());
	this->m_lfilenames = filenames; // element-by-element copy
}

void TrackingDataModel::SetRawImageFileNames(std::vector<std::string> &filenames)
{
	if(m_num_time_points==-1)
		m_num_time_points = filenames.size();
	else
		assert(m_num_time_points == (int)filenames.size());
	this->m_rfilenames = filenames; // element-by-element copy
}

void TrackingDataModel::SetTracksFileName(const char * filename)
{
	this->m_tfilename = filename;
}

void TrackingDataModel::ReadRawImages()
{
	for(unsigned int counter=0; counter<m_rfilenames.size(); counter++)
	{
		m_rimages.push_back(readImage<InputImageType>(m_rfilenames[counter].c_str()));
	}
	//DEBUG3("Read %d raw files\n", m_rfilenames.size());
}

void TrackingDataModel::ReadLabelImages()
{
	for(unsigned int counter=0; counter<m_lfilenames.size(); counter++)
	{
		m_limages.push_back(readImage<LabelImageType>(m_lfilenames[counter].c_str()));
	}
}

void TrackingDataModel::ReadTracks()
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

std::vector<int> TrackingDataModel::GetLabelsAlongZ(double x_arg, double y_arg , double t_arg)
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
	
	std::sort(labels.begin(),labels.end());
	std::vector<int>::iterator new_last = std::unique(labels.begin(),labels.end());
//	printf("About to delete duplicate elements\n");
	labels.erase(new_last,labels.end());
	
	return labels;
}

void TrackingDataModel::GenerateObjects()
{
	for(int counter=0; counter<m_num_time_points; counter++)
	{
		std::vector<ftk::Object> tobjects;
		if(m_limages[counter]->GetLargestPossibleRegion().GetSize()[2] > 1)
		{
			typedef ftk::LabelImageToFeatures<unsigned char, short, 3> FeatureCalcType;
			FeatureCalcType::Pointer labfilter  = FeatureCalcType::New();
			labfilter->SetImageInputs(getRawImagePointer(counter),getLabelImagePointer(counter));
			labfilter->ComputeHistogramOn();
			labfilter->SetLevel(3);
			labfilter->ComputeTexturesOn();
			labfilter->Update();

			std::vector<short> labels = labfilter->GetLabels();
			for(unsigned int lco = 0; lco < labels.size(); lco++)
			{
				ftk::IntrinsicFeatures *features = labfilter->GetFeatures(labels[lco]);
				features->num = labels[lco];
				features->time = counter;
				features->tag = 0;//FIXME
				tobjects.push_back(GetNewTrackPointObject(labfilter->GetFeatures(labels[lco])));
			}
		}
		m_objects.push_back(tobjects);
	}
}

ftk::Object TrackingDataModel::GetNewTrackPointObject(ftk::IntrinsicFeatures * features)
{
	std::string name = "TrackPoint";
	ftk::Object object(name);
	object.SetId(features->time*MAX_TRACKS+features->num);
	object.SetValidity(1);
	object.SetDuplicated(0);
	object.SetClass(-1);

	if(features == NULL)
		return object;

	ftk::Object::Point c;
	c.x = (int)features->Centroid[0];
	c.y = (int)features->Centroid[1];
	c.z = (int)features->Centroid[2];
	c.t = 0;
	object.AddCenter(c);

	ftk::Object::Box b;
	b.min.x = (int)features->BoundingBox[0];
	b.max.x = (int)features->BoundingBox[1];
	b.min.y = (int)features->BoundingBox[2];
	b.max.y = (int)features->BoundingBox[3];
	b.min.z = (int)features->BoundingBox[4];
	b.max.z = (int)features->BoundingBox[5];
	b.min.t = 0;
	b.max.t = 0;
	object.AddBound(b);

	std::vector< float > f(0);
	for (int i=0; i< ftk::IntrinsicFeatures::N; ++i)
	{
		f.push_back( features->ScalarFeatures[i] );
	}

	object.SetFeatures( f );

	return object;
}


QVariant TrackingDataModel::data(const QModelIndex &index, int role) const
{
	return QVariant("1");
}
