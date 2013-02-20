#ifndef _TRACKING_DATA_MODEL_H
#define _TRACKING_DATA_MODEL_H

#include "helpers.h"
#include "ftkObject.h"
#include "ftkLabelImageToFeatures.h"
#include <QStandardItemModel>
#include <QObject>
#include <QMap>


#define MAX_TRACKS 15000

class TrackingDataModel: public QStandardItemModel
{
	Q_OBJECT
	public:
		TrackingDataModel(){
			m_num_time_points = -1;
		}
		void SetLabelImageFileNames(std::vector<std::string> &filenames);
		void SetRawImageFileNames(std::vector<std::string> &filenames);
		void SetTracksFileName(const char * filename);
		void InitializeModel();
		//void mergeTracks(TrackSelection *,int);
		//void deleteTracks(TrackSelection *);
		std::vector<std::vector<ftk::Object> > * getObjects(){return &(this->m_objects);}
		
		unsigned int GetSizeX(){ return m_limages[0]->GetLargestPossibleRegion().GetSize()[0];}
		unsigned int GetSizeY(){ return m_limages[0]->GetLargestPossibleRegion().GetSize()[1];}
		unsigned int GetSizeZ(){ return m_limages[0]->GetLargestPossibleRegion().GetSize()[2];}
		unsigned int GetSizeT(){ return m_limages.size();}
		InputImageType::Pointer getRawImagePointer(int ind)
		{
			return m_rimages[ind];
		}
		LabelImageType::Pointer getLabelImagePointer(int ind)
		{
			return m_limages[ind];
		}
		int getNumTimePoints() const
		{
			return m_num_time_points;
		}
		std::vector<int> GetLabelsAlongZ(double, double,double);
		void saveLabels();
	signals:
		//void labelsChanged(int);
	public slots:
		//flags()
		QVariant data(const QModelIndex &index, int role) const;
		//data()
		//headerData()
		//rowCount()
		//columnCount()
		//void UpdateFeaturesForTimePoint(int);
	private:
		//functions
		ftk::Object GetNewTrackPointObject(ftk::IntrinsicFeatures *);
		void GenerateObjects();
		void ReadLabelImages();
		void ReadRawImages();
		void ReadTracks();

		//variables
		typedef ftk::Object ObjectsType;
		std::vector<std::vector<ftk::Object> > m_objects;
		std::vector<std::string> m_lfilenames;// label image filenames
		std::vector<std::string> m_rfilenames;// raw image filenames
		std::string m_tfilename; // tracks filename
		std::vector<LabelImageType::Pointer> m_limages;// label images
		std::vector<InputImageType::Pointer> m_rimages;// raw images
		std::vector<TrackPoint> m_tracks; // tracks
		int m_num_time_points;
};

#endif
