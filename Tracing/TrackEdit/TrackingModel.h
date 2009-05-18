#ifndef _TRACKING_MODEL_H
#define _TRACKING_MODEL_H
#include "helpers.h"
#include "TrackSelection.h"
#include <QObject>

class TrackingModel: public QObject
{
	Q_OBJECT
	public:
		TrackingModel();
		void SetLabelImageFileNames(std::vector<std::string> &filenames);
		void SetRawImageFileNames(std::vector<std::string> &filenames);
		void SetTracksFileName(const char * filename);
		void ReadLabelImages();
		void ReadRawImages();
		void ReadTracks();
		void mergeTracks(TrackSelection *,int);
		void deleteTracks(TrackSelection *);
		std::vector<std::vector<ftk::LabelImageFeatures>> *getFeatures(){return &features;}
		void GenerateFeatures();
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
			return num_time_points;
		}
		std::vector<int> GetLabelsAlongZ(double, double,double);
		void saveLabels();
	signals:
		void labelsChanged(int);
	public slots:
		void UpdateFeaturesForTimePoint(int);
	private:
		typedef ftk::LabelImageFeatures FeatureType;
		std::vector<std::vector<FeatureType>> features;
		std::vector<std::string> m_lfilenames;// label image filenames
		std::vector<std::string> m_rfilenames;// raw image filenames
		std::string m_tfilename; // tracks filename
		std::vector<LabelImageType::Pointer> m_limages;// label images
		std::vector<InputImageType::Pointer> m_rimages;// raw images
		std::vector<TrackPoint> m_tracks; // tracks
		int num_time_points;
};

#endif