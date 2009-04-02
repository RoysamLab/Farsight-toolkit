#include <stdio.h>
#include "TrackingImageView.h"
#include "TrackingVolumeView.h"
#include "TrackingKymoView.h"
#include "helpers.h"
#include "QApplication.h"

#pragma warning(disable: 4996)

std::vector<std::string> getLabelFileNames(int n)
{
	std::vector<std::string> ret;
	for(int counter=1; counter< n; counter++)
	{
		char buff[1024];
		sprintf(buff,"D:/ucb dataset/output/ena/labeled_tracks/dataset5/labeled_dataset5_%d_w1.tif",counter);
		ret.push_back(buff);
	}
	return ret;
}

std::vector<std::string> getRawFileNames(int n)
{
	std::vector<std::string> ret;
	for(int counter=1; counter< n; counter++)
	{
		char buff[1024];
		sprintf(buff,"D:/ucb dataset/output/ena/cache/dataset5/unmixed_dataset5_%d_1.tif",counter);
		ret.push_back(buff);
	}
	return ret;
}
int main(int argc, char**argv)
{
	QApplication app(argc,argv);
	TrackingModel * trmodel = new TrackingModel();
	int n = 10;
	std::vector<std::string> label_filenames = getLabelFileNames(n);
	std::vector<std::string> raw_filenames = getRawFileNames(n);
	trmodel->SetLabelImageFileNames(label_filenames);
	trmodel->SetRawImageFileNames(raw_filenames);
	//trmodel->SetTracksFileName("D:/ucb dataset/output/ena/new_feature_files/tracks_dataset3_w1.tks");
	DEBUG3("About to start TrackingImageView\n");
	trmodel->ReadRawImages();
	trmodel->ReadLabelImages();
	trmodel->GenerateFeatures();
	TrackingImageView * imview = new TrackingImageView(trmodel);
	//TrackingImageView * imview1 = new TrackingImageView(trmodel);
	TrackingVolumeView * volview = new TrackingVolumeView(trmodel);
	TrackingKymoView * kymoview = new TrackingKymoView(trmodel);
	TrackSelection *select = new TrackSelection(n);
	imview->SetSelection(select);
//	imview1->SetSelection(select);
	volview->SetSelection(select);
	kymoview->SetSelection(select);
	app.exec();
	return 0;
}