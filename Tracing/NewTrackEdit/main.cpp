#include <stdio.h>
#include <vector>
#include "TrackingDataModel.h"
#include <QApplication>
#include <QTableView>

std::vector<std::string> getLabelFileNames(int n)
{
	std::vector<std::string> ret;
	for(int counter=10; counter< 10+n; counter++)
	{
		char buff[1024];
		sprintf(buff,"D:/ucb dataset/output/ena/labeled_tracks_original/dataset1/labeled_dataset1_%d_w1.tif",counter);
		ret.push_back(buff);
	}
	return ret;
}

std::vector<std::string> getRawFileNames(int n)
{
	std::vector<std::string> ret;
	for(int counter=10; counter< 10+n; counter++)
	{
		char buff[1024];
		sprintf(buff,"D:/ucb dataset/output/ena/unmixed_dataset1/unmixed_dataset1_%d_1.tif",counter);
		ret.push_back(buff);
	}
	return ret;
}
int main(int argc, char**argv)
{
	QApplication app(argc,argv);
	TrackingDataModel * trmodel = new TrackingDataModel();
	int n = 5;
	std::vector<std::string> label_filenames = getLabelFileNames(n);
	std::vector<std::string> raw_filenames = getRawFileNames(n);
	trmodel->SetLabelImageFileNames(label_filenames);
	trmodel->SetRawImageFileNames(raw_filenames);
	trmodel->InitializeModel();
	QTableView* tview = new QTableView();
	QItemSelectionModel * selModel = new QItemSelectionModel(trmodel);
	tview->setModel(trmodel);
	tview->setSelectionModel(selModel);
	tview->update();
    trmodel->setColumnCount(1);
    tview->show();
	app.exec();
  delete trmodel;
	return 0;
}
