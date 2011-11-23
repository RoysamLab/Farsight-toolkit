#ifndef SPDMAINWINDOW_H
#define SPDMAINWINDOW_H

#include <QWidget>
#include <QLineEdit>
#include <QLabel>
#include <QPushButton>
#include <QListWidget>
#include "SPDAnalysisModel.h"
#include "ftkGUI/GraphWindow.h"
#include "ftkGUI/HistoWindow.h"
#include "ClusClus/HeatmapWindow.h"

class SPDMainWindow : public QWidget
{
    Q_OBJECT

public:
    explicit SPDMainWindow(QWidget *parent = 0);
	void setModels(vtkSmartPointer<vtkTable> table = NULL, ObjectSelection * sels = NULL, ObjectSelection * sels2 = NULL);
    ~SPDMainWindow();
	void GetProgressionTreeOrder(std::vector<long int> &order);

private slots:
    void browse();
    void load();
	void loadTestData();
    void clusterFunction();
	void generateMST();
	void autoProcess();
	void emdFunction();
	void showPSM();
	void showPSMHist();
	void viewProgression();
	//void saveSelectedFeatures();
	void updateSelMod();
	void editThreshold();
	void editPercentage();

private:
	SPDAnalysisModel *SPDModel;
    QLabel *dataFileLabel;

    QLabel *dataFileName;
    QPushButton *browseButton;
    QPushButton *loadButton;
	QPushButton *loadTestButton;
	
    QLabel *featureNumLabel;
    QLabel *featureNum;
    QLabel *sampleNumLabel;
    QLabel *sampleNum;

    QLabel *clusterCoherenceLabel;
    QLineEdit *clusterCoherenceBox;
    QLabel *clusterMergeLabel;
    QLineEdit *clusterMergeBox;
    QPushButton *clusterButton;
	QListWidget *listWidget;
	QPushButton *generateMSTButton;
	QPushButton *autoProcButton;
	QPushButton *emdButton;
	QLineEdit *emdThresBox;
	QLineEdit *emdPercentageBox;
	QLabel *psmLable;
	QLabel *psmPerLable;
    QPushButton *psmButton;
	QPushButton *psmHisButton;
	QLabel *psdtLable;   // progression sample discovery tree
	QLineEdit *psdModuleSelectBox;  // select similar modules
    QPushButton *psdtButton;
	//QPushButton *saveFeatureButton;  // save the selected features to file

    QString FileName;
	GraphWindow *graph;
	Heatmap *heatmap;
	HistoWindow *histo;

	vnl_vector<int> optimalleaforder;
	vnl_vector<int> selMod;

	vtkSmartPointer<vtkTable> data;
	ObjectSelection *selection;
	ObjectSelection *selection2;

};

#endif // SPDMAINWINDOW_H
