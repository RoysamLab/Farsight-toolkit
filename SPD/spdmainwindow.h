#ifndef SPDMAINWINDOW_H
#define SPDMAINWINDOW_H

#include <QWidget>
#include <QLineEdit>
#include <QLabel>
#include <QPushButton>
#include <QListWidget>
#include <QDoubleSpinBox>
#include "SPDAnalysisModel.h"
#include "ftkGUI/GraphWindow.h"
#include "ftkGUI/HistoWindow.h"
#include "ClusClus/ProgressionHeatmapWindow.h"

class SPDMainWindow : public QWidget
{
    Q_OBJECT

public:
    explicit SPDMainWindow(QWidget *parent = 0);
	void setModels(vtkSmartPointer<vtkTable> table = NULL, ObjectSelection * sels = NULL, ObjectSelection * sels2 = NULL);
    ~SPDMainWindow();
	void GetProgressionTreeOrder(std::vector<long int> &order);

protected:
	virtual void closeEvent(QCloseEvent *event);

private slots:
    void browse();
    void load();
	void loadTestData();
    void clusterFunction();
	void generateMST();
	void clusterCells();
	void emdFunction();
	void showPSM();
	void showPSMHist();
	void viewProgression();
	//void saveSelectedFeatures();
	void updateSelMod();
	void editThreshold();
	void editPercentage();
	void showProgressionHeatmap();

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
    QDoubleSpinBox *clusterCoherenceBox;
	QLabel *sampleCoherenceLabel;
	QDoubleSpinBox *sampleCoherenceBox;
    QLabel *clusterMergeLabel;
    QDoubleSpinBox *clusterMergeBox;
    QPushButton *clusterButton;
	QPushButton *cellClusterButton;
	QListWidget *listWidget;
	QLabel *mstLabel;
	QPushButton *generateMSTButton;
	QLabel *emdLabel;
	QPushButton *emdButton;
	QDoubleSpinBox *emdThresBox;
	QLineEdit *emdPercentageBox;
	QLabel *psmLable;
	QLabel *psmPerLable;
    QPushButton *psmButton;
	QPushButton *psmHisButton;
	QLabel *psdtLable;   // progression sample discovery tree
	QLineEdit *psdModuleSelectBox;  // select similar modules
    QPushButton *psdtButton;
	QLabel *heatmapLabel;
	QPushButton *heatmapButton;  // show progression heatmap

    QString FileName;
	GraphWindow *graph;
	ProgressionHeatmap *simHeatmap;
	HistoWindow *histo;
	ProgressionHeatmap *progressionHeatmap;

	vnl_vector<int> optimalleaforder;
	vnl_vector<int> selMod;

	vtkSmartPointer<vtkTable> data;
	ObjectSelection *selection;
	ObjectSelection *selection2;

};

#endif // SPDMAINWINDOW_H
