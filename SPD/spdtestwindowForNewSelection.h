#ifndef SPDTESTWINDOW_H
#define SPDTESTWINDOW_H

#include <QWidget>
#include <QLineEdit>
#include <QLabel>
#include <QPushButton>
#include <QListWidget>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include "SPDAnalysisModel.h"
//#include "ftkGUI/GraphWindow.h"
#include "ftkGUI/GraphWindowForNewSelection.h"
#include "ftkGUI/HistoWindow.h"
#include "ftkGUI/PlotWindow.h"
#include "ProgressionHeatmapWindow.h"
#include "HeatmapWindowForNewSelection.h"
#include "ftkGUI/SelectiveClustering.h"

class SPDWindowForNewSelection : public QWidget
{
    Q_OBJECT

public:
    explicit SPDWindowForNewSelection(QWidget *parent = 0);
	void setModels(vtkSmartPointer<vtkTable> table = NULL, ObjectSelection * sels = NULL, SelectiveClustering * clusterSelection = NULL, ObjectSelection * sels2 = NULL);
    ~SPDWindowForNewSelection();
	void GetProgressionTreeOrder(std::vector<long int> &order);

protected:
	void split(std::string& s, char delim, std::vector< unsigned int>& indexVec);
	void GetFeatureOrder(std::vector< unsigned int> &selID, std::vector< int> &selIdOrder, std::vector< int> &unselIdOrder);
	bool IsExist(std::vector< unsigned int> vec, unsigned int value);
	virtual void closeEvent(QCloseEvent *event);
	void closeSubWindows();
	vtkSmartPointer<vtkTable> GetSubTableExcludeItems( vtkSmartPointer<vtkTable> table, std::set<long int> &selItems);

protected slots:
    void browse();
    void load();
	void loadContrastData();
    void clusterFunction();
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
	void regenerateProgressionTree();
	void updateProgressionType();
	void ReRunSPDAnlysis();
	void ReColorProgressionTree(int nfeature);

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
    QLabel *kNearestNeighborLabel;
    QDoubleSpinBox *kNearestNeighborBox;
    QPushButton *clusterButton;
	QPushButton *cellClusterButton;
	QListWidget *listWidget;

	QLabel *emdLabel;
	QLabel *progressionOverDistance;
	QCheckBox *bcheckBox;   // progression overall or over distance to device
	QPushButton *emdButton;
	QDoubleSpinBox *emdThresBox;
	QLineEdit *emdPercentageBox;
	QLabel *psmLable;
	QLabel *psmPerLable;
    QPushButton *psmButton;
	QPushButton *psmHisButton;
	QLabel *psdtLable;   // progression sample discovery tree
	QLineEdit *psdModuleSelectBox;  // select similar modules
	QLabel *maxVetexIdLabel;  // max id to seperate the data
	QSpinBox *maxVetexIdEdit; // max id
    QPushButton *psdtButton;
	QLabel *heatmapLabel;
	QPushButton *heatmapButton;  // show progression heatmap  // now shows the progression over distance to device
	QDoubleSpinBox *distanceThres;  // distance threshold for calculating percentage
   
	QString FileName;
	GraphWindowForNewSelection *graph;
	ProgressionHeatmap *simHeatmap;
	HistoWindow *histo;

	HeatmapForNewSelection *progressionHeatmap;
	HeatmapForNewSelection *HeatmapWin;
	PlotWindow *plot;

	vnl_vector<int> optimalleaforder;
	vnl_vector<int> selMod;

	vtkSmartPointer<vtkTable> data;
	std::set<long int> excludedIds;
	ObjectSelection *selection;  /** selection for threshold and IDs */
	ObjectSelection *selection2; 
	SelectiveClustering * ClusterSelections;

	std::vector< double> sampleVec;
	std::vector< double> percentVec;
	std::vector< unsigned int> selFeatureID;
	std::vector< int> selOrder;
	std::vector< int> unselOrder;


};

#endif // SPDTESTWINDOW_H
