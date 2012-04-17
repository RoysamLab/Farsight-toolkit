#ifndef SPDMAINWINDOW_H
#define SPDMAINWINDOW_H

#include <QWidget>
#include <QLineEdit>
#include <QLabel>
#include <QPushButton>
#include <QListWidget>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include "SPDAnalysisModel.h"
#include "ftkGUI/GraphWindow.h"
#include "ftkGUI/HistoWindow.h"
#include "ProgressionHeatmapWindow.h"
#include "ClusClus/HeatmapWindow.h"

class SPDMainWindow : public QWidget
{
    Q_OBJECT

public:
    explicit SPDMainWindow(QWidget *parent = 0);
	void setModels(vtkSmartPointer<vtkTable> table = NULL, ObjectSelection * sels = NULL, ObjectSelection * sels2 = NULL);
    ~SPDMainWindow();
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
	void regenerateProgressionTree();
	void updateProgressionType();
	void ReRunSPDAnlysis();

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
    QPushButton *psdtButton;
	QLabel *heatmapLabel;
	QPushButton *heatmapButton;  // show progression heatmap  // now shows the progression over distance to device

    QString FileName;
	GraphWindow *graph;
	ProgressionHeatmap *simHeatmap;
	HistoWindow *histo;
	Heatmap *progressionHeatmap;
	Heatmap *HeatmapWin;

	vnl_vector<int> optimalleaforder;
	vnl_vector<int> selMod;

	vtkSmartPointer<vtkTable> data;
	std::set<long int> excludedIds;
	ObjectSelection *selection;  /** selection for threshold and IDs */
	ObjectSelection *selection2; 

	std::vector< unsigned int> selFeatureID;
	std::vector< int> selOrder;
	std::vector< int> unselOrder;
};

#endif // SPDMAINWINDOW_H
