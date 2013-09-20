#ifndef SPDMSTMODULEMATCH_H
#define SPDMSTMODULEMATCH_H

#include <QWidget>
#include <QLineEdit>
#include <QLabel>
#include <QPushButton>
#include <QListWidget>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QCheckBox>
#include "SPDAnalysisModel.h"
#include "ftkGUI/GraphWindow.h"
#include "ftkGUI/HistoWindow.h"
#include "ftkGUI/PlotWindow.h"
#include "TrendHeatmapWindow.h"
#include "HeatmapWindow.h"

class SPDMSTModuleMatch : public QWidget
{
    Q_OBJECT

public:
    SPDMSTModuleMatch(QWidget *parent = 0);
	void setModels(vtkSmartPointer<vtkTable> table = NULL, ObjectSelection * sels = NULL, ObjectSelection * sels2 = NULL);
    ~SPDMSTModuleMatch();
	void GetTrendTreeOrder(std::vector<long int> &order);
	vtkSmartPointer<vtkTable> NormalizeTable(vtkSmartPointer<vtkTable> table);

protected:
	void split(std::string& s, char delim, std::vector< unsigned int>& indexVec);
	void GetFeatureOrder(std::vector< unsigned int> &selID, std::vector< int> &selIdOrder, std::vector< int> &unselIdOrder);
	void showHeatmapAfterFeatureClustering();
	bool IsExist(std::vector< unsigned int> vec, unsigned int value);
	virtual void closeEvent(QCloseEvent *event);
	void closeSubWindows();
	vtkSmartPointer<vtkTable> GetSubTableExcludeItems( vtkSmartPointer<vtkTable> table, std::set<long int> &selItems);

protected slots:
    void browse();
    void load();
	void loadContrastData();
	void showOriginalHeatmap();
    void clusterFunction();
	void emdFunction();
	void showPSM();
	void viewTrend();
	void updateSelMod();
	void editThreshold();
	void showTrendHeatmap();
	void regenerateTrendTree();
	void ReRunSPDAnlysis();
	void ReColorTrendTree(int nfeature);

private:
	SPDAnalysisModel *SPDModel;
    QLabel *dataFileLabel;

    QLabel *dataFileName;
    QPushButton *browseButton;
    QPushButton *loadButton;
	//QPushButton *loadTestButton;
	
    QLabel *featureNumLabel;
    QLabel *featureNum;
    QLabel *sampleNumLabel;
    QLabel *sampleNum;

    QLabel *clusterCoherenceLabel;
    QDoubleSpinBox *clusterCoherenceBox;
	QLabel *nBinLabel;
    QSpinBox *nBinBox;
    //QPushButton *clusterButton;

	QLabel *emdLabel;
	//QLabel *progressionOverDistance;
	//QCheckBox *bcheckBox;   // Trend overall or over distance to device
	QPushButton *emdButton;
	QDoubleSpinBox *emdThresBox;
	QLineEdit *emdPercentageBox;
	QLabel *psmLable;
	QLabel *psmPerLable;
    QPushButton *psmButton;
	QLabel *psdtLable;   // Trend sample discovery tree
	QLineEdit *psdModuleSelectBox;  // select similar modules
    QPushButton *psdtButton;
	QPushButton *heatmapButton;  // show Trend heatmap  // now shows the Trend over distance to device
	QLabel *continSelectLabel;   
	QCheckBox *continSelectCheck;
	QString FileName;
	GraphWindow *graph;
	TrendHeatmap *simHeatmap;
	HistoWindow *histo;
	Heatmap *originalHeatmap;
	Heatmap *progressionHeatmap;
	Heatmap *HeatmapWin;
	PlotWindow *plot;

	vnl_vector<int> optimalleaforder;

	vtkSmartPointer<vtkTable> data;
	std::set<long int> excludedIds;
	ObjectSelection *selection;  /** selection for threshold and IDs */
	ObjectSelection *selection2; 

	std::vector< unsigned int> selFeatureID;
	std::vector< int> selOrder;
	std::vector< int> unselOrder;
	vtkSmartPointer<vtkTable> tableAfterCellCluster;
	std::map< int, int> indexMap;
	std::vector<int> connectedComponent;
	int connectedNum;
	bool bconnected;
	std::set<unsigned int> selMod;
};

#endif // SPDMSTMODULEMATCH_H
