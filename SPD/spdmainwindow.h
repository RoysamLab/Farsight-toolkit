#ifndef SPDMAINWINDOW_H
#define SPDMAINWINDOW_H

#include <QWidget>
#include <QLineEdit>
#include <QLabel>
#include <QPushButton>
#include <QListWidget>
#include "SPDAnalysisModel.h"
#include "ftkGUI/GraphWindow.h"
#include "ClusClus/HeatmapWindow.h"

class SPDMainWindow : public QWidget
{
    Q_OBJECT

public:
    explicit SPDMainWindow(QWidget *parent = 0);
    ~SPDMainWindow();

private slots:
    void browse();
    void load();
	void loadTestData();
    void clusterFunction();
	void generateMST();
	void showMST();
	void emdFunction();

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
	QPushButton *showMSTButton;
	QPushButton *emdButton;
    QString FileName;

	GraphWindow *graph;
	Heatmap *heatmap;
};

#endif // SPDMAINWINDOW_H
