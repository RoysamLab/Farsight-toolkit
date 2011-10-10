#ifndef SPDMAINWINDOW_H
#define SPDMAINWINDOW_H

#include <QWidget>
#include <QLineEdit>
#include <QLabel>
#include <QPushButton>
#include <QListWidget>
#include "SPDAnalysisModel.h"

class SPDMainWindow : public QWidget
{
    Q_OBJECT

public:
    explicit SPDMainWindow(QWidget *parent = 0);
    ~SPDMainWindow();

private slots:
    void browse();
    void load();
    void clusterFunction();
	void generateMST();
	void showMST();
	void AddClusterModuleToList();

private:
	SPDAnalysisModel *SPDModel;

private:
    QLabel *dataFileLabel;

    QLabel *dataFileName;
    QPushButton *browseButton;
    QPushButton *loadButton;

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

    QString FileName;
};

#endif // SPDMAINWINDOW_H
