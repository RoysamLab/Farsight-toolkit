#ifndef SPDMAINWINDOW_H
#define SPDMAINWINDOW_H

#include <QWidget>
#include <QLineEdit>
#include <QLabel>
#include <QPushButton>
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
    void showResult();
	void generateMST();
	void showMST();

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
    QPushButton *clusterResultButton;

	QPushButton *generateMSTButton;
	QLabel *mstState;
	QPushButton *showMSTButton;

    QString FileName;
};

#endif // SPDMAINWINDOW_H
