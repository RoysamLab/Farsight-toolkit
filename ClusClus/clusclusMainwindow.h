#ifndef CLUSCLUSMAINWINDOW_H
#define CLUSCLUSMAINWINDOW_H

#include <QWidget>
#include <QLineEdit>
#include <QLabel>
#include <QPushButton>
#include <QButtonGroup>
#include <QRadioButton>
#include "clusclus.h"
#include "clusgap.h"
#include "Dendrogram.h"

class ClusClusMainWindow : public QWidget
{
    Q_OBJECT

public:
	ClusClusMainWindow(QWidget *parent = 0);
	~ClusClusMainWindow();

private:
    QLabel *dataFileLabel;

    QLabel *dataFileName;
    QPushButton *browseButton;
    QPushButton *loadButton;
	
	QLabel *featureNumLabel;
    QLabel *featureNum;
    QLabel *sampleNumLabel;
    QLabel *sampleNum;
	QLabel *rclusterNumLabel;
    QLabel *rclusterNum;
	QLabel *cclusterNumLabel;
    QLabel *cclusterNum;

    QLabel *rtrialsNumLable;
    QLineEdit *rtrialsNumBox;
    QLabel *rgapsNumLable;
    QLineEdit *rgapsNumBox;

	QLabel *ctrialsNumLable;
    QLineEdit *ctrialsNumBox;
    QLabel *cgapsNumLable;
    QLineEdit *cgapsNumBox;


    QPushButton *clusterButton;
    QPushButton *biclusterButton;
	QPushButton *rcomputeGapButton;
	QPushButton *ccomputeGapButton;
	QPushButton *rgenerateDendroButton;
	QPushButton *cgenerateDendroButton;

	QRadioButton *singleradio;
    QRadioButton *averageradio;
	QRadioButton *completeradio;
    QButtonGroup *radiogrp;

	QString FileName;

private slots:
    void browse();
    void load();
    void runcluster();
    void runbicluster();
	void samplecomputegap();
	void featurecomputegap();
	void determinlinkmode(int id);
	void generatedendrogram1();
	void generatedendrogram2();

private:
	clusclus *cc1;
	clusclus *cc2;
	clusgap  *cg1;
	clusgap  *cg2;
};
#endif // CLUSCLUSMAINWINDOW_H