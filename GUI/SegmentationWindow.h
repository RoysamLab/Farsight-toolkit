#ifndef SEGMENTATIONWINDOW_H
#define SEGMENTATIONWINDOW_H

#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QGridLayout>
#include <QtGui/QMainWindow>
#include <QtGui/QSlider>
#include <QtGui/QWidget>
#include <QtGui/QDockWidget>
#include <QtGui/QCheckBox>
#include <QtGui/QGroupBox>
#include <QtGui/QVBoxLayout>
#include <QtGui/QCloseEvent>
#include <QtGui/QPushButton>
#include <QtGui/QSpinBox>
#include <QtGui/QLabel>

#include "SegmentationView.h"


class SegmentationWindow : public QWidget
{
    Q_OBJECT;

public:
	SegmentationWindow(QWidget *parent = 0);
	void SetModels(SegmentationModel *sModel);
	void SetChannelImage(ftk::Image *image);
	void SetLabelImage(ftk::Image *image);
	//void loadImage( QString &fileName );
	//void AddChannelImage( QString &fileName );
	//void AddLabelImage( QString &fileName );

signals:
	void closing(QWidget *widget);
	
protected:
	void showEvent( QShowEvent * event );
	void moveEvent ( QMoveEvent * event );
	void closeEvent(QCloseEvent *event);

private slots:
	void updateChannels(bool){updateChannels();};
    
private:
	//QVector<QImage>* loadImage( QString &fileName );
	void updateVSlider(void);
	void updateHSlider(void);
	void createChannelWindow(ftk::Image *image);
	void closeChannelWindow(void);
	void updateChannels(void);

	SegmentationView *segview;
    QSlider *vSlider;
	QSpinBox *vSpin;
	QLabel *vLabel;
    QSlider *hSlider;
	QSpinBox *hSpin;
	QLabel *hLabel;
	QWidget *channelWidget;
	QCheckBox ** chBoxes;		//pointers to all my checkboxes

	int numZSlices;
	int numTSlices;
	int numChannels;
 };

#endif

