//QDialog for Saving screenshots
#ifndef SCREENSHOT_H
#define SCREENSHOT_H

#include <iostream>
#include <QtGui>
#include <QApplication>

class ScreenShotDialog : public QDialog
{
	Q_OBJECT
	
public:
	ScreenShotDialog(QWidget* parent, QString fileName, QString imageDir);

	QString getfileName();
	QString getDir();
	int getMagnification();
	//QString getSWCfileName();
	//QString getJPGfileName();
	//bool keeporiginalSWCfileName();
	//bool keeporiginalJPGfileName();
	
private slots:	
	void Browse();
	void save();

private:
	QComboBox * directoryComboBox;
	QSpinBox * ZoomSpinBox;
	QPushButton *createButton(const QString &text, const char *member);
	QComboBox *createComboBox(const QString &text = QString());

	QLineEdit * fileNameLine;
//
//	QGroupBox *saveSWCGroupBox;
//	QGroupBox *saveJPGGroupBox;
////set up for swc files
//	QComboBox *swcdirectoryComboBox;
//	QPushButton *swcbrowseButton;
//	QPushButton *swcmoreButton;
//	QWidget *swcextension;
//	QRadioButton *originalswcfileNameButton, *renumberswcfileNameButton, *renameswcfileNameButton;
//	QLineEdit *nameswcfileNameLine;
////set up for jpg files
//	QComboBox *jpgdirectoryComboBox;
//	QPushButton *jpgbrowseButton;
//	QPushButton *jpgmoreButton;
//	QWidget* jpgextension;
//	QRadioButton *originaljpgfileNameButton, *renumberjpgfileNameButton, *renamejpgfileNameButton;
//	QLineEdit *namejpgfileNameLine;\
//
	QPushButton * browseButton;
	QPushButton *OkButton;
	QPushButton *CancelButton;
//
//	//variables to return to the caller
	QString curdirectory;
	QString imageDir;
//	QString curdirectoryjpg;
//	QString swcfileName;
//	QString jpgfileName;
//	bool changeswcfileName;
//	bool changejpgfileName;
	QString fileName;
	int magnifynum;
	
};

#endif