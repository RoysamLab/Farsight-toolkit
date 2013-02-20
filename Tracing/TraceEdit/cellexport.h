//QDialog for Saving Cell Export files (SWC, JPG)
#ifndef CELLEXPORT_H
#define CELLEXPORT_H

#include <iostream>
#include <QtGui>
#include <QApplication>

/**
 * @author Audrey Cheong
 */
class SaveCellExportDialog : public QDialog
{
	Q_OBJECT
	
public:
	SaveCellExportDialog(QWidget* parent, QString curdirectoryswc, QString curdirectoryjpg, QString swcfileName, QString jpgfileName, bool changeswcfileName, bool changejpgfileName);

	QString getSWCDir();
	QString getJPGDir();
	QString getSWCfileName();
	QString getJPGfileName();
	bool differentSWCfileName();
	bool differentJPGfileName();
	bool getSave();
	
private slots:	
	void swcBrowse();
	void jpgBrowse();
	void swcfilenaming();
	void jpgfilenaming();
	void save();

private:
	QPushButton *createButton(const QString &text, const char *member);
	QComboBox *createComboBox(const QString &text = QString());

	QGroupBox *saveSWCGroupBox, *saveJPGGroupBox;
//set up for swc files
	QComboBox *swcdirectoryComboBox;
	QPushButton *swcbrowseButton, *swcmoreButton;
	QWidget *swcextension;
	QRadioButton *originalswcfileNameButton, *renumberswcfileNameButton, *renameswcfileNameButton;
	QLineEdit *nameswcfileNameLine;
//set up for jpg files
	QComboBox *jpgdirectoryComboBox;
	QPushButton *jpgbrowseButton, *jpgmoreButton;
	QWidget* jpgextension;
	QRadioButton *originaljpgfileNameButton, *renumberjpgfileNameButton, *renamejpgfileNameButton;
	QLineEdit *namejpgfileNameLine;

	QPushButton *OkButton, *CancelButton;

	//variables to return to the caller
	QString curdirectoryswc, curdirectoryjpg, swcfileName, jpgfileName;
	bool changeswcfileName, changejpgfileName;
	bool saveclicked;
	
};

#endif