//QDialog for Saving screenshots
#ifndef SCREENSHOT_H
#define SCREENSHOT_H

#include <iostream>
#include <QtGui>
#include <QApplication>

/**
 * @author Audrey Cheong
 */
class ScreenShotDialog : public QDialog
{
	Q_OBJECT
	
public:
	ScreenShotDialog(QWidget* parent, QString fileName, QString imageDir);

	QString getfileName();
	QString getDir();
	int getMagnification();
	bool getSave();
	bool getBaseline();
	
private slots:	
	void Browse();
	void save();

private:
	QComboBox * directoryComboBox;
	QLineEdit * fileNameLine;
	QPushButton * browseButton;
	QSpinBox * ZoomSpinBox;
	QPushButton *OkButton, *CancelButton;
	QCheckBox *BaselineBox;
	
	QPushButton *createButton(const QString &text, const char *member);
	QComboBox *createComboBox(const QString &text = QString());

	bool saveclicked;
	QString imageDir;
//variables to return to the caller
	QString curdirectory;
	QString fileName;
	int magnifynum;
	
};

#endif