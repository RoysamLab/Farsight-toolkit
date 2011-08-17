#ifndef TrackingDialog_H
#define TrackingDialog_H

#include <QtCore/QString>
#include <QtGui/QDialog>
//#include <QtGui/QDir>
#include <QtGui/QFileDialog>
#include <QtGui/QVBoxLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QLabel>
#include <QtGui/QDialogButtonBox>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QPushButton>
#include <QtGui/QLineEdit>
#include <QtGui/QWidget>
#include <QtGui/QTabWidget>
#include <QtGui/QFrame>
#include <QtCore/QFileInfo>
#include <string>
#include <utility>
#include <vector>


class ParametersTab : public QWidget
 {
     Q_OBJECT

 public:
     ParametersTab( QWidget *parent = 0);
	 // Spin Boxes:
	QDoubleSpinBox * xspacingSpinBox;
	QDoubleSpinBox * yspacingSpinBox;
	QDoubleSpinBox * zspacingSpinBox;
	QDoubleSpinBox * distvarSpinBox;
	QDoubleSpinBox * distmeanSpinBox;
	QDoubleSpinBox * timevarSpinBox;
	QDoubleSpinBox * timemeanSpinBox;
	QDoubleSpinBox * overlapvarSpinBox;
	QDoubleSpinBox * overlapmeanSpinBox;
	QDoubleSpinBox * volumevarSpinBox;
	QDoubleSpinBox * mergesplitpriorSpinBox;
	QDoubleSpinBox * appeardisappearpriorSpinBox;
	QDoubleSpinBox * translationpriorSpinBox;
	QDoubleSpinBox * boundarydistmeanSpinBox;
	QDoubleSpinBox * boundarydistvarSpinBox;


 private:
	 // Labels:
	QLabel * xspacingLabel;
	QLabel * yspacingLabel;
	QLabel * zspacingLabel;
	QLabel * distvarLabel;
	QLabel * distmeanLabel;
	QLabel * timevarLabel;
	QLabel * timemeanLabel;
	QLabel * overlapvarLabel;
	QLabel * overlapmeanLabel;
	QLabel * volumevarLabel;
	QLabel * mergesplitpriorLabel;
	QLabel * appeardisappearpriorLabel;
	QLabel * translationpriorLabel;
	QLabel * boundarydistmeanLabel;
	QLabel * boundarydistvarLabel;

 };

 class FoldersTab : public QWidget
 {
     Q_OBJECT

 public:
	 FoldersTab( QWidget *parent = 0);
	// Line Edits:
	QLineEdit * numbersdirectoryLineEdit;
	QLineEdit * entropyfileLineEdit;
	QLineEdit * entropydirectoryLineEdit;
	QLineEdit * debugfileLineEdit;
	QLineEdit * debugdirectoryLineEdit;
	QLineEdit * debugdprefixLineEdit;
	QLineEdit * resultprefixLineEdit;
	QLineEdit * resultdirectoryLineEdit;

 private slots:
	void browseNumbersDirectory();
	void browseEntropyDirectory();
	void browseDebugDirectory();
	void browseResultDirectory();


 private:
	//Functions:
	 QPushButton * createButton(const QString &text, const char *member);
	 QLineEdit * createLineEdit(const QString &text = QString());
	// Labels
	QLabel * loadLabel;
	QLabel * saveLabel;
	QLabel * numbersdirectoryLabel;
	QLabel * entropyfileLabel;
	QLabel * entropydirectoryLabel;
	QLabel * debugfileLabel;
	QLabel * debugdirectoryLabel;
	QLabel * debugdprefixLabel;
	QLabel * resultprefixLabel;
	QLabel * resultdirectoryLabel;
	
	// PushButtons:
	 QPushButton * browseNumbersDirectoryButton;
	 QPushButton * browseEntropyDirectoryButton;
	 QPushButton * browseDebugDirectoryButton;
	 QPushButton * browseResultDirectoryButton;
	 QString lastdirpath;

 };



class TrackingDialog : public QDialog
{
	Q_OBJECT
public:
	// Functions:
	TrackingDialog(QWidget *parent = 0);
	std::vector<std::pair<std::string,float> > getParameters(void);
	std::vector<std::pair<std::string,std::string> > getFolders(void);

private:
	ParametersTab * paramtertab;
	FoldersTab * foldertab;

private slots:

private:
	 QTabWidget *tabWidget;
	// Push Buttons:
	QDialogButtonBox  *trackdlgButton;

};
#endif
		
