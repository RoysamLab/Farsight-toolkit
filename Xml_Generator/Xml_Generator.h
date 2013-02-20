#ifndef XMLGENERATOR_H
#define XMLGENERATOR_H



//QT INCLUDES
#include <QtGui/QLabel>
#include <QtGui/QDialog>
#include <QtGui/QPushButton>
#include <QtGui/QWidget>
#include <QtGui/QButtonGroup>
#include <QtGui/QVBoxLayout>
#include <QtGui/QMessageBox>
#include <QtGui/QLineEdit>
#include <QtGui/QFileDialog>
#include <QtGui/QMainWindow>
#include <QIODevice>
#include <QVector>
#include <QXmlStreamWriter>
#include <QFile>
#include<QtGui/QWizard>
#include <QtGui/QLabel>
#include <QtGui/QDialog>
#include <QtGui/QMenuBar>
#include <QtGui/QMenu>
#include <iostream>
#include<string>
#include <QIODevice>


// Class to generate the main Window

class Xml_Generator : public QMainWindow
{
    Q_OBJECT;


	
public:
	Xml_Generator(QWidget * parent = 0, Qt::WindowFlags flags = 0);
	~Xml_Generator();
	
    


protected:
	void createMenu();


protected slots:
	void getProjectName();
	
	void addChild();
	void nowSave();
	

protected:
	//Project Definition menu
	QMenu *Project_Definition;
	QAction *ProjectName;
	QAction *addchildAction;
	QAction *saveAction;
	QAction *exitAction;
	
	


};

//Class to Generate Dialog box for taking the Project Definition Name

class Enter_Name : public QDialog
{ 
Q_OBJECT;


	
public:
	Enter_Name(QWidget *parent=0);
	QString pName;
	void write(void);

	
private:
	QLineEdit *projectName;
	QPushButton * okButton;
	QPushButton * cancelButton;
	
};

//Class to Write an Xml File

class Generate_Xml : public QMainWindow
{
Q_OBJECT;

public:
	Generate_Xml(QWidget * parent = 0);
    ~Generate_Xml();
	//call this function to create xml
	void CreateXMLFile(QString s);
	// Call this to append child to the xml file
	void appendXMLFile(QLineEdit *p,QVector<QLineEdit *>  s,int track_subfield,int track_tabchild,int track_tabsibling);

	//Call this to append siblings to the xml file
	int rewriteXMLFile(QVector<QLineEdit *> p,int a,int track_tabsibling);
	
	
};

// Class to create a Add Child Dialog

class Add_Child : public QDialog
{
Q_OBJECT;

public slots:
	void addsubfield();
	void addfeature();
	void save();
	void remSubField();
	void addchild();//Not Yet Implemented
	;
public:
	Add_Child(QWidget *parent=0);
	~Add_Child();
	void write_subfield(void);
private:
	QXmlStreamWriter *xmlWriter;
	QVBoxLayout * inputsLayout;
	QLineEdit *ChildName;
	QVector<QLineEdit *> childName;
	QVector<QHBoxLayout *> iLayouts;
	QVector<QLabel *> inputLabels;
	QPushButton *addSubField;
	QPushButton *addFeature;
	QPushButton *okButton;
	QPushButton *cancelButton;
	QPushButton *saveButton;
	QPushButton *deleteButton;
	QPushButton *addChild;
	int track_subfield;
	int track_tabchild;
	int track_tabsibling;
	
};
#endif