#include "Xml_Generator.h"

Xml_Generator::Xml_Generator(QWidget * parent, Qt::WindowFlags flags)
: QMainWindow(parent,flags)
{


	createMenu();

	

	
	
	  activateWindow();

}


//destructor for the Xml_Generator
Xml_Generator::~Xml_Generator()
{
}



//Constructor for the class Enter_Name


Enter_Name::Enter_Name(QWidget *parent)
: QDialog(parent){

this->setWindowTitle(tr("Project Definition Name"));
this->setModal(false);
QLabel *Label = new QLabel(tr("Please enter the name for the project definition: "));
QLabel *Label1 = new QLabel(tr("Enter Name:"));
QVBoxLayout *vbox = new QVBoxLayout;
QHBoxLayout *hbox = new QHBoxLayout;
projectName = new QLineEdit();
projectName->setMinimumWidth(100);
projectName->setFocusPolicy(Qt::StrongFocus);



okButton = new QPushButton(tr("Ok"));
connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
okButton->setDefault(true);
okButton->setAutoDefault(true);
	
cancelButton = new QPushButton(tr("Cancel"));
connect(cancelButton, SIGNAL(clicked()), this, SLOT(reject()));
cancelButton->setDefault(false);
cancelButton->setAutoDefault(false);

hbox->addWidget(okButton);
hbox->addWidget(cancelButton);


vbox->addWidget(Label);
vbox->addStretch(1);

vbox->addWidget(Label1);
vbox->addWidget(projectName);
vbox->addStretch(1);
vbox->addLayout(hbox);



this->setLayout(vbox);


}

//Constructor for Generate_Xml

Generate_Xml::Generate_Xml(QWidget * parent)
: QMainWindow(parent)
{

}

//Destructor for Generate_Xml

Generate_Xml :: ~Generate_Xml()
{
}


Add_Child::Add_Child(QWidget *parent)
:QDialog(parent){
this->track_subfield=0;
this->track_tabchild=0;
this->track_tabchild++;
this->track_tabsibling=0;
this->setWindowTitle(tr("Add Child"));
this->setModal(false);
QLabel *Label = new QLabel(tr("Please Enter the name for the child:"));
QLabel *Label1 = new QLabel(tr("Enter Name:"));
inputsLayout = new QVBoxLayout;
QHBoxLayout *hbox = new QHBoxLayout;
QHBoxLayout *hbox1 = new QHBoxLayout;
QHBoxLayout *hbox2 = new QHBoxLayout;
ChildName = new QLineEdit();
ChildName->setMaximumWidth(100);
ChildName->setMinimumWidth(100);
ChildName->setFocusPolicy(Qt::StrongFocus);

///Ok Button Initialization

okButton = new QPushButton(tr("Ok"));
connect(okButton,SIGNAL(clicked()),this,SLOT(accept()));
okButton->setDefault(false);
okButton->setAutoDefault(false);

///cancelButton Initialization

cancelButton = new QPushButton(tr("Cancel"));
connect(cancelButton, SIGNAL(clicked()), this, SLOT(reject()));
cancelButton->setDefault(false);
cancelButton->setAutoDefault(false);


/// Addsubfield Button Initialization


addSubField = new QPushButton(tr("Add SubFields"));
connect(addSubField, SIGNAL(clicked()), this, SLOT(addsubfield()));
addSubField->setDefault(false);
addSubField->setAutoDefault(false);



///  AddFeature Button Initialization


addFeature = new QPushButton(tr("Add Features"));
connect(addFeature, SIGNAL(clicked()), this, SLOT(addfeature()));
addFeature->setDefault(false);
addFeature->setAutoDefault(false);

///SaveButton Initialization

saveButton = new QPushButton(tr("Save"));
connect(saveButton, SIGNAL(clicked()), this, SLOT(save()));
saveButton->setDefault(false);
saveButton->setAutoDefault(false);



//Delete Button Initialized


deleteButton = new QPushButton(tr("Delete Subfield"));
connect(deleteButton, SIGNAL(clicked()), this, SLOT(remSubField()));
deleteButton->setDefault(false);
deleteButton->setAutoDefault(false);

// Add Child Button Initialized

addChild = new QPushButton(tr("Add Child"));
connect(addChild,SIGNAL(clicked()),this,SLOT(addchild()));
addChild->setDefault(false);
addChild->setAutoDefault(false);


 inputsLayout->addWidget(Label);
hbox->addWidget(Label1);
hbox->addWidget(ChildName);
 inputsLayout->addLayout(hbox);
hbox1->addWidget(addSubField);
hbox1->addWidget(addFeature);
hbox1->addWidget(saveButton);
hbox1->addWidget(deleteButton);
hbox1->addWidget(addChild);
hbox2->addWidget(okButton);
hbox2->addWidget(cancelButton);
 inputsLayout->addLayout(hbox1);
 inputsLayout->addLayout(hbox2);
this->setLayout( inputsLayout);









}

//Destructor for the Add_Child Class

Add_Child::~Add_Child()
{
}

void Add_Child::addsubfield()
{

	QLineEdit *null=0;
	childName.push_back(null);
  this->track_subfield++;
  this->track_tabsibling++;
	QLabel *label = new QLabel(tr("Name ")); 
	inputLabels.push_back( label );
QLineEdit *in = new QLineEdit();
	in->setMinimumWidth(200);
	in->setFocusPolicy(Qt::StrongFocus);
	childName.push_back( in );

	QHBoxLayout *layout = new QHBoxLayout;
	layout->addWidget( inputLabels.back() );
	layout->addWidget( childName.back() );
	layout->addStretch(1);
	iLayouts.push_back(layout);
	inputsLayout->addLayout( iLayouts.back() );
	
	
	childName.back()->setFocus();


}
void Add_Child::addfeature()

{

QLabel *label = new QLabel(tr("Feature Name ")); 
	inputLabels.push_back( label );
QLineEdit *in = new QLineEdit();
	in->setMinimumWidth(200);
	in->setFocusPolicy(Qt::StrongFocus);
	childName.push_back( in );

	QHBoxLayout *layout = new QHBoxLayout;
	layout->addWidget( inputLabels.back() );
	layout->addWidget( childName.back() );
	layout->addStretch(1);
	iLayouts.push_back(layout);
	inputsLayout->addLayout( iLayouts.back() );
	childName.back()->setFocus();
	
	QLabel *label1 = new QLabel(tr("Value ")); 
	inputLabels.push_back( label1 );
QLineEdit *in1 = new QLineEdit();
	in1->setMinimumWidth(200);
	in1->setFocusPolicy(Qt::StrongFocus);
	childName.push_back( in1 );

	QHBoxLayout *layout1 = new QHBoxLayout;
	layout1->addWidget( inputLabels.back() );
	layout1->addWidget( childName.back() );
	layout1->addStretch(1);
	iLayouts.push_back(layout1);
	inputsLayout->addLayout( iLayouts.back() );
	childName.back()->setFocus();

}
void Add_Child::save()
{
this->write_subfield();

  }

void Add_Child::remSubField()
{
if(inputLabels.size()==0)
return;

delete inputLabels.back();
delete childName.back();
inputLabels.remove(inputLabels.size()-1);
	childName.remove(childName.size()-1);

	delete iLayouts.back();
	iLayouts.remove(iLayouts.size()-1);


}

void Add_Child::addchild()
{
}


/*----------------------------------------------------------------------------*/
/* function: createMenu()                                                    */
/*                                                                            */
/* This function is used to create Menus that are associated with various     */
/* functions in the application. In order to add a menu, we need to do the    */
/* following:                                                                 */
/* 1.) Define a QMenu type (e.g., QMenu *fileMenu) and add it to menuBar()    */
/* 2.) Define QAction elements (e.g., QAction *openAction) associated with    */
/*     each QMenu                                                             */
/*																			  */
/*In order to create an Action, we need to do the							  */
/* following:                                                                 */
/* 1.) Define a QAction (e.g., QAction *openAction)                           */
/* 2.) Label the QAction element (e.g., openAction = new QAction(QIcon(":src/ */
/*     images/open.png"), tr("&Open..."), this). The QIcon argumenet is       */
/*     optional.                                                              */
/* 3.) Add optional "setShortcut" and "setStatusTip".                         */
/* 4.) Finally, bind this item with a "connect" that essentially calls the    */
/*     module to implement the operation (e.g.,                               */
/*     connect(openAction, SIGNAL(triggered()), this, SLOT(loadImage())). In  */
/*     this example, "loadImage()" is the module that is being called. These  */
/*     modules should be defined as "private" operators in the main class.    */
/*     The actual routines performing the operations (e.g., an image          */
/*     thresholding operation) must be accessed from within the called module.*/
/*	   Finally, after all these action's, we bind them to a "QActionGroup".       */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Xml_Generator::createMenu()
{

	//For handling the Project_Definition Menu
	Project_Definition= menuBar()->addMenu(tr("Project Definition"));
	ProjectName = new QAction(tr("Project Name..."), this);
	connect(ProjectName, SIGNAL(triggered()), this, SLOT(getProjectName()));
	ProjectName->setShortcut(tr("Ctrl+P"));
	Project_Definition->addAction(ProjectName);



	addchildAction = new QAction(tr("Add Child"), this);
	addchildAction->setShortcut(tr("Ctrl+C"));
	connect(addchildAction, SIGNAL(triggered()), this, SLOT(addChild()));
	Project_Definition->addAction(addchildAction);



	saveAction = new QAction(tr("Save"), this);
	saveAction->setShortcut(tr("Ctrl+S"));
	connect(saveAction, SIGNAL(triggered()), this, SLOT(nowSave()));
	Project_Definition->addAction(saveAction);


	exitAction= new QAction(tr("Exit"), this);
	exitAction->setShortcut(tr("Ctrl+X"));
	connect(exitAction, SIGNAL(triggered()), this, SLOT(close()));
	Project_Definition->addAction(exitAction);

}



//Function to open a Dialog box to get the Project Definition Name


void Xml_Generator::getProjectName()
{    

	
	Enter_Name *a=new Enter_Name(this);
	if(a->exec())
	{
	a->write();
	}

	a->show();
}



void Xml_Generator::addChild()
{

	Add_Child *a = new Add_Child(this);

	if(a->exec())
	{   
	
	}
	a->show();

}








////Function for creating Xml_File


void Generate_Xml::CreateXMLFile(QString s)
{


	QFile Default("default1.xml");
	QXmlStreamWriter *xmlWriter = new QXmlStreamWriter();
	Default.remove();
	if (!Default.open(QIODevice::Append))
	{
	/* show wrror message if not able to open file */
	QMessageBox::warning(0, "Read only", "The file is in read only mode");
	}	
	else
	{
		std::string p= s.toStdString();
		std::cout << p;
	 std::cout << "You are now writing to the file";
     
	 xmlWriter->setAutoFormatting (true);
	 xmlWriter->setDevice(&Default);
	 xmlWriter->writeStartDocument();
	 xmlWriter->writeStartElement("ProjectDefinition");
	 xmlWriter->writeAttribute("name",s);//This is only useful
	
	 
	 xmlWriter->writeCharacters(" ");
	 Default.close();
	 //xmlWriter->writeEndElement();
	 /*xmlWriter->writeEndDocument();*/


	}


	delete xmlWriter;

}


void Generate_Xml::appendXMLFile(QLineEdit *s,QVector<QLineEdit *> p,int track_subfield,int track_tabchild,int track_tabsibling)
{
	
	int j=0;
	QFile Default("default1.xml");
	std::cout << "You are now inside appendXmlFile";
    QXmlStreamWriter *xmlWriter = new QXmlStreamWriter();
    QXmlStreamWriter *xmlWriter2 = new QXmlStreamWriter();
	if (!Default.open(QIODevice::Append))
	{
	/* show wrror message if not able to open file */
	QMessageBox::warning(0, "Read only", "The file is in read only mode");
	}	
	else
	{  
		QString tem;
		QString temp1;
		
		 xmlWriter->setAutoFormatting (true);
		 xmlWriter->setDevice(&Default);
		 xmlWriter2->setDevice(&Default);
		 xmlWriter->autoFormatting();
		 tem=s->displayText();
		 xmlWriter2->writeEndDocument();
		 for(int i=0;i<track_tabchild;i++)
		 {
			 xmlWriter->writeCharacters("\t");
		 }
		 xmlWriter->writeStartElement(tem);
		 /* std::cout << "You are now appending to the file";*/
		
		  
	 
	 /*xmlWriter1->setDevice(&Default);*/
	 
	 /*xmlWriter1->autoFormatting();*/
	  /*xmlWriter1->setAutoFormatting (true);*/
			
	 
	
		
		
		
	   
		
		//while(i<t)
		//{
		//	if(p.at(i)==0 && i != (t-1))
		//	{
		//		i++;
		//		
		//	tem=(p.at(i))->displayText();
		//	/*temp=tem.toStdString();*/
		//	xmlWriter->writeStartElement(tem);
		//	    i++;
		//		
		//		
		//	}
		//	if(p.at(i-1)!= 0 && p.at(i)== 0)
		//	{
		//		//xmlWriter1->writeCharacters(" ");
		//		xmlWriter->writeEndDocument();
		//		i++;
		//		
		//		if(i=t)
		//			break;
		//		
		//	}
		//	if(p.at(i)!= 0)
		//	{
		//		tem=(p.at(i))->displayText();
		//		i++;
		//		
		//		temp1=(p.at(i))->displayText();
		//		i++;
		//		
		//		xmlWriter->writeAttribute(tem,temp1);
		//		
		//	}
		//}
		
			 
			 xmlWriter->writeCharacters(" ");
		Default.close();
	
		for(int i=0;i<track_subfield;++i)
			j=rewriteXMLFile(p,j+1,track_tabsibling);
		
	
    
	 
	 //xmlWriter->writeStartDocument();
	 
	 //xmlWriter->writeAttribute("name"," ");//This is only useful
	
	
	 //xmlWriter->writeEndDocument();


	}
	delete xmlWriter2;
	bool a;
	Default.open(QIODevice::Append);
	QXmlStreamWriter *xmlWriter1 = new QXmlStreamWriter();
	xmlWriter1->setDevice(&Default);
	xmlWriter1->writeEndDocument();
	a=xmlWriter->autoFormatting();
	xmlWriter->setAutoFormatting(a);
	//xmlWriter->writeCharacters(" ");
	for(int i=0;i<track_tabchild;i++)
		 {
			 xmlWriter->writeCharacters("\t");
		 }
	xmlWriter->writeEndDocument();
	delete xmlWriter;
	delete xmlWriter1;

}

int Generate_Xml::rewriteXMLFile(QVector<QLineEdit *> p,int j,int track_tabsibling)
{
	bool tab;
	int tab_get;
	QFile Default("default1.xml");
	
	QXmlStreamWriter *xmlWriter = new QXmlStreamWriter();
	QXmlStreamWriter *newline = new QXmlStreamWriter();
	xmlWriter->setDevice(&Default);
	newline->setDevice(&Default);
	xmlWriter->setAutoFormatting (true);
	xmlWriter->setAutoFormatting(true);
	tab=xmlWriter->autoFormatting();
	
	xmlWriter->setAutoFormattingIndent(-10);
	tab_get=xmlWriter->autoFormattingIndent();
	
		
	
	if (!Default.open(QIODevice::Append))
	{
	/* show wrror message if not able to open file */
	QMessageBox::warning(0, "Read only", "The file is in read only mode");
	}	
	QString temp1;
	QString temp2;
		 
		 
	temp1=(p.at(j))->displayText();
		
	newline->writeEndDocument();

	for(int i=0;i<track_tabsibling;i++)//changed tab here actaully multiplied by two
	xmlWriter->writeCharacters("\t");

	xmlWriter->writeStartElement(temp1);
		
	j++;
				
	while(p.at(j)!=0)
	{
	temp1=(p.at(j))->displayText();
	j++;
	temp2=(p.at(j))->displayText();
	j++;
	xmlWriter->writeAttribute(temp1,temp2);
	}

	xmlWriter->writeEndElement();

	Default.close();
	delete xmlWriter;
	delete newline;
	return j;
}

void Enter_Name::write()
{
	QString temp;
	temp=projectName->displayText();
	Generate_Xml *file_pointer = new Generate_Xml(this);
    file_pointer->CreateXMLFile(temp);
	
	return;


}

void Add_Child::write_subfield()
{
	
	QLineEdit *null=0;
	childName.push_back(null);
    QString s = ChildName->displayText();
	std::string p = s.toStdString();
	std:: cout << p;
	Generate_Xml *file_pointer = new Generate_Xml(this);
    file_pointer->appendXMLFile(this->ChildName,this->childName,this->track_subfield,this->track_tabchild,this->track_tabsibling);
	QMessageBox::warning(0, "Saving", "If want to add more child go to menu and press Add Child");

}

void Xml_Generator::nowSave()
{
	QXmlStreamWriter *xmlWriter = new QXmlStreamWriter();
	QFile Default("default.xml");
	
	Default.open(QIODevice::Append);
	xmlWriter->setDevice(&Default);
	
	xmlWriter->writeStartElement("/ProjectDefinition");
	xmlWriter->writeCharacters(" ");
	
	xmlWriter->setAutoFormatting (true);
	xmlWriter->setAutoFormatting(true);
	Default.close();
	delete xmlWriter;
	QMessageBox::warning(0, "Saved Status", "Your file has been saved by the name default1.xml");


}