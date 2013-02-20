#include "ActiveValidationWizard.h"

Active_Validation_Wizard::Active_Validation_Wizard(QWidget *parent) :
    QWidget(parent)
{
	QGridLayout *mainLayout = new QGridLayout;

    int frameStyledata = QFrame::Sunken | QFrame::Panel;
    dataFileNamedata = new QLabel;
    dataFileNamedata->setFrameStyle(frameStyledata);

    int frameStylelable = QFrame::Sunken | QFrame::Panel;
    dataFileNamelable = new QLabel;
    dataFileNamelable->setFrameStyle(frameStylelable);

    browseButtondata = new QPushButton(tr("Data Browse"));
    loadButtondata = new QPushButton(tr("Load"));
	browseButtonlable = new QPushButton(tr("Label Browse"));
    loadButtonlable = new QPushButton(tr("Load"));

    featureNumLabel = new QLabel(tr("Feature size:"));
    featureNum = new QLabel;
    featureNum->setFrameStyle(frameStyledata);
    sampleNumLabel = new QLabel(tr("Sample size:"));
    sampleNum = new QLabel;
    sampleNum->setFrameStyle(frameStyledata);

    numbinLable = new QLabel(tr(" Bin Number:"));
    numbinBox = new QLineEdit;
    deltaLable = new QLabel(tr("Delta:"));
    deltaBox = new QLineEdit;
	runtimeLable = new QLabel(tr(" Run times:"));
    runtimeBox = new QLineEdit;

    runButton = new QPushButton(tr("Run"));

    connect(browseButtondata, SIGNAL(clicked()), this, SLOT(browsedata()));
    connect(loadButtondata, SIGNAL(clicked()), this, SLOT(loaddata()));
	connect(browseButtonlable, SIGNAL(clicked()), this, SLOT(browselable()));
    connect(loadButtonlable, SIGNAL(clicked()), this, SLOT(loadlable()));
	connect(runButton, SIGNAL(clicked()), this, SLOT(Run()));


    for ( int col = 0; col<= 3; col++)
    {
        mainLayout->setColumnMinimumWidth(col,100);
        mainLayout->setColumnStretch(col, 1);
    }
    for ( int row = 1; row <= 5; row++)
    {
        mainLayout->setRowMinimumHeight(row,20);
        mainLayout->setRowStretch(row, 1);
    }

	mainLayout->addWidget(browseButtondata, 0, 0);
	mainLayout->addWidget(dataFileNamedata, 0, 1, 1, 2);
	mainLayout->addWidget(loadButtondata, 0, 3);

	mainLayout->addWidget(browseButtonlable, 1, 0);
	mainLayout->addWidget(dataFileNamelable, 1, 1, 1, 2);
	mainLayout->addWidget(loadButtonlable, 1, 3);

	mainLayout->addWidget(sampleNumLabel, 2, 0);
    mainLayout->addWidget(sampleNum, 2, 1);
    mainLayout->addWidget(featureNumLabel, 2, 2);
    mainLayout->addWidget(featureNum, 2, 3);

    mainLayout->addWidget(numbinLable, 3, 0);
    mainLayout->addWidget(numbinBox, 3, 1);
	mainLayout->addWidget(deltaLable, 3, 2);
    mainLayout->addWidget(deltaBox, 3, 3);

	mainLayout->addWidget(runtimeLable, 4, 0);
    mainLayout->addWidget(runtimeBox, 4, 1);
	mainLayout->addWidget(runButton, 4, 2, 1, 2);

    setLayout(mainLayout);
}

Active_Validation_Wizard::~Active_Validation_Wizard()
{
}

void Active_Validation_Wizard::browsedata()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose Data Files"),
                                                    QDir::currentPath(), tr("Text files (*.txt)"));
    if (!fileName.isEmpty())
    {
        dataFileNamedata->setText(fileName);
        this->FileNamedata = fileName;
    }
}

void Active_Validation_Wizard::loaddata()
{
	std::string file = this->FileNamedata.toStdString();
	this->ReadFiledata(file.c_str());
	this->sampleNum->setText( QString::number(this->numsamp));
	this->featureNum->setText( QString::number(this->numfeat));
}

void Active_Validation_Wizard::browselable()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose Label Files"),
                                                    QDir::currentPath(), tr("Text files (*.txt)"));
    if (!fileName.isEmpty())
    {
        dataFileNamelable->setText(fileName);
        this->FileNamelable = fileName;
    }
}

void Active_Validation_Wizard::loadlable()
{
	std::string file = this->FileNamelable.toStdString();
	this->ReadFilelable(file.c_str());
}

void Active_Validation_Wizard::Run()
{
	int numbin;
	int runtime;
	double delta;
		
	std::string snumbin = this->numbinBox->text().toStdString();
	if(snumbin.length()>0)
	{
		numbin = atoi(snumbin.c_str());
	}
	else
	{
		std::cout <<"Please input valid number of bin !!!"<<std::endl;
	}

	std::string sdelta = this->deltaBox->text().toStdString();
	if(sdelta.length()>0)
	{
		delta = atof(sdelta.c_str());
	}
	else
	{
		std::cout <<"Please input valid delta !!!"<<std::endl;
	}

	std::string sruntime = this->runtimeBox->text().toStdString();
	if(sruntime.length()>0)
	{
		runtime = atoi(sruntime.c_str());
	}
	else
	{
		std::cout <<"Please input valid run time !!!"<<std::endl;
	}

	ActiveValidation *active_validation = new ActiveValidation();
	this->sumnumsam = 0;
	this->sumnumit = 0;
	for(int i = 0; i < runtime; i++ )
	{
		active_validation->Initializing(this->data, this->label, numbin, delta,(unsigned long)i);
		active_validation->Sampling();
		active_validation->TrueAccuracy();
		this->phat.push_back(active_validation->phat);
		this->varphat.push_back(sqrt( (double)active_validation->varphat));
		this->numiteration.push_back(active_validation->numiteration);
		this->numsampled.push_back(active_validation->numsampled);
		this->varphattrue.push_back( (active_validation->phat - active_validation->trueaccuracy) / active_validation->phat );
	
		sumnumsam += active_validation->numsampled;
		sumnumit += active_validation->numiteration;



		if (i == 1)
		{
			std::cout <<"The true accuracy is "<<active_validation->trueaccuracy<<std::endl;
		}
		if(i%100 == 0)
			std::cout <<"run time left" <<runtime - i<<std::endl;
	}
	this->sumnumsam = this->sumnumsam/(double)runtime;
	this->sumnumit = this->sumnumit/(double)runtime;
	std::cout <<"Averrage sample number for "<<runtime<<" is "<<this->sumnumsam<<std::endl;
	std::cout <<"Averrage iteration number for "<<runtime<<" is "<<this->sumnumit<<std::endl;
	this->WriteFile();
	delete active_validation;
}

void Active_Validation_Wizard::ReadFiledata(const char *filename)
{
	FILE *fp = fopen(filename,"r");
	int num_samples = 0;
	int num_features = 0;
	if(fp == NULL)
	{
		fprintf(stderr,"can't open input file %s\n",filename);
		exit(1);
	}

	int n=0;
	while(1)
	{
		int c = fgetc(fp);
		switch(c)
		{
			case '\n':				
				(num_samples)++;	
				n++;
				if(num_features == 0)num_features = n;
				break;	
			case '\t':
				n++;
				break;
			case EOF:
				goto out;
			default:
				;
		}
	}
	out:
	rewind(fp);

	std::vector<std::vector<double> > tempdata;
	for(int i=0; i<num_samples; i++)
	{
		std::vector<double > tempvector;
		double temp;

		for(int j=0; j<num_features; j++)
		{
			fscanf(fp, "%lf", &temp);
			tempvector.push_back(temp);
		}
		tempdata.push_back(tempvector);
	}

	this->numsamp = num_samples;
	this->numfeat = num_features;
	this->data = tempdata;
	std::cout<<"Data file loaded !!!"<<std::endl;
}

void Active_Validation_Wizard::ReadFilelable(const char *filename)
{
	FILE *fp = fopen(filename,"r");
	int num_samples = 0;
	int num_features = 0;
	if(fp == NULL)
	{
		fprintf(stderr,"can't open input file %s\n",filename);
		exit(1);
	}

	int n=0;
	while(1)
	{
		int c = fgetc(fp);
		switch(c)
		{
			case '\n':				
				(num_samples)++;	
				n++;
				if(num_features == 0)num_features = n;
				break;	
			case '\t':
				n++;
				break;
			case EOF:
				goto out;
			default:
				;
		}
	}
	out:
	rewind(fp);

	std::vector<int > templable;
	for(int i=0; i<num_samples; i++)
	{
		double temp;
		fscanf(fp, "%lf", &temp);
		templable.push_back((int)temp);
	}

	this->label = templable;

	std::cout<<"Label file loaded !!!"<<std::endl;
}

void Active_Validation_Wizard::WriteFile()
{
	const char* filename = "Act_Val_result.txt";
	FILE *fp = fopen(filename,"w");
	// head names
	fprintf(fp,"%s", "phat");
	fprintf(fp,"\t");
	fprintf(fp,"%s", "varphat");
	fprintf(fp,"\t");
	fprintf(fp,"%s", "numiteration");
	fprintf(fp,"\t");
	fprintf(fp,"%s", "numsample");
	fprintf(fp,"\t");
	fprintf(fp,"%s", "difference_phat_trueaccuracy");
	fprintf(fp,"\t");
	fprintf(fp,"\n");

	//runtime result
	int counter = 0;
	for(int i = 0; i < this->phat.size(); i++)
	{
		fprintf(fp,"%f", this->phat[i]);
		fprintf(fp,"\t");
		fprintf(fp,"%f", this->varphat[i]);
		fprintf(fp,"\t");
		fprintf(fp,"%d", this->numiteration[i]);
		fprintf(fp,"\t");
		fprintf(fp,"%d", this->numsampled[i]);
		fprintf(fp,"\t");
		fprintf(fp,"%f", this->varphattrue[i]);
		fprintf(fp,"\n");
		if(this->varphattrue[i] < 0.05 && this->varphattrue[i] > -0.05 )
		{
			counter++;
		}
	}
	std::cout<< "there are  "<<counter << " times that the phat is within  0.05 of true accuracy !"<<std::endl;
	fprintf(fp,"%s", "average numsample");
	fprintf(fp,"\t");
	fprintf(fp,"%f", this->sumnumsam);
	fprintf(fp,"\t");
	fprintf(fp,"%s", "average iteration");
	fprintf(fp,"\t");
	fprintf(fp,"%f", this->sumnumit);
	fprintf(fp,"\t");
	fprintf(fp,"\n");
	fclose(fp);
}