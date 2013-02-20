/*Open a window that ask the user to click on a button called "Neighbor images".  
Another dialogbox pops up for the user to select image files to be montaged. Set
the number of columns and rows.  The neighboring image files are matched and saved
in a text file (neighborimage.txt)
*/

#include "window.h"

MainWindow::MainWindow()
{
	status = statusBar();
	status->showMessage(tr("Ready"));

	selectButton = new QPushButton(tr("Select image files"));
	selectButton->setDefault(true);
	connect(selectButton, SIGNAL(clicked()), this, SLOT(selectfiles()));

	numofiles = new QLabel(tr("No files selected"));

	columnTitle = new QLabel(tr("Column:"));
	columnBox = new QSpinBox;
	rowTitle = new QLabel(tr("Row:"));
	rowBox = new QSpinBox;
	
	testRowColTitle = new QLabel(tr("Test if the selected files form a complete rectangle"));
	testRowColButton = new QPushButton(tr("Test"));
	testRowColButton->setEnabled(false);
	connect(testRowColButton, SIGNAL(clicked()),this, SLOT(testforrectangle()));

	twidget = new QTableWidget(this);

	neighborTitle = new QLabel(tr("<font size=\"1000\"><b>Neighboring images</b></font>"));

	nwidget = new QTableWidget(this);

	listofileCheckBox = new QCheckBox(tr("Save list of &selected files"), this);
	listofileCheckBox->setChecked(true);
	neighborCheckBox = new QCheckBox(tr("Save list of &image pairs"), this);
	neighborCheckBox->setChecked(true);

	autoSaveButton = new QPushButton(tr("AutoSave"));
	autoSaveButton->setEnabled(false);
	connect(autoSaveButton, SIGNAL(clicked()), this, SLOT(autoSave()));

	saveButton = new QPushButton(tr("Save text file(s)"));
	saveButton->setDefault(true);
	saveButton->setEnabled(false);
	connect(saveButton, SIGNAL(clicked()), this, SLOT(save()));

	saveSelectionButton = new QPushButton(tr("Save selection"));
	saveSelectionButton->setEnabled(false);
	connect(saveSelectionButton, SIGNAL(clicked()), this, SLOT(saveSelection()));

	qcancel = new QPushButton(tr("Close"));
	connect(qcancel, SIGNAL(clicked()), qApp, SLOT(quit())); //quit is a function inside the class qApp

	//Window Layout
	QHBoxLayout *rowcolLayout = new QHBoxLayout;
	rowcolLayout->addWidget(columnTitle);
	rowcolLayout->addWidget(columnBox);
	rowcolLayout->addWidget(rowTitle);
	rowcolLayout->addWidget(rowBox);

	QVBoxLayout *rowcoltestLayout = new QVBoxLayout;
	rowcoltestLayout->addLayout(rowcolLayout);
	rowcoltestLayout->addWidget(testRowColTitle);
	rowcoltestLayout->addWidget(testRowColButton);
	rowcoltestLayout->addStretch();

	QHBoxLayout *tableLayout = new QHBoxLayout;
	tableLayout->addWidget(twidget);
	tableLayout->addLayout(rowcoltestLayout);

	QHBoxLayout *bottomLayout = new QHBoxLayout;
	bottomLayout->addWidget(autoSaveButton);
	bottomLayout->addWidget(saveButton);
	bottomLayout->addWidget(saveSelectionButton);
	bottomLayout->addWidget(qcancel);

	QVBoxLayout *finalLayout = new QVBoxLayout;
	finalLayout->addWidget(numofiles);
	finalLayout->addWidget(selectButton);
	finalLayout->addLayout(tableLayout);
	finalLayout->addWidget(neighborTitle);
	finalLayout->addWidget(nwidget);
	finalLayout->addWidget(listofileCheckBox);
	finalLayout->addWidget(neighborCheckBox);
	finalLayout->addLayout(bottomLayout);

	QGroupBox *MainGroupBox = new QGroupBox(tr("Montage Images"));
	MainGroupBox->setLayout(finalLayout);
	this->setCentralWidget(MainGroupBox);
}

void MainWindow::selectfiles()
{
	////disable save buttons if not yet tested for rectangularity
	//saveButton->setEnabled(false);
	//autoSaveButton->setEnabled(false);
	//saveSelectionButton->setEnabled(false);

	QStringList files; 
	files.clear();
	LastDirectory = settings.value("neighborDir", " ").toString();

	files = QFileDialog::getOpenFileNames(	this,				//no parent class so it's null
											tr("Select file(s)"),	//title of window
											LastDirectory,		//directory
											tr("Images (*.tif)"));	//show only tif files
	
	if (!files.isEmpty())
	{
		//get the path from the first element of "files" and set it as the LastDirectory
		QFileInfo fileinfo = QFileInfo(files.at(0));
		this->LastDirectory = fileinfo.path();
		this->settings.setValue("neighborDir", LastDirectory);
		settings.sync ();
		testRowColButton->setEnabled(true);
	}
	else
	{
		if (vect.isEmpty())
		{
			testRowColButton->setEnabled(false);
		}
		return;
	}
	
	//set up table to show the list of selected files in the GUI
	twidget->setColumnCount(1);
	QStringList tabletitle;
	tabletitle<<tr("List of selected files");
	twidget->setHorizontalHeaderLabels(tabletitle);

	QStringList list = files;	//put files in a list
	twidget->setRowCount((int) list.size());	//set number of rows in the table equal to the number of selected files
	//QString cutlist;

	for (int i = 0; i < (int) list.size(); i++)
	{
		list[i] = list[i].section('\\', -1);	//remove pathname and show filename only
		//cutlist = list[i];
		QTableWidgetItem *newItem = new QTableWidgetItem(tr("%1").arg(list[i]));
		twidget->setItem(0,i,newItem);
	}

	//QVector containing data from QStringList
	vect = QVector<QString>::fromList(list);
	numofiles->setText(QString(tr("number of files: %1")).arg((int) vect.size()));
	status->showMessage(QString(tr("%1 file have been selected")).arg((int) vect.size()));
	
	//limit columnBox and rowBox to the maximum number of files selected
	columnBox->setRange(1,(int) vect.size()); 
	rowBox->setRange(1,(int) vect.size());
}

//The column and row must form a complete rectangle to make the matrix
void MainWindow::testforrectangle()
{
	int total = (columnBox->value())*(rowBox->value());
	if (total == (int) vect.size())
	{
		//Enable save buttons if row and column values form a complete rectangle
		selectButton->setDefault(false);
		saveButton->setEnabled(true);
		autoSaveButton->setEnabled(true);
		saveSelectionButton->setEnabled(true);

		//setup table based upon the given column and row value
		WIDTH = columnBox->value();
		HEIGHT = rowBox->value();
		nwidget->setColumnCount(WIDTH);
		nwidget->setRowCount(HEIGHT);
		
		int location=0;
		int col = 0;

		QVector<QString> row;
		QString fileName;
		
		//set up matrix named mat
		for (int j = 0; j < HEIGHT; j++)
		{
			for  (int i = 0; i < WIDTH; i++)
			{
					fileName = vect[location].section('\\', -1);	//remove pathname and show filename only for the text file
					row.push_back( fileName );	//append files to row vector
					location++;

					QTableWidgetItem *newItem = new QTableWidgetItem(tr("%1").arg(fileName));
					nwidget->setItem(j,i,newItem);
			}
			matrix.push_back( row );	//make a matrix vector of row vectors
			col++;
			row.clear();
		}
		statusBar()->showMessage(tr("SUCCESS"));
	}
	else
	{
		//clear neighbor table and disable save buttons if row and column value are incorrect
		nwidget->setColumnCount(0);
		nwidget->setRowCount(0);
		saveButton->setEnabled(false);
		autoSaveButton->setEnabled(false);
		saveSelectionButton->setEnabled(false);
		statusBar()->showMessage(tr("FIX the row and/or column setting."));
	}
}

void MainWindow::autoSave()
{
	QList<QTableWidgetSelectionRange> ranges = nwidget->selectedRanges();
	
	//check if listofileCheckBox is checked to save the list of selected files in a text file
	if (listofileCheckBox->checkState())
	{
		QString path;
		path = LastDirectory + "\\listofimages.txt";
		std::ofstream outfileselect;
		outfileselect.open( path.toStdString().c_str());

		for (int i = 0; i < (int) vect.size(); i++)
		{
			outfileselect << vect[i].toStdString() << "\n";
		}
		outfileselect.close();
	}
	//check if neighborCheckBox is checked to save the list of IMAGE PAIRS
	if (neighborCheckBox->checkState())
	{
		QString path;
		QString pairstring;
		QVector < QString > pairlist;
		path = LastDirectory + "\\neighborimagelist.txt";
		std::ofstream outfileneigh;
		outfileneigh.open( path.toStdString().c_str());
		//std::cout << path.toStdString()<<std::endl;

		//paralisting column-wise
		for (int j = 0; j < HEIGHT; j++)
		{
			for  (int i = 1; i < WIDTH; i++)
			{
				pairstring = matrix[j][i-1] + " " + matrix[j][i];	//Append two QStrings into one QString
				pairlist.push_back( pairstring );						//Append pairstring to the pairlist
				outfileneigh << pairstring.toStdString() << "\n";
			}
		}
		//paralisting row-wise
		for (int j = 1; j < HEIGHT; j++)
		{
			for  (int i = 0; i < WIDTH; i++)
			{
				pairstring = matrix[j-1][i] + " " + matrix[j][i];
				pairlist.push_back( pairstring );
				outfileneigh << pairstring.toStdString() << "\n";
			}
		}
		outfileneigh.close(); 
	}

	//shows user the status after pressing the autoSave button
	if (listofileCheckBox->checkState() || neighborCheckBox->checkState())
	{
		if (listofileCheckBox->checkState() && neighborCheckBox->checkState())
		{
			statusBar()->showMessage(tr("TWO FILES AUTOSAVED"));
		}
		else
		{
			statusBar()->showMessage(tr("ONE FILE AUTOSAVED"));
		}
	}
	else
	{
		statusBar()->showMessage(tr("Please click a checkbox"));
	}
}

void MainWindow::save()
{
	//check if listofileCheckBox is checked to save the list of selected files
	if (listofileCheckBox->checkState())
	{
		//write list of selected files to a text file named "imagelist.txt"
		std::ofstream flist;
		
		QString newfile = QFileDialog::getSaveFileName( this,
														tr("Save list of selected files"),
														QString(LastDirectory + "/imagelist.txt"),
														tr("Text files (*.txt)"));

		//std::cout << newfile.toStdString()<<std::endl;
		flist.open(newfile.toStdString().c_str() );

		for (int i = 0; i < (int) vect.size(); i++)
		{
			flist << vect[i].toStdString() << "\n";
		}
		flist.close();
	}
	//check if neighborCheckBox is checked to save the list of IMAGE PAIRS in a text file "imagepairs.txt"
	if (neighborCheckBox->checkState())
	{
		QString pairstring;
		QVector < QString > pairlist;
		
		QString neighfile = QFileDialog::getSaveFileName(	this,
															tr("Save List of image pairs"),
															QString(LastDirectory + "/imagepairs.txt"),
															tr("Text files (*.txt)"));
		
		std::ofstream oFile;

		//std::cout << neighfile.toStdString()<<std::endl;
		oFile.open(neighfile.toStdString().c_str() );  //convert neighfile (QString) to QStdString to C String

		//paralisting column-wise
		for (int j = 0; j < HEIGHT; j++)
		{
			for  (int i = 1; i < WIDTH; i++)
			{
				pairstring = matrix[j][i-1] + " " + matrix[j][i];	//Append two QStrings into one QString
				pairlist.push_back( pairstring );						//Append pairstring to the pairlist
				oFile << pairstring.toStdString() << "\n";
			}
		}
		//paralisting row-wise
		for (int j = 1; j < HEIGHT; j++)
		{
			for  (int i = 0; i < WIDTH; i++)
			{
				pairstring = matrix[j-1][i] + " " + matrix[j][i];
				pairlist.push_back( pairstring );
				oFile << pairstring.toStdString() << "\n";
			}
		}
		oFile.close();
	}
	//show status to user after pressing the save button
	if (listofileCheckBox->checkState() || neighborCheckBox->checkState())
	{
		statusBar()->showMessage(tr("DONE"));
	}
	else
	{
		statusBar()->showMessage(tr("Please click a checkbox"));
	}
}

void MainWindow::saveSelection()
{
	//get indexes of the selected cells
	QList<QTableWidgetSelectionRange> ranges = nwidget->selectedRanges();

	if (!ranges.isEmpty())
	{
		//define variables
		int initialcolumn = ranges.first().leftColumn();
		int finalcolumn = ranges.first().rightColumn();
		int initialrow = ranges.first().topRow();
		int finalrow = ranges.first().bottomRow();

		QString path;

		/*std::cout << "row count: " << ranges.first().rowCount() << "\n";
		std::cout << "column count: " << ranges.first().columnCount() << "\n";
		std::cout << "top row: " << ranges.first().topRow() << "\n";
		std::cout << "bottom row: " << ranges.first().bottomRow() << std::endl;
		std::cout << "left column: " << ranges.first().leftColumn() << "\n";
		std::cout << "right column: " << ranges.first().rightColumn() << std::endl;*/

		//save list of selected files from the table selection in a text file
		if (listofileCheckBox->checkState())
		{
			std::ofstream outfileselect;
			
			QString newfile = QFileDialog::getSaveFileName( this,
															tr("Save list of selected files"),
															QString(LastDirectory + "/imagelistselected.txt"),
															tr("Text files (*.txt)"));

			//std::cout << newfile.toStdString()<<std::endl;
			outfileselect.open(newfile.toStdString().c_str() );

			/*path = LastDirectory + "\\listofimages.txt";
			outfileselect.open( path.toStdString().c_str());*/

			for (int j = initialrow; j <= finalrow; j++)
			{
				for (int i = initialcolumn; i <= finalcolumn; i++)
				{
					outfileselect << matrix[j][i].toStdString() << "\n";
					std::cout << matrix[j][i].toStdString() << "\n";
				}
			}
			outfileselect.close();
		}
		//save image pairs from the table selection in a text file
		if (neighborCheckBox->checkState())
		{
			std::ofstream outfileneigh;
			
			QString neighfile = QFileDialog::getSaveFileName( this,
												tr("Save List of image pairs"),
												QString(LastDirectory + "/imagepairsselected.txt"),
												tr("Text files (*.txt)"));

			QString pairstring;
			QVector < QString > pairlist;
			//path = LastDirectory + "\\neighborimagelist.txt";
			outfileneigh.open( neighfile.toStdString().c_str());

			//paralisting column-wise
			for (int j = initialrow; j <= finalrow; j++)
			{
				for  (int i = initialcolumn+1; i <= finalcolumn; i++)
				{
					//std::cout << "Matrix[" << j << "][" << (i-1) << "]" << std::endl;
					pairstring = matrix[j][i-1] + " " + matrix[j][i];	//Append two QStrings into one QString
					pairlist.push_back( pairstring );					//Append pairstring to the pairlist
					outfileneigh << pairstring.toStdString() << "\n";
				}
			}
			
			//paralisting row-wise
			for (int j = initialrow+1; j <= finalrow; j++)
			{
				for  (int i = initialcolumn; i <= finalcolumn; i++)
				{
					pairstring = matrix[j-1][i] + " " + matrix[j][i];
					pairlist.push_back( pairstring );
					outfileneigh << pairstring.toStdString() << "\n";
				}
			}
			outfileneigh.close(); 
			
			}
		if (listofileCheckBox->checkState() || neighborCheckBox->checkState())
		{
			statusBar()->showMessage(tr("Selection files SAVED."));
		}
		else
		{
			statusBar()->showMessage(tr("Please click a checkbox"));
		}
	}
	else
	{
		statusBar()->showMessage(tr("No selection has been made."));
	}
}