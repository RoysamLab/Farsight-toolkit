//#include <QDialog>
//#include <QtGui>

#include "window.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
	app.setOrganizationName("FARSIGHT Toolkit");
	app.setOrganizationDomain("www.farsight-toolkit.org");
	app.setApplicationName("neighborfile");
	app.setApplicationVersion("V1.0");
	MainWindow mainWin;
	mainWin.show();
	return app.exec();
}

	//make a window
	//QWidget *window = new QWidget;
	//window->setWindowTitle("Neighbor Files...");

	/*QFileDialog *dialog = new QFileDialog;
	dialog->setNameFilter("Images (*.tif)");
	dialog->setViewMode(QFileDialog::Detail);
	
	QObject::connect(dialog, SIGNAL(click()),*/



	/*fileName = QFileDialog::getOpenFileName(this,
     "Open Image", "/directory.......................", "Image Files (*.tif)");*/

	//dialog.setNameFilter("Images (*.tif)");
	//dialog.setViewMode(QFileDialog::Detail);

	/*QStringList fileNames;
	if (dialog.exec())
	{
		fileNames = dialog.selectedFiles();
	}*/

	//bool myswitch = false;

	//QPushButton *neighborbutton = new QPushButton("Neighbor Images");
 //   QObject::connect(neighborbutton, SIGNAL(clicked()),
 //                     NULL, SLOT(myswitchtrue(&myswitch)));
 //   neighborbutton->show();


	//if (myswitch)
	//{



//		QStringList files; 
//		while (files.isEmpty())	//keeps asking user to select files, window reopens if "cancel" is clicked
//		{
//			files= QFileDialog::getOpenFileNames(NULL //no parent class so it's null
//			,
//			"Select file(s)",	//title of window
//			"",					//directory
//			"Images (*.tif)");	//show only tif files
//
//		}
//
//		QStringList list = files;	//put files in a list
//		QStringList::Iterator it = list.begin();	//goes through the list starting from the first file
//		while(it != list.end())						//keeps printing out each file until the last file is shown
//		{
//			 std::cout<< (*it).toStdString() << "\n";
//			 ++it;
//		}
//
//	//	FindDialog *dialog = new FindDialog;
//		//window->show();			//print results
//
//		//QVector containing data from QStringList
//		QVector<QString> vect = QVector<QString>::fromList(list);
//		
//		//Form matrix and get neighbors
//		int WIDTH=0;
//		int HEIGHT=0;
//		
//		int location=0;
//		int col = 0;
//
//		QVector<QString> row;
//		QVector< QVector< QString > > matrix;
//		QString fileName;
//
//		std::cout << "\nEnter number of columns: ";
//		std::cin >> WIDTH;
//		
//		std::cout << "Enter number of rows: ";
//		std::cin >> HEIGHT;
//		
//		//set up matrix named mat
//		for (int j = 0; j < HEIGHT; j++)
//		{
//			for  (int i = 0; i < WIDTH; i++)
//			{
//				fileName = vect[location].section('\\', -1);
//				row.push_back( fileName );	//append files to row vector
//				location++;
//				
//			}
//			matrix.push_back( row );	//make a matrix vector of row vectors
//		//	std::cout << " " << matrix[0][col];
//			col++;
//			row.clear();
//		}
//			
//		//print results
//		std::cout << "\n" << "Matrix of files" << "\n";
//		for (int j = 0; j < HEIGHT; j++)
//		{
//			for  (int i = 0; i < WIDTH; i++)
//			{
//				std::cout << " " << matrix[j][i].toStdString();	//convert the QVector matrix to a StdString for outputing
//			}
//			std::cout << "\n";
//		}
//
//		QString pairstring;
//		QVector < QString > pairlist;
//
//		std::cout << "\n" << "My neighbors" << "\n";
//
//		//write to text file
//		std::ofstream oFile;
//		oFile.open("neighborimage.txt");
//
//
//		//paralisting column-wise
//		for (int j = 0; j < HEIGHT; j++)
//		{
//			for  (int i = 1; i < WIDTH; i++)
//			{
//				pairstring = matrix[j][i-1] + "\t\t" + matrix[j][i];	//Append two QStrings into one QString
//				pairlist.push_back( pairstring );						//Append pairstring to the pairlist
//				std::cout << pairstring.toStdString() << "\n";
//
//				oFile << pairstring.toStdString() << "\n";
//
//			}
//		}
//		
//		//paralisting row-wise
//		for (int j = 1; j < HEIGHT; j++)
//		{
//			for  (int i = 0; i < WIDTH; i++)
//			{
//				pairstring = matrix[j-1][i] + "\t\t" + matrix[j][i];
//				pairlist.push_back( pairstring );
//				std::cout << pairstring.toStdString() << "\n";
//
//				oFile << pairstring.toStdString() << "\n";
//			}
//		}
//
//		oFile.close();
//	//}
//
//    return app.exec();
//}
//
//void myswitchtrue(bool *myswitch)
//{
//	*myswitch = true;
//}