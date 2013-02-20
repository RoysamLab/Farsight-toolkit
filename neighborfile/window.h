#include <QtGui>
#include <QApplication>

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

class MainWindow : public QMainWindow
{
	Q_OBJECT

	public:
		MainWindow();
		QStatusBar *status;
		QSpinBox *columnBox;
		QSpinBox *rowBox;
		QPushButton	*selectButton;
		QPushButton	*autoSaveButton;
		QPushButton	*saveButton;
		QPushButton	*saveSelectionButton;
		QPushButton	*testRowColButton;
		QPushButton	*qcancel;
		QLabel *numofiles;
		QLabel *neighborTitle;
		QLabel *testRowColTitle;
		QLabel *columnTitle;
		QLabel *rowTitle;
		QCheckBox *listofileCheckBox;
		QCheckBox *neighborCheckBox;
		QTableWidget *twidget;
		QTableWidget *nwidget;

		int WIDTH;
		int HEIGHT;
		QVector<QString> vect;
		QVector< QVector< QString > > matrix;
		QString LastDirectory;
		

	private:
		QSettings settings;

	public slots:
		void selectfiles();
		void testforrectangle();
		void autoSave();
		void save();
		void saveSelection();

};

#endif