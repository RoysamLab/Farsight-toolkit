/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/
#ifndef GalleryDialog_H
#define GalleryDialog_H

#include <QtGui/QAction>
#include <QtGui/QWidget>
#include <QtGui/QVBoxLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QGridLayout>
#include <QtGui/QDialog>
#include <QtGui/QColorDialog>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtCore/QStringList>
#include <QtCore/QSignalMapper>
#include <QtGui/QTableWidget>
#include <QtGui/QTableWidgetItem>
#include <QtGui/QLineEdit>
#include <QtGui/QPlainTextEdit>
#include <QHeaderView>
#include <QRadioButton>
#include <QtGui/QPainter>

#include <vtkTable.h>
#include <vtkSmartPointer.h>
#include <vtkVariant.h>

#define MAX(a,b) (((a) > (b))?(a):(b))
#define MIN(a,b) (((a) < (b))?(a):(b))

class GalleryDialog : public QDialog
{
	Q_OBJECT
public:
	GalleryDialog(std::vector<std::pair<QImage,int> > gallery ,QVector<QColor>,QWidget *parent=0);

private slots:

	
private:
	//QMap<QString, QColor> m_snapshot;
	//QSignalMapper *signalMapper;

};


#endif
