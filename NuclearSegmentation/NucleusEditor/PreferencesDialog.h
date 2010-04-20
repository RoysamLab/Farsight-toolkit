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
#ifndef PREFERENCESDIALOG_H
#define PREFERENCESDIALOG_H

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
#include <QtCore/QMap>

class PreferencesDialog : public QDialog
{
	Q_OBJECT
public:
	PreferencesDialog(QMap<QString, QColor> colorItemsMap, QWidget *parent = 0);
	QMap<QString, QColor> GetColorItemsMap(){ return m_colorItemsMap; };

private slots:
	void chooseColor(const QString & colorItem);
	
private:
	QMap<QString, QColor> m_colorItemsMap;
	QMap<QString, QLabel *> colorLabelsMap;

	QSignalMapper *signalMapper;
	
};


#endif
