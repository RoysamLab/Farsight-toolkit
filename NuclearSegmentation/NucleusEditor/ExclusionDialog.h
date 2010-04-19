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
#ifndef EXCLUSIONDIALOG_H
#define EXCLUSIONDIALOG_H

#include <QtGui/QAction>
#include <QtGui/QWidget>
#include <QtGui/QVBoxLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QDialog>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QFileDialog>
#include <QtGui/QInputDialog>
#include <QtGui/QMessageBox>
#include <QtGui/QSpinBox>
#include <QtGui/QPainter>
#include <QtCore/QFileInfo>


class ExclusionDialog : public QDialog
{
	Q_OBJECT
public:
	ExclusionDialog(QImage * pImage = NULL, QWidget *parent = 0);
	int getMargin(int x);

private slots:
	void updatePreview();
private:
	void addSpin(QString label, int min, int deflt, int max, QString units);

	QVBoxLayout * masterLayout;
	QVBoxLayout * spinLayout;
	QVector<QSpinBox *> spins;
	QLabel *preview;
	QPushButton *okButton;

	QImage * baseImage;
};


#endif
