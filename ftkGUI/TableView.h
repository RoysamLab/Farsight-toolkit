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

#ifndef _TABLEVIEW_H
#define _TABLEVIEW_H

#include "AbstractModel.h"
#include <QtGui/QItemSelection>
#include <QtGui/QItemSelectionModel>
#include <QtGui/QTableView>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include <vtkQtTableModelAdapter.h>

class TableView : public QWidget
{
	Q_OBJECT;

public:
	TableView(AbstractModel * model, QWidget *parent = 0);
	~TableView();

private slots:
	void slotQtSelectionChanged(const QItemSelection&, const QItemSelection&);
	void slotObjSelectionChanged(int, int);

private:
	bool m_Selecting;
	AbstractModel *m_Model;
	QVBoxLayout *m_Layout;
	QTableView *m_Table;
	vtkQtTableModelAdapter *m_TableAdapter;
};

#endif