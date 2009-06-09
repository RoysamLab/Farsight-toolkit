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

#include "TableView.h"
#include <QtCore/QModelIndexList>

//Constructor
TableView::TableView(AbstractModel * model, QWidget *parent) 
 : QWidget(parent)
{
	this->m_Table = new QTableView();
	this->m_Model = model;
	this->m_TableAdapter = new vtkQtTableModelAdapter( const_cast<vtkTable*>(model->GetTable(0)), 0 );
	this->m_Table->setModel(this->m_TableAdapter);
	this->m_Table->setEditTriggers(QAbstractItemView::NoEditTriggers);
	this->m_Table->setSelectionMode(QAbstractItemView::ExtendedSelection);
	this->m_Table->setSelectionBehavior(QAbstractItemView::SelectRows);
	this->m_Table->setAlternatingRowColors(true);
	this->m_Table->resizeRowsToContents();
	this->m_Table->resizeColumnsToContents();

	this->m_Layout = new QVBoxLayout();
	this->m_Layout->addWidget(m_Table);
	this->m_Layout->setContentsMargins(0,0,0,0);

	this->setLayout(m_Layout);
	this->setWindowTitle(tr("Table"));
	this->setAttribute ( Qt::WA_DeleteOnClose );

	m_Selecting = false;	//Used to make sure that we don't keep updating view when we change the selections.
	QObject::connect(this->m_Table->selectionModel(), SIGNAL(selectionChanged(const QItemSelection&,const QItemSelection&)), \
					 this,                            SLOT(slotQtSelectionChanged(const QItemSelection&, const QItemSelection&)));

	QObject::connect(this->m_Model->GetObjectSelection(), SIGNAL(changed(int, int)), \
					 this,								  SLOT(slotObjSelectionChanged(int, int)));

}

TableView::~TableView()
{
	delete this->m_Table;
	delete this->m_TableAdapter;
}

void TableView::slotQtSelectionChanged(const QItemSelection& s1, const QItemSelection& s2)
{ 
	m_Selecting = true;

	// Convert from a QModelIndexList to rows in table
	const QModelIndexList qmil = this->m_Table->selectionModel()->selectedRows();

	this->m_Model->GetObjectSelection()->clear();
	for (int i = 0; i < qmil.size(); ++i) 
	{
		this->m_Model->GetObjectSelection()->add( qmil.at(i).row(), 0);
    }

	m_Selecting = false;

}

void TableView::slotObjSelectionChanged(int id, int t)
{
	if(m_Selecting)			//I'm the one doing selections, so do not respond
		return;

	std::vector<int> rows = this->m_Model->GetObjectSelection()->getObjects();

	this->m_Table->clearSelection();
	for (int i = 0; i < (int)rows.size(); ++i) 
	{
		this->m_Table->selectRow( rows.at(i) );
    }
}