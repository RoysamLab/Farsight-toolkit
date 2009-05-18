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