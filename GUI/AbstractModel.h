#ifndef _ABSTRACTMODEL_H
#define _ABSTRACTMODEL_H

#include <QObject>

#include <vtkVariantArray.h>
#include <vtkTable.h>

#include <string>
#include <vector>

#include "ObjectSelection.h"

//****************************************************************************
// NOTES:
//  1. MUST emit changed when data changes!!
//****************************************************************************
class AbstractModel : public QObject
{
	Q_OBJECT;

public:
	AbstractModel(QObject * parent = 0) : QObject(parent){m_selection = NULL};
	~AbstractModel(){ if(m_selection) delete m_selection; };

	virtual int NumberOfColumns() = 0;						//Number of Features
	virtual int NumberOfRows() = 0;							//Number of Objects
	virtual int NumberOfTables() = 0;						//Number of Time Points
	virtual std::string ColumnName(int col) = 0;			//Feature name

	virtual bool IsColumnEditable(int col) = 0;				//Is this column editable?
	virtual bool IsColumnEditable(std::string colName) = 0;
					
	virtual const vtkTable * GetTable(int t) = 0;
	virtual const vtkVariantArray * GetRow(int t, int row) = 0;
	virtual const vtkAbstractArray * GetColumn(int t, int col) = 0;
	virtual const vtkAbstractArray * GetColumn(int t, std::string colName) = 0;
	virtual vtkVariant GetValue(int t, int row, int col) = 0;
	virtual vtkVariant GetValue(int t, int row, std::string colName) = 0;

	virtual void SetValue(int t, int row, int col, vtkVariant value) = 0;
	virtual void SetValue(int t, int row, std::string colName, vtkVariant value) = 0;
	virtual void AddColumn(std::string colName) = 0;

	virtual void RemoveColumn(int col) = 0;
	virtual void RemoveColumn(std::string colName) = 0;

	ObjectSelection * GetObjectSelection(){ return m_selection; };
	
signals:
	void Changed(int t, int row);							//use -1 for all tables or all rows

protected:
	ObjectSelection * m_selection;

private:

};

#endif