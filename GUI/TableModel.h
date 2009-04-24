#ifndef _TABLEMODEL_H
#define _TABLEMODEL_H

#include "AbstractModel.h"

//****************************************************************************
// NOTES:
//  1. MUST emit changed when data changes!!
//**************************************************************************
class TableModel : public AbstractModel
{
	Q_OBJECT;

public:
	TableModel(QObject * parent = 0);
	~TableModel();

	virtual void LoadTable(std::string filename);					//Load from a file
	virtual void SaveTable(int t, std::string filename = ""){};		//Save to file

	int NumberOfColumns(){ return m_NumColumns; };			//Number of Features
	int NumberOfRows(){ return m_NumRows; };				//Number of Objects
	int NumberOfTables() { return (int)m_Tables.size(); };	//Number of Time Points
	std::string ColumnName(int col);						//Feature name

	virtual bool IsColumnEditable(int col) {return false;};			//Is this column editable?
	virtual bool IsColumnEditable(std::string colName) {return false;};
					
	const vtkTable * GetTable(int t);
	const vtkVariantArray * GetRow(int t, int row);
	const vtkAbstractArray * GetColumn(int t, int col);
	const vtkAbstractArray * GetColumn(int t, std::string colName);
	vtkVariant GetValue(int t, int row, int col);
	vtkVariant GetValue(int t, int row, std::string colName);

	void SetValue(int t, int row, int col, vtkVariant value);
	void SetValue(int t, int row, std::string colName, vtkVariant value);
	
	virtual void AddColumn(std::string colName, bool editable){};
	virtual void RemoveColumn(int col){};
	virtual void RemoveColumn(std::string colName){};
	
signals:
	//IN ABSTRACT CLASS:
	//void Changed(int t, int row);							//use -1 for all tables or all rows

protected:
	std::vector< vtkTable * > m_Tables;
	int m_NumColumns;
	int m_NumRows;

private:
	
};

#endif