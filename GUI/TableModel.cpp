#include "TableModel.h"
#include <vtkStringArray.h>
#include <fstream>

//Constructor
TableModel::TableModel(QObject * parent) 
: AbstractModel(parent)
{
	m_Tables.clear();
	int m_numColumns = 0;
	int m_numRows = 0;
}

TableModel::~TableModel()
{
	for(int t=0; t<(int)m_Tables.size(); ++t)
	{
		m_Tables.at(t)->Delete();
	}
	m_Tables.clear();
}


std::string TableModel::ColumnName(int col)
{
	if(m_Tables.size() > 0)
		return m_Tables.at(0)->GetColumnName(col);
	else
		return "";
}

					
const vtkTable * TableModel::GetTable(int t)
{
	if((int)m_Tables.size() > t)
		return m_Tables.at(t);
	else
		return NULL;
}

const vtkVariantArray * TableModel::GetRow(int t, int row)
{
	if((int)m_Tables.size() > t && m_NumRows > row)
		return m_Tables.at(t)->GetRow(row);
	else
		return NULL;
}

const vtkAbstractArray * TableModel::GetColumn(int t, int col)
{
	if((int)m_Tables.size() > t && m_NumColumns > col)
		return m_Tables.at(t)->GetColumn(col);
	else
		return NULL;
}

const vtkAbstractArray * TableModel::GetColumn(int t, std::string colName)
{
	if((int)m_Tables.size() > t)
		return m_Tables.at(t)->GetColumnByName(colName.c_str());
	else
		return NULL;
}

vtkVariant TableModel::GetValue(int t, int row, int col)
{
	if((int)m_Tables.size() > t && m_NumRows > row && m_NumColumns > col)
		return m_Tables.at(t)->GetValue(row,col);
	else
		return vtkVariant(0);
}

vtkVariant TableModel::GetValue(int t, int row, std::string colName)
{
	if((int)m_Tables.size() > t && m_NumRows > row)
		return m_Tables.at(t)->GetValueByName(row,colName.c_str());
	else
		return vtkVariant(0);
}

void TableModel::SetValue(int t, int row, int col, vtkVariant value)
{
	if((int)m_Tables.size() > t && m_NumRows > row && m_NumColumns > col)
		if( IsColumnEditable(col) )
		{
			m_Tables.at(t)->SetValue(row, col, value);
			emit Changed(t,row);
		}
}

void TableModel::SetValue(int t, int row, std::string colName, vtkVariant value)
{
	if((int)m_Tables.size() > t && m_NumRows > row)
		if( IsColumnEditable(colName) )
		{
			m_Tables.at(t)->SetValueByName(row, colName.c_str(), value);
			emit Changed(t,row);
		}
}

//Load from a file
//File format is simple white space delimited.
//first row should be feature names.
//Assumed that all features are text
void TableModel::LoadTable(std::string filename)
{
	const int MAXLINESIZE = 1024;
	char line[MAXLINESIZE];

	//NOW LOAD THE DATA FROM FILE:
	ifstream file; 
	file.open( filename.c_str() );
	if ( !file.is_open() )
	{
		std::cerr << "Failed to Load Document: " << file << std::endl;
		return;
	}
	
	vtkTable *table = vtkTable::New();

	//Get the field names:
	file.getline(line, MAXLINESIZE, '\n');
	char * pch = strtok (line," \t");
	while (pch != NULL)
	{
		vtkStringArray * arr = vtkStringArray::New();
		arr->SetName(pch);
		table->AddColumn(arr);
		arr->Delete();
		pch = strtok (NULL, " \t");
	}
	
	//Get the data:
	file.getline(line, MAXLINESIZE, '\n');
	while ( !file.eof() )
	{
		vtkVariantArray * arr = vtkVariantArray::New();
		pch = strtok (line," \t");
		while (pch != NULL)
		{
			arr->InsertNextValue( vtkVariant(pch) );
			pch = strtok (NULL, " \t");
		}
		table->InsertNextRow(arr);
		arr->Delete();
		file.getline(line, MAXLINESIZE, '\n');
	}
	file.close();

	//Check to be sure that this new table matches size of previous tables
	if( m_Tables.size() > 0 )
	{
		if( m_NumColumns == table->GetNumberOfColumns() &&
			m_NumRows == table->GetNumberOfRows() )
		{
			m_Tables.push_back(table);
			m_selection->IncrementTMax();
			emit Changed(m_Tables.size()-1, -1 );
		}
	}
	else
	{
		m_Tables.push_back(table);
		m_selection->IncrementTMax();
		m_NumColumns = table->GetNumberOfColumns();
		m_NumRows = table->GetNumberOfRows();
		emit Changed(m_Tables.size()-1, -1 );
	}
}