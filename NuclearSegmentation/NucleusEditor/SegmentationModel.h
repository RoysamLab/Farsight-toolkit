#ifndef SEGMENTATIOINMODEL_H
#define SEGMENTATIOINMODEL_H

#include "FarsightConfig.h"

//QT INCLUDES
#include <QtCore/QObject>
#include <QtGui/QStandardItemModel>
#include <QtGui/QItemSelectionModel>
#include <QtGui/QColor>

//OUTSIDE INCLUDES
//#include <SegmentationCommon/ftkSegmentationResult.h>
#include <NuclearSegmentation/ftkNuclearSegmentation.h>


//**************************************************************************
// This class constructs the selectionModel, and model used 
// in the QT model/view framework to view a segmentationResult.
// 
//**************************************************************************
class SegmentationModel : public QObject
{
	Q_OBJECT

public:
	
	SegmentationModel(ftk::NuclearSegmentation *segresult);
	~SegmentationModel();

	void SetOutliers( vector<int> o );
	void ShowOutliers(bool show);
	bool HasOutliers(void);

	QStandardItemModel *GetModel(){ return model; };
	QItemSelectionModel *GetSelectionModel(){ return selectionModel; };

	int ColumnForID(){ return columnForID; };
	int RowForID(int id);
	int NumFeatures(){ return numFeatures; };
	int NumObjects(){ return numObjects; };
	int ColumnForColor(){ return columnForColor; };
	QMap<int, QColor> ColorMap(){ return colorMap; };
	ftk::NuclearSegmentation *SegResult(void){ return segResult; };

signals:
	void modelChanged(void);

public slots:
	void deleteTrigger(void);
	void mergeTrigger(void);

private:
	int columnForID;
	int numFeatures;
	int numObjects;
	int columnForClass;
	int columnForOutliers;
	int columnForColor;
	QMap<int,QColor> colorMap;
	QMap<int, int> LabelToRowMap;		//A label to a row in the model

	ftk::NuclearSegmentation *segResult;

	//I control these to make sure they line up with the information in segmentation Result.
	QStandardItemModel *model;
	QItemSelectionModel *selectionModel;

	//Functions:
	void SyncModel();
	void UpdateColors();
	void updateMapping();
	bool neighbor(ftk::Object *obj1, ftk::Object *obj2);
};

#endif
