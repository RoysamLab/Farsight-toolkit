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
	void s_modelChanged(const QModelIndex & topLeft, const QModelIndex & bottomRight);
	void modelChanged();

public slots:
	void updateMapping(void);
	void deleteTrigger(void);
	void mergeTrigger(void);
	void splitTrigger(void);
	void addTrigger(void);
	void startSplitTrigger(void);
	void endSplitTrigger(void);

	bool isSplittingMode() { return SplittingMode; };
	void addPointToSplitList(int x, int y, int z);
	int getSizeOfSplittingList() { return (int)pointsForSplitting.size(); };
	std::vector<int> getPointFromSplittingList(int id);

	bool isAddMode() { return addMode; };
	int getSizeOfAddList(){ return (int)pointsForAdd.size(); };
	void addPointToAddList(int x, int y, int z);
	std::vector<int> getStartAddPoint();
	void abortAdd();




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
	void MostDataInModelChanged(void);
	void SyncModel();
	void UpdateColors();
	bool neighbor(ftk::Object *obj1, ftk::Object *obj2);

	//added by Yousef 7-30-2009
	bool SplittingMode; //to be used to indicate that we are in splitting mode
	bool addMode;  //indicates we are currently attempting to add an object (in add mode)
	std::vector<ftk::Object::Point> pointsForSplitting;
	std::vector<ftk::Object::Point> pointsForAdd;
};

#endif
