/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkConvertMDLTracesToTrees.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkConvertMDLTracesToTrees - 
// .SECTION Description


#ifndef __vtkConvertMDLTracesToTrees_h
#define __vtkConvertMDLTracesToTrees_h

#include "itkImage.h"
#include "itkShapeLabelObject.h"

#include "vtkObject.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include <vtkstd/list>
#include <vtkstd/utility>
#include <vtkstd/vector>

class vtkIdTypeArray;
class vtkMutableDirectedGraph;
class vtkPolyLine;
class vtkTree;
typedef itk::Image<unsigned char, 3> ImageType;
typedef unsigned long LabelType;
typedef itk::ShapeLabelObject< LabelType, 3 > LabelObjectType;
typedef itk::LabelMap< LabelObjectType > LabelMapType;

class VTK_EXPORT vtkConvertMDLTracesToTrees : public vtkObject
{
public:
  vtkTypeRevisionMacro(vtkConvertMDLTracesToTrees, vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkConvertMDLTracesToTrees *New();

  vtkGetObjectMacro(TracesData, vtkPolyData);
  vtkSetObjectMacro(TracesData, vtkPolyData);

  void LoadSomas(const char *filename);
  void LoadTracesData(const char *filename);
  int CreateTrees();
  void WriteToSWC(const char *filename);
  void WriteTrees(const char *basefilename);
  void FindTraceRoots();
  void WriteCellsToSWCFiles(const char *filename);

protected:
  vtkConvertMDLTracesToTrees();
  virtual ~vtkConvertMDLTracesToTrees(); 
  void GenerateTrace(vtkMutableDirectedGraph *trace, vtkIdType parent,
                     int cellID, double *end_point,
                     vtkIdTypeArray *nodeIDsToLineIDs);
  void CheckForOrphans();
  void WriteTraceToSWC(vtkTree *trace, vtkIdType currentNode,
                       long *currentSWCID, long parentSWCID,
                       double *parentEndPoint);
  void GetTraceEndPoint(int cellID, bool getFirstPoint, double *point);
  unsigned long GetSomaAtPoint(double point[3]);
  bool VerifyRootTrace(int cellID);
  bool TraceContainsPoint(int cellID, double point[3]);

private:
  vtkConvertMDLTracesToTrees(const vtkConvertMDLTracesToTrees&);  // Not implemented.
  void operator=(const vtkConvertMDLTracesToTrees&);  // Not implemented.
  vtkstd::vector<vtkstd::pair<vtkTree *, int> > Traces;
  vtkSmartPointer<vtkIdTypeArray> Orphans;
  /**
   * a list of trace lines that are rooted in the somas
   * the vtkIdType (the key) is the cell id of the trace line
   * The first element of the pair (the double[3]) is the end point of the
   * trace line that is NOT in the soma.
   * The second element of the pair (the int) is the ID of the soma that the
   * trace connects to.
   **/
  vtkstd::map<vtkIdType, vtkstd::pair<double *, int> > TraceRoots;
  vtkstd::map<LabelType, bool> SomaHasTraces;
  vtkPolyData *TracesData;
  bool WriteGraphs;
  ofstream FileWriter;
  LabelMapType::Pointer Somas;
  ImageType::Pointer SomaImage;
};

#endif
