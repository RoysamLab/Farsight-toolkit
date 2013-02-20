#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkRGBPixel.h"

#include "vtkAdjacentVertexIterator.h"
#include "vtkDataSetAttributes.h"
#include "vtkGraphReader.h"
#include "vtkGraphWriter.h"
#include "vtkIdTypeArray.h"
#include "vtkImageData.h"
#include "vtkPoints.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkObjectFactory.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyLine.h"
#include "vtkSmartPointer.h"
#include "vtkTIFFReader.h"
#include "vtkTree.h"

#include <vtksys/SystemTools.hxx>

#include "vtkPlotEdges.h"
#include "vtkConvertMDLTracesToTrees.h"

#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkShapeLabelMapFilter.h" 

#include "itkLabelMapToBinaryImageFilter.h"

typedef itk::RGBPixel<unsigned char> RGBPixelType;
typedef itk::Image<RGBPixelType, 3> RGBImageType;

vtkStandardNewMacro(vtkConvertMDLTracesToTrees);
vtkCxxRevisionMacro(vtkConvertMDLTracesToTrees, "$Revision: 1107 $");

////////////////////////////////////////////////////////////////////////////////
vtkConvertMDLTracesToTrees::vtkConvertMDLTracesToTrees()
{
  this->TracesData = NULL;
  this->Somas = NULL;
}

////////////////////////////////////////////////////////////////////////////////
vtkConvertMDLTracesToTrees::~vtkConvertMDLTracesToTrees()
{
  for(vtkstd::vector<vtkstd::pair<vtkTree *, int> >::iterator
      traceItr = this->Traces.begin();
      traceItr != this->Traces.end();
      ++traceItr)
    {
    (*traceItr).first->Delete();
    }
  vtkstd::map<vtkIdType, vtkstd::pair<double *, int> >::iterator rootItr;
  for(rootItr = this->TraceRoots.begin();
      rootItr != this->TraceRoots.end(); ++rootItr)
    {
    delete rootItr->second.first;
    }
  this->TracesData->Delete();
}

////////////////////////////////////////////////////////////////////////////////
void vtkConvertMDLTracesToTrees::PrintSelf(ostream& os, vtkIndent indent)
{
}

////////////////////////////////////////////////////////////////////////////////
void vtkConvertMDLTracesToTrees::FindTraceRoots()
{
  this->TraceRoots.clear();
  vtkstd::map<vtkIdType, vtkstd::pair<double *, int> > candidateRoots;
  //step 1: find all the traces who have an end point in the soma
  for(int cellID = 0; cellID < this->TracesData->GetNumberOfCells();
      cellID++)
    {
    double *first_point = new double[3];
    double *last_point = new double[3];
    this->GetTraceEndPoint(cellID, true, first_point);
    this->GetTraceEndPoint(cellID, false, last_point);
    unsigned long somaAtPoint = this->GetSomaAtPoint(first_point);
    if(somaAtPoint != 0)
      {
      this->SomaHasTraces[somaAtPoint] = true;
      candidateRoots[cellID] = vtkstd::make_pair(last_point, somaAtPoint);
      delete [] first_point;
      }
    else
      {
      somaAtPoint = this->GetSomaAtPoint(last_point);
      if(somaAtPoint != 0)
        {
        this->SomaHasTraces[somaAtPoint] = true;
        candidateRoots[cellID] = vtkstd::make_pair(first_point, somaAtPoint);
        delete [] last_point;
        }
      else
        {
        delete [] first_point;
        delete [] last_point;
        }
      }
    }

  //step 2: loop over the candidate roots to determine which are true trace
  //roots.
  //
  //Note to self: this whole process seems inefficient because we're basically
  //doing the same search twice (once here and once in GenerateTrace), but it is
  //important that we start with the right list of roots to build our trees...
  //vtkstd::list<vtkstd::pair<vtkIdType, double *> >::iterator candItr;
  vtkstd::map<vtkIdType, vtkstd::pair<double *, int> >::iterator candItr;
  for(candItr = candidateRoots.begin();
      candItr != candidateRoots.end(); ++candItr)
    {
    if(this->VerifyRootTrace( candItr->first ))
      {
      this->TraceRoots[candItr->first] = candItr->second;
      this->Orphans->SetValue(candItr->first, -1);
      }
    else
      {
      delete [] candItr->second.first;
      }
    }
  cout << "number of root traces found: " << this->TraceRoots.size() << endl;
}

////////////////////////////////////////////////////////////////////////////////
void vtkConvertMDLTracesToTrees::GetTraceEndPoint(int cellID,
                                                  bool getFirstPoint,
                                                  double *point)
{
  vtkPolyLine *polyLine = reinterpret_cast<vtkPolyLine *>
    (this->GetTracesData()->GetCell(cellID));
  vtkSmartPointer<vtkPoints> points = polyLine->GetPoints();
  if(getFirstPoint)
    {
    points->GetPoint(0, point);
    }
  else
    {
    points->GetPoint(points->GetNumberOfPoints() - 1, point);
    }
}

////////////////////////////////////////////////////////////////////////////////
unsigned long vtkConvertMDLTracesToTrees::GetSomaAtPoint(double point[3])
{
  LabelMapType::IndexType index;
  LabelMapType::PointType position;
  LabelMapType::PixelType pixel;
  position[0] = point[0];
  position[1] = point[1];
  if(position[1] < 0)
    {
    position[1] *= -1.0;
    }
  position[2] = point[2];
  this->Somas->TransformPhysicalPointToIndex(position, index);

  //check if the soma image is non-zero at this point
  pixel = this->Somas->GetPixel(index);
  if(pixel != 0)
    {
    //If it is you're at the root of a trace.  Return the soma ID at this point
    unsigned long matchedSoma = this->Somas->GetLabelObject(index)->GetLabel();
    this->SomaHasTraces[matchedSoma] = true;
    return matchedSoma;
    }
  return 0;
}

//Determine whether or not a trace is connected to other traces on both
//ends.  True trace roots do not have neighbors on both ends.
////////////////////////////////////////////////////////////////////////////////
bool vtkConvertMDLTracesToTrees::VerifyRootTrace(int rootID)
{
  bool headNeighborFound = false;
  bool tailNeighborFound = false;
  vtkPolyLine *polyLine = reinterpret_cast<vtkPolyLine *>
    (this->GetTracesData()->GetCell(rootID));
  vtkSmartPointer<vtkPoints> points = polyLine->GetPoints();
  double first_point[3];
  double last_point[3];
  points->GetPoint(0, first_point);
  points->GetPoint(points->GetNumberOfPoints() - 1, last_point);
  for(int cellID = 0; cellID < this->TracesData->GetNumberOfCells();
      cellID++)
    {
    if(cellID == rootID)
      {
      continue;
      }
    if(!headNeighborFound)
      {
      if(this->TraceContainsPoint(cellID, first_point))
        {
        headNeighborFound = true;
        }
      }
    if(!tailNeighborFound)
      {
      if(this->TraceContainsPoint(cellID, last_point))
        {
        tailNeighborFound = true;
        }
      }
    if(headNeighborFound && tailNeighborFound)
      {
      return false;
      }
    }
  if(!(headNeighborFound && tailNeighborFound))
    {
    return true;
    }
  cerr << "WARNING: logic dictates that we should never go here..." << endl;
  return false;
}

//returns true if point is one of the end points of the trace identified by
//cellID
////////////////////////////////////////////////////////////////////////////////
bool vtkConvertMDLTracesToTrees::TraceContainsPoint(int cellID, double point[3])
{
  vtkSmartPointer<vtkPolyLine> trace = reinterpret_cast<vtkPolyLine *>
      (this->GetTracesData()->GetCell(cellID));
  vtkSmartPointer<vtkPoints> points = trace->GetPoints();
  double start[3];
  double end[3];
  points->GetPoint(0, start);
  points->GetPoint(points->GetNumberOfPoints() - 1, end);
  if(start[0] == point[0] && start[1] == point[1] &&
     start[2] == point[2])
    {
    return true;
    }
  if(end[0] == point[0] && end[1] == point[1] &&
     end[2] == point[2])
    {
    return true;
    }
  return false;
}

////////////////////////////////////////////////////////////////////////////////
int vtkConvertMDLTracesToTrees::CreateTrees()
{
  if(this->TracesData == NULL)
    {
    cerr << "ERROR: Call LoadTracesData before CreateTrees." << endl;
    return 1;
    }

  if(this->Somas.IsNull())
    {
    cerr << "ERROR: Call LoadSomas before CreateTrees." << endl;
    return 1;
    }

  //initialize "orphans": a list of lines whose parents have not been found yet.
  this->Orphans = vtkSmartPointer<vtkIdTypeArray>::New();
  this->Orphans->SetNumberOfValues(this->TracesData->GetNumberOfCells());

  for(int cellID = 0; cellID < this->GetTracesData()->GetNumberOfCells();
      cellID++)
    {
    this->Orphans->SetValue(cellID, cellID);
    }

  //find all the lines that are rooted in somas
  cout << "searching for trace roots" << endl;
  this->FindTraceRoots();

  cout << "generating traces" << endl;
  vtkstd::map<vtkIdType, vtkstd::pair<double *, int> >::iterator rootItr;
  for(rootItr = this->TraceRoots.begin();
      rootItr != this->TraceRoots.end(); ++rootItr)
    {
      //create a new tree for this root trace
      vtkSmartPointer<vtkMutableDirectedGraph> trace =
        vtkSmartPointer<vtkMutableDirectedGraph>::New();
      vtkIdType root = trace->AddVertex();

      vtkIdTypeArray *nodeIDsToLineIDs = vtkIdTypeArray::New(); 
      nodeIDsToLineIDs->InsertValue(root, -1);
      nodeIDsToLineIDs->SetName("idmap");

      //call GenerateTrace to flesh out this new trace
      this->GenerateTrace(trace, root, rootItr->first, rootItr->second.first,
                        nodeIDsToLineIDs); 

      //associate a node to line mapping with this tree and add it to the
      //vector of traces
      trace->GetVertexData()->AddArray(nodeIDsToLineIDs);
      vtkTree *tree = vtkTree::New();
      tree->CheckedShallowCopy(trace);
      nodeIDsToLineIDs->Delete();
      this->Traces.push_back(vtkstd::make_pair(tree, rootItr->second.second));
      }
  this->CheckForOrphans();
  return 1;
}

//this function recursively searches for lines that are connected to each
//other's end points, growing a tree in the process.
////////////////////////////////////////////////////////////////////////////////
void vtkConvertMDLTracesToTrees::GenerateTrace(vtkMutableDirectedGraph *trace,
                                             vtkIdType parent, int cellID,
                                             double *end_point,
                                             vtkIdTypeArray *nodeIDsToLineIDs)
{
  vtkIdType newNode = trace->AddChild(parent);
  nodeIDsToLineIDs->InsertValue(newNode, cellID);

  double first_point[3];
  double last_point[3];
  vtkPolyLine *polyLine;
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  //for each that hasn't found its parent yet
  for(int lineID = 0; lineID < this->GetTracesData()->GetNumberOfCells();
      lineID++)
    {
    if( (cellID == lineID) || (this->Orphans->GetValue(lineID) == -1) )
      {
      continue;
      }

    polyLine = reinterpret_cast<vtkPolyLine *>(this->GetTracesData()->GetCell(lineID));
    points = polyLine->GetPoints();
    points->GetPoint(0, first_point);
    points->GetPoint(points->GetNumberOfPoints() - 1, last_point);
    //if either of its end points match the end point of this new node
    //create a new node associated with that line, with this node as its parent
    if(first_point[0] == end_point[0] && first_point[1] == end_point[1] &&
       first_point[2] == end_point[2])
      {
      if(lineID == 452)
        {
        cout << "I should never go here 1" << endl;
        }
      this->Orphans->SetValue(lineID, -1);
      this->GenerateTrace(trace, newNode, lineID, last_point, nodeIDsToLineIDs); 
      }
    if(last_point[0] == end_point[0] && last_point[1] == end_point[1] &&
       last_point[2] == end_point[2])
      {
      if(lineID == 452)
        {
        cout << "I should never go here 2" << endl;
        }
      this->Orphans->SetValue(lineID, -1);
      this->GenerateTrace(trace, newNode, lineID, first_point, nodeIDsToLineIDs); 
      }
    }
}

////////////////////////////////////////////////////////////////////////////////
void vtkConvertMDLTracesToTrees::CheckForOrphans()
{
  int numOrphans = 0;
  for(int i = 0; i < this->Orphans->GetSize(); i++)
    {
    if(this->Orphans->GetValue(i) != -1)
      {
      numOrphans++;
      }
    }
  if(numOrphans > 0)
    {
    cout << numOrphans << " traces could not be connected to a soma." << endl;
    }
}

//filename should point to a .vtk file that contains the output of a tracing
//routine
////////////////////////////////////////////////////////////////////////////////
void vtkConvertMDLTracesToTrees::LoadTracesData(const char *filename)
{
  vtkSmartPointer<vtkPolyDataReader> tracesReader = 
    vtkSmartPointer<vtkPolyDataReader>::New();
  tracesReader->SetFileName(filename);

  vtkSmartPointer<vtkPlotEdges> plotEdges =
    vtkSmartPointer<vtkPlotEdges>::New();
  plotEdges->SetInput(tracesReader->GetOutput());
  cout << "Running vtkPlotEdges.  Depending on the size of your traces data "
       << "this may take a while." << endl;
  plotEdges->Update();

  this->SetTracesData(plotEdges->GetOutput());
}

//filename should point to a .tif file that contains the output of a soma
//segementation routine
////////////////////////////////////////////////////////////////////////////////
void vtkConvertMDLTracesToTrees::LoadSomas(const char *filename)
{
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(filename);

  //we use a BinaryImageToShapeLabelMapFilter to convert the binary
  //image into a collection of objects
  typedef itk::BinaryImageToShapeLabelMapFilter< ImageType, LabelMapType >
    ConverterType;
  ConverterType::Pointer converter = ConverterType::New();
  converter->SetInput(reader->GetOutput());
  this->Somas = converter->GetOutput();
  converter->Update();

  //initialize the label -> bool SomaHasTraces map to all false
  //these values will be switched to true in FindTraceRoots
  this->SomaHasTraces.clear();
  unsigned int numLabels = this->Somas->GetNumberOfLabelObjects();
  for(unsigned int label=1; label<= numLabels; label++)
    {
    this->SomaHasTraces[label] = false;
    }
}

//write out the trees generated by this class to a series of .vtk files
////////////////////////////////////////////////////////////////////////////////
void vtkConvertMDLTracesToTrees::WriteTrees(const char *basefilename)
{
  vtkSmartPointer<vtkGraphWriter> graphWriter =
    vtkSmartPointer<vtkGraphWriter>::New();
  vtkstd::string outFileName = vtkstd::string(basefilename);
  vtkstd::string basename = outFileName.substr(0, outFileName.rfind(".vtk"));
  vtkstd::string currentFileName;
  int itr = 0;
  for(vtkstd::vector<vtkstd::pair<vtkTree *, int> >::iterator
      traceItr = this->Traces.begin();
      traceItr != this->Traces.end();
      ++traceItr)
    {
    graphWriter->SetInput( (*traceItr).first );
    currentFileName = basename;
    vtkstd::ostringstream os;
    os << itr;
    currentFileName += os.str();
    currentFileName += ".vtk";
    graphWriter->SetFileName(currentFileName.c_str());
    graphWriter->Update();
    itr++;
    }
}

//write out a separate .swc file for each cell
////////////////////////////////////////////////////////////////////////////////
void vtkConvertMDLTracesToTrees::WriteCellsToSWCFiles(const char *basefilename)
{
  if(this->Traces.size() < 1)
    {
    cout << "No traces to write, aborting. " << endl;
    return;
    }
  vtkstd::string outFileName = vtkstd::string(basefilename);
  vtkstd::string basename = outFileName.substr(0, outFileName.rfind(".swc"));
  vtkstd::string currentFileName;
  int itr = 0;
  unsigned int numLabels = this->Somas->GetNumberOfLabelObjects();
  for(unsigned int label=1; label<= numLabels; label++)
    {

    //skip somas that don't have any traces coming out of them.
    if(this->SomaHasTraces[label] == false)
      {
      continue;
      }

    //generate a new filename and open it for writing
    currentFileName = basename;
    vtkstd::ostringstream os;
    os << itr;
    currentFileName += "-";
    currentFileName += os.str();
    currentFileName += ".swc";
    this->FileWriter.open(currentFileName.c_str());
    if(!this->FileWriter.is_open())
      {
      cerr << "unable to open " << currentFileName.c_str() << endl;
      return;
      }

    //reinitialize currentSWCID to zero for each new file
    long int currentSWCID = 0;

    //write out this soma's centroid as the first point
    const LabelObjectType * labelObject = this->Somas->GetLabelObject(label);
    LabelObjectType::CentroidType centroid = labelObject->GetCentroid(); 
    this->FileWriter << ++currentSWCID << " 1 " << centroid[0] << " "
                     << centroid[1] << " " << centroid[2] << " 0 " << "-1"
                     << endl;
    long int somaSWCID = currentSWCID;

    //iterate over the list of traces, searching for those who connect to this
    //soma
    for(vtkstd::vector<vtkstd::pair<vtkTree *, int> >::iterator
        traceItr = this->Traces.begin();
        traceItr != this->Traces.end();
        ++traceItr)
      {
      if( (*traceItr).second != (int)label)
        {
        continue;
        }

    //find the xyz point that connects the trace to the soma
    vtkTree *tree = (*traceItr).first;
    vtkIdTypeArray *vertexArray = reinterpret_cast<vtkIdTypeArray *>
      (tree->GetVertexData()->GetArray(0));
    vtkIdType cellID = vertexArray->GetValue(1);
    double *nonSomaEndPoint = this->TraceRoots[cellID].first;
    vtkPolyLine *traceLine = reinterpret_cast<vtkPolyLine *>
      (this->GetTracesData()->GetCell(cellID));
    vtkPoints *tracePoints = traceLine->GetPoints();
    double xyz[3];
    tracePoints->GetPoint(0, xyz);
    if(xyz[0] == nonSomaEndPoint[0] && xyz[1] == nonSomaEndPoint[1] &&
       xyz[2] == nonSomaEndPoint[2])
      {
      tracePoints->GetPoint(tracePoints->GetNumberOfPoints() - 1, xyz);
      }

    //don't write out the root of the tree, as it isn't associated with
    //a vtkPolyLine
    this->WriteTraceToSWC(tree, 1, &currentSWCID,
                          somaSWCID, xyz);

      }
    this->FileWriter.close();
    itr++;
    }
}

//write out the tree structure generated by this class to .swc format
////////////////////////////////////////////////////////////////////////////////
void vtkConvertMDLTracesToTrees::WriteToSWC(const char *filename)
{
  this->FileWriter.open(filename);
  if(!this->FileWriter.is_open())
    {
    cerr << "unable to open " << filename << endl;
    return;
    }

  cout << "Writing " << this->Traces.size() << " trace(s) to " << filename
       << endl;

  long int currentSWCID = 0;

  //write out the somas first.
  unsigned int numLabels = this->Somas->GetNumberOfLabelObjects();
  vtkstd::vector<int> somaSWCIDs(numLabels + 1);
  for(unsigned int label=1; label<= numLabels; label++)
    {
    //skip somas that don't have any traces coming out of them.
    if(this->SomaHasTraces[label] == false)
      {
      continue;
      }
    // we don't need a SmartPointer of the label object here, because the
    // reference is kept in the label map.
    const LabelObjectType * labelObject = this->Somas->GetLabelObject(label);
    LabelObjectType::CentroidType centroid = labelObject->GetCentroid(); 
    this->FileWriter << ++currentSWCID << " 1 " << centroid[0] << " "
                     << centroid[1] << " " << centroid[2] << " 0 " << "-1"
                     << endl;
    //keep track of this mapping between somas and SWC IDs.
    somaSWCIDs[label] = currentSWCID;
    }

  //write out all the traces
  vtkTree *tree;
  int somaID = -1;
  int somaSWCID = -1;
  for(vtkstd::vector<vtkstd::pair<vtkTree *, int> >::iterator
      traceItr = this->Traces.begin();
      traceItr != this->Traces.end();
      ++traceItr)
    {
    tree = (*traceItr).first;
    somaID = (*traceItr).second;
    somaSWCID = somaSWCIDs[somaID];

    //find the xyz point that connects the trace to the soma
    vtkIdTypeArray *vertexArray = reinterpret_cast<vtkIdTypeArray *>
      (tree->GetVertexData()->GetArray(0));
    vtkIdType cellID = vertexArray->GetValue(1);
    double *nonSomaEndPoint = this->TraceRoots[cellID].first;
    vtkPolyLine *traceLine = reinterpret_cast<vtkPolyLine *>
      (this->GetTracesData()->GetCell(cellID));
    vtkPoints *tracePoints = traceLine->GetPoints();
    double xyz[3];
    tracePoints->GetPoint(0, xyz);
    if(xyz[0] == nonSomaEndPoint[0] && xyz[1] == nonSomaEndPoint[1] &&
       xyz[2] == nonSomaEndPoint[2])
      {
      tracePoints->GetPoint(tracePoints->GetNumberOfPoints() - 1, xyz);
      }

    //don't write out the root of the tree, as it isn't associated with
    //a vtkPolyLine
    this->WriteTraceToSWC(tree, 1, &currentSWCID,
                          somaSWCID, xyz);
    }
}

//recursive helper function for WriteToSWC.  Note that this is currently
//hardcoding radius to 1 and class to 3.
////////////////////////////////////////////////////////////////////////////////
void vtkConvertMDLTracesToTrees::WriteTraceToSWC(vtkTree *trace,
                                                vtkIdType currentNode,
                                                long *currentSWCID,
                                                long parentSWCID,
                                                double *parentEndPoint)
{
  //get the line associated with this node in the tree
  vtkIdTypeArray *vertexArray = reinterpret_cast<vtkIdTypeArray *>
    (trace->GetVertexData()->GetArray(0));
  vtkIdType cellID = vertexArray->GetValue(currentNode);
  //special case: the trace's root isn't associated with a line of its own.
  vtkPolyLine *currentLine = reinterpret_cast<vtkPolyLine *>
    (this->GetTracesData()->GetCell(cellID));
  vtkPoints *currentPoints = currentLine->GetPoints();

  //figure out which end of this trace connects to the parent
  bool firstPointConnectsToParent;
  double xyz[3];
  currentPoints->GetPoint(0, xyz);
  if(xyz[0] == parentEndPoint[0] && xyz[1] == parentEndPoint[1] &&
     xyz[2] == parentEndPoint[2])
    {
    if(cellID == 452)
      {
      cout << "my first point connects to my parent: " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << endl;
      }
    firstPointConnectsToParent = true;
    }
  else
    {
    currentPoints->GetPoint(currentPoints->GetNumberOfPoints() - 1, xyz);
    if(xyz[0] == parentEndPoint[0] && xyz[1] == parentEndPoint[1] &&
       xyz[2] == parentEndPoint[2])
      {
      firstPointConnectsToParent = false;
      }
    else
      {
      cout << "BUG: neither my first nor last point connects to my parent!" << endl;
      cout << "parent SWC == " << parentSWCID << ", current ID == "
           << *currentSWCID << ", cell ID == " << cellID << endl;
      }
    }

  //connect the first point to the parent's SWC ID
  if(firstPointConnectsToParent)
    {
    currentPoints->GetPoint(1, xyz);
    }
  else
    {
    currentPoints->GetPoint(currentPoints->GetNumberOfPoints() - 2, xyz);
    }
  this->FileWriter << ++(*currentSWCID) << " 3 " << xyz[0] << " " << xyz[1]
                   << " " << xyz[2] << " 1 " << parentSWCID << endl;

  //now iterate over the trace's points, writing out each one, connecting it to
  //the previous point.
  if(firstPointConnectsToParent)
    {
    for(vtkIdType pointID = 2; pointID < currentPoints->GetNumberOfPoints();
        pointID++)
      {
      currentPoints->GetPoint(pointID, xyz);
      this->FileWriter << ++(*currentSWCID) << " 3 " << xyz[0] << " " << xyz[1]
                       << " " << xyz[2] << " 1 " << (*currentSWCID) << endl;
      }
    }
  else
    {
    for(vtkIdType pointID = currentPoints->GetNumberOfPoints() - 3; pointID > -1;
        pointID--)
      {
      currentPoints->GetPoint(pointID, xyz);
      this->FileWriter << ++(*currentSWCID) << " 3 " << xyz[0] << " " << xyz[1]
                       << " " << xyz[2] << " 1 " << (*currentSWCID) << endl;
      }
    }
 
  //then call this function on the children
  long branchPointID = *currentSWCID; 
  vtkAdjacentVertexIterator *childItr = vtkAdjacentVertexIterator::New();
  trace->GetChildren(currentNode, childItr);
  int numChildren = 0;
  while(childItr->HasNext())
    {
    numChildren++;
    this->WriteTraceToSWC(trace, childItr->Next(), currentSWCID, 
                          branchPointID, xyz);
    }
  childItr->Delete();
}

