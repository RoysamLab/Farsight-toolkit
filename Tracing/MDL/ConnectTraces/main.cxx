#include "vtkConvertMDLTracesToTrees.h"

int main(int argc, char **argv)
  {
  if(argc < 4)
    {
    cerr << argv[0] << " <traces .vtk file> <soma segmentation file> "
         << "<output .swc file>" << endl;
    return 0;
    }
  vtkSmartPointer<vtkConvertMDLTracesToTrees> treeBuilder =
    vtkSmartPointer<vtkConvertMDLTracesToTrees>::New();
  treeBuilder->LoadTracesData(argv[1]);
  treeBuilder->LoadSomas(argv[2]);
  treeBuilder->CreateTrees();
  treeBuilder->WriteCellsToSWCFiles(argv[3]);
  //treeBuilder->WriteToSWC(argv[3]);
  return 1;
  }

