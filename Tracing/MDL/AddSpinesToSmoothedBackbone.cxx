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

#include "vtkAppendPolyData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkSmartPointer.h"

#include <iostream>
using std::cerr;
using std::endl;

#define VTK_CREATE(type, var) vtkSmartPointer<type> var = vtkSmartPointer<type>::New()

int main(int argc, char **argv)
  {
  if(argc < 4)
    {
    cerr << argv[0] << " <backbone .vtk file> <spines .vtk file> "
         << "<output .vtk file>" << endl; 
    return 0;
    }

	VTK_CREATE(vtkPolyData, backbone);
	VTK_CREATE(vtkPolyData, spines);
	VTK_CREATE(vtkPolyData, output);
	VTK_CREATE(vtkPolyDataReader, reader);
	VTK_CREATE(vtkPolyDataWriter, writer);
  VTK_CREATE(vtkAppendPolyData, append);

  reader->SetFileName(argv[1]);
  reader->Update();
  backbone->DeepCopy(reader->GetOutput());

  reader->SetFileName(argv[2]);
  reader->Update();
  spines->DeepCopy(reader->GetOutput());

  append->AddInput(backbone);
  append->AddInput(spines);
  output = append->GetOutput();

  writer->SetInput(output);
  writer->SetFileName(argv[3]);
  writer->Update();

  return 1;
  }

