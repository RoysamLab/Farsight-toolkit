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

#include "MDLWizard.h"
#include <QtGui/QApplication>
//#include <QtCore/QObject>

int main (int argc, char* argv[])
  {
  QApplication app(argc, argv);
	MDLWizard *wizard = new MDLWizard();
  if(argc < 2)
    {
    wizard->show();
    }
  else
    {
    //run MDL non-interactively.  Check that we have all the required parameters.
    if(argc < 4)
      {
      cerr << argv[0] << " <input image> <backbone output> <spines output> "
           << "[connected components size] [vector magnitude] [morph strength] "
           << "[edge range] [weight factor] [BSpline order] [number of levels] "
           << "[graph prune size]" << endl;
      delete wizard;
      return 1;
      }
    else
      {
      wizard->SetInteractiveExecution(false);
      wizard->SetInputImage(argv[1]);
      wizard->SetBackboneFile(argv[2]);
      wizard->SetSpinesFile(argv[3]);
      if(argc > 4)
        {
        wizard->SetConnectedComponentsSize(argv[4]);
        }
      if(argc > 5)
        {
        wizard->SetVectorMagnitude(argv[5]);
        }
      if(argc > 6)
        {
        wizard->SetMorphStrength(argv[6]);
        }
      if(argc > 7)
        {
        wizard->SetEdgeRange(argv[7]);
        }
      if(argc > 8)
        {
        wizard->SetWeightFactor(argv[8]);
        }
      if(argc > 9)
        {
        wizard->SetBSplineOrder(argv[9]);
        }
      if(argc > 10)
        {
        wizard->SetBSplineLevels(argv[10]);
        }
      if(argc > 11)
        {
        wizard->SetGraphPruneSize(argv[11]);
        }
      wizard->show();
      wizard->next();
      wizard->RunVolumeProcess();
      }
    }
  int retval = app.exec();
  return retval;
  }
