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

#include "MDL2Wizard.h"
#include "MDL2WizardHelper.h"
#include "tinyxml/tinyxml.h"
#include <QtGui/QApplication>

int main (int argc, char* argv[])
  {
  if(argc < 2)
    {
    QApplication app(argc, argv);
    MDL2Wizard *wizard = new MDL2Wizard();
    wizard->show();
    return app.exec();
    }

  //otherwise run MDL non-interactively.  Check that we have all the required
  //parameters.
  if(argc < 3)
    {
    cerr << argv[0] << " <parameters.xml> <image1> [image2 image3...]"
         << endl;
    return 1;
    }

  //make sure the .xml file exists
  const char *xmlFileName = argv[1];
  QFileInfo xmlInfo(xmlFileName);
  if(!xmlInfo.exists())
    {
    cerr << "Error opening " << xmlFileName << endl;
    return 1;
    }

  //parse parameter values from xml
  TiXmlDocument doc(xmlFileName);
  doc.LoadFile();
  TiXmlHandle docHandle( &doc );

  TiXmlElement* parametersElement = docHandle.FirstChild("parameters").Element();
  if (!parametersElement)
    {
    cerr << "Aborting, no <parameters> tag found in " << xmlFileName << endl;
    return 1;
    }
  TiXmlElement* parameterElement =
    parametersElement->FirstChildElement("parameter");
  if(!parameterElement)
    {
    cerr << "Aborting, no <parameter> tag found in " << xmlFileName << endl;
    return 1;
    }
  //if we get this far, assume its a legit MDL parameters xml file
  QString ComponentsSize = QString(parameterElement->Attribute("value"));
  parameterElement = parameterElement->NextSiblingElement("parameter");
  QString VectorMagnitude = QString(parameterElement->Attribute("value"));
  parameterElement = parameterElement->NextSiblingElement("parameter");
  QString EdgeRange1 = QString(parameterElement->Attribute("value"));
  parameterElement = parameterElement->NextSiblingElement("parameter");
  QString MorphStrength1 = QString(parameterElement->Attribute("value"));
  parameterElement = parameterElement->NextSiblingElement("parameter");
  QString Order = QString(parameterElement->Attribute("value"));
  parameterElement = parameterElement->NextSiblingElement("parameter");
  QString Levels = QString(parameterElement->Attribute("value"));
  parameterElement = parameterElement->NextSiblingElement("parameter");
  QString EdgeRange2 = QString(parameterElement->Attribute("value"));
  parameterElement = parameterElement->NextSiblingElement("parameter");
  QString MorphStrength2 = QString(parameterElement->Attribute("value"));

  //make sure that this input file exists
  const char *inputFileName = argv[2];
  QFileInfo inputInfo(inputFileName);
  if(!inputInfo.exists())
    {
    cerr << "Error opening " << inputFileName << endl;
    return 1;
    }

  //initialize the QWizard
  QApplication app(argc, argv);
  MDL2Wizard *wizard = new MDL2Wizard();
  wizard->SetInteractiveExecution(false);
  wizard->ComponentsSizeInput->setText(ComponentsSize);
  wizard->VectorMagnitudeInput->setText(VectorMagnitude);
  wizard->EdgeRangeInput1->setText(EdgeRange1);
  wizard->MorphStrengthInput1->setText(MorphStrength1);
  wizard->OrderInput->setText(Order);
  wizard->LevelsInput->setText(Levels);
  wizard->EdgeRangeInput2->setText(EdgeRange2);
  wizard->MorphStrengthInput2->setText(MorphStrength2);

  //use the wizard to process this file
  cout << "Processing " << inputFileName << endl;
  wizard->restart();
  wizard->SetInputImage(inputFileName);
  QString backboneFileName = inputInfo.absolutePath() + "/" + inputInfo.completeBaseName();
  backboneFileName += "_backbone.vtk";
  wizard->SetBackboneFile(backboneFileName.toStdString().c_str());
  QString skelFileName = inputInfo.absolutePath() + "/" + inputInfo.completeBaseName();
  skelFileName += "_skel.vtk";
  wizard->SetSkeletonFile(skelFileName.toStdString().c_str());
  wizard->show();
  wizard->next();
  wizard->MaskUsingGraphCuts();
  app.exec();
  }

