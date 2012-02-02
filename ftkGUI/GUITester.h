/*=========================================================================
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

#ifndef __GUITESTER_H
#define __GUITESTER_H

#include <QWidget>
#include "vtkSmartPointer.h"

class pqTestUtility;
class vtkImageData;
class vtkRenderWindow;
class vtkTesting;

class GUITester : public QWidget
{
  Q_OBJECT
public:
  GUITester(QWidget *parent = 0);
  ~GUITester();
  bool playTestAndCompareResults( QString filename );
  void playTestFile( QString filename );
  bool compareResults();
  bool compareResults( QString testImgFileName );
  void SetBaselineImage(const char *fn);
  void SetRenderWindow(vtkRenderWindow *rw);
  void SetThreshold(double t);

public slots:
  void record();
  void play();

private:
  Q_DISABLE_COPY(GUITester)
  pqTestUtility *TestUtility;
  vtkSmartPointer<vtkTesting> Testing;
  double Threshold;
  bool BaselineSet;
  QString BaselineImage;
};

#endif
