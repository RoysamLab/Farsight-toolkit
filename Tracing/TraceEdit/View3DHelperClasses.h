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

class vtkSlider2DCallbackBrightness : public vtkCommand
{
public:
  static vtkSlider2DCallbackBrightness *New() 
    { return new vtkSlider2DCallbackBrightness; }
  virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
      vtkSliderWidget *sliderWidget = 
        reinterpret_cast<vtkSliderWidget*>(caller);
      double value = static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue();

      vtkSmartPointer<vtkPiecewiseFunction> opacityTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();
    opacityTransferFunction->AddPoint(2,0.0);
    //opacityTransferFunction->AddPoint(2+256*(1-value),0.2);
    opacityTransferFunction->AddPoint(50,value);
     this->volume->GetProperty()->SetScalarOpacity(opacityTransferFunction);
    }
  vtkSlider2DCallbackBrightness() {

  }

  vtkSmartPointer<vtkVolume> volume;
};

class vtkSlider2DCallbackContrast : public vtkCommand
{
public:
  static vtkSlider2DCallbackContrast *New() 
    { return new vtkSlider2DCallbackContrast; }
  virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
      vtkSliderWidget *sliderWidget = 
        reinterpret_cast<vtkSliderWidget*>(caller);
      double value = static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue();

      /*vtkSmartPointer<vtkPiecewiseFunction> colorTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();
    colorTransferFunction->AddPoint(1.0,0.0f,0.0f,0.0f);
    colorTransferFunction->AddPoint(255*(1-value),0.5f,0.5f,0.0f);*/
    vtkColorTransferFunction *colorTransferFunction =
      vtkSmartPointer<vtkColorTransferFunction>::New();
    colorTransferFunction->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
    colorTransferFunction->AddRGBPoint(255*(1-value),1,0,0);
     this->volume->GetProperty()->SetColor(colorTransferFunction);
    }
  vtkSlider2DCallbackContrast() {

  }

  vtkSmartPointer<vtkVolume> volume;
};

class vtkSlider2DCallbackContourThreshold : public vtkCommand
{
public:
  static vtkSlider2DCallbackContourThreshold *New() 
    { return new vtkSlider2DCallbackContourThreshold; }
  virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
      vtkSliderWidget *sliderWidget = 
        reinterpret_cast<vtkSliderWidget*>(caller);
      double value = static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue();
    printf("value of threshold = %d\n",int(value*255.0));
    cfilter->Print(std::cout);
    cfilter->GetInput()->Print(std::cout);
    double curr_value = cfilter->GetValue(0);
    if(fabs(curr_value-value*255)>10)
    {
    cfilter->SetValue(0,255*value);
    //ping me if this is really necessary and I'll show you a cross-platform
    //way to do it.  -Zack
    //Sleep(1000);
    }

    }
  vtkSlider2DCallbackContourThreshold() {

  }

  vtkSmartPointer<vtkContourFilter> cfilter;
};

class vtkSubRep: public vtkPlaybackRepresentation
{
public:
  static vtkSubRep * New(){ return new vtkSubRep;}
  virtual void Play() 
  {
    slice_counter++; 
    slice_counter = slice_counter % (im_pointer.size());
    vmapper->SetInput((im_pointer)[slice_counter]);
  }
  virtual void Stop() { 
    cout<< "stopped"; 
    slice_counter=0;
    vmapper->SetInput((im_pointer)[slice_counter]);
  }
  virtual void ForwardOneFrame() 
  {
    slice_counter++; 
    slice_counter = slice_counter % (im_pointer.size());
    vmapper->SetInput((im_pointer)[slice_counter]);
    cout << "forward one frame\n";
  }
  virtual void BackwardOneFrame() 
  {
    slice_counter--; 
    slice_counter = slice_counter+im_pointer.size();
    slice_counter = slice_counter % (im_pointer.size());
    vmapper->SetInput((im_pointer)[slice_counter]);
    cout << "backward one frame\n";
  }
  virtual void JumpToBeginning() 
  {
    slice_counter=0;
    vmapper->SetInput((im_pointer)[slice_counter]);
    cout << "jump to beginning\n";
  }
  virtual void JumpToEnd() 
  {
    slice_counter= im_pointer.size()-1;
    vmapper->SetInput((im_pointer)[slice_counter]);
    cout << "jump to end\n";
  }
  int slice_counter;
  std::vector<vtkSmartPointer<vtkImageData> > im_pointer;
  vtkSmartPointer<vtkVolumeMapper> vmapper;
};
