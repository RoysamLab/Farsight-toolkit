/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include "vtkDataSetMapper.h"


class vtkSlider2DCallbackBrightness : public vtkCommand
{
public:
  static vtkSlider2DCallbackBrightness *New() 
    { return new vtkSlider2DCallbackBrightness; }
  virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
      vtkSliderWidget *sliderWidget = 
        reinterpret_cast<vtkSliderWidget*>(caller);
      double value = static_cast<vtkSliderRepresentation*>(sliderWidget->GetRepresentation())->GetValue();

      vtkSmartPointer<vtkPiecewiseFunction> opacityTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();
    opacityTransferFunction->AddPoint(2,0.0);
    opacityTransferFunction->AddPoint(2+256*(1-value),0.2);
    //opacityTransferFunction->AddPoint(255,value);
     
      this->volume->GetProperty()->SetScalarOpacity(opacityTransferFunction);
    // this->imActor->GetProperty()->SetScalarOpacity(opacityTransferFunction);
    // this->imActor1->GetProperty()->SetScalarOpacity(opacityTransferFunction);		
    }
  vtkSlider2DCallbackBrightness() {

  }

  vtkSmartPointer<vtkVolume> volume;
 vtkSmartPointer<vtkVolume> imActor;	
 vtkSmartPointer<vtkVolume> imActor1;
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
    vtkColorTransferFunction *colorTransferFunction = vtkColorTransferFunction::New();
    colorTransferFunction->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
    colorTransferFunction->AddRGBPoint(255*(1-value),0.5,0.5,0.5);
     this->volume->GetProperty()->SetColor(colorTransferFunction);
     //this->imActor->GetProperty()->SetColor(colorTransferFunction);
     //this->imActor1->GetProperty()->SetColor(colorTransferFunction);

    }
  vtkSlider2DCallbackContrast() {

  }

  vtkSmartPointer<vtkVolume> volume;
  vtkSmartPointer<vtkVolume> imActor;	
 vtkSmartPointer<vtkVolume> imActor1;
};


class vtkSlider2DCallbackSeedSize : public vtkCommand
{
public:
  static vtkSlider2DCallbackSeedSize *New() 
    { return new vtkSlider2DCallbackSeedSize; }
  virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
      vtkSliderWidget *sliderWidget = 
        reinterpret_cast<vtkSliderWidget*>(caller);
      this->Origglyph->SetScaleFactor(static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue());
      this->Addglyph->SetScaleFactor(static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue());
      this->delOrigglyph->SetScaleFactor(static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue());
	  this->delAddedglyph->SetScaleFactor(static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue());
	  this->CPglyph->SetScaleFactor(static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue());
	  this->FPglyph->SetScaleFactor(static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue());
	  this->GTglyph->SetScaleFactor(static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue());
	  
	  //this->handleRep->SetHandleSize(static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue());
    }
  vtkSlider2DCallbackSeedSize():Origglyph(0),Addglyph(0),delOrigglyph(0), delAddedglyph(0), CPglyph(0), FPglyph(0), GTglyph(0) {}
  vtkGlyph3D *Origglyph;
  vtkGlyph3D *Addglyph;
  vtkGlyph3D *delOrigglyph; 
  vtkGlyph3D *delAddedglyph; 
  vtkGlyph3D *CPglyph; 
  vtkGlyph3D *FPglyph; 
  vtkGlyph3D *GTglyph; 
  
  
  
  //vtkSphereHandleRepresentation *handleRep;
};

/*
class vtkSlider2DCallbackZSlice : public vtkCommand
{
public:
  static vtkSlider2DCallbackZSlice *New() 
    { return new vtkSlider2DCallbackZSlice; }
  virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
      vtkSliderWidget *sliderWidget = 
        reinterpret_cast<vtkSliderWidget*>(caller);
            
      float p = floor((static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue()));
      this->Reslice->SetResliceAxesOrigin(0,0,p);
      this->im_mapper->SetInput(this->Reslice->GetOutput());
      this->imActor->SetMapper(im_mapper);
      	
    }
  
  //vtkSlider2DCallbackZSlice():Reslice(0) {}
  vtkSmartPointer<vtkImageReslice> Reslice;
  vtkSmartPointer<vtkDataSetMapper> im_mapper;
  vtkSmartPointer<vtkActor> imActor; 
};
*/




// This callback is responsible for setting the seed label.
class vtkSeedCallback : public vtkCommand
{
public:
  static vtkSeedCallback *New() 
    { return new vtkSeedCallback; }
  virtual void Execute(vtkObject*, unsigned long, void*)
    {

        this->x = this->handle->GetWorldPosition();
    	this->bounds = this->vol->GetBounds();
       	//int x = 0;
     	//int y = 0;


   //double* p1; 
   //double* p2; 
   double placePoint2[2] = {this->x[0]-(this->bounds[1]/2.0),this->x[1]-(this->bounds[3]/2.0)};
   double placePoint1[2] = {(this->bounds[3]/2.0)-this->x[1],this->x[2]-(this->bounds[5]/2.0)};
    
    if((x[0]>=bounds[0] && x[0]<=bounds[1]) && (x[1]>=bounds[2] && x[1]<=bounds[3]) && (x[2]>=bounds[4] && x[2]<=bounds[5]))
       {
        this->reslice->SetResliceAxesOrigin(x[0],(bounds[3]-bounds[2])/2.0,(bounds[5]-bounds[4])/2.0);
        cout<<x[0]<<endl;
		cout<<(bounds[3]-bounds[2])/2.0<<endl;
		cout<<(bounds[5]-bounds[4])/2.0<<endl;
		this->im_mapper = vtkDataSetMapper::New();
        this->im_mapper->SetInput(this->reslice->GetOutput());
        this->reslice1->SetResliceAxesOrigin((bounds[1]-bounds[0])/2.0,(bounds[3]-bounds[2])/2.0,x[2]); 
        this->im_mapper1 = vtkDataSetMapper::New();
        this->im_mapper1->SetInput(this->reslice1->GetOutput());
		this->imActor->SetPosition(0.0,0.0,0.0);
		this->imActor1->SetPosition(0.0,0.0,0.0);
        this->handle1->SetWorldPosition(placePoint1);	
		this->handle2->SetWorldPosition(placePoint2);
 		this->QVTK1->GetRenderWindow()->Render();
        this->QVTK2->GetRenderWindow()->Render();
       }
	
}
  vtkSeedCallback() : handle(0) {} 
  double* x;
  vtkPointHandleRepresentation3D *handle;
  QVTKWidget *QVTK1;		
  QVTKWidget *QVTK2;
  vtkSmartPointer<vtkImageReslice> reslice; 	
  vtkSmartPointer<vtkImageReslice> reslice1;
  vtkSmartPointer<vtkDataSetMapper> im_mapper;
  vtkSmartPointer<vtkDataSetMapper> im_mapper1;
  vtkSmartPointer<vtkActor> imActor;
  vtkSmartPointer<vtkActor> imActor1;
  vtkSmartPointer<vtkRenderer> Renderer1;
  vtkSmartPointer<vtkRenderer> Renderer2;
  vtkSmartPointer<vtkVolume> vol;
vtkPointHandleRepresentation2D *handle1;
  vtkPointHandleRepresentation2D *handle2;
  double* bounds;
};



// This callback is responsible for setting the seed label.
class vtkSeedCallback1 : public vtkCommand
{
public:
  static vtkSeedCallback1 *New() 
    { return new vtkSeedCallback1; }
  virtual void Execute(vtkObject*, unsigned long, void*)
    {
       
	this->handle1pos = this->handle1->GetWorldPosition();
	this->handle2pos = this->handle2->GetWorldPosition();
	this->bounds = this->vol->GetBounds();
	this->bounds[0] = this->handle2pos[0];
	this->bounds[1] = this->handle1pos[0]*-1.0;
    this->handle2->SetWorldPosition(this->bounds);
	this->bounds1 = this->vol->GetBounds();
    this->placePoint[0] = (this->bounds1[1]/2.0)+(this->handle2pos[0]);
    this->placePoint[1] =(bounds1[3]/2.0)+this->handle2pos[1];
    this->placePoint[2] = (bounds1[5]/2.0)+this->handle1pos[1];
    this->handle->SetWorldPosition(this->placePoint);
	
	this->reslice1->SetResliceAxesOrigin((bounds1[1]-bounds1[0])/2.0,(bounds1[3]-bounds1[2])/2.0,this->placePoint[2]); 
    this->im_mapper1 = vtkDataSetMapper::New();
    this->im_mapper1->SetInput(this->reslice1->GetOutput());
	//this->imActor1->SetPosition(0.0,0.0,0.0);
	
	this->QVTK2->GetRenderWindow()->Render();
	this->QVTK1->GetRenderWindow()->Render();
	this->QVTK->GetRenderWindow()->Render();
//	cout<<"I am moved"<<endl;	
  }

  vtkSeedCallback1() : handle1(0),handle2(0){}
  double* handle1pos;
  double* bounds;
  double* bounds1;
  double* handle2pos;
  vtkPointHandleRepresentation2D *handle1;
  vtkPointHandleRepresentation2D *handle2;
  vtkPointHandleRepresentation3D *handle;	
  QVTKWidget *QVTK2;
  vtkSmartPointer<vtkVolume> vol;
  double placePoint[3];
  QVTKWidget *QVTK;	
  vtkSmartPointer<vtkDataSetMapper> im_mapper;
  vtkSmartPointer<vtkDataSetMapper> im_mapper1;
  vtkSmartPointer<vtkActor> imActor;
  vtkSmartPointer<vtkActor> imActor1;
  vtkSmartPointer<vtkImageReslice> reslice1;
  QVTKWidget *QVTK1;

  };




// This callback is responsible for setting the seed label.
class vtkSeedCallback2 : public vtkCommand
{
public:
  static vtkSeedCallback2 *New() 
    { return new vtkSeedCallback2; }
  virtual void Execute(vtkObject*, unsigned long, void*)
    {
     
	this->handle2pos = this->handle2->GetWorldPosition();
	this->handle1pos = this->handle1->GetWorldPosition();
	this->bounds = this->vol->GetBounds();
	this->bounds[0] = this->handle2pos[1]*-1.0;
	this->bounds[1] = this->handle1pos[1];
    this->handle1->SetWorldPosition(bounds);
    this->bounds1 = this->vol->GetBounds();
    this->placePoint[0] = (this->bounds1[1]/2.0)+(this->handle2pos[0]);
    this->placePoint[1] =(bounds1[3]/2.0)+this->handle2pos[1];
    this->placePoint[2] = (bounds1[5]/2.0)+this->handle1pos[1];
    this->handle->SetWorldPosition(this->placePoint);
    
	this->reslice->SetResliceAxesOrigin(this->placePoint[0],(bounds1[3]-bounds1[2])/2.0,(bounds1[5]-bounds1[4])/2.0);
	cout<<this->placePoint[0]<<endl;
	cout<<(bounds1[3]-bounds1[2])/2.0<<endl;
	cout<<(bounds1[5]-bounds1[4])/2.0<<endl;
	
	this->im_mapper = vtkDataSetMapper::New();
    this->im_mapper->SetInput(this->reslice->GetOutput());
	this->imActor->SetPosition(0.0,0.0,0.0);
	
	this->QVTK2->GetRenderWindow()->Render();
	this->QVTK1->GetRenderWindow()->Render();
	this->QVTK->GetRenderWindow()->Render();
	
}
  vtkSeedCallback2() : handle1(0),handle2(0){}
  double* handle1pos;
  double* bounds;
  double* handle2pos;
  vtkPointHandleRepresentation2D *handle1;
  vtkPointHandleRepresentation2D *handle2;	
  vtkPointHandleRepresentation3D *handle;	
  QVTKWidget *QVTK2;
  QVTKWidget *QVTK1;
  QVTKWidget *QVTK;
  vtkSmartPointer<vtkVolume> vol;
  double placePoint[3];						
  double* bounds1;	

  vtkSmartPointer<vtkDataSetMapper> im_mapper;
  vtkSmartPointer<vtkActor> imActor;
  vtkSmartPointer<vtkImageReslice> reslice;
  	};


