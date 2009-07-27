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

////////////////////////////////////////////////////////////////////////////
// File: template.h
//	3-24-98
//
// This file contains the declaration of a a rectangular template class 
// as well as a vector class. 
// The value of a template (vector) cell represents absolute shifts from 
// the base of the template (vector) to the given cell. Hence, such shifts
// depends on the image being operated on. All templates are composed of 
// five rows, with the middle row containing no shift values, hence need 
// not be created. The width of a vector, on the other hand, is one.
// The length and direction of a template (vector) is given at creation time.

#ifndef template_h
#define template_h


//////////////////////////////////////////////////////////////////////////
//////////////// base template class  ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

// base class template
// contains common fucntionality such as that needed to position, shift,
// and apply the template. Such operations are the same regardless whether
// the template is a left or right one. Notice that this class is declared 
// abstract
// Note: Remember to always set the image pointer to the image where
// the template is to be applied.
class CTemplate
{
public:

	// CTOR. needs the template's length, direction, and the image it is
	// created on.
	CTemplate(int row, int col, int slice, float Hdir, float Vdir, int length);
	// DTOR
	virtual ~CTemplate();

	// construct the template rows
	void ConstructTemplate();

	// set the base of the template to the given position
	void Position(int x, int y, int z = -1);
	// set the base of the template to the given point
	void Position(CPoint &);

	// draw the template in the given image
	void DisplayTemplate(C3DImage &);

	// Apply the template at the given location and return the response
	int Apply(CPoint *, int templateWidth = 0);

	// Apply the templage at the given location and return the response
	// this function is the same like the above but it takes a pointer to
	// the image data at the base of the template.
	int Apply(unsigned char *dataPtr);

	// Shift the template by the given amount in a direction perpendicular
	// to its direction
	int ShiftAndApply(int index);

	// Calculate the template responses at all points along the
	// perpendicular vector based on the given point and return
	// the maximum of such responses. The argument "pPoint" need to be passed if the location of the
	// maximum response is needed. 
	int CalculateMaxResponse(CPoint *, CPoint *atPoint = 0, CVessel *aVessel = 0, CPoint *pPoints = 0);

	// compute the template responses for different shift values and different
	// lengths
	int CalculateResponses(CPoint *, CVessel *,
								  CPoint aPoints[][MaxTemplateLength],
								  int aiResponses[][MaxTemplateLength]);

	// function that calculates just individual responses
//	int CalculateIndividualResponses(CPoint *aPoint, CPoint aPoints[][MaxTemplateLength],
//											 short aiResponses[][MaxTemplateLength]);
	int CalculateKernelResponses(const CPoint &aPoint, int prev_shift,
											CPoint aPoints[][MaxTemplateLength],
											short aiResponses[][MaxTemplateLength]);

	// all subsequent operations will be applied on this image
	void AssociateWithImage(CImage *anImage) { m_pMyImage = anImage; }

	// all subsequent operations will be applied on this image
	void AssociateWithImage(C3DImage &anImage) { m_pMy3DImage = &anImage; }


	// return the direction orthogonal to my direction (i.e. direction of shift)
	// this depends on whether the template is left or right template, hence it
	// it is defined as pure virtual and the derived classes are required to 
	// provide the proper definition
	virtual int GetInPlaneShiftDir() = 0;
	virtual int GetOrthogonalShiftDir() = 0;
	
	// return a pointer to a vector pointing in the direction that
	// will shift this template a way from a dendrite's centerline
	virtual CVector *GetShiftFromCenterVector() = 0;

	// return a pointer to a vector pointing in the direction that
	// will shift this template towards a dendrite's centerline
	virtual CVector *GetShiftToCenterVector() = 0;

	///////////////////////////////////////////////////////////////////////////
	// Compute the shifts for the given row. Before shifting, the row lies
	// in a direction parallel to the x-axis and passes through the point (0, y)
	void ComputeShifts(int *row, int y);


	////////////////
	// data members
	////////////////
	
	// the direction of the template (in degrees)
	float m_fHdir;
	float m_fVdir;

	// the length of the template
	int m_iLength;
	// the number of rows, columns, and slices in the corresponding image
	int m_iRows;
	int m_iCols;
	int m_iSlices;
	// a pointer to the address of the template center in an image
	// notice that this pointer is set by the method "Apply" and is 
	// updated by the method "Shift"
	unsigned char *m_puchBase;
	// the template vectors. These contain absolute shifts from the "base"
	// for an image of a given image
	int *m_piPlus2;
	int *m_piPlus1;
	int *m_piMinus1;
	int *m_piMinus2;

	// the responses of the template vectors. These make shifting and
	// calculating the responses computationally cheaper.
	int m_iPlus2Resp;
	int m_iPlus1Resp;
	int m_iMinus1Resp;
	int m_iMinus2Resp;

	// the best estimate of the forground and background pixel
	// values after applying the template to the image
	int m_iForegroundPixelValue;
	int m_iBackgroundPixelValue;
	// a pointer to the image where this templage will be applied
	CImage   *m_pMyImage;
	C3DImage *m_pMy3DImage;
};


////////////////////////////////////////////////////////////////
// class: CHLeftTemplate
//
// represents a horizontal left template
class CHLeftTemplate : public CTemplate
{
public:

	// CTOR. needs the template's length, direction, and the image it is
	// created on.
	CHLeftTemplate(int row, int col, int slice, float Hdir, float Vdir, int length);
	// DTOR
	virtual ~CHLeftTemplate();
	
	
	int GetInPlaneShiftDir();
	int GetOrthogonalShiftDir();

	// return a pointer to a vector pointing in the direction that
	// will shift this template a way from a dendrite's centerline
	CVector *GetShiftFromCenterVector();

	// return a pointer to a vector pointing in the direction that
	// will shift this template towards a dendrite's centerline
	CVector *GetShiftToCenterVector();
};

/////////////////////////////////////////////////////////////////
// class: CHTemplateRight
// 
// Represents a horizontal right template
class CHRightTemplate : public CTemplate
{
public:

	// CTOR. needs the template's length, direction, and the image it is
	// created on.
	CHRightTemplate(int row, int col, int slice, float Hdir, float Vdir, int length);
	// DTOR
	virtual ~CHRightTemplate();

	int GetInPlaneShiftDir();
	int GetOrthogonalShiftDir();

	// return a pointer to a vector pointing in the direction that
	// will shift this template a way from a dendrite's centerline
	CVector *GetShiftFromCenterVector();

	// return a pointer to a vector pointing in the direction that
	// will shift this template towards a dendrite's centerline
	CVector *GetShiftToCenterVector();
};

/////////////////////////////////////////////////////////////////
// class: CVTemplateLeft
// 
// Represents a vertical left template
class CVLeftTemplate : public CTemplate
{
public:

	// CTOR. needs the template's length, direction, and the image it is
	// created on.
	CVLeftTemplate(int row, int col, int slice, float Hdir, float Vdir, int length);
	// DTOR
	virtual ~CVLeftTemplate();

	int GetInPlaneShiftDir();
	int GetOrthogonalShiftDir();

	// return a pointer to a vector pointing in the direction that
	// will shift this template a way from a dendrite's centerline
	CVector *GetShiftFromCenterVector();

	// return a pointer to a vector pointing in the direction that
	// will shift this template towards a dendrite's centerline
	CVector *GetShiftToCenterVector();
};

/////////////////////////////////////////////////////////////////
// class: CVTemplateRigth
// 
// Represents a vertical right template
class CVRightTemplate : public CTemplate
{
public:

	// CTOR. needs the template's length, direction, and the image it is
	// created on.
	CVRightTemplate(int row, int col, int slice, float Hdir, float Vdir, int length);
	// DTOR
	virtual ~CVRightTemplate();

	
	int GetInPlaneShiftDir();
	int GetOrthogonalShiftDir();

	// return a pointer to a vector pointing in the direction that
	// will shift this template a way from a dendrite's centerline
	CVector *GetShiftFromCenterVector();

	// return a pointer to a vector pointing in the direction that
	// will shift this template towards a dendrite's centerline
	CVector *GetShiftToCenterVector();
};

#endif

