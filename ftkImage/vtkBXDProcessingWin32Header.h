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

/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkBXDProcessingWin32Header.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkBXDProcessingWin32Header - manage Windows system differences
// .SECTION Description
// The vtkBXDProcessingWin32Header captures some system differences between Unix
// and Windows operating systems. 

#ifndef __vtkBXDProcessingWin32Header_h
#define __vtkBXDProcessingWin32Header_h

//#include <VTKBXDConfigure.h>

//#if defined(WIN32) && !defined(VTKBXD_STATIC)
//#if defined(vtkBXDProcessing_EXPORTS)
//#define VTK_BXD_PROCESSING_EXPORT __declspec( dllexport )
//#else
//#define VTK_BXD_PROCESSING_EXPORT __declspec( dllimport )
//#endif
//#else
#define VTK_BXD_PROCESSING_EXPORT
//#endif

#endif
