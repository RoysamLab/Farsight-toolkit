/*=========================================================================

  Program:   BioImageXD
  Module:    $RCSfile: vtkLSMReader.h,v $
  Language:  C++
  Date:      $Date: 2003/11/04 21:26:04 $
  Version:   $Revision: 1.28 $


 Copyright (C) 2005  BioImageXD Project
 See CREDITS.txt for details

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

=========================================================================*/

#define TIF_NEWSUBFILETYPE 254
#define TIF_IMAGEWIDTH 256
#define TIF_IMAGELENGTH 257
#define TIF_BITSPERSAMPLE 258
#define TIF_COMPRESSION 259
#define TIF_PHOTOMETRICINTERPRETATION 262
#define TIF_STRIPOFFSETS 273
#define TIF_SAMPLESPERPIXEL 277
#define TIF_STRIPBYTECOUNTS 279
#define TIF_PLANARCONFIGURATION 284
#define TIF_PREDICTOR 317
#define TIF_COLORMAP 320
#define TIF_CZ_LSMINFO 34412

// .NAME vtkLSMReader - read LSM files
// .SECTION Description
// vtkLSMReader is a source object that reads LSM files.
// It should be able to read most any LSM file
//
// .SECTION Thanks
// This class was developed as a part of the BioImageXD Project.
// The BioImageXD project includes the following people:
//
// Dan White <dan@chalkie.org.uk>
// Kalle Pahajoki <kalpaha@st.jyu.fi>
// Pasi Kankaanp‰‰ <ppkank@bytl.jyu.fi>
// 


#ifndef __vtkLSMReader_h
#define __vtkLSMReader_h

//#include "vtkImageSource.h"
#include "vtkImageAlgorithm.h"
#include "vtkIntArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkUnsignedCharArray.h"

#define TIFF_BYTE 1
#define TIFF_ASCII 2
#define TIFF_SHORT 3
#define TIFF_LONG 4
#define TIFF_RATIONAL 5

#define LSM_MAGIC_NUMBER 42

#define LSM_COMPRESSED 5

#define VTK_FILE_BYTE_ORDER_BIG_ENDIAN 0
#define VTK_FILE_BYTE_ORDER_LITTLE_ENDIAN 1

#include "vtkBXDProcessingWin32Header.h"                                            

class VTK_BXD_PROCESSING_EXPORT vtkLSMReader : public vtkImageAlgorithm
{
public:
 
  static vtkLSMReader *New();
  vtkTypeMacro(vtkLSMReader,vtkImageAlgorithm);
  virtual void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Get the file extensions for this format.
  // Returns a string with a space separated list of extensions in 
  // the format .extension
  const char* GetFileExtensions()
    {
    return ".lsm .LSM";
    }

  int GetHeaderIdentifier();
  int IsValidLSMFile();
  int IsCompressed();
  int GetNumberOfTimePoints();
  int GetNumberOfChannels();
  int OpenFile();

  int GetChannelColorComponent(int,int);
  char* GetChannelName(int);
  void SetFileName(const char *);
  //void ExecuteInformation();
  int RequestInformation (
  vtkInformation       * vtkNotUsed( request ),
  vtkInformationVector** vtkNotUsed( inputVector ),
  vtkInformationVector * outputVector);    
  void SetUpdateTimePoint(int);
  void SetUpdateChannel(int);

  void SetDataByteOrderToBigEndian();
  void SetDataByteOrderToLittleEndian();
  void SetDataByteOrder(int);
  int GetDataByteOrder();
  const char *GetDataByteOrderAsString();

  // Description:
  // Set/Get the byte swapping to explicitly swap the bytes of a file.
  vtkSetMacro(SwapBytes,int);
  virtual int GetSwapBytes() {return this->SwapBytes;}
  vtkBooleanMacro(SwapBytes,int);

  vtkGetStringMacro(FileName);
  vtkGetVector3Macro(VoxelSizes,double);
  vtkGetVectorMacro(Dimensions,int,5);
  vtkGetVectorMacro(NumberOfIntensityValues,int,4);
  vtkGetVectorMacro(DataSpacing,double,3);
  vtkGetMacro(Identifier,unsigned short);
  vtkGetMacro(NewSubFileType,unsigned int);
  vtkGetMacro(Compression,unsigned int);
  vtkGetMacro(SamplesPerPixel,unsigned int);
  vtkGetMacro(ScanType,unsigned short);
  vtkGetMacro(DataType,int);
  vtkGetObjectMacro(TimeStampInformation,vtkDoubleArray);
  vtkGetObjectMacro(ChannelColors,vtkIntArray);

  vtkImageData* GetTimePointOutput(int,int);

protected:

  vtkLSMReader();
  ~vtkLSMReader();

  int TIFF_BYTES(unsigned short);
  int BYTES_BY_DATA_TYPE(int);
  void ClearFileName();
  void Clean();
  unsigned long ReadImageDirectory(ifstream *,unsigned long);
  int AllocateChannelNames(int);
  int SetChannelName(const char *,int);
  int ClearChannelNames();
  int FindChannelNameStart(const char *, int);
  int ReadChannelName(const char *, int, char *);
  int ReadChannelColorsAndNames(ifstream *,unsigned long);
  int ReadTimeStampInformation(ifstream *,unsigned long);
  int ReadLSMSpecificInfo(ifstream *,unsigned long);
  int AnalyzeTag(ifstream *,unsigned long);
  int NeedToReadHeaderInformation();
  void NeedToReadHeaderInformationOn();
  void NeedToReadHeaderInformationOff();
  unsigned long SeekFile(int);
  unsigned long GetOffsetToImage(int, int);
  ifstream *GetFile();

  int RequestUpdateExtent (
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector);
  
int RequestData(

  vtkInformation *vtkNotUsed(request),

  vtkInformationVector **vtkNotUsed(inputVector),

  vtkInformationVector *outputVector);



  
  //void ExecuteData(vtkDataObject *out);
  void CalculateExtentAndSpacing(int extent[6],double spacing[3]);
  //  void DecodeHorizontalDifferencing(unsigned char *,int);
  //  void DecodeLZWCompression(unsigned  char *,int);
  //  unsigned char *decode_string(unsigned char *buffer,unsigned int code);
  //  unsigned int input_code(char *input);

  int SwapBytes;

  int IntUpdateExtent[6];
  unsigned long OffsetToLastAccessedImage;
  int NumberOfLastAccessedImage;
  int FileNameChanged;
  ifstream *File;
  char *FileName;
  double VoxelSizes[3];
  int Dimensions[5];// x,y,z,time,channels
  int NumberOfIntensityValues[4];
  unsigned short Identifier;
  unsigned int NewSubFileType;
  vtkUnsignedShortArray *BitsPerSample;
  unsigned int Compression;
  vtkUnsignedIntArray *StripOffset;
  unsigned int SamplesPerPixel;
  vtkUnsignedIntArray *StripByteCount;
  unsigned int LSMSpecificInfoOffset;
  unsigned short PhotometricInterpretation;
  unsigned long ColorMapOffset;
  unsigned short PlanarConfiguration;
  unsigned short Predictor;
  unsigned short ScanType;
  int DataScalarType;
  
  double DataSpacing[3];
  int DataExtent[6];
  int NumberOfScalarComponents;
  int DataType;
  unsigned long ChannelInfoOffset;
  vtkIntArray *ChannelColors;
  char **ChannelNames;
  vtkDoubleArray *TimeStampInformation;

  unsigned char CharPointerToUnsignedChar(char *);
  int CharPointerToInt(char *);
  unsigned int CharPointerToUnsignedInt(char *);
  short CharPointerToShort(char *);
  unsigned short CharPointerToUnsignedShort(char *);
  double CharPointerToDouble(char *);

  int ReadInt(ifstream *,unsigned long *);
  unsigned int ReadUnsignedInt(ifstream *,unsigned long *);
  short ReadShort(ifstream *,unsigned long *);
  unsigned short ReadUnsignedShort(ifstream *,unsigned long *);
  double ReadDouble(ifstream *,unsigned long *);
  int ReadFile(ifstream *,unsigned long *,int,char *,bool swap=0);
  int ReadData(ifstream *,unsigned long *,int,char *);


private:
  vtkLSMReader(const vtkLSMReader&);  // Not implemented.
  void operator=(const vtkLSMReader&);  // Not implemented.
};
#endif
