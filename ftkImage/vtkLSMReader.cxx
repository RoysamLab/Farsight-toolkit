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

  Program:   BioImageXD
  Module:    $RCSfile: vtkLSMReader.cxx,v $
  Language:  C++
  Date:      $Date: 2003/08/22 14:46:02 $
  Version:   $Revision: 1.39 $


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
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif
#include "vtkLSMReader.h"
#include "vtkObjectFactory.h"
#include "vtkImageData.h"
#include "vtkSource.h"
#include "vtkPointData.h"
#include "vtkByteSwap.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include <ctime>

#define PRT_EXT(ext) ext[0],ext[1],ext[2],ext[3],ext[4],ext[5]
#define PRT_EXT2(ext) ext[0]<<","<<ext[1]<<","<<ext[2]<<","<<ext[3]<<","<<ext[4]<<","<<ext[5]

vtkStandardNewMacro(vtkLSMReader);

vtkLSMReader::vtkLSMReader()
{
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);      
  this->Clean();
}

vtkLSMReader::~vtkLSMReader()
{
  this->ClearFileName();
  this->ClearChannelNames();
  this->ChannelColors->Delete();
  this->BitsPerSample->Delete();
  this->StripOffset->Delete();
  this->StripByteCount->Delete();
}

void vtkLSMReader::ClearFileName()
{
  if (this->File)
    {
    this->File->close();
    delete this->File;
    this->File = NULL;
    }
  
  if (this->FileName)
    {
    delete [] this->FileName;
    this->FileName = NULL;
    }
}
//----------------------------------------------------------------------------
// This function sets the name of the file. 
void vtkLSMReader::SetFileName(const char *name)
{
  if ( this->FileName && name && (!strcmp(this->FileName,name)))
    {
    return;
    }
  if (!name && !this->FileName)
    {
    return;
    }
  if (this->FileName)
    {
    delete [] this->FileName;
    }
  if (name)
    {
    this->FileName = new char[strlen(name) + 1];
    strcpy(this->FileName, name);
    }
  else
    {
    this->FileName = NULL;
    }
  this->NeedToReadHeaderInformationOn();
  this->Modified();
}


void vtkLSMReader::Clean()
{
  
  this->IntUpdateExtent[0] = this->IntUpdateExtent[1] = this->IntUpdateExtent[2] = this->IntUpdateExtent[4] = 0;
  this->IntUpdateExtent[3] = this->IntUpdateExtent[5] = 0;
    
  this->DataExtent[0] = this->DataExtent[1] = this->DataExtent[2] = this->DataExtent[4] = 0;
  this->DataExtent[3] = this->DataExtent[5] = 0;    
  this->OffsetToLastAccessedImage = 0;
  this->NumberOfLastAccessedImage = 0;
  this->FileNameChanged = 0;
  this->FileName = NULL;
  this->File = NULL;
  this->VoxelSizes[0] = this->VoxelSizes[1] = this->VoxelSizes[2] = 0.0;
  this->Identifier = 0;
    
  this->DataSpacing[0] = this->DataSpacing[1] = this->DataSpacing[2] =  1.0f;
  this->Dimensions[0] = this->Dimensions[1] = this->Dimensions[2] = this->Dimensions[3] = this->Dimensions[4] = 0;
  this->NewSubFileType = 0;
  this->BitsPerSample = vtkUnsignedShortArray::New();
  this->BitsPerSample->SetNumberOfTuples(4);
  this->BitsPerSample->SetNumberOfComponents(1);  
  this->Compression = 0;
  this->StripOffset = vtkUnsignedIntArray::New();
  this->StripOffset->SetNumberOfTuples(4);
  this->StripOffset->SetNumberOfComponents(1);  
  this->SamplesPerPixel = 0;
  this->StripByteCount = vtkUnsignedIntArray::New();
  this->StripByteCount->SetNumberOfTuples(4);
  this->StripByteCount->SetNumberOfComponents(1);  
  this->Predictor = 0;
  this->PhotometricInterpretation = 0;
  this->PlanarConfiguration = 0;
  this->ColorMapOffset = 0;
  this->LSMSpecificInfoOffset = 0;
  this->NumberOfIntensityValues[0] = this->NumberOfIntensityValues[1] = 
  this->NumberOfIntensityValues[2] = this->NumberOfIntensityValues[3] = 0;
  this->ScanType = 0;
  this->DataType = 0;
  this->ChannelColors = vtkIntArray::New();
  this->ChannelNames = NULL;
  this->TimeStampInformation = vtkDoubleArray::New();
}

int vtkLSMReader::OpenFile()
{
  if (!this->FileName)
    {
    vtkErrorMacro(<<"FileName must be specified.");
    return 0;
    }

  // Close file from any previous image
  if (this->File)
    {
    this->File->close();
    delete this->File;
    this->File = NULL;
    }
  
  // Open the new file
#ifdef _WIN32
  this->File = new ifstream(this->FileName, ios::in | ios::binary);
#else
  this->File = new ifstream(this->FileName, ios::in);
#endif
  if (! this->File || this->File->fail())
    {
    vtkErrorMacro(<< "OpenFile: Could not open file " <<this->FileName);
    return 0;
    }
  return 1;
}


void vtkLSMReader::SetDataByteOrderToBigEndian()
{
#ifndef VTK_WORDS_BIGENDIAN
  this->SwapBytesOn();
#else
  this->SwapBytesOff();
#endif
}

void vtkLSMReader::SetDataByteOrderToLittleEndian()
{
#ifdef VTK_WORDS_BIGENDIAN
  this->SwapBytesOn();
#else
  this->SwapBytesOff();
#endif
}

void vtkLSMReader::SetDataByteOrder(int byteOrder)
{
  if ( byteOrder == VTK_FILE_BYTE_ORDER_BIG_ENDIAN )
    {
    this->SetDataByteOrderToBigEndian();
    }
  else
    {
    this->SetDataByteOrderToLittleEndian();
    }
}

int vtkLSMReader::GetDataByteOrder()
{
#ifdef VTK_WORDS_BIGENDIAN
  if ( this->SwapBytes )
    {
    return VTK_FILE_BYTE_ORDER_LITTLE_ENDIAN;
    }
  else
    {
    return VTK_FILE_BYTE_ORDER_BIG_ENDIAN;
    }
#else
  if ( this->SwapBytes )
    {
    return VTK_FILE_BYTE_ORDER_BIG_ENDIAN;
    }
  else
    {
    return VTK_FILE_BYTE_ORDER_LITTLE_ENDIAN;
    }
#endif
}

const char *vtkLSMReader::GetDataByteOrderAsString()
{
#ifdef VTK_WORDS_BIGENDIAN
  if ( this->SwapBytes )
    {
    return "LittleEndian";
    }
  else
    {
    return "BigEndian";
    }
#else
  if ( this->SwapBytes )
    {
    return "BigEndian";
    }
  else
    {
    return "LittleEndian";
    }
#endif
}


char* vtkLSMReader::GetChannelName(int chNum)
{
  if(!this->ChannelNames || chNum < 0 || chNum > this->GetNumberOfChannels()-1)
    {
    vtkDebugMacro(<<"GetChannelName: Illegal channel index!");
    return (char *)"";
    }
  return this->ChannelNames[chNum];
}

int vtkLSMReader::ClearChannelNames()
{
  vtkDebugMacro(<<"clearing " << this->GetNumberOfChannels()<<"channel names");
   if(!this->ChannelNames || this->GetNumberOfChannels() < 1)
    {
    return 0;
    }
  
  for(int i=0;i<this->GetNumberOfChannels();i++)
    {
    delete [] this->ChannelNames[i];
    }
  delete [] this->ChannelNames;

  return 0;
}

int vtkLSMReader::AllocateChannelNames(int chNum)
{
  this->ClearChannelNames();
  vtkDebugMacro(<<"allocating space for "<<chNum<<"channel names");
  this->ChannelNames = new char*[chNum];
  if(!this->ChannelNames)
    {
    vtkErrorMacro(<<"Could not allocate memory for channel name table!");
    return 1;
    }
  for(int i=0;i<chNum;i++)
    {
    this->ChannelNames[i] = NULL;
    }
  return 0;
}

int vtkLSMReader::SetChannelName(const char *chName, int chNum)
{
  char *name;
  int length;
  if(!chName || chNum > this->GetNumberOfChannels())
    {
    return 0;
    }
  if(!this->ChannelNames)
    {
    this->AllocateChannelNames(this->GetNumberOfChannels());
    }
  
  length = (int)strlen(chName);
  vtkDebugMacro(<<"length="<<length);    
  name = new char[length+1];
  if(!name)
    {
    vtkErrorMacro(<<"Could not allocate memory for channel name");
    return 1;
    }
  strncpy(name,chName,length);
  name[length] = 0;
  this->ChannelNames[chNum] = name;
  return 0;
}

int vtkLSMReader::FindChannelNameStart(const char *nameBuff, int length)
{
  int i;
  char ch;
  for(i=0;i<length;i++)
    {
    ch = *(nameBuff+i);
    if(ch > 32)
      {
      break;
      }
    }
  if(i >= length)
    {
    vtkWarningMacro(<<"Start of the channel name may not be found!");
    }
  return i;
}

int vtkLSMReader::ReadChannelName(const char *nameBuff, int length, char *buffer)
{
  int i;
  char component;
  for(i=0;i<length;i++)
    {
    component = *(nameBuff+i);
    *(buffer+i) = component;
    if(component == 0)
      {
      break;
      }
    }
  return i;
}

int vtkLSMReader::ReadChannelColorsAndNames(ifstream *f,unsigned long start)
{
  int colNum,nameNum,sizeOfStructure,sizeOfNames,nameLength, nameSkip;
  unsigned long colorOffset,nameOffset,pos;
  char *nameBuff,*colorBuff,*name,*tempBuff;
  unsigned char component;

  colorBuff = new char[5];

  pos = start;
  // Read size of structure

  sizeOfStructure = this->ReadInt(f,&pos);
  //vtkDebugMacro(<<"size of structure = "<<sizeOfStructure<<"\n");
  // Read number of colors
  colNum = this->ReadInt(f,&pos);
  // Read number of names
  nameNum = this->ReadInt(f,&pos);
  //vtkDebugMacro(<<"nameNum="<<nameNum);
  sizeOfNames = sizeOfStructure - ( (10*4) + (colNum*4) );
  //vtkDebugMacro(<<"sizeofNames="<<sizeOfNames<<"\n");
  nameBuff = new char[sizeOfNames+1];
  name = new char[sizeOfNames+1];

  if(colNum != this->GetNumberOfChannels())
    {
    vtkWarningMacro(<<"Number of channel colors is not same as number of channels!");
    }
  if(nameNum != this->GetNumberOfChannels())
    {
    vtkWarningMacro(<<"Number of channel names is not same as number of channels!");
    }

  // Read offset to color info
  colorOffset = this->ReadInt(f,&pos) + start;
  // Read offset to name info
  nameOffset = this->ReadInt(f,&pos) + start;

  //vtkDebugMacro(<<"colorOffset="<<colorOffset);
  //vtkDebugMacro(<<"nameOffset="<<nameOffset);
  //vtkDebugMacro(<<"number of colors"<< colNum);
  this->ChannelColors->Reset();
  this->ChannelColors->SetNumberOfValues(3*(colNum+1));
  this->ChannelColors->SetNumberOfComponents(3);

    
  // Read the colors
  for(int j=0;j<colNum;j++)
    {
    this->ReadFile(f,&colorOffset,4,colorBuff,1);
    
    for(int i=0;i<3;i++)
      {
        component = this->CharPointerToUnsignedChar(colorBuff+i);        
        //printf("Setting value %d to %d\n",i+((colNum+1)*j),component);
        //this->ChannelColors->SetValue(i+((colNum+1)*j),component);
        this->ChannelColors->SetValue(i+(3*j),component);
      }
    }
  //vtkDebugMacro(<<"read colors, now reading"<<sizeOfNames<<"from "<<nameOffset);
  
/*  for(int j=0;j<colNum;j++) {
      this->ChannelColors->SetValue(0+((colNum+1)*j),255);
      this->ChannelColors->SetValue(1+((colNum+1)*j),255);
      this->ChannelColors->SetValue(2+((colNum+1)*j),255);
  }*/
  
  this->ReadFile(f,&nameOffset,sizeOfNames,nameBuff,1);
  //printf("namebuf=\n");
  //for(int k=0;k<sizeOfNames;k++) {
  //   printf("%c",nameBuff[k]);
  //}
  //printf("\n");
  nameLength = nameSkip = 0;
  tempBuff = nameBuff;
  for(int i=0;i<nameNum;i++)
    {
    nameSkip = this->FindChannelNameStart(tempBuff,sizeOfNames-nameSkip);
    nameLength = this->ReadChannelName(tempBuff+nameSkip,sizeOfNames-nameSkip,name);
    
    tempBuff += nameSkip + nameLength;
    //vtkDebugMacro(<<"Setting channel "<<i<<"name");
    this->SetChannelName(name,i);
    }
  
  delete [] nameBuff;
  delete [] name;
  delete [] colorBuff;
  return 0;
}

int vtkLSMReader::ReadTimeStampInformation(ifstream *f,unsigned long offset)
{
  int numOffStamps = 0;
  if( offset == 0 ) // position is 0 for non-timeseries files!
    {
    return 0;
    }
  offset += 4;
  numOffStamps = this->ReadInt(f,&offset);
  if(numOffStamps != this->GetNumberOfTimePoints())
    {
//    vtkWarningMacro(<<"Number of time stamps does not correspond to the number off time points!");
    }
  this->TimeStampInformation->Reset();
  this->TimeStampInformation->SetNumberOfTuples(numOffStamps);
  this->TimeStampInformation->SetNumberOfComponents(1);
  for(int i=0;i<numOffStamps;i++)
    {
    this->TimeStampInformation->SetValue(i,this->ReadDouble(f,&offset));
    }
  return 0;
}

/* Read the TIF_CZ_LSMINFO entry described in Table 17 of the LSM file format specification
 *
 *
 */
int vtkLSMReader::ReadLSMSpecificInfo(ifstream *f,unsigned long pos)
{
  unsigned long offset;
  vtkDebugMacro("ReadLSMSpecificInfo(stream,"<<pos<<")\n");
  pos += 2 * 4; // skip over the start of the LSMInfo
                // first 4 byte entry if magic number
                // second is number of bytes in this structure

  // Then we read X
  this->NumberOfIntensityValues[0] = this->ReadInt(f,&pos); 
  // vtkByteSwap::Swap4LE((int*)&this->NumberOfIntensityValues[0]);
  this->Dimensions[0] = this->NumberOfIntensityValues[0];
  // Y
  this->NumberOfIntensityValues[1] = this->ReadInt(f,&pos); 
  this->Dimensions[1] = this->NumberOfIntensityValues[1];
  // and Z dimension
  this->NumberOfIntensityValues[2] = this->ReadInt(f,&pos); 
  this->Dimensions[2] = this->NumberOfIntensityValues[2];
  //printf("Got dimensions from LSM file=%d,%d,%d\n",Dimensions[0],Dimensions[1],Dimensions[2]);
  vtkDebugMacro(<<"Dimensions =" << Dimensions[0]<<","<<Dimensions[1]<<","<<Dimensions[2]<<"\n");
  // Read number of channels
  this->Dimensions[4] = this->ReadInt(f,&pos); 

  // Read number of timepoints
  this->NumberOfIntensityValues[3] = this->ReadInt(f,&pos);
  this->Dimensions[3] = this->NumberOfIntensityValues[3];

  // Read datatype, 1 for 8-bit unsigned int
  //                2 for 12-bit unsigned int
  //                5 for 32-bit float (timeseries mean of ROIs)
  //                0 if the channels have different types
  //                In that case, u32OffsetChannelDataTypes
  //                has further info
  this->DataType = this->ReadInt(f,&pos);

  // Skip the width and height of thumbnails
  pos += 2 * 4;

  // Read voxel sizes
  this->VoxelSizes[0] = this->ReadDouble(f,&pos);
  this->VoxelSizes[1] = this->ReadDouble(f,&pos);
  this->VoxelSizes[2] = this->ReadDouble(f,&pos);
  vtkDebugMacro("Voxel size="<<VoxelSizes[0]<<","<<VoxelSizes[1]<<","<<VoxelSizes[2]<<"\n");

  // Skip over OriginX,OriginY,OriginZ which are not used
  pos += 3*8;
  
  // Read scan type which is 
  // 0 for normal x-y-z scan
  // 1 for z-scan (x-z plane)
  // 2 for line scan
  // 3 for time series x-y
  // 4 for time series x-z
  // 5 time series mean of ROIs
  // 6 time series x y z
  // 7 spline scan
  // 8 spline plane x-z
  // 9 time series spline plane
  // 10 point mode
  this->ScanType = this->ReadShort(f,&pos);

  // skip over SpectralScan flag
  // if 0, no spectral scan
  // if 1, image has been acquired with spectral scan mode with a "meta" detector
  // skip over DataType, Offset to vector overlay, Offset to input LUT
  pos += 1*2 + 4*4;// + 1*8 + 3*4;
  
  // Read OffsetChannelColors, which is an offset to channel colors and names
  this->ChannelInfoOffset = this->ReadUnsignedInt(f,&pos);
  vtkDebugMacro(<<"Channel info offset (from addr"<<pos<<")="<<this->ChannelInfoOffset<<"\n");
  this->ReadChannelColorsAndNames(f,this->ChannelInfoOffset);

  // Skip time interval in seconds (8 bytes)
  // Skip offset to channel datatypes, if they differ from above
  // Skip scan information (device settings)
  // SKip Zeiss Vision KS-3D speific data
  pos += 1*8 + 3*4;

  // Read timestamp information
  offset = this->ReadUnsignedInt(f,&pos);
  this->ReadTimeStampInformation(f,offset);

  
  return 1;
}

int vtkLSMReader::AnalyzeTag(ifstream *f,unsigned long startPos)
{
  unsigned short type,length,tag;
  unsigned long readSize;
  int value, dataSize,i;
  char tempValue[4],tempValue2[4];
  char *actualValue = NULL;
    //vtkDebugMacro(<<"Analyze tag start pos="<<startPos<<"\n");
  tag = this->ReadUnsignedShort(f,&startPos);
  type = this->ReadUnsignedShort(f,&startPos);
  length = this->ReadUnsignedInt(f,&startPos);
   
  this->ReadFile(f,&startPos,4,tempValue);

  for(i=0;i<4;i++)tempValue2[i]=tempValue[i];
#ifdef VTK_WORDS_BIGENDIAN
  //vtkDebugMacro(<<"Swapping byte order...");
  vtkByteSwap::Swap4LE((unsigned int*)tempValue2);
#endif
  value = this->CharPointerToUnsignedInt(tempValue2);
    //vtkDebugMacro(<<"value="<<value<<"\n");
  // if there is more than 4 bytes in value, 
  // value is an offset to the actual data
  dataSize = this->TIFF_BYTES(type);
  readSize = dataSize*length;
  if(readSize > 4 && tag != TIF_CZ_LSMINFO)
    {
    actualValue = new char[readSize];
    startPos = value;
   if(tag == TIF_STRIPOFFSETS ||tag == TIF_STRIPBYTECOUNTS) 
      // vtkDebugMacro(<<"Reading actual value from "<<startPos<<"to " << startPos+readSize);
    if( !this->ReadFile(f,&startPos,readSize,actualValue) ) return 0;
    }
  else
    {
      actualValue = new char[4];
      //strcpy(actualValue,tempValue);
      // stupid..
      for(int o=0;o<4;o++)actualValue[o] = tempValue[o];
    }
    //vtkDebugMacro(<<"Analuzing tag"<<tag);
  switch(tag)
    {
    case TIF_NEWSUBFILETYPE: 
      //vtkDebugMacro(<<"New subfile type="<<value);
    this->NewSubFileType = value;
       
      
      /*
      vtkByteSwap::Swap4LE((unsigned int*)actualValue);
      {
    unsigned int subfileType = this->CharPointerToUnsignedInt(actualValue);
    vtkDebugMacro(<<"Subfiletype="<<subfileType<<"value="<<value);
    this->NewSubFileType = subfileType;
    }*/
      break;
    
    case TIF_IMAGEWIDTH: 
#ifdef VTK_WORDS_BIGENDIAN
      vtkByteSwap::Swap4LE((unsigned int*)actualValue);
      //vtkDebugMacro(<<"Image width="<<value);
#endif
      //this->Dimensions[0] = this->CharPointerToUnsignedInt(actualValue);
      //this->Dimensions[0] = value;
      break;
    
    case TIF_IMAGELENGTH:
#ifdef VTK_WORDS_BIGENDIAN
      vtkByteSwap::Swap4LE((unsigned int*)actualValue);
      //this->Dimensions[1] = this->CharPointerToUnsignedInt(actualValue);
      vtkDebugMacro(<<"Image length="<<value);
#endif
      //this->Dimensions[1] = value;
      break;
    
    case TIF_BITSPERSAMPLE:
#ifdef VTK_WORDS_BIGENDIAN
      vtkByteSwap::Swap2LE((unsigned short*)actualValue);
#endif
      this->BitsPerSample->SetNumberOfValues(length);
      unsigned short bitsPerSample;
      for(i=0;i<length;i++)
    {
      bitsPerSample = this->CharPointerToUnsignedShort(actualValue + (this->TIFF_BYTES(TIFF_SHORT)*i));
      //      vtkDebugMacro(<<"Bits per sample " << i<<"="<<bitsPerSample<<"\n");
      this->BitsPerSample->SetValue(i,bitsPerSample);
    }
    break;
    
    case TIF_COMPRESSION:
#ifdef VTK_WORDS_BIGENDIAN
      vtkByteSwap::Swap2LE((unsigned short*)actualValue);
#endif
      this->Compression = this->CharPointerToUnsignedShort(actualValue);
      break;
    
    case TIF_PHOTOMETRICINTERPRETATION:
#ifdef VTK_WORDS_BIGENDIAN
      vtkByteSwap::Swap2LE((unsigned short*)actualValue);
#endif
      this->PhotometricInterpretation = this->CharPointerToUnsignedShort(actualValue);
      break;
    
    case TIF_STRIPOFFSETS:
      //      vtkDebugMacro(<<"Number of values="<<length);
      this->StripOffset->SetNumberOfValues(length);
#ifdef VTK_WORDS_BIGENDIAN
      vtkByteSwap::Swap4LERange((unsigned int*)actualValue,length);
#endif
    if(length>1) {
          for(i=0;i<length;i++)
        {
          //          unsigned int stripOffset = this->CharPointerToUnsignedInt(actualValue + (this->TIFF_BYTES(TIFF_LONG)*i));
          unsigned int* offsets = (unsigned int*)actualValue;
          //    this->StripOffset->SetValue(i,this->CharPointerToUnsignedInt(actualValue + (this->TIFF_BYTES(TIFF_LONG)*i)));
          unsigned int stripOffset=offsets[i];
          //vtkDebugMacro(<<"Strip offset to "<<i<<"="<<stripOffset);   
          this->StripOffset->SetValue(i,stripOffset);
        }
    } else {
        //vtkDebugMacro(<<"Strip offset to only channel="<<value);
        this->StripOffset->SetValue(0,value);
    }
      break;
    
    case TIF_SAMPLESPERPIXEL:
#ifdef VTK_WORDS_BIGENDIAN
      vtkByteSwap::Swap4LE((unsigned int*)actualValue);
#endif
      this->SamplesPerPixel = this->CharPointerToUnsignedInt(actualValue);
      //      vtkDebugMacro(<<"Samples per pixel="<<SamplesPerPixel<<"\n");
      break;
    
    case TIF_STRIPBYTECOUNTS:
#ifdef VTK_WORDS_BIGENDIAN
      vtkByteSwap::Swap4LERange((unsigned int*)actualValue,length);
#endif      
      this->StripByteCount->SetNumberOfValues(length);

    if(length>1) {
          for(i=0;i<length;i++)
        {
          //unsigned int* counts = (unsigned int*)actualValue;
          unsigned int bytecount = this->CharPointerToUnsignedInt(actualValue + (this->TIFF_BYTES(TIFF_LONG)*i));
          
            this->StripByteCount->SetValue(i,bytecount);
        //  this->StripByteCount->SetValue(i,counts[i]);
            //vtkDebugMacro(<<"Strip byte count of " << i <<"="<<counts[i] <<"("<<bytecount<<")");
        }
    } else {
         //vtkDebugMacro(<<"Bytecount of only strip="<<value);
         this->StripByteCount->SetValue(0,value);
    }
      break;
    case TIF_PLANARCONFIGURATION:
#ifdef VTK_WORDS_BIGENDIAN
      vtkByteSwap::Swap2LE((unsigned short*)actualValue);
#endif
      this->PlanarConfiguration = this->CharPointerToUnsignedShort(actualValue);
      break;
    case TIF_PREDICTOR:
#ifdef VTK_WORDS_BIGENDIAN
      vtkByteSwap::Swap2LE((unsigned short*)actualValue);
#endif
      this->Predictor = this->CharPointerToUnsignedShort(actualValue);
      break;
    case TIF_COLORMAP:
#ifdef VTK_WORDS_BIGENDIAN
      vtkByteSwap::Swap4LE((unsigned int*)actualValue);
#endif
      //this->ColorMapOffset = this->CharPointerToUnsignedInt(actualValue);
      break;
    case TIF_CZ_LSMINFO:

      this->LSMSpecificInfoOffset = value;
      break;
    }

  if(actualValue)    
    {
    vtkDebugMacro(<<"Deleting actual value...");
    delete [] actualValue;
    }
    vtkDebugMacro(<<"done\n");
  return 0;
}


/*------------------------------------------------------------------------------------------*/

int vtkLSMReader::GetHeaderIdentifier()
{  
  return this->Identifier;
}

int vtkLSMReader::IsValidLSMFile()
{
  if(this->GetHeaderIdentifier() == LSM_MAGIC_NUMBER) return 1;
  return 0;
}

int vtkLSMReader::IsCompressed()
{
  return (this->Compression == LSM_COMPRESSED ? 1 : 0);
}

int vtkLSMReader::GetNumberOfTimePoints()
{
  return this->Dimensions[3];
}

int vtkLSMReader::GetNumberOfChannels()
{
  return this->Dimensions[4];
}

unsigned long vtkLSMReader::GetOffsetToImage(int slice, int timepoint)
{
  return this->SeekFile(slice+(timepoint*this->Dimensions[2]));
}

unsigned long vtkLSMReader::SeekFile(int image)
{
  unsigned long offset = 4, finalOffset;
  int i = 0;
  //int readSize = 4;
  //unsigned short numberOfTags = 0;  
  int imageCount = image+1;

  if(this->OffsetToLastAccessedImage && (this->NumberOfLastAccessedImage < image))
    {
    offset = this->OffsetToLastAccessedImage;
    imageCount = image - this->NumberOfLastAccessedImage;
    }
  else
    {
    offset = (unsigned long)this->ReadInt(this->GetFile(),&offset);
    }

  offset = this->ReadImageDirectory(this->GetFile(),offset);
  do
    {
    // we count only image directories and not thumbnail images
    // subfiletype 0 = images
    // subfiletype 1 = thumbnails
    if(this->NewSubFileType == 0) 
      {
      i++;
      }
    finalOffset = offset;
    offset = this->ReadImageDirectory(this->GetFile(),offset);
    }while(i<imageCount && offset != 0);

  this->OffsetToLastAccessedImage = finalOffset;
  this->NumberOfLastAccessedImage = image;

  return finalOffset;
}

unsigned long vtkLSMReader::ReadImageDirectory(ifstream *f,unsigned long offset)
{
  unsigned short numberOfTags=0;
  unsigned long nextOffset = offset;
  
  //vtkDebugMacro(<<"Reading unsigned short from "<<offset<<"\n");
  numberOfTags = this->ReadUnsignedShort(f,&offset);
    //vtkDebugMacro(<<"Number of tags="<<numberOfTags<<"\n");
  for(int i=0;i<numberOfTags;i++)
    {   
    this->AnalyzeTag(f,offset);
    //vtkDebugMacro(<<"Tag analyed...\n");
    if(this->NewSubFileType == 1) {
    //vtkDebugMacro(<<"Found thumbnail, get next");
    break; //thumbnail image
    }
    offset = offset + 12;
    //vtkDebugMacro(<<"New offset="<<offset);
    }
  nextOffset += 2 + numberOfTags * 12;
    //vtkDebugMacro(<<"Next offset is "<<nextOffset);
  return this->ReadUnsignedInt(f,&nextOffset);
}

/*
void vtkLSMReader::DecodeHorizontalDifferencing(unsigned char *buffer, int size)
{
  for(int i=1;i<size;i++)
    {
      //printf("%d",*(buffer+i));
      *(buffer+i) = *(buffer+i) + *(buffer+i-1);
    }
//  printf("\n");
}
*/


//----------------------------------------------------------------------------
// Convert to Imaging API
int vtkLSMReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  unsigned long offset, imageOffset;;
  unsigned char *buf, *tempBuf;
  int size,readSize,numberOfPix,timepoint,channel;
  time_t start, end;
  int outExtent[6];
  //vtkDebugMacro(<<"Executing data. Foo.");

  
  // get the info object
  //printf("vtkLSMReader RequestData\n");
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  //vtkImageData *data = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkImageData *data;
  data = this->AllocateOutputData(outInfo->Get(vtkDataObject::DATA_OBJECT()));

   data->GetPointData()->GetScalars()->SetName("LSM Scalars");
    
    //data->SetExtent(outExtent);
    //data->AllocateScalars();
  //data->SetUpdateExtent(data->GetWholeExtent());

    
    data->GetExtent(outExtent);  
    //printf("Extent for LSM reader=%d,%d,%d,%d,%d,%d\n",PRT_EXT(outExtent));
  if(!this->Identifier)
    {
    vtkDebugMacro(<<"Can not execute data since information has not been executed!");
    return 0;
    }

  // if given time point or channel index is bigger than maximum,
  // we use maximum
  timepoint = (this->IntUpdateExtent[3]>this->GetNumberOfTimePoints()-1?this->GetNumberOfTimePoints()-1:this->IntUpdateExtent[3]);
  channel = (this->IntUpdateExtent[4]>this->GetNumberOfChannels()-1?this->GetNumberOfChannels()-1:this->IntUpdateExtent[4]);
  //printf("Timepoint=%d, channel=%d\n",timepoint,channel);
  //int nSlices = (outExtent[5]-outExtent[4])+1;
  //printf("Allocating memory for %d slices\n",nSlices);
  numberOfPix = this->Dimensions[0]*this->Dimensions[1]*this->Dimensions[2];
  //printf("Allocating %d bytes\n",numberOfPix);
  vtkDebugMacro(<<"Channel = "<<channel<<",numberOfPix="<<numberOfPix);
  size = numberOfPix*this->BYTES_BY_DATA_TYPE(this->DataType);

  // this buffer will be deleted by the vtkXXXArray when the array is destroyed.
  buf = new unsigned char[size];
  tempBuf = buf;

  start = time (NULL);
  for(int i=outExtent[4];i<=outExtent[5];i++)

  //for(int i=0;i<1;i++)
    {
    //UpdateProgress(i/float(this->Dimensions[2]));
       
//  UpdateProgress((i- outExtent[4])/
//                         (outExtent[5] - outExtent[4] + 1.0));        
    imageOffset = this->GetOffsetToImage(i,timepoint);
    vtkDebugMacro(<<"Offset to image "<<i<<"="<<imageOffset);
    this->ReadImageDirectory(this->GetFile(),imageOffset);
    offset = this->StripOffset->GetValue(channel);
    vtkDebugMacro(<<"Offset to channel"<<offset);
    readSize = this->StripByteCount->GetValue(channel);

    vtkDebugMacro(<<"Strip byte count="<<readSize);
    this->ReadFile(this->GetFile(),&offset,readSize,(char *)tempBuf,1);

    /*
    if(this->IsCompressed())
      {
      this->DecodeLZWCompression(tempBuf,readSize);
      //this->DecodeHorizontalDifferencing(tempBuf,readSize);
      }
    */
    tempBuf += readSize;
    }
  end = time (NULL);

  vtkDebugMacro(<<"Dataset generation time: "<<end-start);


  vtkUnsignedCharArray *uscarray;
  vtkUnsignedShortArray *ussarray;
  if(this->BYTES_BY_DATA_TYPE(this->DataType) > 1)
    {
    ussarray = vtkUnsignedShortArray::New();
    ussarray->SetNumberOfComponents(1);
    ussarray->SetNumberOfValues(numberOfPix);
      
    ussarray->SetArray((unsigned short *)buf, numberOfPix, 0);
    data->GetPointData()->SetScalars(ussarray);
    
    ussarray->Delete();     
    }
  else
    {
    uscarray = vtkUnsignedCharArray::New();
    uscarray->SetNumberOfComponents(1);
    uscarray->SetNumberOfValues(numberOfPix);
    
    uscarray->SetArray(buf, numberOfPix, 0);
    data->GetPointData()->SetScalars(uscarray);
    
    uscarray->Delete();     
    }
 
   //printf("Data executed\n"); 
    return 1;
  
}

int vtkLSMReader::RequestUpdateExtent (
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  int uext[6], ext[6];
    
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  //vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

  outInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),ext);
  // Get the requested update extent from the output.
  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), uext);
  
  // If they request an update extent that doesn't cover the whole slice
  // then modify the uextent 
  if(uext[1] < ext[1] ) uext[1] = ext[1];
  if(uext[3] < ext[3] ) uext[3] = ext[3];
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), uext,6);
  //request->Set(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT(), uext,6);
  return 1;    
}

int vtkLSMReader::RequestInformation (

  vtkInformation       * vtkNotUsed( request ),

  vtkInformationVector** vtkNotUsed( inputVector ),

  vtkInformationVector * outputVector)
{
  unsigned long startPos;
  unsigned int imageDirOffset;

  //char buf[12];
  
   //printf("RequestInformation\n");
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
   
  this->SetDataByteOrderToLittleEndian();

  if(!this->NeedToReadHeaderInformation())
    {
         //printf("Won't read info\n");
    vtkDebugMacro(<<"Don't need to read header information");
    return 1;
    }
    
  vtkDebugMacro(<<"Executing information.");

  if(!this->OpenFile())
    {
    this->Identifier = 0;
    return 0;
    }

  startPos = 2;  // header identifier

  this->Identifier = this->ReadUnsignedShort(this->GetFile(),&startPos);
  //printf("Identifier=%d\n",Identifier);
  if(!this->IsValidLSMFile())
    {
    vtkErrorMacro("Given file is not a valid LSM-file.");
    return 0;
    }
  
  imageDirOffset = this->ReadUnsignedInt(this->GetFile(),&startPos);
  //vtkDebugMacro(<<"Image dir offset="<<imageDirOffset<<"\n");
  // get information from the first image directory
  this->ReadImageDirectory(this->GetFile(),imageDirOffset);
  //vtkDebugMacro(<<"Read image directory\n");
  if(this->LSMSpecificInfoOffset)
    {        
      ReadLSMSpecificInfo(this->GetFile(),(unsigned long)this->LSMSpecificInfoOffset);
        //printf("Got LSM specific info\n");
    }
  else
    {
      vtkErrorMacro("Did not found LSM specific info!");
        //printf("Failed to read info\n");
      return 0;
    }
  if( !(this->ScanType == 6 || this->ScanType == 0 || this->ScanType == 3) )
    {
      vtkErrorMacro("Sorry! Your LSM-file must be of type 6 LSM-file (time series x-y-z) or type 0 (normal x-y-z) or type 3 (2D + time). Type of this File is " <<this->ScanType);
      return 0;
    }
  
  vtkDebugMacro(<<"Executing information: first directory has been read.");

  if(this->IsCompressed())
    {
      vtkDebugMacro("Can't handle compressed data!");
      return 0;
    }

    
  this->CalculateExtentAndSpacing(this->DataExtent,this->DataSpacing);
  //printf("%x Calculated extent of data, it is %d,%d,%d,%d,%d,%d\n",this,PRT_EXT(this->DataExtent));
  //printf("Spacing=%f,%f,%f\n",DataSpacing[0],DataSpacing[1],DataSpacing[2]);
    outInfo->Set(vtkDataObject::SPACING(), this->DataSpacing, 3);
  
    
//  this->GetOutput()->SetUpdateExtent(this->GetOutput()->GetWholeExtent());
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),

               this->DataExtent, 6);    
  
//  this->GetOutput()->SetNumberOfScalarComponents(1);
    this->NumberOfScalarComponents = 1;
  if(this->DataType > 1)
    {
      this->DataScalarType = VTK_UNSIGNED_SHORT;
      //this->GetOutput()->SetScalarType(VTK_UNSIGNED_SHORT);
    }
  else
    {
        this->DataScalarType = VTK_UNSIGNED_CHAR;
      //this->GetOutput()->SetScalarType(VTK_UNSIGNED_CHAR);
    }
  vtkDataObject::SetPointDataActiveScalarInfo(outInfo, this->DataScalarType,

    this->NumberOfScalarComponents);
    
  this->NeedToReadHeaderInformationOff();
  vtkDebugMacro(<<"Executing information: executed.");
    //printf("Read dimensions, %d,%d,%d,%d,%d\n",Dimensions[0],Dimensions[1],Dimensions[2],Dimensions[3],Dimensions[4]);
    return 1;
}

void vtkLSMReader::CalculateExtentAndSpacing(int extent[6],double spacing[3])
{
  
  extent[0] = extent[2] = extent[4] = 0;
  extent[1] = this->Dimensions[0]-1;
  extent[3] = this->Dimensions[1]-1;
  extent[5] = this->Dimensions[2]-1;

  // Instead of normalized spacing, return the real voxel sizes, only scaled to micrometer range.
  // This is so that the physical information is carried to the ITK side as well, when doing 
  // calculations
  spacing[0] = int(this->VoxelSizes[0]*1000000);
  if(spacing[0]<1)spacing[0]=1;
  spacing[1] = this->VoxelSizes[1]/this->VoxelSizes[0];
  spacing[2] = this->VoxelSizes[2]/this->VoxelSizes[0];
//  spacing[0] = this->VoxelSizes[0]*1000000;
//  spacing[1] = this->VoxelSizes[1]*1000000;
//  spacing[2] = this->VoxelSizes[2]*1000000;    
//  printf("Spacing is now %.3f, %.3f, %.3f\n",spacing[0],spacing[1],spacing[2]);
}

//----------------------------------------------------------------------------

int vtkLSMReader::GetChannelColorComponent(int ch, int component)
{
  if(ch < 0 || ch > this->GetNumberOfChannels()-1 || component < 0 || component > 2)
    {
        //printf("ch%d not in limits (%d)\n",ch,this->GetNumberOfChannels());
    return 0;
    }
    // printf("Returning component %d= %d\n",ch*3+component,*this->ChannelColors->GetPointer((ch*3)+component));
  
  //i+((colNum+1)*j  
  return *(this->ChannelColors->GetPointer((ch*3) + component));
}

vtkImageData* vtkLSMReader:: GetTimePointOutput(int timepoint, int channel)
{
  this->SetUpdateTimePoint(timepoint);
  this->SetUpdateChannel(channel);
  return this->GetOutput();
}

void vtkLSMReader::SetUpdateTimePoint(int timepoint)
{
  if(timepoint < 0 || timepoint == this->IntUpdateExtent[3]) 
    {
    return;
    }
  this->IntUpdateExtent[3] = timepoint;
  this->Modified();
}

void vtkLSMReader::SetUpdateChannel(int ch)
{
  if(ch < 0 || ch == this->IntUpdateExtent[4])
    {
    return;
    }
  this->IntUpdateExtent[4] = ch;
  this->Modified();
}

void vtkLSMReader::NeedToReadHeaderInformationOn()
{
  this->FileNameChanged = 1;
}

void vtkLSMReader::NeedToReadHeaderInformationOff()
{
  this->FileNameChanged = 0;
}

int vtkLSMReader::NeedToReadHeaderInformation()
{
  return this->FileNameChanged;
}

ifstream *vtkLSMReader::GetFile()
{
  return this->File;
}

int vtkLSMReader::BYTES_BY_DATA_TYPE(int type)
{
  int bytes = 1;
  switch(type)
    {
    case(1): 
      return 1;
    case(2):
      return 2;
    case(5):
      return 4;
    }
  return bytes;
}

int vtkLSMReader::TIFF_BYTES(unsigned short type)
{
  int bytes = 1;
  switch(type)
    {
    case(TIFF_BYTE): 
      return 1;
    case(TIFF_ASCII):
    case(TIFF_SHORT): 
      return 2;
    case(TIFF_LONG):
    case(TIFF_RATIONAL):
      return 4;
    }
  return bytes;
}

unsigned char vtkLSMReader::CharPointerToUnsignedChar(char *buf)
{
  return *((unsigned char*)(buf));
}

int vtkLSMReader::CharPointerToInt(char *buf)
{
  return *((int*)(buf));
}

unsigned int vtkLSMReader::CharPointerToUnsignedInt(char *buf)
{
  return *((unsigned int*)(buf));
}

short vtkLSMReader::CharPointerToShort(char *buf)
{
  return *((short*)(buf));
}

unsigned short vtkLSMReader::CharPointerToUnsignedShort(char *buf)
{
  return *((unsigned short*)(buf));
}

double vtkLSMReader::CharPointerToDouble(char *buf)
{
  return *((double*)(buf));
}

int vtkLSMReader::ReadInt(ifstream *f,unsigned long *pos)
{
  char buff[4];
  this->ReadFile(f,pos,4,buff);
#ifdef VTK_WORDS_BIGENDIAN
  vtkByteSwap::Swap4LE((int*)buff);
#endif
  return CharPointerToInt(buff);
}

unsigned int vtkLSMReader::ReadUnsignedInt(ifstream *f,unsigned long *pos)
{
  char buff[4];
  this->ReadFile(f,pos,4,buff);
#ifdef VTK_WORDS_BIGENDIAN
  vtkByteSwap::Swap4LE((unsigned int*)buff);
#endif
  return this->CharPointerToUnsignedInt(buff);
}

short vtkLSMReader::ReadShort(ifstream *f,unsigned long *pos)
{
  char buff[2];
  this->ReadFile(f,pos,2,buff);
#ifdef VTK_WORDS_BIGENDIAN
  vtkByteSwap::Swap2LE((short*)buff);
#endif  
  return this->CharPointerToShort(buff);
}

unsigned short vtkLSMReader::ReadUnsignedShort(ifstream *f,unsigned long *pos)
{
  char buff[2];
  this->ReadFile(f,pos,2,buff);
#ifdef VTK_WORDS_BIGENDIAN
  vtkByteSwap::Swap2LE((unsigned short*)buff);
#endif
  return this->CharPointerToUnsignedShort(buff);
}

double vtkLSMReader::ReadDouble(ifstream *f,unsigned long *pos)
{
  char buff[8];
  this->ReadFile(f,pos,8,buff);
#ifdef VTK_WORDS_BIGENDIAN
  vtkByteSwap::Swap8LE((double*)buff);
#endif  
  return this->CharPointerToDouble(buff);
}

int vtkLSMReader::ReadData(ifstream *f,unsigned long *pos,int size,char *buf)
{
  return this->ReadFile(f,pos,size,buf,1);
}

int vtkLSMReader::ReadFile(ifstream *f,unsigned long *pos,int size,char *buf,bool swap)
{
  f->seekg(*pos,ios::beg);
  f->read(buf,size);
#ifdef VTK_WORDS_BIGENDIAN
  if(swap) {
    vtkByteSwap::SwapLERange(buf,size);
  }      
#endif
  
  if( !f ) return 0;
  *pos = *pos + size;
  return 1;
}

void vtkLSMReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "Identifier: " << this->Identifier <<"\r\n";
  os << indent << "Dimensions: " << this->Dimensions[0] << "," << this->Dimensions[1] << ","<<this->Dimensions[2] << "\r\n";
  os << indent << "Time points: " << this->Dimensions[3] << "\r\n";
  os << "Number of channels: " << this->Dimensions[4] << "\r\n";
  os << "\r\n";
  os << indent << "Number of intensity values X: " << this->NumberOfIntensityValues[0] << "\r\n";
  os << indent << "Number of intensity values Y: " << this->NumberOfIntensityValues[1] << "\r\n";
  os << indent << "Number of intensity values Z: " << this->NumberOfIntensityValues[2] << "\r\n";
  os << indent << "Number of intensity values Time: " << this->NumberOfIntensityValues[3] << "\r\n";
  os << indent << "Voxel size X: " << this->VoxelSizes[0] << "\r\n";
  os << indent << "Voxel size Y: " << this->VoxelSizes[1] << "\r\n";
  os << indent << "Voxel size Z: " << this->VoxelSizes[2] << "\r\n";
  os << "\r\n";
  os << indent << "Scan type: " << this->ScanType << "\r\n";
  os << indent << "Data type: " << this->DataType << "\r\n";
  os << indent << "Compression: " << this->Compression << "\r\n";
  os << "\r\n";
  os << indent << "Planar configuration: " << this->PlanarConfiguration << "\r\n";
  os << indent << "Photometric interpretation: " << this->PhotometricInterpretation << "\r\n";
  os << indent << "Predictor: " << this->Predictor << "\r\n";
  os << indent << "Channel info:\n";

  for(int i=0;i<this->Dimensions[4];i++)
    {
        os << indent << indent << this->GetChannelName(i)<<",("<<this->GetChannelColorComponent(i,0)<<","<<this->GetChannelColorComponent(i,1)<<","<<this->GetChannelColorComponent(i,2)<<")\n";
    }
  os << indent << "Strip byte counts:\n";

  for(int i=0;i<this->Dimensions[4];i++)
    {
      os << indent << indent << this->StripByteCount->GetValue(i) << std::endl;
    }
}

