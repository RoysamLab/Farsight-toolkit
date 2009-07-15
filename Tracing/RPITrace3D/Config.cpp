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

///////////////////////////////////////////////////////////////////////////////
// File: Config.cpp
// By: Khalid Al-Kofahi
// Date: 1-10-2002
//
//
//      ****** Description ******
//
// This file loads and saves the parameters, paths, etc, set in the 
// configuration file specified at the command prompt. See the configuration
// file for the required format.
// 
// At the time the original tracing program was written STL was not part of the
// STD library. Hence, I did not use it. This class could've been implemented
// using a simple multimap, but I did not want to introduce STL at this point!
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <string>
//#include "getopt.h"

#include "Config.h"

//By Yousef: 10-25-2007
//The header file needed for xml parsing
#include "tinyxml/tinyxml.h"

//#pragma comment( lib, "libxml2.lib" ) 
/////////////////////////////////////////////

using namespace std;

CConfig::CConfig()
{
  m_aData = NULL;
  m_iNumOfElements = 0;
  m_bReadConfigFromFile = false;
  m_iImagePadding = 0;

  //load default values

  m_Tracing.m_iGridSpacing            = 5;  //15 //ch00:15
  m_Tracing.m_iMinimumTemplateLength        = 10;  //10 //ch00:10
  m_Tracing.m_iMaximumTemplateLength        = 18;  //18 //ch00:18 
  m_Tracing.m_iRelativeShiftDistance        = 2;//2     //ch00:2
  m_Tracing.m_iDirectionalDegreeOfFreedom     = 7;//7     //ch00:7
  m_Tracing.m_iMaximumStepSize          = 10;//10   //ch00:10
  m_Tracing.m_iContrastThresholdMultiplier    = 4;  //3  //ch00:10
  m_Tracing.m_iMaximumAllowedStoppingViolations = 3;  //2   //ch00:3
  m_Tracing.m_iMinimumShiftDistance       = 2;//2     //ch00:2
  m_Tracing.m_iMaximumShiftDistance       = 20;  //10 //ch00:10
  m_Tracing.m_bMedianTracing            = true;

  m_Tracing.m_bDoNotMergeTraces         = false;    
  //Yousef
  m_bDetectBranches = false;
  ////////

  m_Soma.m_bDetectSoma = false;
  m_Soma.m_iMinSomaArea = 50; //100
  m_Soma.m_fSomaSensitivity = 0.55f;
  m_Soma.m_ReadSomaImage = false;

  m_Output.m_bSaveVesselLengthInFile = false;
  m_Output.m_bLabelVesselsOnImage = false;
  m_Output.m_bFlipY = true;
  m_Output.m_bHeatedObject = false;
  m_Output.m_bDisplayVerifiedSeeds = false;
  m_Output.m_bColorTrace = false;
  m_Output.m_bWriteSeeds = false;
  m_Output.m_bSaveTracedPoints = true;  
  m_Output.m_bWriteResultOutputFiles = true;

  m_QA.m_bComputeQA = false;
  m_QA.alpha      = 0.5f;
  m_QA.w        = 3.1f;
  m_QA.epsilon    = 0.1f; 

  m_bDirectionality = false;
  m_bNeurolucidaFillSoma = false;m_bImagewideSeedResponse = false;

  m_bReadSeedCandidatesFile = true;
  m_bReadVerifiedSeedsFile = false;

  m_fScale = 1.0;
  m_b10XFlag = true;
  m_ProcessType = std::string("dendrite");

  m_bDisableVesselMerging = false;

  m_bWriteOutputFiles = true;

  
}

// Utility Function for comparing two data entries  
int CompareConfigEntryName(const void* A, const void* B)
{
  return (strcmp(((TConfigEntry *) A)->m_achName,
    ((TConfigEntry *) B)->m_achName));
}



///////////////////////////////////////////////////////////////////////////////
// Method: ReadConfigurationFile
// 
// Read the configuration file and set the parameters
void CConfig::ReadConfigurationFile(const char* achFName)
{
  char achBuffer[1024];
  char achBuffer2[1024];
  char* pchStr = NULL;
  int iCounter = 0;

  ifstream inFile(achFName);
  if (! inFile)
  {
    cout << "Fatal Error: Could not open configuration file " << achFName
      << " .... Terminating Program." << endl;
    exit(0);
  }

  // count the number of non-comment lines in the file    
  while (inFile)
  {
    inFile.getline(achBuffer, 1024, '\n');
    if (achBuffer[0] != '\0' && achBuffer[0] != '!')
      iCounter++;
  }
  if (iCounter == 0)
  {
    cout << "Fatal Error: Empty Configuration File " << achFName
      << " .... Terminating Program." << endl;
    exit(0);
  }
  // create The data, and read it
  m_aData = new TConfigEntry[iCounter];
  inFile.close();
  inFile.open(achFName);
  ifstream inFile2(achFName);
  while (inFile2)
  {
    inFile2.getline(achBuffer, 1024, '\n');
    if (achBuffer[0] != '\0' && achBuffer[0] != '!')
    {
      strcpy(achBuffer2, achBuffer);

      pchStr = strtok(achBuffer, "\t ");

      if (pchStr)
      {
        strcpy(m_aData[m_iNumOfElements].m_achName, pchStr);
      }

      pchStr = strtok(NULL, "\t ");

      if (pchStr)
      {
        if (*pchStr == ':')
          pchStr = strtok(NULL, "\t " );

        strcpy(m_aData[m_iNumOfElements].m_achValue, pchStr);
      }
      m_iNumOfElements++;
    }
  }

  inFile2.close();  

  // sort the data so that we can search it
  qsort(m_aData, m_iNumOfElements, sizeof(TConfigEntry), CompareConfigEntryName);
  SetValuesFromConfigurationFile();
}

void CConfig::SetValuesFromConfigurationFile()
{
  char* pchValue = NULL;

  // Trace Parameters
  if ((pchValue = GetStringValue("GridSpacing")) != NULL)
    SetGridSpacing(atoi(pchValue));
  if ((pchValue = GetStringValue("MinimumTemplateLength")) != NULL)
    SetMinimumTemplateLength(atoi(pchValue));
  if ((pchValue = GetStringValue("MaximumTemplateLength")) != NULL)
    SetMaximumTemplateLength(atoi(pchValue));
  if ((pchValue = GetStringValue("RelativeShiftDistance")) != NULL)
    SetRelativeShiftDistance(atoi(pchValue));
  if ((pchValue = GetStringValue("DirectionalDegreeOfFreedom")) != NULL)
    SetDirectionalDegreeOfFreedom(atoi(pchValue));
  if ((pchValue = GetStringValue("MaximumStepSize")) != NULL)
    SetMaximumStepSize(atoi(pchValue));
  if ((pchValue = GetStringValue("ContrastThresholdMultiplier")) != NULL)
    SetContrastThresholdMultiplier(atoi(pchValue));
  if ((pchValue = GetStringValue("MaximumAllowedStoppingViolations")) != NULL)
    SetMaximumAllowedStoppingViolations(atoi(pchValue));
  if ((pchValue = GetStringValue("MinimumShiftDistance")) != NULL)
    SetMinimumShiftDistance(atoi(pchValue));
  if ((pchValue = GetStringValue("MaximumShiftDistance")) != NULL)
    SetMaximumShiftDistance(atoi(pchValue));
  if ((pchValue = GetStringValue("MedianTracing")) != NULL) 
    SetMedianTracing((strcmp(pchValue, "yes") == 0 || strcmp(pchValue, "y") == 0));
  std::cout << (strcmp(pchValue, "yes") == 0 || strcmp(pchValue, "y") == 0) << std::endl;
  std::cout << "t" << pchValue << "      t" << std::endl;

}

///////////////////////////////////////////////////////////////////////////////
// Method: GetStringValue(char achString)
//
// Get the value of the given string. Return NULL if not found
// The algorithm uses binary search
char* CConfig::GetStringValue(const char* achString)
{
  int iHi = m_iNumOfElements - 1;
  int iLow = 0;
  int iMid = 0;
  int iComp = 0;

  while (iLow <= iHi)
  {
    iMid = (iLow + iHi) / 2;
    iComp = strcmp(m_aData[iMid].m_achName, achString);
    if (iComp < 0)
      iLow = iMid + 1;
    else if (iComp > 0)
      iHi = iMid - 1;
    else
      return m_aData[iMid].m_achValue;
  }

  // string not found
  return NULL;
}

//convenience method to handle error checking when retrieving XML attributes
const char * CConfig::GetXMLAttribute(TiXmlElement *element, const char *attributeName)
{
  const char *attribute = element->Attribute(attributeName);
  if(!attribute)
    {
    cerr << "ERROR: No XML attribute called '" << attributeName <<
            "' in element '" << element->Value() << "'" << endl;
    exit(0);
    }
  return attribute;
}

//By Yousef 10/25/2007
//This function reads the tracing settings from an xml document.
//The name of the xml file is sent as an input argument to the tracing program
//Updated by Zack 2/18/09 to use TinyXML
void CConfig::XmlReadSettings (std::string filename)
{
  TiXmlDocument doc(filename.c_str());
  doc.LoadFile();
  TiXmlHandle docHandle( &doc );
  TiXmlElement* parameterElement = docHandle.FirstChild("parameters").FirstChild("parameter").Element();
  const char *parameterName;
  const char *parameterValue;
  while(parameterElement)
    {
    parameterName = this->GetXMLAttribute(parameterElement, "name");
    parameterValue = this->GetXMLAttribute(parameterElement, "value");
    if(strcmp(parameterName, "input_image") == 0)
      {
      if(strstr(parameterValue, ".pic") || strstr(parameterValue, ".PIC"))
        {
        this->m_ImageType= "pic";
        }
      else
        {
        if(strstr(parameterValue, ".txt") || strstr(parameterValue, ".TXT") )
          {
          this->m_ImageType = "pgm";
          }
        else
          {
          cout << "Unrecognizable image type: " << parameterValue << endl;
          exit(0);
          }
        }
      string imageNameStr = parameterValue;
      this->m_ImageName = imageNameStr.substr(0,imageNameStr.rfind("."));
      }
    else if(strcmp(parameterName, "SeedsSamplingDensity") == 0)
      {
      this->SetGridSpacing(atoi(parameterValue));
      }
    else if(strcmp(parameterName, "MinTemplateLength") == 0)
      {
      this->SetMinimumTemplateLength(atoi(parameterValue));
      }
    else if(strcmp(parameterName, "MaxTemplateLength") == 0)
      {
      this->SetMaximumTemplateLength(atoi(parameterValue));
      }
    else if(strcmp(parameterName, "MinShiftDistance") == 0)
      {
      this->SetMinimumShiftDistance(atoi(parameterValue));
      }
    else if(strcmp(parameterName, "MaxShiftDistance") == 0)
      {
      this->SetMaximumShiftDistance(atoi(parameterValue));
      }
    else if(strcmp(parameterName, "DirectionFreedom") == 0)
      {
      this->SetDirectionalDegreeOfFreedom(atoi(parameterValue));
      }
    else if(strcmp(parameterName, "StepSize") == 0)
      {
      this->SetMaximumStepSize(atoi(parameterValue));
      }
    else if(strcmp(parameterName, "TracingSensitivity") == 0)
      {
      this->SetContrastThresholdMultiplier(atoi(parameterValue));
      }
    else if(strcmp(parameterName, "MaxAllowedStoppingViolations") == 0)
      {
      this->SetMaximumAllowedStoppingViolations(atoi(parameterValue));
      }
    else if(strcmp(parameterName, "RelativeShiftDistance") == 0)
      {
      this->SetRelativeShiftDistance(atoi(parameterValue));
      }
    else if(strcmp(parameterName, "DetectSoma") == 0)
      {
      if( *parameterValue == 'Y' || *parameterValue == 'y')
        {
        this->m_Soma.m_bDetectSoma = true;
        }
      else
        {
        this->m_Soma.m_bDetectSoma = false;
        }
      }
    else if(strcmp(parameterName, "MergeTraceEndPoints") == 0)
      {
      if(*parameterValue == 'Y' || *parameterValue == 'y')
        {
        this->m_bDisableVesselMerging = false;          
        }
      else
        {
        this->m_bDisableVesselMerging = true;
        }
      }
    else if(strcmp(parameterName, "DetectBranchPoints") == 0)
      {
      if( *parameterValue == 'Y' || *parameterValue == 'y')
        {
        this->m_bDetectBranches = true;
        }
      else
        {
        this->m_bDetectBranches = false;
        }
      }
    //we don't currently use the "MaxBranchSearchDistance" attribute
    parameterElement = parameterElement->NextSiblingElement();
    }
}

void CConfig::SetGridSpacing (int value)
{
  m_Tracing.m_iGridSpacing = value;
  std::cout << "Option: Grid spacing = " << m_Tracing.m_iGridSpacing << std::endl;
}

void CConfig::SetMinimumTemplateLength(int value)
{
  m_Tracing.m_iMinimumTemplateLength = value;
  std::cout << "Option: Minimum template length = " << m_Tracing.m_iMinimumTemplateLength << std::endl;
}

void CConfig::SetMaximumTemplateLength(int value)
{
  m_Tracing.m_iMaximumTemplateLength = value;
  std::cout << "Option: Maximum template length = " << m_Tracing.m_iMaximumTemplateLength << std::endl;
}

void CConfig::SetRelativeShiftDistance(int value)
{
  m_Tracing.m_iRelativeShiftDistance = value;
  std::cout << "Option: Relative shift distance = " << m_Tracing.m_iRelativeShiftDistance << std::endl;
}

void CConfig::SetDirectionalDegreeOfFreedom(int value)
{
  m_Tracing.m_iDirectionalDegreeOfFreedom = value;
  std::cout << "Option: Directional degree of freedom = " << m_Tracing.m_iDirectionalDegreeOfFreedom << std::endl;
}

void CConfig::SetMaximumStepSize(int value)
{
  m_Tracing.m_iMaximumStepSize = value;
  std::cout << "Option: Maximum step size = " << m_Tracing.m_iMaximumStepSize << std::endl;
}

void CConfig::SetContrastThresholdMultiplier(int value)
{
  m_Tracing.m_iContrastThresholdMultiplier = value;
  std::cout << "Option: Contrast threshold multiplier = " << m_Tracing.m_iContrastThresholdMultiplier << std::endl;
}

void CConfig::SetMaximumAllowedStoppingViolations(int value)
{
  m_Tracing.m_iMaximumAllowedStoppingViolations = value;
  std::cout << "Option: Maximum allowed stopping violations = " << m_Tracing.m_iMaximumAllowedStoppingViolations << std::endl;
} 

void CConfig::SetMinimumShiftDistance(int value)
{
  m_Tracing.m_iMinimumShiftDistance = value;
  std::cout << "Option: Minimum shift distance = " << m_Tracing.m_iMinimumShiftDistance << std::endl;
}

void CConfig::SetMaximumShiftDistance(int value)
{
  m_Tracing.m_iMaximumShiftDistance = value;
  std::cout << "Option: Maximum shift distance = " << m_Tracing.m_iMaximumShiftDistance << std::endl;
}

void CConfig::SetMedianTracing(bool value)
{
  m_Tracing.m_bMedianTracing = value;
  std::cout << "Option: Median tracing = " << m_Tracing.m_bMedianTracing << std::endl;
}

void CConfig::SetImagePadding(int value)
{
  m_iImagePadding = value;
  std::cout << "Option: Image padding = " << m_iImagePadding << std::endl;
}

void CConfig::Options()
{
  cout << "Options: \n"  
    << "\tSome options can be in a longer form\n"
    << "\te.g., -s 3 and --step=3 sets the step size to 3\n\n"
    << "-G <ConfigFileName>\n"
    << "\tRead configurations from <ConfigFileName>\n"
    << "\n--g=val or -g val\n"
    << "\tGrid spacing=val\n"
    << "\n--L_min=val or -a val\n"
    << "\tMinimum template length=val\n"
    << "\n--L_max=val or -b val\n"
    << "\tMaximum template length=val\n"
    << "\n--n_shift=val or -r val\n"
    << "\tRelative shift distance=val\n"
    << "\n--n_rotate=val or -o val\n"
    << "\tDirectional degree of freedom=val\n"
    << "\n--step=val or -s val\n"
    << "\tMaximum step size=val\n"
    << "\n--t_c=val or -c val\n"
    << "\tContrast threshold multiplier=val\n"
    << "\n--v=val or -v val\n"
    << "\tMaximum allowed stopping violations=val\n"
    << "\n--input_path=val\n"
    << "\tSet path for input images to val\n"
    << "\n--output_path=val\n"
    << "\tSet path for output images and files to val\n"
    << "\n--image_padding=val\n"
    << "\tSet image padding in z-axis to val\n"
    << "\n--detectsoma\n"
    << "\tEnable soma detection\n"
    << "\n--nomerge\n"
    << "\tDisable trace merging\n"
    << "\n--options\n"
    << "\tPrints this message\n"    
    << endl;
  exit(0);
}

void CConfig::Usage()
{
  cout << "\n\nUsage:\n" 
    << "\t<Executable> <ImageFileName> <Options>\n"
    << "\nWhere,\n"
    << "\t<ImageFileName> is the name of the 3D image in BioRad Pic Format"
    << "\t                or a txt file with list of PGM slices, 1 per line\n"
    << "\n  use --options for help on available options\n"
    << "\n==========> ALL COMMANDS ARE CASE SENSITIVE <==========\n"
    << endl;
  exit(0);
}
