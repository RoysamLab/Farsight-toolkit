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
// File: Config.h
//
// The declaration for the configuration file
//
#ifndef _CONFIG_H_
#define _CONFIG_H_

class TiXmlElement;

typedef struct _ConfigEntry
{
	char m_achName[128];
	char m_achValue[256];
} TConfigEntry;

class CConfig
{
public:

	CConfig();

	~CConfig()
	{
		if (m_aData)
			delete [] m_aData;
		m_iNumOfElements = 0;
	}

	// read the configuration file
	void ReadConfigurationFile(const char* achFName);

	// get the value associated with the given string
	char* GetStringValue(const char* pchStringName);

	/////////////////
	// Data
	//
	// an array of configuration entries.
	// value
	TConfigEntry* m_aData;
	// the number of elements
	int m_iNumOfElements;
	void ProcessCommandLine (int argc, char *argv[]);
  const char* GetXMLAttribute(TiXmlElement *element, const char *attributeName);
	//By Yousef 10-25-2007/////////////////////
	void XmlReadSettings(std::string filename);
	///////////////////////////////////////////

	// Get configuration items
	int GetGridSpacing () { return m_Tracing.m_iGridSpacing; }
	int GetMinimumTemplateLength() { return m_Tracing.m_iMinimumTemplateLength; }
	int GetMaximumTemplateLength() { return m_Tracing.m_iMaximumTemplateLength; }
	int GetRelativeShiftDistance() { return m_Tracing.m_iRelativeShiftDistance; }
	int GetDirectionalDegreeOfFreedom() { return m_Tracing.m_iDirectionalDegreeOfFreedom; }
	int GetMaximumStepSize() { return m_Tracing.m_iMaximumStepSize; }
	int GetContrastThresholdMultiplier() { return m_Tracing.m_iContrastThresholdMultiplier; }
	int GetMaximumAllowedStoppingViolations() { return m_Tracing.m_iMaximumAllowedStoppingViolations; }
	int GetMinimumShiftDistance() { return m_Tracing.m_iMinimumShiftDistance; }
	int GetMaximumShiftDistance() { return m_Tracing.m_iMaximumShiftDistance; }
	bool GetMedianTracing() { return m_Tracing.m_bMedianTracing; }

	int GetImagePadding() { return m_iImagePadding; }
	std::string GetOutputPath() { return m_OutputPath; }
	std::string GetInputPath() { return m_InputPath; }
	std::string GetImageName() { return m_ImageName; }
	std::string GetImageType() { return m_ImageType; }
	std::string GetVesselnessImage() { return m_VesselnessImage; }

	bool GetDetectSoma() { return m_Soma.m_bDetectSoma; }
	
	bool GetDisableVesselMerging() { return m_bDisableVesselMerging; }

	bool GetWriteOutputFiles() { return m_bWriteOutputFiles; }

	//Yousef
	bool GetDetectBranches() { return m_bDetectBranches; }

	bool GetQA() { return m_QA.m_bComputeQA; }

private:
	
	void SetValuesFromConfigurationFile();
	void Usage();
	void Options();

	int ConvertToHex (char chChar, int iOrder);

	class TracingParameters {
	public:
		int		m_iGridSpacing;  
		int		m_iMinimumTemplateLength; 
		int		m_iMaximumTemplateLength;
		int		m_iRelativeShiftDistance;
		int		m_iDirectionalDegreeOfFreedom;
		int		m_iMaximumStepSize;
		int		m_iContrastThresholdMultiplier;
		int		m_iMaximumAllowedStoppingViolations;  
		bool	m_bDoNotMergeTraces;
		int		m_iMinimumShiftDistance;
		int		m_iMaximumShiftDistance;
		bool	m_bMedianTracing;
	} m_Tracing;
	void SetGridSpacing (int value);
	void SetMinimumTemplateLength(int value);
	void SetMaximumTemplateLength(int value);
	void SetRelativeShiftDistance(int value);
	void SetDirectionalDegreeOfFreedom(int value);
	void SetMaximumStepSize(int value);
	void SetContrastThresholdMultiplier(int value);
	void SetMaximumAllowedStoppingViolations(int value);
	void SetMinimumShiftDistance(int value);
	void SetMaximumShiftDistance(int value);
	void SetMedianTracing(bool value);  // "yes" or "no"

	class SomaParameters {
	public:
		bool	m_bDetectSoma;
		int		m_iMinSomaArea;
		float	m_fSomaSensitivity;
		bool	m_ReadSomaImage;
		std::string m_SomaImageFN;
	} m_Soma;

	class OutputOptions {
	public:
		bool	m_bSaveVesselLengthInFile; 
		bool	m_bSaveTracedPoints;
		bool	m_bColorTrace;
		bool	m_bWriteSeeds;
		bool	m_bLabelVesselsOnImage;
		bool	m_bFlipY;
		bool	m_bHeatedObject;
		bool	m_bDisplayVerifiedSeeds;
		bool	m_bWriteResultOutputFiles;
	} m_Output;

	class QualityAssessmentOptions {
	public:
		bool	m_bComputeQA;
		float	alpha;
		float	w;			// representation efficiency coefficient
		float	epsilon;
	} m_QA;
	
	bool			m_bReadConfigFromFile;
	std::string		m_ConfigFileName;
	std::string		m_ImageName;
	std::string		m_ImageType;
	int				m_iImagePadding;
	void			SetImagePadding(int value);
	std::string		m_InputPath;
	std::string		m_OutputPath;
	std::string		m_VesselnessImage;

	bool	m_bDisableVesselMerging;

	bool	m_bWriteOutputFiles;

	bool	m_bDirectionality;
	
	bool	m_bNeurolucidaFillSoma;
	bool	m_bImagewideSeedResponse;
	bool	m_bReadSeedCandidatesFile;
	bool	m_bReadVerifiedSeedsFile;

	//Try this..Yousef
	bool	m_bDetectBranches;
	
	std::string m_SeedCandidatesFN;
	std::string m_VerifiedSeedsFN;

	std::string	m_ProcessType;
	
	float	m_fScale;

	std::string m_InputImageFullPath;

	// The image is 10X magnification. This is used to decide which detect soma function to be used.
	bool	m_b10XFlag;
};

#endif
