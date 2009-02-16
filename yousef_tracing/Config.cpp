///////////////////////////////////////////////////////////////////////////////
// File: Config.cpp
// By: Khalid Al-Kofahi
// Date: 1-10-2002
//
//
//			****** Description ******
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
//#include <getopt.h>

#include "Config.h"

using namespace std;

CConfig::CConfig()
{
	m_aData = NULL;
	m_iNumOfElements = 0;
	m_bReadConfigFromFile = false;
	m_iImagePadding = 0;

	//load default values

	m_Tracing.m_iGridSpacing						= 15;  
	m_Tracing.m_iMinimumTemplateLength				= 10; 
	m_Tracing.m_iMaximumTemplateLength				= 18;		
	m_Tracing.m_iRelativeShiftDistance				= 2;
	m_Tracing.m_iDirectionalDegreeOfFreedom			= 7;
	m_Tracing.m_iMaximumStepSize					= 3;
	m_Tracing.m_iContrastThresholdMultiplier		= 3;
	m_Tracing.m_iMaximumAllowedStoppingViolations	= 1;  
	m_Tracing.m_bDoNotMergeTraces					= false;
	m_Tracing.m_iMinimumShiftDistance				= 1;
	m_Tracing.m_iMaximumShiftDistance				= 30;
	m_Tracing.m_bMedianTracing						= true;

	m_Soma.m_bDetectSoma = true;
	m_Soma.m_iMinSomaArea = 100;
	m_Soma.m_fSomaSensitivity = 0.55f;
	m_Soma.m_ReadSomaImage = false;

	m_Output.m_bSaveVesselLengthInFile = false;
	m_Output.m_bLabelVesselsOnImage = false;
	m_Output.m_bFlipY = true;
	m_Output.m_bHeatedObject = false;
	m_Output.m_bDisplayVerifiedSeeds = false;
	m_Output.m_bColorTrace = true;
	m_Output.m_bWriteSeeds = false;
	m_Output.m_bSaveTracedPoints = true;	
	m_Output.m_bWriteResultOutputFiles = true;

	m_QA.m_bComputeQA	= false;
	m_QA.alpha			= 0.5f;
	m_QA.w				= 3.1f;
	m_QA.epsilon		= 0.1f;	

	m_bDirectionality = false;
	m_bNeurolucidaFillSoma = false;
	m_bImagewideSeedResponse = false;

	m_bReadSeedCandidatesFile = false;
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
char* CConfig::GetStringValue(char* achString)
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

void CConfig::ProcessCommandLine (int argc, char *argv[])
{

	if (argc < 2)
	{
		Usage();
		exit(0);
	}

	string ImageFileName = argv[1];
	if (strstr(ImageFileName.c_str(), ".pic") || strstr(ImageFileName.c_str(), ".PIC"))
		m_ImageType= "pic";
	else
	{
		if( strstr(ImageFileName.c_str(), ".txt") || strstr(ImageFileName.c_str(), ".TXT") )
		{
			m_ImageType = "pgm";
		}
		else
		{
			cout << "Unrecognizable image type." << endl;
			exit(0);
		}
	}
	
	m_ImageName = ImageFileName.substr(0,ImageFileName.rfind("."));

	if (argc == 2) // leave all options at default value
		return;

	//Comment the rest of the function for now until we add a new parameters reader that does not use getopt library
//	int c;
//	while (1)
//	{
//		static struct option long_options[] =
//		{
//			/* These options set a flag. */
////			{"verbose", no_argument,       &verbose_flag, 1},
//			
//			// Tracing Options
//			{"g",		required_argument, 0, 'g'}, // Grid spacing
//			{"L_min",	required_argument, 0, 'a'}, // Minimum template length
//			{"L_max",	required_argument, 0, 'b'}, // Maximum template length
//			{"n_shift", required_argument, 0, 'r'}, // Relative shift distance
//			{"n_rotate",required_argument, 0, 'o'}, // Directional degree of freedom
//			{"step",	required_argument, 0, 's'}, // Maximum step size
//			{"t_c",		required_argument, 0, 'c'}, // Contrast threshold multiplier
//			{"v",		required_argument, 0, 'v'}, // Maximum allowed stopping violations
//			{"median_trace",	required_argument, 0, 2}, // Trace using median responses
//
//			{"minsd",	required_argument, 0, 'i'}, // Minimum shift distance
//			{"maxsd",	required_argument, 0, 'j'}, // Maximum shift distance
//
//			
//
//			{"help",	no_argument, 0,	0}, // Display usage message
//			{"options",	no_argument, 0,	1}, // Display available options
//
//			{"input_path",	required_argument, 0,	3}, // input path
//			{"output_path", required_argument, 0,	4}, // input path
//
//			{"image_padding", required_argument, 0,	5}, // input path
//			{"vesselness_image", required_argument, 0, 6}, // filename of the vesselness image
//
//			{"qa", no_argument, 0, 'q'},
//			{"detectsoma",	no_argument, 0,	7}, // Enable soma detection
//			{"nomerge",	no_argument, 0,	8}, // Disable vessel merging
//			{"noresults",	no_argument, 0,	9}, // Don't output files 
//	
//			// IO options
//			{"G",		required_argument, 0, 'G'}, // Read configuration from file
//
//			{0, 0, 0, 0}
//		};
//		/* getopt_long stores the option index here. */
//		int option_index = 0;
//
//		c = getopt_long (argc, argv, "g:a:b:r:o:s:c:v:G:i:j:2:3:4:5:6:q",
//			long_options, &option_index);
//
//		/* Detect the end of the options. */
//		if (c == -1)
//			break;
//
//		switch (c)
//		{
//		case 0:
//			Usage();
//			break;
//		case 1:
//			Options();
//			break;
//		case 2:
//			SetMedianTracing((strcmp(optarg, "yes") == 0 || strcmp(optarg, "y") == 0));
//			break;
//		case 3:
//			m_InputPath = std::string(optarg);
//			break;
//		case 4:
//			m_OutputPath = std::string(optarg);
//			break;
//		case 5:
//			SetImagePadding(atoi(optarg));			
//			break;
//		case 6:
//			m_VesselnessImage = std::string(optarg);
//			break;
//		case 7:
//			std::cout << "Option: Soma Detection = Enabled" << std::endl;
//			m_Soma.m_bDetectSoma = true;
//			break;
//		case 8:
//			std::cout << "Option: Vessel Merging = Disabled" << std::endl;
//			m_bDisableVesselMerging = true;
//			break;
//		case 9:
//			std::cout << "Option: Write Output Files = Disabled" << std::endl;
//			m_bWriteOutputFiles = false;
//			break;
//
//		case 'q':
//			m_QA.m_bComputeQA = true;
//			break;
//
//		case 'g':		// Grid spacing
//			SetGridSpacing (atoi(optarg));
//			break;
//		case 'a':		// Minimum template length
//			SetMinimumTemplateLength (atoi(optarg));
//			break;
//		case 'b':		// Maximum template length
//			SetMaximumTemplateLength (atoi(optarg));
//			break;
//		case 'r':		// Relative shift distance
//			SetRelativeShiftDistance (atoi(optarg));
//			break;
//		case 'o':		// Directional degree of freedom
//			SetDirectionalDegreeOfFreedom (atoi(optarg));
//			break;
//		case 's':		// Maximum step size
//			SetMaximumStepSize (atoi(optarg));
//			break;
//		case 'c':		// Contrast threshold multiplier
//			SetContrastThresholdMultiplier (atoi(optarg));
//			break;
//		case 'v':		// Maximum allowed stopping violations
//			SetMaximumAllowedStoppingViolations (atoi(optarg));
//			break;
//		case 'i':		// Maximum allowed stopping violations
//			SetMinimumShiftDistance (atoi(optarg));
//			break;
//		case 'j':		// Maximum allowed stopping violations
//			SetMaximumShiftDistance (atoi(optarg));
//			break;
//
//		case 'G':
//			//m_bReadConfigFromFile = true;
//			//m_ConfigFileName = std::string(optarg);
//			std::cout << "Reading configfile is not supported yet" << std::endl;
//			break;
//
//		case '?':
//			/* getopt_long already printed an error message. */
//			break;
//
//		default:
//			abort ();
//		}
//	}
//
//	if(m_bReadConfigFromFile)
//		ReadConfigurationFile(m_ConfigFileName.c_str());
//
//	std::cout << "Input image: " << m_ImageName << " (type: " << m_ImageType << " )" << std::endl;
//	if(m_QA.m_bComputeQA)
//		std::cout << "Vesselness image: " << m_VesselnessImage << std::endl;
//	std::cout << "Input path: " << m_InputPath << std::endl;
//	std::cout << "Output path: " << m_OutputPath << std::endl;
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
		<< "-G <ConfigFileName>       Read configurations from <ConfigFileName>\n"
		<< "--g=val or -g val         Grid spacing=val\n"
		<< "--L_min=val or -a val     Minimum template length=val\n"
		<< "--L_max=val or -b val     Maximum template length=val\n"
		<< "--n_shift=val or -r val   Relative shift distance=val\n"
		<< "--n_rotate=val or -o val  Directional degree of freedom=val\n"
		<< "--step=val or -s val      Maximum step size=val\n"
		<< "--t_c=val or -c val       Contrast threshold multiplier=val\n"
		<< "--v=val or -v val         Maximum allowed stopping violations=val\n"
		<< "--options                 Prints this message\n"
		<< "\n==========> ALL OPTIONS ARE CASE SENSITIVE <==========\n"
		<< endl;
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
}
