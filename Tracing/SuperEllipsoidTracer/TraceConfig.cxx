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

//tracing_config.cpp
#include "TraceConfig.h"

#include <string>


#define MY_ENCODING "ISO-8859-1"

TraceConfig::TraceConfig()  {
     GridSpacing = 15;
     StepRatio = 0.5;
     AspectRatio = 1.35;
     THRESHOLD = 0.5;
     minContrast = 3.0;
     MaximumVesselWidth = 20.0;
     MinimumVesselLength = 2.0;
	 FitIterations = 100.0;
     MinimumVesselWidth = 1.5;
	 StartTHRESHOLD = 0.3;
	 Spacing[0] = 0.6;
	 Spacing[1] = 0.6;
	 Spacing[2] = 1.2;
	 numDataFiles = 0;
	 UseMultiscaleHessianFilter = 0;
}

TraceConfig::~TraceConfig() { }

void TraceConfig::SetFileNames(char* fname) 	{
	this->InputFileNames.push_back(std::string(fname));
	this->OutputFileNames.push_back(std::string(fname) + std::string("_SEs.txt"));
}

void TraceConfig::SetGridSpacing(char * gs) 	{
	this->GridSpacing = atoi(gs);
}

void TraceConfig::SetAspectRatio(char * ar)	{
	this->AspectRatio = atof(ar);
}

/*
bool TraceConfig::LoadParameters(char * ParameterFileName)  {

	TiXmlDocument doc(ParameterFileName);
	if (!doc.LoadFile()) {
		return false;
	}

	TiXmlElement* root = doc.FirstChildElement( "SuperElliposoidTracing" );
	if ( root )		{
		TiXmlAttribute* rAttrib = root->FirstAttribute();
		while (rAttrib)	{
			if (!strcmp(rAttrib->Name(),"InputFileName"))	{
				this->InputFileName = std::string(rAttrib->Value());
				std::cout << "Input file name: " << this->InputFileName << std::endl;
			}
			else if (!strcmp(rAttrib->Name(),"OutputFileName"))	{
				this->OutputFileName = std::string(rAttrib->Value());
				std::cout << "Output file name: " << this->OutputFileName << std::endl;
			}
			rAttrib=rAttrib->Next();
		}
		TiXmlElement* param = root->FirstChildElement( "Parameter" );
		while ( param )	{
			// TODO: look for data as well as Parameter tag HERE
			TiXmlAttribute* pAttrib = param->FirstAttribute();
			if (pAttrib){
				if (!strcmp(pAttrib->Name(),"GridSpacing"))	{
					double temp = -1.0;
					if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS) {
						if (!ParseDoubleInput(temp, 15 , 100, "GridSpacing") )	{
							return false;
						}
						GridSpacing = temp;
					}
					std::cout << "GridSpacing:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
				}
				else if (!strcmp(pAttrib->Name(),"StepRatio"))	{
					double temp = -1.0;
					if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)	{
						if (!ParseDoubleInput(temp, 0.2 , 1.0, "StepRatio") )	{
							return false;
						}
						StepRatio = temp;
					}
					std::cout << "StepRatio:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
				}
				else if (!strcmp(pAttrib->Name(),"XYSpacing"))	{
					double temp = -1.0;
					if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)	{
						if (!ParseDoubleInput(temp, 0.0 , 100.0, "XYSpacing") )	{
							return false;
						}
						Spacing[0] = temp;
						Spacing[1] = temp;
					}
					std::cout << "XYspacing:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
				}

				else if (!strcmp(pAttrib->Name(),"ZSpacing"))	{
					double temp = -1.0;
					if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)	{
						if (!ParseDoubleInput(temp, 0.0 , 100.0, "ZSpacing") )	{
							return false;
						}
						Spacing[2] = temp;
					}
					std::cout << "ZSpacing:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
				}

				else if (!strcmp(pAttrib->Name(),"MaxModelAspectRatio"))	{
					double temp = -1.0;
					if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)	{
						if (!ParseDoubleInput(temp, 1.5 , 3.0, "MaxModelAspectRatio") )	{
							return false;
						}
						AspectRatio = temp;
					}
					std::cout << "MaxModelAspectRatio:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
				}
				else if (!strcmp(pAttrib->Name(),"THRESHOLD"))	{
					double temp = -1.0;
					if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS) {
						if (!ParseDoubleInput(temp, 0.1 , 1.0, "THRESHOLD") )	{
							return false;
						}
						THRESHOLD = temp;
					}
					std::cout << "THRESHOLD:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
				}
				else if(!strcmp(pAttrib->Name(),"minContrast"))	{
					double temp = -1.0;
					if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)	{
						if (!ParseDoubleInput(temp, 3.0 , 20.0, "minContrast") )	{
							return false;
						}
						minContrast = temp;
					}
					std::cout << "minContrast:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
				}
				else if(!strcmp(pAttrib->Name(),"MaximumVesselWidth"))	{
					double temp = -1.0;
					if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)	{
						if (!ParseDoubleInput(temp, 10.0 , 100.0, "MaximumVesselWidth") )	{
							return false;
						}
						MaximumVesselWidth = temp/2.0;
					}
					std::cout << "MaximumVesselWidth:  " << 2*MaximumVesselWidth << " ("<< pAttrib->Value()  <<")" << std::endl;
				}
				else if(!strcmp(pAttrib->Name(),"MinimumVesselWidth"))	{
					double temp = -1.0;
					if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)	{
						if (!ParseDoubleInput(temp, 3.0 , 10.0, "MinimumVesselWidth") )	{
							return false;
						}
						MinimumVesselWidth = temp/2.0;
					}
					std::cout << "MinimumVesselWidth:  " << 2*MaximumVesselWidth << " ("<< pAttrib->Value()  <<")" << std::endl;
				}

				else if(!strcmp(pAttrib->Name(),"MinimumVesselLength"))	{
					double temp = -1.0;
					if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)  {
						if (!ParseDoubleInput(temp, 3.0 , 20.0, "MinimumVesselWidth") )	{
							return false;
						}
						MinimumVesselLength = temp;
					}
					std::cout << "MinimumVesselLength:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
				}

				else if(!strcmp(pAttrib->Name(),"StartTHRESH"))	{
					double temp = -1.0;
					if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)	{
						if (!ParseDoubleInput(temp, 0.1 , 1.0, "StartTHRESH") )	{
							return false;
						}
						StartTHRESHOLD = temp;
					}
					std::cout << "StartTHRESH:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
				}

				else if(!strcmp(pAttrib->Name(),"StartIterations"))	{
					int temp = -1;
					if (pAttrib->QueryIntValue(&temp)==TIXML_SUCCESS)
						FitIterations = temp;
					std::cout << "StartIterations:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
				}
				else {
					std::cout << "UNRECOGNIZED TAG:  " << pAttrib->Value() << std::endl;
				}
			}
			param = param->NextSiblingElement();
		}
	}
	return true;
}
*/
bool TraceConfig::ParseDoubleInput(double val, double minV, double maxV, const char* errmsg)	{
	if((val > maxV) || (val < minV))	{
		std::cout << "Input " << errmsg << " exceeds range. Value read" << val << " whereas the range is [min:" << minV << 
			" , max: " << maxV << "]" << std::endl;
		return false;
	}
	return true;
}


bool TraceConfig::LoadParameters(char * ParameterFileName)  {

	TiXmlDocument doc(ParameterFileName);
	if (!doc.LoadFile()) {
		return false;
	}

	TiXmlElement* root = doc.FirstChildElement( "SuperElliposoidTracing" );
	if ( root )		{
		TiXmlNode* node = root->FirstChild();
		while ( node )	{
			// Check if this is NODE
			if (node->Type() != TiXmlNode::ELEMENT) {
				continue;
			}
			//decide data tag or parameter tag	
			if (!strncmp(node->Value(), "Data", 4)) {
				TiXmlAttribute* rAttrib = node->ToElement()->FirstAttribute();
				while (rAttrib)	{
					if (!strcmp(rAttrib->Name(),"InputFileName"))	{
						this->InputFileNames.push_back(std::string(rAttrib->Value()));
						std::cout << "Input file name: " << this->InputFileNames[numDataFiles] << std::endl;
					}
					else if (!strcmp(rAttrib->Name(),"OutputFileName"))	{
						this->OutputFileNames.push_back(std::string(rAttrib->Value()));
						std::cout << "Output file name: " << this->OutputFileNames[numDataFiles] << std::endl;
					}
					rAttrib=rAttrib->Next();
				}
				numDataFiles++;
			}
			else if (!strncmp(node->Value(), "Parameter", 9)) {
				TiXmlAttribute* pAttrib = node->ToElement()->FirstAttribute();
				if (pAttrib){
					if (!strcmp(pAttrib->Name(),"GridSpacing"))	{
						double temp = -1.0;
						if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS) {
							if (!ParseDoubleInput(temp, 15 , 100, "GridSpacing") )	{
								return false;
							}
							GridSpacing = temp;
						}
						std::cout << "GridSpacing:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
					}
					else if (!strcmp(pAttrib->Name(),"StepRatio"))	{
						double temp = -1.0;
						if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)	{
							if (!ParseDoubleInput(temp, 0.2 , 1.0, "StepRatio") )	{
								return false;
							}
							StepRatio = temp;
						}
						std::cout << "StepRatio:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
					}
					else if (!strcmp(pAttrib->Name(),"UseMultiscaleHessianFilter"))	{
						int temp = 0;
						if (pAttrib->QueryIntValue(&temp)==TIXML_SUCCESS)	{
							UseMultiscaleHessianFilter = temp;
						}
						std::cout << "UseMultiscaleHessianFilter:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
					}
					else if (!strcmp(pAttrib->Name(),"XYSpacing"))	{
						double temp = -1.0;
						if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)	{
							if (!ParseDoubleInput(temp, 0.0 , 100.0, "XYSpacing") )	{
								return false;
							}
							Spacing[0] = temp;
							Spacing[1] = temp;
						}
						std::cout << "XYspacing:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
					}

					else if (!strcmp(pAttrib->Name(),"ZSpacing"))	{
						double temp = -1.0;
						if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)	{
							if (!ParseDoubleInput(temp, 0.0 , 100.0, "ZSpacing") )	{
								return false;
							}
							Spacing[2] = temp;
						}
						std::cout << "ZSpacing:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
					}

					else if (!strcmp(pAttrib->Name(),"MaxModelAspectRatio"))	{
						double temp = -1.0;
						if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)	{
							if (!ParseDoubleInput(temp, 1.35 , 3.0, "MaxModelAspectRatio") )	{
								return false;
							}
							AspectRatio = temp;
						}
						std::cout << "MaxModelAspectRatio:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
					}
					else if (!strcmp(pAttrib->Name(),"THRESHOLD"))	{
						double temp = -1.0;
						if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS) {
							if (!ParseDoubleInput(temp, 0.1 , 1.0, "THRESHOLD") )	{
								return false;
							}
							THRESHOLD = temp;
						}
						std::cout << "THRESHOLD:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
					}
					else if(!strcmp(pAttrib->Name(),"minContrast"))	{
						double temp = -1.0;
						if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)	{
							if (!ParseDoubleInput(temp, 3.0 , 20.0, "minContrast") )	{
								return false;
							}
							minContrast = temp;
						}
						std::cout << "minContrast:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
					}
					else if(!strcmp(pAttrib->Name(),"MaximumVesselWidth"))	{
						double temp = -1.0;
						if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)	{
							if (!ParseDoubleInput(temp, 10.0 , 100.0, "MaximumVesselWidth") )	{
								return false;
							}
							MaximumVesselWidth = temp/2.0;
						}
						std::cout << "MaximumVesselWidth:  " << 2*MaximumVesselWidth << " ("<< pAttrib->Value()  <<")" << std::endl;
					}
					else if(!strcmp(pAttrib->Name(),"MinimumVesselWidth"))	{
						double temp = -1.0;
						if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)	{
							if (!ParseDoubleInput(temp, 3.0 , 10.0, "MinimumVesselWidth") )	{
								return false;
							}
							MinimumVesselWidth = temp/2.0;
						}
						std::cout << "MinimumVesselWidth:  " << 2*MaximumVesselWidth << " ("<< pAttrib->Value()  <<")" << std::endl;
					}

					else if(!strcmp(pAttrib->Name(),"MinimumVesselLength"))	{
						double temp = -1.0;
						if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)  {
							if (!ParseDoubleInput(temp, 3.0 , 20.0, "MinimumVesselWidth") )	{
								return false;
							}
							MinimumVesselLength = temp;
						}
						std::cout << "MinimumVesselLength:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
					}

					else if(!strcmp(pAttrib->Name(),"StartTHRESH"))	{
						double temp = -1.0;
						if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)	{
							if (!ParseDoubleInput(temp, 0.1 , 1.0, "StartTHRESH") )	{
								return false;
							}
							StartTHRESHOLD = temp;
						}
						std::cout << "StartTHRESH:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
					}

					else if(!strcmp(pAttrib->Name(),"StartIterations"))	{
						int temp = -1;
						if (pAttrib->QueryIntValue(&temp)==TIXML_SUCCESS)
							FitIterations = temp;
						std::cout << "StartIterations:  " << temp << " ("<< pAttrib->Value()  <<")" << std::endl;
					}
					else {
						std::cout << "UNRECOGNIZED TAG:  " << pAttrib->Value() << std::endl;
					}
				}
			}
			node = node->NextSibling();
		}
	}
	return true;
}




