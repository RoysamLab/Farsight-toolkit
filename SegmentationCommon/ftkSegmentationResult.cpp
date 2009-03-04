/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
#include "ftkSegmentationResult.h"
#include <iostream>
#include <fstream>
#include <iomanip>

namespace ftk 
{

//SegmentationResult Constructor
SegmentationResult::SegmentationResult(string projpath, string projname)
{
	projectPath = projpath;			//Location of all files for this project
	projectName = projname;			//A string that is common to all the beginning of all filenames for this project

	programName.clear();
	dataFilenames.clear();		
	resultFilenames.clear(); 
	myParameters.clear();
	myObjects.clear();
	featureNames.clear();
	maxID = 0;
	errorMessage.clear();
}

Object* SegmentationResult::GetObjectPtr(int id)
{
	int index = IdToIndexMap[id];
	return &myObjects.at(index);
}

bool SegmentationResult::RestoreFromXML()
{
	programName.clear();
	dataFilenames.clear();		
	resultFilenames.clear(); 
	myParameters.clear();
	myObjects.clear();
	featureNames.clear();

	TiXmlDocument doc;
	if ( !doc.LoadFile( XmlFilename(1).c_str() ) )
	{
		errorMessage = "Unable to load XML File";
		return 0;
	}

	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "SegmentationResult" ) != 0 )
	{
		errorMessage = "Incorrect XML root Element: ";
		errorMessage.append(rootElement->Value());
		return 0;
	}
		
	programName = rootElement->Attribute("program");

	//Parents we know of: datafilename,resultfilename,object,parameter
	TiXmlElement* parentElement = rootElement->FirstChildElement();
	while (parentElement)
	{
		const char * parent = parentElement->Value();

		if ( strcmp( parent, "datafile" ) == 0 )
		{
			dataFilenames.push_back( parentElement->GetText() );
		}
		else if ( strcmp( parent, "resultfile" ) == 0 )
		{
			resultFilenames.push_back( parentElement->GetText() );
		}
		else if ( strcmp( parent, "parameter" ) == 0 )
		{
			Parameter p;
			p.name = parentElement->Attribute("name");
			p.value = atoi( parentElement->GetText() );
			myParameters.push_back( p );
		}
		else if ( strcmp( parent, "object" ) == 0 )
		{
			Object o = parseObject(parentElement);
			if ( o.GetType() != "null" )
			{
				myObjects.push_back( o );
				int id = o.GetId();
				IdToIndexMap[id] = (int)myObjects.size() - 1;
				if(id > maxID) maxID = id;
			}
		}
		else
		{
			errorMessage = "Unrecognized parent element: ";
			errorMessage.append(parent);
			return 0;
		}
		parentElement = parentElement->NextSiblingElement();
	} // end while(parentElement)

	//doc.close();

	return 1;
}

Object SegmentationResult::parseObject(TiXmlElement *objectElement)
{
	if ( strcmp( objectElement->Value(), "object" ) != 0 )
	{
		errorMessage = "This is not an object: ";
		errorMessage.append(objectElement->Value());
		return Object("null");
	}

	Object object( objectElement->Attribute("type") );

	TiXmlElement *member = objectElement->FirstChildElement();
	while(member)
	{
		const char* memberName = member->Value();
		if ( strcmp( memberName, "id" ) == 0 )
		{
			object.SetId( atoi( member->GetText() ) );
		}
		else if ( strcmp( memberName, "validity" ) == 0 )
		{
			object.SetValidity( atoi( member->GetText() ) );
		}
		else if ( strcmp( memberName, "duplicated" ) == 0 )
		{
			object.SetDuplicated( atoi( member->GetText() ) );
		}
		else if ( strcmp( memberName, "class" ) == 0 )
		{
			object.SetClass( atoi( member->GetText() ) );
		}
		else if ( strcmp( memberName, "center") == 0 )
		{
			object.AddCenter( parseCenter(member) );
		}
		else if ( strcmp( memberName, "bound") == 0 )
		{
			object.AddBound( parseBound(member) );
		}
		else if ( strcmp( memberName, "features") == 0 )
		{
			object.SetFeatures( parseFeatures(member) );
		}
		else if ( strcmp( memberName, "EditHistory") == 0 )
		{
			TiXmlElement *record = member->FirstChildElement();
			while(record)
			{
				if ( strcmp( record->Value(), "record" ) == 0 )
				{
					Object::EditRecord r;
					r.date = record->Attribute("date");
					r.description = record->GetText();
					object.AddEditRecord( r );
				}
				record = record->NextSiblingElement();
			}
		}
		else
		{
			errorMessage = "Unrecognized object member: ";
			errorMessage.append(memberName);
			return Object("null");
		}
		member = member->NextSiblingElement();
	}

	return object;
}

Object::Point SegmentationResult::parseCenter(TiXmlElement *centerElement)
{
	Object::Point center;

	TiXmlElement *coord = centerElement->FirstChildElement();
	while(coord)
	{
		const char* coordName = coord->Value();
		if ( strcmp( coordName, "x" ) == 0 )
		{
			center.x = atoi( coord->GetText() );
		}
		else if ( strcmp( coordName, "y" ) == 0 )
		{
			center.y = atoi( coord->GetText() );
		}
		else if ( strcmp( coordName, "z" ) == 0 )
		{
			center.z = atoi( coord->GetText() );
		}
		else if ( strcmp( coordName, "t" ) == 0 )
		{
			center.t = atoi( coord->GetText() );
		}
		coord = coord->NextSiblingElement();
	}
	return center;
}

Object::Box SegmentationResult::parseBound(TiXmlElement *boundElement)
{
	Object::Box bound;

	TiXmlElement *coord = boundElement->FirstChildElement();
	while(coord)
	{
		const char* coordName = coord->Value();
		if ( strcmp( coordName, "xmin" ) == 0 )
		{
			bound.min.x = atoi( coord->GetText() );
		}
		else if ( strcmp( coordName, "ymin" ) == 0 )
		{
			bound.min.y = atoi( coord->GetText() );
		}
		else if ( strcmp( coordName, "zmin" ) == 0 )
		{
			bound.min.z = atoi( coord->GetText() );
		}
		else if ( strcmp( coordName, "tmin" ) == 0 )
		{
			bound.min.t = atoi( coord->GetText() );
		}
		else if ( strcmp( coordName, "xmax" ) == 0 )
		{
			bound.max.x = atoi( coord->GetText() );
		}
		else if ( strcmp( coordName, "ymax" ) == 0 )
		{
			bound.max.y = atoi( coord->GetText() );
		}
		else if ( strcmp( coordName, "zmax" ) == 0 )
		{
			bound.max.z = atoi( coord->GetText() );
		}
		else if ( strcmp( coordName, "tmax" ) == 0 )
		{
			bound.max.t = atoi( coord->GetText() );
		}
		coord = coord->NextSiblingElement();
	}
	return bound;
}

vector< float > SegmentationResult::parseFeatures(TiXmlElement *featureElement)
{
	vector< float > tempFeatures(0);
	vector< string > tempNames(0); 

	//Extract all of the features for this object
	TiXmlElement *feature = featureElement->FirstChildElement();
	while(feature)
	{
		tempNames.push_back( feature->Value() );
		tempFeatures.push_back( float( atof( feature->GetText() ) ) );
		
		feature = feature->NextSiblingElement();
	}

	//Now test for existing featureNames
	if (featureNames.size() == 0) //Meaning hasn't been set yet
	{
		featureNames = tempNames; //Set my feature names to these names
		return tempFeatures;
	}

	//Have already set the names, so need to make sure we return the features in the correct order
	if( featureNames.size() != tempNames.size() )
	{
		errorMessage = "Features do not match";
		return vector< float >(0);
	}
		
	vector<float> retFeatures( tempFeatures.size() );
	unsigned int numMatches = 0;
	for (unsigned int i=0; i<tempFeatures.size(); ++i)
	{
		string featName = tempNames[i];
		for (unsigned int j=0; j<featureNames.size(); ++j)
		{
			if( featName == featureNames[j] )
			{
				retFeatures[j] = tempFeatures[i];
				++numMatches;
			}
		}
	}

	if ( numMatches != tempFeatures.size() )
	{
		errorMessage = "Features do not match";
		return vector< float >(0);
	}

	return retFeatures;
}


bool SegmentationResult::WriteToXML()
{
	TiXmlDocument doc;   
 
	TiXmlElement * root = new TiXmlElement( "SegmentationResult" );  
	doc.LinkEndChild( root );  
	root->SetAttribute("program", "Yousef_Nucleus_Seg");

	TiXmlComment * comment = new TiXmlComment();
	comment->SetValue(" Segmentation Results/Parameters/Features/Edits " );  
	root->LinkEndChild( comment );  
 
	//Attach data filenames:
	for (unsigned int fnum = 0; fnum < dataFilenames.size(); ++fnum)
	{
		TiXmlElement *element = new TiXmlElement("datafile");
		element->LinkEndChild( new TiXmlText( dataFilenames[fnum].c_str() ) );
		root->LinkEndChild(element);
	}

	//Attach result filenames:
	for (unsigned int fnum = 0; fnum < resultFilenames.size(); ++fnum)
	{
		TiXmlElement *element = new TiXmlElement("resultfile");
		element->LinkEndChild( new TiXmlText( resultFilenames[fnum].c_str() ) );
		root->LinkEndChild(element);
	}

	//Attach parameters
	for (unsigned int pnum = 0; pnum < myParameters.size(); ++pnum)
	{
		TiXmlElement *element = new TiXmlElement("parameter");
		element->SetAttribute( "name", myParameters[pnum].name );
		element->LinkEndChild( new TiXmlText( NumToString(myParameters[pnum].value) ) );
		root->LinkEndChild(element);

	}

	//Attach Objects
	for (unsigned int onum=0; onum<myObjects.size(); ++onum)
	{
		root->LinkEndChild( GetObjectElement(myObjects[onum]) );
	}

	doc.SaveFile( XmlFilename(1).c_str() );

	return true;
}

TiXmlElement* SegmentationResult::GetObjectElement(Object object)
{
	TiXmlElement *objectElement = new TiXmlElement("object");
	objectElement->SetAttribute("type", object.GetType());

	if ( object.GetId() != -1 )
	{
		TiXmlElement *idElement = new TiXmlElement("id");
		idElement->LinkEndChild( new TiXmlText( NumToString(object.GetId()) ) );
		objectElement->LinkEndChild(idElement);
	}
	if ( object.GetValidity() != -1 )
	{
		TiXmlElement *vElement = new TiXmlElement("validity");
		vElement->LinkEndChild( new TiXmlText( NumToString(object.GetValidity()) ) );
		objectElement->LinkEndChild(vElement);
	}
	if ( object.GetDuplicated() != -1 )
	{
		TiXmlElement *dElement = new TiXmlElement("duplicated");
		dElement->LinkEndChild( new TiXmlText( NumToString(object.GetDuplicated()) ) );
		objectElement->LinkEndChild(dElement);
	}
	if ( object.GetClass() != -1 )
	{
		TiXmlElement *cElement = new TiXmlElement("class");
		cElement->LinkEndChild( new TiXmlText( NumToString(object.GetClass()) ) );
		objectElement->LinkEndChild(cElement);
	}
	vector<Object::Point> centers = object.GetCenters();
	for (unsigned int c=0; c<centers.size(); ++c)
	{
		TiXmlElement *oElement = new TiXmlElement("center");
		TiXmlElement *xElement = new TiXmlElement("x");
		xElement->LinkEndChild( new TiXmlText( NumToString(centers[c].x) ) );
		oElement->LinkEndChild(xElement);
		TiXmlElement *yElement = new TiXmlElement("y");
		yElement->LinkEndChild( new TiXmlText( NumToString(centers[c].y) ) );
		oElement->LinkEndChild(yElement);
		TiXmlElement *zElement = new TiXmlElement("z");
		zElement->LinkEndChild( new TiXmlText( NumToString(centers[c].z) ) );
		oElement->LinkEndChild(zElement);
		TiXmlElement *tElement = new TiXmlElement("t");
		tElement->LinkEndChild( new TiXmlText( NumToString(centers[c].t) ) );
		oElement->LinkEndChild(tElement);
		objectElement->LinkEndChild(oElement);
	}
	vector<Object::Box> bounds = object.GetBounds();
	for (unsigned int b=0; b<bounds.size(); ++b)
	{
		TiXmlElement *bElement = new TiXmlElement("bound");
		TiXmlElement *xminElement = new TiXmlElement("xmin");
		xminElement->LinkEndChild( new TiXmlText( NumToString(bounds[b].min.x) ) );
		bElement->LinkEndChild(xminElement);
		TiXmlElement *yminElement = new TiXmlElement("ymin");
		yminElement->LinkEndChild( new TiXmlText( NumToString(bounds[b].min.y) ) );
		bElement->LinkEndChild(yminElement);
		TiXmlElement *zminElement = new TiXmlElement("zmin");
		zminElement->LinkEndChild( new TiXmlText( NumToString(bounds[b].min.z) ) );
		bElement->LinkEndChild(zminElement);
		TiXmlElement *tminElement = new TiXmlElement("tmin");
		tminElement->LinkEndChild( new TiXmlText( NumToString(bounds[b].min.t) ) );
		bElement->LinkEndChild(tminElement);
		TiXmlElement *xmaxElement = new TiXmlElement("xmax");
		xmaxElement->LinkEndChild( new TiXmlText( NumToString(bounds[b].max.x) ) );
		bElement->LinkEndChild(xmaxElement);
		TiXmlElement *ymaxElement = new TiXmlElement("ymax");
		ymaxElement->LinkEndChild( new TiXmlText( NumToString(bounds[b].max.y) ) );
		bElement->LinkEndChild(ymaxElement);
		TiXmlElement *zmaxElement = new TiXmlElement("zmax");
		zmaxElement->LinkEndChild( new TiXmlText( NumToString(bounds[b].max.z) ) );
		bElement->LinkEndChild(zmaxElement);
		TiXmlElement *tmaxElement = new TiXmlElement("tmax");
		tmaxElement->LinkEndChild( new TiXmlText( NumToString(bounds[b].max.t) ) );
		bElement->LinkEndChild(tmaxElement);
		objectElement->LinkEndChild(bElement);
	}
	if(featureNames.size() > 0)
	{
		TiXmlElement *fsElement = new TiXmlElement("features");
		vector<float> features = object.GetFeatures();
		for(unsigned int f=0; f<features.size(); ++f)
		{
			TiXmlElement *fElement = new TiXmlElement( featureNames[f].c_str() );
			fElement->LinkEndChild( new TiXmlText( NumToString(features[f]) ) );
			fsElement->LinkEndChild(fElement);
		}
		objectElement->LinkEndChild(fsElement);
	}
	vector<Object::EditRecord> records = object.getHistory();
	if(records.size() > 0)
	{
		TiXmlElement *histElement = new TiXmlElement("EditHistory");
		for (unsigned int r=0; r<records.size(); ++r)
		{
			TiXmlElement *recElement = new TiXmlElement("record");
			recElement->SetAttribute("date",records[r].date);
			recElement->LinkEndChild( new TiXmlText( records[r].description ) );
			histElement->LinkEndChild(recElement);
		}
		objectElement->LinkEndChild(histElement);
	}
	return objectElement;
}

bool SegmentationResult::WriteToMETA()
{
	//This function writes the features to a text file that can be read be MetaNeural program
	string fname = projectName;
	fname.append("_META.txt");

	ofstream outFile; 
	outFile.open(PrependProjectPath(fname).c_str(), ios::out | ios::trunc );
	if ( !outFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << outFile << std::endl;
		return false;
	}
	//Now write out the features
	for(unsigned int obj = 0; obj < myObjects.size(); ++obj)
	{
		if( myObjects.at(obj).GetValidity() == false )
			continue;

		vector<float> feats = myObjects.at(obj).GetFeatures();
		for(unsigned int f = 0; f < feats.size(); ++f)
		{
			outFile << NumToString(feats.at(f)) << "\t";
		}
		outFile << (int)myObjects.at(obj).GetClass() << "\t";
		outFile << myObjects.at(obj).GetId() << endl;
	}
	outFile.close();

	//Now print header file:
	fname.clear();
	fname = projectName;
	fname.append("_META_Header.txt");
	outFile.open(PrependProjectPath(fname).c_str(), std::ios::out | std::ios::trunc );
	if ( !outFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << outFile << std::endl;
		return false;
	}

	for(int i=0; i<featureNames.size(); ++i)
	{
		outFile << featureNames.at(i) << "\n";
	}

	outFile << "ID" << "\n";
	outFile << "CLASS" << std::endl;

	outFile.close();

	return true;
}

string SegmentationResult::XmlFilename(bool full = false)
{
	string f = "";
	if(full) 
		f = projectPath;
	f.append("/");
	f.append(projectName);
	f.append(".xml");
	return f;
}

string SegmentationResult::PrependProjectPath(string file)
{
	string returnName = projectPath;
	returnName.append("/");
	returnName.append(file);
	return returnName;
}

string SegmentationResult::NumToString(double d)
{
	stringstream out;
	out << setprecision(2) << fixed << d;	//Default is to use 2 decimal places
	return out.str();
}

string SegmentationResult::NumToString(int i)
{
	stringstream out;
	out << i ;	 
	return out.str();
}

string SegmentationResult::NumToString(double d, int p)
{
	stringstream out;
	out << setprecision(p) << fixed << d;	
	return out.str();
}

} //end namespace ftk
