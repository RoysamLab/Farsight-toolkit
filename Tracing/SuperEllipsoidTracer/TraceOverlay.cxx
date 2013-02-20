
# include <iostream>
#include "itkCommand.h"
#include "itkVector.h"
#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkRGBPixel.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "TraceNode.h"
#include "tinyxml.h"

typedef itk::Image< unsigned short, 3 > InputImageType;
typedef itk::Image< unsigned char, 3 > ImageType;

bool ReadNodeXMLFile(const char*, std::vector<TraceNode*>& );

int main(const int argc, char** argv)	{
	if (argc != 4)	{
		std::cout << "Usage: "<<argv[0] << " [ImageFile] [TraceFile] [OverlayFile]" <<std::endl;
		return EXIT_FAILURE;
	}
	
	std::cout << "Reading " << argv[1] << " ..." ;
    typedef itk::ImageFileReader< InputImageType > ReaderType;
    ReaderType::Pointer reader  = ReaderType::New();
    reader->SetFileName( argv[1] );

	typedef itk::RescaleIntensityImageFilter<InputImageType, ImageType> RescalerType;
	RescalerType::Pointer rescaler = RescalerType::New();
	rescaler->SetOutputMaximum(255);
	rescaler->SetOutputMinimum(0);
	rescaler->SetInput(reader->GetOutput());
    ImageType::Pointer im2 = rescaler->GetOutput();
	try {
		im2->Update();
	}
	catch (itk::ExceptionObject &e)	{
		std::cerr << "Could not read file "<< std::endl << e << std::endl;
	}

	std::cout << "done! " <<std::endl << "Reading " << argv[2] << " ...";

	std::vector<TraceNode*> NodeContainer;
	if (ReadNodeXMLFile(argv[2], NodeContainer))	{
		std::cout << "done! " <<std::endl << "Overlaying trace .... wait!! " << std::endl;
		

		// create an empty image
		ImageType::Pointer binImg = ImageType::New();
		binImg->SetRegions(im2->GetBufferedRegion());
		binImg->Allocate();
		binImg->FillBuffer(0);
		
		// create the nodes
		std::vector<TraceNode*>::iterator itr;
		std::map<long, itk::Vector<double,3> > IDmap;
		for(itr = NodeContainer.begin(); itr != NodeContainer.end(); ++itr)	{
			IDmap[(*itr)->ID] = (*itr)->loc;
		}

		for(itr = NodeContainer.begin(); itr != NodeContainer.end(); ++itr)	{
			itk::Vector<double,3> n1 = (*itr)->loc;
			for (int i = 0; i < (*itr)->nbrID.size() ; i++) {
				itk::Vector<double,3> n2 = IDmap[(*itr)->nbrID[i]];
				itk::Vector<double,3> d = n2 - n1;
				itk::Index<3> n;
				double dn = d.GetNorm();
				d.Normalize();
				double k = 0.0f;
				while(k < dn)	{
					for (int o1=-1; o1<1; o1++)	{
						for (int o2=-1; o2<1; o2++)	{
							for (int o3=-1; o3<1; o3++)	{
								n[0] = static_cast<long>(n1[0] + k*d[0]) + o1;
								n[1] = static_cast<long>(n1[1] + k*d[1]) + o2;
								n[2] = static_cast<long>(n1[2] + k*d[2]) + o3;
								if (binImg->GetBufferedRegion().IsInside(n))	{
									binImg->SetPixel(n,255);
								}
							}
						}
					}
					k += 1.0;
				}
				n[0] = static_cast<long>(n2[0]);
				n[1] = static_cast<long>(n2[1]);
				n[2] = static_cast<long>(n2[2]);
				binImg->SetPixel(n,255);
			}
		}

		std::string OverlayFileName(argv[3]);
		OverlayFileName.replace(OverlayFileName.length()-4,8,"_bin.tif");
		itk::ImageFileWriter<ImageType>::Pointer w = itk::ImageFileWriter<ImageType>::New();
		w->SetFileName(OverlayFileName);
		w->SetInput(binImg);
		w->Update();

		//for each tile make an RGB image
		typedef itk::RGBPixel<unsigned char> RGBType;
		typedef itk::Image<RGBType,2> RGBImageType;
		itk::ImageFileWriter<RGBImageType>::Pointer w2 = itk::ImageFileWriter<RGBImageType>::New();

		itk::Size<3> sz = im2->GetBufferedRegion().GetSize();
		itk::Size<2> sz2 = {{sz[0], sz[1]}};
		RGBImageType::Pointer rgb = RGBImageType::New();
		rgb->SetRegions(sz2);
		rgb->Allocate();
		for (unsigned int z=0; z<sz[2]; z++) {
			for (unsigned int y=0; y<sz[1]; y++) {
				for (unsigned int x=0; x<sz[0]; x++) {
					itk::Index<3> ndx3 = {x, y, z};
					itk::Index<2> ndx2 = {x, y};
					if (binImg->GetPixel(ndx3) == 0)	{
						ImageType::PixelType val = im2->GetPixel(ndx3);
						RGBImageType::PixelType rval(val);			
						rgb->SetPixel(ndx2, rval);
					}
					else {
						RGBImageType::PixelType rval;
						rval.Fill(0);
						rval.SetGreen(255);
						rgb->SetPixel(ndx2, rval);
					}
				}
			}
			
			char fname[8];
			sprintf(fname,"%03d.png",z);
			OverlayFileName.replace(OverlayFileName.length()-7,7,fname);
			w2->SetFileName(OverlayFileName);
			w2->SetInput(rgb);
			w2->Update();
		}
	}
}


bool ReadNodeXMLFile(const char* xmlfname, std::vector<TraceNode*>& NodeContainer) {
	NodeContainer.reserve(10000);
	TiXmlDocument doc(xmlfname);
	if (!doc.LoadFile()) {
		return false;
	}

	//scan each Superellipse
	TiXmlNode* xmlnode; 

	for ( xmlnode = doc.FirstChild(); xmlnode != 0; xmlnode = xmlnode->NextSibling()) 	{

		//verify if the xmlnode is a type element
		if (xmlnode->Type()!=TiXmlNode::ELEMENT)	{
			continue;
		}

		//verify if the xmlnode is a superellipse, if not 
		if (strcmp(xmlnode->Value(),"Superellipse"))	{
			continue;
		}

		TraceNode *n = new TraceNode();
		TiXmlAttribute* pAttrib = xmlnode->ToElement()->FirstAttribute();
		while (pAttrib)	{
			if (!strcmp(pAttrib->Name(),"ID"))	{
				int temp = -1;
				if (pAttrib->QueryIntValue(&temp)==TIXML_SUCCESS)
					n->ID = temp;
			}
			else if (!strcmp(pAttrib->Name(),"TraceID"))	{
				int temp = -1;
				if (pAttrib->QueryIntValue(&temp)==TIXML_SUCCESS)
					n->TraceID = temp;
			}
			
			else if (!strcmp(pAttrib->Name(),"x"))	{
				double temp = -1.0;
				if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)
					n->loc[0] = temp;
			}

			else if (!strcmp(pAttrib->Name(),"y"))	{
				double temp = -1.0;
				if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)
					n->loc[1] = temp;
			}
			
			else if (!strcmp(pAttrib->Name(),"z"))	{
				double temp = -1.0;
				if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)
					n->loc[2] = temp;
			}
			
			pAttrib=pAttrib->Next();
		}

		
		TiXmlNode* nbr; 
		for ( nbr = xmlnode->FirstChild(); nbr != 0; nbr = nbr->NextSibling())		{
			TiXmlAttribute* nAttr = nbr->ToElement()->FirstAttribute();
			if (!strcmp(nAttr->Name(),"ID"))	{
				int temp = -1;
				if (nAttr->QueryIntValue(&temp)==TIXML_SUCCESS)
					n->nbrID.push_back(temp);
			}
		}
		
		//store in container
		NodeContainer.push_back(n);
	}
	return true;

}

