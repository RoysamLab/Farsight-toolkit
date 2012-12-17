#include "RenderObject.h"
#include "TraceView3D.h"
#include <stdio.h>
#include <stdlib.h>


RenderObject::RenderObject(QObject *parent)
: QObject(parent)
{
	//this->computedImageFeatureThreshold.clear();
	//this->QVTK = 0;
	//m_abort = false;
}

RenderObject::~RenderObject()
{
	//mutex.lock();
	//m_abort = true;
	//mutex.unlock();

	//wait();
	//if(this->QVTK)
	//{
	//	delete this->QVTK;
	//}

	//while(!this->computedImageFeatureThreshold.empty())
	//   {
	//   delete this->computedImageFeatureThreshold.back();
	//   this->computedImageFeatureThreshold.pop_back();
	//   }

}

void RenderObject::runMNTTracer(){

	//for(int i=0;i<100;i++){
	//	Sleep(100);
	//}

	std::cout<<"in RunMNT Tracer"<<std::endl;
	char buffer [33];
	itoa (cost_threshold,buffer,10);
	clock_t start_time = clock();
	std::string InputFilename = std::string(imageFile);
	std::string SWCFilename = InputFilename;
	SWCFilename.erase(SWCFilename.length()-4,SWCFilename.length());
	//SWCFilename.append(buffer);
	SWCFilename.append("_ANT.swc");
	try 
	{	
		MultipleNeuronTracer * MNT = new MultipleNeuronTracer();
		clock_t LoadCurvImage_start_time = clock();
		/*MNT->LoadCurvImage(InputFilename, 1);*/
		itk::CastImageFilter<ImageType3D,UnsignedCharImageType >::Pointer caster2 = itk::CastImageFilter<ImageType3D,UnsignedCharImageType>::New();
		caster2->SetInput(inputImage);
		MNT->intensity_threshold = _featureThreshold->getintensity_threshold();//intensity_threshold;
		MNT->contrast_threshold = _featureThreshold->getContrastThreshold();//contrast_threshold;
		MNT->device = _featureThreshold->getdevice();//this->device;
		MNT->offshoot = _featureThreshold->getoffshoot();//this->offshoot;
		MNT->debris_threshold = _featureThreshold->getDebrisThreshold();//this->debris_threshold;
		MNT->LoadCurvImage_1(inputImage, 1);
		std::cout << "LoadCurvImage took: " << (clock() - LoadCurvImage_start_time)/(float) CLOCKS_PER_SEC << std::endl;
		clock_t ReadStartPoints_start_time = clock();
		//
		std::vector< itk::Index<3> > centroid_list = getCentroidList();		
		std::vector< itk::Index<3> > soma_Table = getSomaTable(centroid_list, startx, starty, startz, 600,600,300,widthx,widthy,widthz);
		MNT->ReadStartPoints_1(soma_Table, 0);
		//MNT->ReadStartPoints(seedsFile,1);
		std::cout << "ReadStartPoints took: " << (clock() - ReadStartPoints_start_time)/(float) CLOCKS_PER_SEC << std::endl;
		clock_t SetCostThreshold_start_time = clock();
		std::cout<<"cost_threshold before"<<cost_threshold<<std::endl;
		MNT->SetCostThreshold(cost_threshold);
		clock_t LoadSomaImage_start_time = clock();
		MNT->LoadSomaImage_1(somaImage);
		std::cout << "LoadSomaImage took: " << (clock() - LoadSomaImage_start_time)/(float) CLOCKS_PER_SEC << std::endl;

		clock_t RunTracing_start_time = clock();
		bool gvfTracing = _featureThreshold->isGVFTracer();
		if(gvfTracing)
		{
			MNT->RunGVFTracing(false);
		}else{
			MNT->RunTracing();
		}
		
		MNT->RemoveSoma( somaImage );
		std::cout << "RunTracing took: " << (clock() - RunTracing_start_time)/(float) CLOCKS_PER_SEC << std::endl;
		
		vtkSmartPointer< vtkTable > swcTable = MNT->GetSWCTable(0);
		_featureThreshold->saveTraces(swcTable);
		std::stringstream ssx, ssy, ssz;
		ssx << startx; ssy << starty; ssz << startz;
		std::string outPath = _featureThreshold->getDir();
		std::string swcFilename = outPath + "/Trace_" + ssx.str() + "_" + ssy.str() + "_" + ssz.str() + "_ANT.swc";
		_featureThreshold->setSWCFilePath(swcFilename);
		std::cout<<"swc file name"<<swcFilename<<std::endl;
		
		
		int xTile= 600;
		int yTile= 600;
		int zTile= 300;
		int x = (xTile/2<startx) ? xTile/2:startx; 
		int y = (yTile/2<starty) ? yTile/2:starty; 
		int z = (zTile/2<startz) ? zTile/2:startz; 
		std::cout<<"	"<<x<<"	"<<y<<"	"<<z<<std::endl;
		this->WriteCenterTrace(swcTable,x,y,z,swcFilename);
		//MNT->WriteSWCFile(swcFilename, 1);
		delete MNT;
	} catch (std::exception& e)
	{
		std::cout << e.what() << std::endl;
	}


	// emit the signal once the tracing is over and this is the signal to kill the thread and merge the traces
	emit finished();
	emit finished(_featureThreshold);
	return;

}
std::vector< itk::Index<3> > RenderObject::getCentroidList()
{	
	std::cout << std::endl << seedsFile;
	vtkSmartPointer<vtkTable> somaCentroidsTable = ftk::LoadTable(seedsFile);

	std::vector< itk::Index<3> > centroid_list;
	for(int r=0; r<(int)somaCentroidsTable->GetNumberOfRows(); ++r)
	{
		int cx = somaCentroidsTable->GetValue(r, 0).ToInt();
		int cy = somaCentroidsTable->GetValue(r, 1).ToInt();
		int cz = somaCentroidsTable->GetValue(r, 2).ToInt();
		
		itk::Index<3> cen;
		cen[0] = cx; cen[1] = cy; cen[2] = cz; 
		centroid_list.push_back(cen);
	}
// 	std::cout << std::endl << "-----------------END";
// 	std::cout << std::endl << "-----------------END";
	return centroid_list;
}


std::vector< itk::Index<3> > RenderObject::getSomaTable( std::vector< itk::Index<3> > centroid_list, int x, int y, int z,int xTile, int yTile, int zTile,int xSize,int ySize, int zSize)
{
	itk::Index<3> centroid;
	centroid[0] = ((x - xTile/2)>0) ? xTile/2:x; 
	centroid[1] = ((y - yTile/2)>0) ? yTile/2:y;
	centroid[2] = ((z - zTile/2)>0) ? zTile/2:z;
	
	itk::Index<3> start;
	start[0] = ((x - xTile/2)>0) ? (x - xTile/2):0;
	start[1] = ((y - yTile/2)>0) ? (y - yTile/2):0;
	start[2] = ((z - zTile/2)>0) ? (z - zTile/2):0;

	itk::Size<3> size;
	size[0] = ((x+xTile/2)<xSize) ? xTile : (xTile/2+xSize-x-1); 
	size[1] = ((y+yTile/2)<ySize) ? yTile : (yTile/2+ySize-y-1);
	size[2] = ((z+zTile/2)<zSize) ? zTile : (zTile/2+zSize-z-1);
	
	std::vector< itk::Index<3> > soma_Table;      
	for(int ctr =0; ctr<centroid_list.size() ; ++ctr)
	{
		itk::Index<3> cen = centroid_list[ctr];
		if( (cen[0]>=start[0]) && (cen[0]<(start[0]+size[0])) && (cen[1]>=start[1]) && (cen[1]<(start[1]+size[1])) && (cen[2]>=start[2]) && (cen[2]<(start[2]+size[2])) )
		{
			itk::Index<3> centroid2;
			centroid2[0] = centroid[0] + cen[0] - x;
			centroid2[1] = centroid[1] + cen[1] - y;
			centroid2[2] = centroid[2] + cen[2] - z;
			
			std::cout<<"selected cen	"<<centroid2[0]<<" "<<centroid2[1]<<" "<<centroid2[2]<<std::endl;
			soma_Table.push_back(centroid2);
		}
	}
	return soma_Table;
}




void RenderObject::WriteCenterTrace(vtkSmartPointer< vtkTable > swcNodes, int x, int y, int z, std::string filename)
{
	std::cout << "Writing SWCImage file " << filename << " with " << swcNodes->GetNumberOfRows() << " nodes...";

	std::vector<int> soma_ids;
	std::vector<int> del_ids;

	for(int r=0; r<(int)swcNodes->GetNumberOfRows(); ++r)
	{
		if(swcNodes->GetValue(r,6).ToInt() == -1)
			soma_ids.push_back(swcNodes->GetValue(r,0).ToInt());
		else
			break;
	}

	for(int i=0; i<soma_ids.size(); ++i)
	{
		if( (swcNodes->GetValue(i,2).ToInt() != x) || (swcNodes->GetValue(i,3).ToInt() != y) || (swcNodes->GetValue(i,4).ToInt() != z) )//if( ((int)swcNodes[i][2] + x > 861) && ((int)swcNodes[i][2] + x < 3861) && ((int)swcNodes[i][3] + y < 8610) )
		{
			del_ids.push_back(soma_ids[i]);
		}
	}

	for(int r=0 ; r<(int)swcNodes->GetNumberOfRows(); ++r)
	{
		if(r+1 > (int)swcNodes->GetNumberOfRows()) break;
		std::vector<int>::iterator posn1 = find(del_ids.begin(), del_ids.end(), swcNodes->GetValue(r,6).ToInt());
		if(posn1 != del_ids.end())
		{
			del_ids.push_back(swcNodes->GetValue(r,0).ToInt());
			swcNodes->RemoveRow(r);
			--r;
		}
	}

	for(int i=0; i<soma_ids.size(); ++i)
	{
		if(i+1 > (int)swcNodes->GetNumberOfRows()) break;
		std::vector<int>::iterator posn1 = find(del_ids.begin(), del_ids.end(), swcNodes->GetValue(i,0).ToInt());
		if(posn1 != del_ids.end())
		{
			swcNodes->RemoveRow(i);
			--i;
		}
	}
	std::cout << "Writing text File" << filename << " with " << swcNodes->GetNumberOfRows() << " nodes...";
	std::ofstream outfile(filename.c_str());

	for (int row = 0; row < (int)swcNodes->GetNumberOfRows(); ++row) 
	{
		for (int col = 0; col < (int)swcNodes->GetNumberOfColumns(); ++col) 
		{
			outfile << swcNodes->GetValue(row,col) << " ";
		}
		outfile << "\n";
	}
	outfile.close();
}

void RenderObject::displayImage(){
	cout<<"write code to display the window"<<endl;
	View3D *View = new View3D();
	View->show();		
	// 
}


void RenderObject::sampleTraces(std::string thresholdVal)
{
	std::cout<<"sample Traces"<<std::endl;
	ImageType3D::RegionType region;
	ImageType3D::IndexType start;
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;

	//ImageType3D::SizeType size;
	/*size[0] = 596;
	size[1] = 600;
	size[2] = 281;*/

	region.SetSize(this->inputImage->GetBufferedRegion().GetSize());
	//region.SetSize(size);
	region.SetIndex(start);
	ImageType3D::SizeType size = region.GetSize();

	ImageType3D::Pointer image = ImageType3D::New();
	image->SetRegions(region);
	image->Allocate();
	image->FillBuffer(0);


	itk::Size<2> maxSize;
	maxSize[0] = inputImage->GetBufferedRegion().GetSize()[0];
	maxSize[1] = inputImage->GetBufferedRegion().GetSize()[1];
	maxSize[2] = inputImage->GetBufferedRegion().GetSize()[2];
	//std::cout<<"maxSize****************"<<maxSize[0]<<"	"<<maxSize[1]<<"	"<<maxSize[2]<<std::endl;
	// get xyz location swc file
	FILE * fp = fopen(thresholdVal.c_str(), "r");
	if(fp==NULL)
	{
		printf("Couldn't open file %s for parsing\n",thresholdVal);
		//	return false;
	}
	char buff[1024];
	int id, type,parent;
	double i,j,k,r;
	int numPoints = 1;
	float total_fg_objectness = 0.0;
	while(!feof(fp))
	{
		if(fgets(buff,1024,fp)==NULL)
		{
			break;
		}
		int pc = 0;
		while(buff[pc]==' '&&(pc<1023))
		{
			pc++;
		}
		if(buff[pc]=='#') // ignoring comment lines for now. We have to read scale from it! TODO
		{
			continue;
		}
		sscanf(buff,"%d %d %lf %lf %lf %lf %d",&id,&type,&i,&j,&k,&r,&parent);
		double radius = r;
		ImageType3D::IndexType pixelIndex;
		pixelIndex[0] = i;
		pixelIndex[1] = j;
		pixelIndex[2] = k;
		this->SWCPoints.push_back(pixelIndex);

		//std::cout<<"i j k r"<<i<<" "<<j<<" "<<k<<"	"<<r<<std::endl;
		try{
			
			if(pixelIndex[0]<maxSize[0]-r && pixelIndex[1]<maxSize[1]-r && pixelIndex[2]<maxSize[2]-r){
				image->SetPixel(pixelIndex, 255);
				//std::cout<<"in if	"<<maxSize[2]-r<<std::endl;
				// set the pixel surrounding the 
				double radius2 = r*r;
				double c1 = i;
				double c2 = j;
				double c3 = k;
				for(unsigned int x = i-r; x < i+r; x++)
				{
					const double dx = static_cast<double>( x -c1 );
					for(unsigned int y = j-r; y < j+r; y++)
					{
						const double dy = static_cast<double>( y -c2 );
						for(unsigned int z = k-r; z< k+r; z++)
						{		
							const double dz = static_cast<double>( z-c3 );
							const double d2 = dx*dx + dy*dy + dz*dz;
							//std::cout<<"d2	"<<d2<<std::endl;
							//std::cout<<"radius	"<<radius2<<std::endl;
							if( d2 < radius2 ){
								ImageType3D::IndexType pixelIndex;
								pixelIndex[0] = x;
								pixelIndex[1] = y;
								pixelIndex[2] = z;
								//std::cout<<"before setPixel"<<std::endl;
									image->SetPixel(pixelIndex, 1);
							}
						}
					}
				}
			}
		}
		catch (itk:: ExceptionObject & err)
		{
			std:: cout << " ExceptionObject caught !" << std :: endl;
			std:: cout << err << std :: endl;
		}

		//float test = this->objectnessImage->GetPixel(pixelIndex);
		////std:: cout << " Vessleness Value at forground" << this->objectnessImage->GetPixel(pixelIndex) <<std :: endl;
		////std:: cout << "Log of Vessleness Value at forground" << test <<std :: endl;
		//if(test != 0){
		//	test = vcl_log(test);
		//	//std:: cout << "in if Log of Vessleness Value at forground" << test <<std :: endl;
		//	total_fg_objectness += test; 
		//}else{
		//	total_fg_objectness += -100000;
		//}
	}
	this->sampledImage = image;
	//sampledFileName.append(thresholdVal);
	//sampledFileName.append(".tif");
	//writeImage(sampledImage,sampledFileName);
	//std::cout<<"total fg objectness		"<<total_fg_objectness<<std::endl;
	fclose(fp);
	//std::cout<<"sample Traces End"<<std::endl;
}

void RenderObject::getSWCPoints(std::string filePath, std::vector<IndexType> &swcPoints)
{
	ImageType3D::RegionType region;

	// get xyz location swc file
	std::string fileName = filePath;

	FILE * fp = fopen(fileName.c_str(), "r");
	//this->ParseFileName(fileName.c_str());
	if(fp==NULL)
	{
		printf("Couldn't open file %s for parsing\n",fileName);
		//	return false;
	}
	char buff[1024];
	int id, type,parent;
	double i,j,k,r;
	int numPoints = 1;

	while(!feof(fp))
	{
		if(fgets(buff,1024,fp)==NULL)
		{
			break;
		}
		int pc = 0;
		while(buff[pc]==' '&&(pc<1023))
		{
			pc++;
		}
		if(buff[pc]=='#') // ignoring comment lines for now. We have to read scale from it! TODO
		{
			continue;
		}
		sscanf(buff,"%d %d %lf %lf %lf %lf %d",&id,&type,&i,&j,&k,&r,&parent);
		double radius = r;
		ImageType3D::IndexType pixelIndex;
		pixelIndex[0] = i;
		pixelIndex[1] = j;
		pixelIndex[2] = k;
		swcPoints.push_back(pixelIndex);
		
	}
	fclose(fp);
}

void RenderObject::writeImage(ImageType3D::Pointer &image,std::string fileName){
	try
	{
		std::cout<<"write Image"<<std::endl;

		//itk::CastImageFilter<ImageType3D,UnsignedCharImageType >::Pointer caster2 = itk::CastImageFilter<ImageType3D,UnsignedCharImageType>::New();
		//caster2->SetInput(image);


		itk::ImageFileWriter<ImageType3D>::Pointer writer = itk::ImageFileWriter<ImageType3D>::New();
		writer->SetInput(image);
		writer->SetFileName(fileName );
		writer->Update();
	}
	catch (itk:: ExceptionObject & err)
	{
		std:: cout << " ExceptionObject caught !" << std :: endl;
		std:: cout << err << std :: endl;
	}

}

void RenderObject::getSiftFeatures(ImageType3D::Pointer &image)
{

	//typedef itk::ScaleInvariantFeatureImageFilter<ImageType3D, Dimension> SiftFilterType;
	//typedef itk::ImageSource< ImageType3D > ImageSourceType;
	//
	//SiftFilterType siftFilter1;
	//SiftFilterType::PointSetTypePointer keypoints1;
	//keypoints1 = siftFilter1.getSiftFeatures(image);


}





//float RenderObject::getCalcThreshold(std::vector<float> &features, std::string type){
//	std::cout<<"in getCalcThreshold"<<std::endl;
//	std::cout<<features[0]<<std::endl;
//	float threshold = 0.00;
//	if(type == "intensity"){
//		//Const = 0
//		//		var	PReg
//		//		2	-0.2046
//		//		4	25.9873
//		//		6	-60.3184
//		//		11	-640.2467
//		//		12	363.7265
//		//		13	987.3842
//		//		14	-406.3363
//		threshold = features[1]*(-0.2046)
//			+features[3]*(25.9873)
//			+features[5]*(-60.3184)
//			+features[10]*(-640.2467)
//			+features[11]*(363.7265)
//			+features[12]*(987.3842)
//			+features[13]*(-406.3363);
//
//	}else if(type=="contrast"){
//		//Const = 0
//		//var	PReg
//		//1		-0.0109
//		//4		-12.4757
//		//6		70.7617
//		//8		-188.6227
//		//9		-13.0697
//		//10	330.6868
//		//12	-397.2717
//		//14	243.9354
//		threshold = features[0]*(-0.0109)
//			+features[3]*(-12.4757)
//			+features[5]*(70.7617)
//			+features[7]*(-188.6227)
//			+features[8]*(-13.0697)
//			+features[9]*(330.6868)
//			+features[11]*(-397.2717)
//			+features[13]*(243.9354);
//	}else{
//		threshold = 0.003;//Return default
//	}
//	return threshold;
//}
std::vector<float> RenderObject::loadThresholds(float minVal,float maxVal,float count)
{
	std::vector<float> vals;
	float increment = (maxVal - minVal)/count;
	for(int i =0;i<=count;i++){
		vals.push_back(minVal);
		minVal = minVal + increment;
    }
	return vals;
}
