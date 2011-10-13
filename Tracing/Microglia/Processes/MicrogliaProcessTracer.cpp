#include "MicrogliaProcessTracer.h"
  
bool CompareNeighbors(std::pair< double, Node *> n1,
                      std::pair< double, Node *> n2)
{
  if(n1.first < n2.first)
    {
    return true;
    }
  return false;
}

Node::Node()
  {
    this->parent = NULL;
    this->IsOpen = true;
  }

MicrogliaProcessTracer::MicrogliaProcessTracer()
{
  this->NodeCounter = 1;
  this->Padding = 1;

  //no default distance threshold
  //this means all nodes will be added to a tree
  this->MaxDistance = std::numeric_limits<double>::max();
}

MicrogliaProcessTracer::~MicrogliaProcessTracer()
{
  std::vector< Node * >::iterator nodeItr;
  for(nodeItr = this->Open.begin(); nodeItr != this->Open.end(); ++nodeItr)
    {
    delete (*nodeItr);
    }
  for(nodeItr = this->Closed.begin(); nodeItr != this->Closed.end(); ++nodeItr)
    {
    delete (*nodeItr);
    }
} 

void MicrogliaProcessTracer::LoadInputImage(std::string fname)
{
  std::cout << "Reading input file "<< fname << std::endl;
  ReaderType::GlobalWarningDisplayOff();
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(fname);
  ImageType3D::Pointer image = reader->GetOutput();
  image->Update();

  this->LoadInputImage(image);
}

void MicrogliaProcessTracer::LoadInputImage(ImageType3D::Pointer &image)
{
  RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetOutputMinimum(0.0);
  rescaler->SetOutputMaximum(1.0);
  rescaler->SetInput(image);
  
  //Median filter 
  MedianFilterType::Pointer medfilt = MedianFilterType::New();
  medfilt->SetNumberOfThreads(16);
  medfilt->SetInput(rescaler->GetOutput());
  ImageType3D::SizeType rad = { {1, 1, 1} };
  medfilt->SetRadius(rad);
  medfilt->Update();
  //this->InputImage = medfilt->GetOutput();
  this->InputImage = rescaler->GetOutput();

  //pad z slices
  itk::Size<3> isz = this->InputImage->GetBufferedRegion().GetSize();
  itk::Size<3> osz = isz;
  osz[2] += 2*this->Padding;
  itk::Index<3> indx, ondx;
  
  this->PaddedInputImage = ImageType3D::New();
  this->PaddedInputImage->SetRegions(osz);
  this->PaddedInputImage->Allocate();
  this->PaddedInputImage->SetSpacing(this->InputImage->GetSpacing());
  
  for(ondx[2] = 0; ondx[2] < (unsigned int)osz[2]; ++ondx[2]) 
  {
    indx[2] = (ondx[2] < this->Padding) ? 0 : ondx[2] - this->Padding;
    indx[2] = (ondx[2] >= (unsigned int)osz[2]-this->Padding) ? isz[2]-1 : indx[2];
    for(ondx[1] = 0; ondx[1] < (unsigned int)osz[1]; ++ondx[1]) 
    {
      indx[1] = ondx[1];
      for(ondx[0] = 0; ondx[0] < (unsigned int)osz[0]; ++ondx[0]) 
      {
        indx[0] = ondx[0];
        this->PaddedInputImage->SetPixel(ondx, this->InputImage->GetPixel(indx));
      }
    }
  }

  std::cout << "Input file size (after zero padding) is " << this->PaddedInputImage->GetBufferedRegion().GetSize() << std::endl;
}

///////////////////////////////////////////////////////////////////////
std::vector< Node * > MicrogliaProcessTracer::ReadListOfPoints(std::string fname)
{
  std::vector< Node * > listOfPoints;
  std::string temp, num;
  std::ifstream infile;
  infile.open(fname.c_str());
  size_t x1, x2;

  while(!infile.eof()) 
  {
    std::getline(infile,temp);
    if (temp.length() < 1)
      continue;
    
    x1 = temp.find_first_of("0123456789.");
    x2 = temp.find_first_not_of("0123456789.",x1);
    if ((x2 - x1) > 10)
      continue;
    
    num = temp.substr(x1,x2-x1);
    int x = atoi(num.c_str());

    x1 = temp.find_first_of("0123456789.",x2+1);
    x2 = temp.find_first_not_of("0123456789.",x1);
    if ((x2 - x1) > 10)
      continue;
    
    num = temp.substr(x1,x2-x1);
    int y = atoi(num.c_str());

    x1 = temp.find_first_of("0123456789.",x2+1);
    x2 = temp.find_first_not_of("0123456789.",x1);
    if (x2 > temp.length())
      x2 = temp.length();
    
    if ((x2 - x1) > 10)
      continue;
    
    num = temp.substr(x1,x2-x1);
    int z = atoi(num.c_str());

    itk::Index<3> idx;
    idx[0] = x;
    idx[1] = y;
    idx[2] = z;
    Node *n = new Node();
    n->ID = this->NodeCounter++;
    n->index = idx;
    listOfPoints.push_back(n);
  }
  infile.close();

  return listOfPoints;
}
	
////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::LoadSeedPoints(std::string fname)
{
  this->Closed = this->ReadListOfPoints(fname);

  //nodes are created open by default, so close all these ones.
  std::vector< Node * >::iterator nodeItr;
  for(nodeItr = this->Closed.begin(); nodeItr != this->Closed.end(); ++nodeItr)
    {
    (*nodeItr)->IsOpen = false;
    }
}

////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::LoadNodes(std::string fname)
{
  this->Open = this->ReadListOfPoints(fname);
}

////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::LoadSomaImage(std::string somaFileName)
{
  typedef itk::ImageFileReader<CharImageType3D> SomaReaderType;
  SomaReaderType::Pointer somaReader = SomaReaderType::New();
  somaReader->SetFileName(somaFileName);
  this->SomaImage = somaReader->GetOutput();
  somaReader->Update();
  this->SomaImage->SetSpacing(this->InputImage->GetSpacing());
}

///////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::RunTracing()
{
  if(!this->InputImage)
    {
    std::cerr << "You must call LoadImage() before RunTracing()" << std::endl;
    return;
    }
  if(this->Closed.empty())
    {
    std::cerr << "You must call LoadSeedPoints() before RunTracing()" << std::endl;
    return;
    }

  if(this->SomaImage)
    {
    this->MaskAwaySomas();
    }

  this->CalculateCriticalPoints();

  itk::ImageRegionConstIterator<ImageType3D> Nit(this->NDXImage, this->NDXImage->GetBufferedRegion());
  for (Nit.GoToBegin(); !Nit.IsAtEnd(); ++Nit) 
  {
    if (Nit.Get() > 0) 
    {
      Node *n = new Node();
      n->ID = this->NodeCounter++;
      n->index = Nit.GetIndex();
      this->Open.push_back(n);
    }
  }
  
  std::cout << "open list size: " << this->Open.size() << std::endl;
  std::cout << "closed list size: " << this->Closed.size() << std::endl;

  //use RATS to threshold the input image.
  //this is used in GetDistanceBetweenPoints 
  typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< ImageType3D, ImageType3D > GradientType;
  GradientType::Pointer gradient = GradientType::New();
  gradient->SetInput( this->InputImage );
  gradient->SetSigma( 0.25 );
  gradient->Update();
  typedef itk::RobustAutomaticThresholdImageFilter< ImageType3D, ImageType3D > ThresholdFilterType;
  ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
  thresholdFilter->SetInput( this->InputImage );
  thresholdFilter->SetGradientImage( gradient->GetOutput() );
  thresholdFilter->SetPow( 1.0 );
  this->ThresholdedImage = thresholdFilter->GetOutput();
  thresholdFilter->Update();
 
  std::cout << "computing adjacency matrix for seed points" << std::endl;
  this->ComputeAdjacencies(this->Closed);
  
  std::cout << "computing adjacency matrix for critical points" << std::endl;
  this->ComputeAdjacencies(this->Open);
  
  this->BuildTrees();
}

////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::CalculateCriticalPoints(void)
{
  std::cout << std::endl<< "Searching for critical points" << std::endl;
  this->NDXImage = ImageType3D::New();
  this->NDXImage->SetRegions(this->PaddedInputImage->GetBufferedRegion());
  this->NDXImage->Allocate();
  this->NDXImage->FillBuffer(0.0f);
  
  //float power = 0.0;
  float power = -0.25;

  for (unsigned int i = 0; i < 6; ++i)
  {
    power = power + 0.25;
    float sigma = vcl_pow(2, power);
    std::cout << "Analysis at sigma " << sigma << std::endl;
    CalculateCriticalPointsAtScale( sigma );
  }
}

////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::CalculateCriticalPointsAtScale( float sigma ) 
{
  typedef itk::LaplacianRecursiveGaussianImageFilter< ImageType3D , ImageType3D> GFilterType;
  GFilterType::Pointer gauss = GFilterType::New();
  gauss->SetInput( this->PaddedInputImage );
  gauss->SetSigma( sigma );
  gauss->SetNormalizeAcrossScale(false);
  gauss->GetOutput()->Update();
  
  itk::ImageRegionIterator<ImageType3D> ittemp(gauss->GetOutput(), gauss->GetOutput()->GetBufferedRegion());
  float tnorm = vcl_pow(sigma, 1.6f);
  for(ittemp.GoToBegin(); !ittemp.IsAtEnd(); ++ittemp) 
  {
    float q = ittemp.Get()*tnorm;
    ittemp.Set(-1.0f*q);
  }

  // set the diagonal terms in neighborhood iterator
  itk::Offset<3>
    xp =  {{2 ,  0 ,   0}},
    xn =  {{-2,  0,    0}},
    yp =  {{0,   2,   0}},
    yn =  {{0,  -2,    0}},
    zp =  {{0,   0,    2}},
    zn =  {{0,   0,   -2}};

  itk::Size<3> rad = {{1,1,1}};
  itk::NeighborhoodIterator<ImageType3D> nit(rad , gauss->GetOutput(), gauss->GetOutput()->GetBufferedRegion());
  itk::ImageRegionIterator<ImageType3D> it(gauss->GetOutput(), gauss->GetOutput()->GetBufferedRegion());

  unsigned int
    xy1 =  17, //{ 1 ,   1 ,  0 },
    xy2 =  9,  //{ -1,  -1 ,  0 },
    xy3 =  15, //{ -1,   1 ,  0 },
    xy4 =  11, //{ 1 ,  -1 ,  0 },

    yz1 =  25, //{ 0 ,   1 ,  1 },
    yz2 =  1,  //{ 0 ,  -1 , -1 },
    yz3 =  19, //{ 0 ,  -1 ,  1 },
    yz4 =  7,  //{ 0 ,   1 , -1 },

    xz1 =  23, //{ 1 ,   0 ,  1 },
    xz2 =  3,  //{-1 ,   0 , -1 },
    xz3 =  21, //{-1 ,   0 ,  1 },
    xz4 =  5;  //{ 1 ,   0 , -1 };

  typedef itk::FixedArray< double, 3 > EigenValuesArrayType;
  typedef itk::Matrix< double, 3, 3 > EigenVectorMatrixType;
  typedef itk::SymmetricSecondRankTensor<double,3> TensorType;

  itk::Size<3> sz = this->PaddedInputImage->GetBufferedRegion().GetSize();
  sz[0] = sz[0] - 3;
  sz[1] = sz[1] - 3; 
  sz[2] = sz[2] - 3;

  it.GoToBegin();
  nit.GoToBegin();
  itk::Vector<float,3> sp = this->PaddedInputImage->GetSpacing();

  long win = long(sigma)/2;
  if (win <2) 
  {
    win = 2;
  }
  
  long ctCnt = 0;
  while(!nit.IsAtEnd()) 
  {
    itk::Index<3> ndx = it.GetIndex();
    if ( (ndx[0] < 2) || (ndx[1] < 2) || (ndx[2] < 2) ||
      (ndx[0] > (unsigned int)sz[0]) || (ndx[1] > (unsigned int)sz[1]) ||
      (ndx[2] > (unsigned int)sz[2]) )
    {
      ++it;
      ++nit;
      continue;
    }

    float a1 = 0.0;
    for (unsigned int i=0; i < 13; ++i)
    {
      a1 += vnl_math_max(nit.GetPixel(i), nit.GetPixel(26 - i));
    }
    
    float val = nit.GetPixel(13) ;

    const float thresh1 = 0.03;   // 3% of maximum theshold from Lowe 2004
    const float thresh2 = 0.001;  // -0.1 percent of range

    if ( ((val - a1/13.0f) > thresh2 ) && ( val > thresh1 ))  
    {
      TensorType h;
      h[0] = gauss->GetOutput()->GetPixel( ndx + xp ) + gauss->GetOutput()->GetPixel( ndx + xn ) - 2*nit.GetPixel( 13 );
      h[3] = gauss->GetOutput()->GetPixel( ndx + yp ) + gauss->GetOutput()->GetPixel( ndx + yn ) - 2*nit.GetPixel( 13 );
      h[5] = gauss->GetOutput()->GetPixel( ndx + zp ) + gauss->GetOutput()->GetPixel( ndx + zn ) - 2*nit.GetPixel( 13 );
      h[1] = nit.GetPixel(xy1) + nit.GetPixel(xy2) - nit.GetPixel(xy3) - nit.GetPixel(xy4);
      h[2] = nit.GetPixel(xz1) + nit.GetPixel(xz2) - nit.GetPixel(xz3) - nit.GetPixel(xz4);
      h[4] = nit.GetPixel(yz1) + nit.GetPixel(yz2) - nit.GetPixel(yz3) - nit.GetPixel(yz4);

      EigenValuesArrayType ev;
      EigenVectorMatrixType em;
      h.ComputeEigenAnalysis (ev, em);

      unsigned int w = ShapeAnalysis(ev);
      float value = vnl_math_abs(ev[0]) + vnl_math_abs(ev[1]) + vnl_math_abs(ev[2]) - vnl_math_abs(ev[w]);
      if (RegisterIndex(value, ndx, sz, win)) 
      {
        this->NDXImage->SetPixel(ndx,value);
        ctCnt++;
      }
    }
    ++it;
    ++nit;
  }
  std::cout <<"Number of CTs at this stage: " << ctCnt <<std::endl;
  
}

unsigned int MicrogliaProcessTracer::ShapeAnalysis(const itk::FixedArray<float, 3> &ev)
{
  unsigned int w;
  float L1, L2, L;
  if ( (ev[0] > ev[1]) && (ev[0] > ev[2]) ) 
  {
    w = 0;
    L = ev[0];
    L1 = ev[1]; 
    L2 = ev[2];
    if (ev[1] > ev[2])
    {
      L1 = ev[2]; 
      L2 = ev[1];   
    }
  }

  else if( (ev[1] > ev[0]) && (ev[1] > ev[2]) ) 
  {
    w = 1;
    L = ev[1];
    L1 = ev[0];
    L2 = ev[2];
    if (ev[0] > ev[2]) 
    {
      L1 = ev[2];
      L2 = ev[0];   
    }
  }

  else  
  {
    w = 2;
    L = ev[2];
    L1 = ev[0];
    L2 = ev[1];
    if (ev[0] > ev[1]) 
    {
      L1 = ev[1];
      L2 = ev[0];   
    }
  }

  return w;
}

bool MicrogliaProcessTracer::RegisterIndex(const float value, itk::Index<3> &ndx, itk::Size<3>& sz, long h = 2) 
{
  itk::Index<3> n;
  bool higherPresent = false;
  for (n[0] = ndx[0]-h; n[0] <= ndx[0]+h; ++n[0]) 
  {
    for (n[1] = ndx[1]-h; n[1] <= ndx[1]+h; ++n[1]) 
    {
      for (n[2] = ndx[2]-h; n[2] <= ndx[2]+h; ++n[2]) 
      {
        if ( (n[0] < 2) || (n[1] < 2) || (n[2] < 2) || (n[0] > (unsigned int)sz[0]) ||
          (n[1] > (unsigned int)sz[1]) || (n[2] > (unsigned int)sz[2]) )
        {
          continue;
        }

        float curval = this->NDXImage->GetPixel(n);
        if (value > curval) 
        {
          this->NDXImage->SetPixel(n,0.0f);
        }
        else if (value < curval) 
        {
          higherPresent = true;       
        }
      }
    }
  }
  if (higherPresent == true) 
  {
    return false;
  }
  
  return true;
}

////////////////////////////////////////////////////////////////////////////////
float MicrogliaProcessTracer::GetRadius(itk::Index<3> idx) 
{
  float r = 2.0f;
  itk::Vector<float,3> m1, m2, m, pos;
  itk::Index<3> ndx;

  pos[0] = idx[0];
  pos[1] = idx[1];
  pos[2] = idx[2];
      
  for (int iter = 0; iter < 20; ++iter) 
  {
    for( int i = 0; i<3; i++) 
    {
      m1[i] = pos[i] - vnl_math_max(2.0f*r, 5.0f);
      m2[i] = pos[i] + vnl_math_max(2.0f*r, 5.0f);
    }

    std::vector<float> c;
    c.reserve(4*4*int(r*r));
    float i1 = 0.0f, i2 = 0.0f, i1s = 0.0f, i2s = 0.0f;
    for (m[2] = m1[2]; m[2] <= m2[2]; m[2]++) 
    {
      for (m[1] = m1[1]; m[1] <= m2[1]; m[1]++) 
      {
        for (m[0] = m1[0]; m[0] <= m2[0]; m[0]++) 
        {
          ndx.CopyWithRound(m);
          itk::Vector<float,3> mm = pos - m;
          float d = mm.GetNorm();
          if (this->PaddedInputImage->GetBufferedRegion().IsInside(ndx)) 
          {
            float val = this->PaddedInputImage->GetPixel(ndx);
            if (d < r) 
            {
              i1 += val;
              ++i1s;
            }
            else 
            {
              i2 += val;
              ++i2s;
            }

            if (vnl_math_abs(d - r) < 0.7f) 
              c.push_back(val);           
          }
        }
      }
    }
    i1 /= i1s;
    i2 /= i2s;
    float dr = 0.0f;
    for (std::vector<float>::iterator it = c.begin(); it < c.end(); ++it) 
      dr += vnl_math_abs((*it) - i1) - vnl_math_abs((*it) - i2);
    
    dr *= 1.0f / float(c.size()); //rate
    dr = vnl_math_max(dr , -1.0f);
    dr = vnl_math_min(dr , 1.0f);
    r -= dr;
    r = vnl_math_max(r , 1.0f) ; 
  }
  return r;
}

////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::BuildTrees()
{
  std::pair< Node *, Node * > parentAndChild;
  Node *parent = NULL;
  Node *child = NULL;

  while( !( this->Open.empty() ) )
    {
    parentAndChild = this->FindClosestOpenNode();
    if(parentAndChild.first == NULL)
      {
      //We've added all the nodes within our distance threshold.
      //Delete those that remain in the open list.
      std::vector< Node * >::iterator nodeItr;
      std::cout << this->Open.size() << " nodes were too far away to get added to a cell model" << std::endl;
      for(nodeItr = this->Open.begin(); nodeItr != this->Open.end(); ++nodeItr)
        {
        delete (*nodeItr);
        }
      this->Open.clear();
      break;
      }
    
    parent = parentAndChild.first;
    child = parentAndChild.second;

    //setup the parent/child relationship between these two nodes
    parent->children.push_back(child);
    
    //remove the child node from the Open list
    bool itsThere = false;
    if(std::find(this->Open.begin(), this->Open.end(), child) != this->Open.end())
      {
      itsThere = true;
      }
    if(!itsThere)
      {
      std::cout << "(Debug) PROBLEM: Before deletion, child is NOT in Open." << std::endl;
      }
    std::vector< Node* >::iterator newEnd =
      std::remove(this->Open.begin(), this->Open.end(), child);
    this->Open.erase(newEnd, this->Open.end());
    if(std::find(this->Open.begin(), this->Open.end(), child) != this->Open.end())
      {
      itsThere = true;
      }
    else
      {
      itsThere = false;
      }
    if(itsThere)
      {
      std::cout << "(Debug) PROBLEM: After deletion, child is still in Open." << std::endl;
      }
    child->parent = parent;
  
    //and add it to the Closed list
    child->IsOpen = false;
    this->Closed.push_back(child);
    }
}

////////////////////////////////////////////////////////////////////////////////
// This function returns (NULL, NULL) when the closest open node is further away
// than our distance threshold
std::pair< Node *, Node * > MicrogliaProcessTracer::FindClosestOpenNode()
{
  std::pair< Node *, Node * > parentAndChild;
  std::pair< double, Node * > distanceAndNeighbor;
  std::list< std::pair< double, Node * > > *neighborList;
  double minDistance = std::numeric_limits<double>::max();
  std::vector< Node * >::iterator nodeItr, nbrItr;
  //std::vector< Node * >::iterator nodeItr, nbrItr;
 
  for(nodeItr = this->Closed.begin(); nodeItr != this->Closed.end(); ++nodeItr)
    {
    neighborList = &(this->AdjacencyMap[(*nodeItr)]);
    bool openNeighborFound = false;
    while(!openNeighborFound && !neighborList->empty())
      {
      distanceAndNeighbor = neighborList->front();
      //std::cout << "front: " << distanceAndNeighbor.first << ", back: " << neighborList->back().first << std::endl;
      if(!distanceAndNeighbor.second->IsOpen)
        {
        neighborList->pop_front();
        }
      else
        {
        openNeighborFound = true;
        if(distanceAndNeighbor.first < minDistance)
          {
          minDistance = distanceAndNeighbor.first;
          parentAndChild.first = (*nodeItr);
          parentAndChild.second = distanceAndNeighbor.second;
          }
        }
      }
    }

 //the child node is closed in the caller function
 return parentAndChild;

/*
  for(nodeItr = this->Open.begin(); nodeItr != this->Open.end(); ++nodeItr)
    {
    itk::Point<double,3> start;
    this->InputImage->TransformIndexToPhysicalPoint( (*nodeItr)->index, start ); 

    for(nbrItr = this->Closed.begin(); nbrItr != this->Closed.end(); ++nbrItr)
      {
      if( (*nodeItr)->ID == (*nbrItr)->ID )
        {
        std::cout << "this should never happen" << std::endl;
        continue;
        }

      itk::Point<double,3> end;
      this->InputImage->TransformIndexToPhysicalPoint( (*nbrItr)->index, end ); 
      
      double d = start.EuclideanDistanceTo(end);
      if( d < minDistance )
        {
        minDistance = d;
        parentAndChild.first = (*nbrItr);
        parentAndChild.second = (*nodeItr);
        }
      }
    }
  //std::cout << "(Debug) " << minDistance << std::endl;
  
  //check distance threshold
  if(minDistance > this->MaxDistance)
    {
    parentAndChild.first = NULL;
    parentAndChild.second = NULL;
    }

  return parentAndChild;
*/
}

////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::WriteSWC( std::string fname )
{
  if (this->Open.size() != 0 )
    {
    std::cout << "ERROR: Open list is not empty.  Run BuildTrees first." << std::endl;
    return;
    }

  std::ofstream output( fname.c_str() );
  
  std::vector< Node * >::iterator nodeItr;
	
  ImageType3D::PointType origin = this->InputImage->GetOrigin();
  for(nodeItr = this->Closed.begin(); nodeItr != this->Closed.end(); ++nodeItr)
    {
    itk::Index<3> idx = (*nodeItr)->index;
    float x = idx[0] + origin[0];
    float y = idx[1] + origin[1];
    float z = idx[2] + origin[2];
    long parentID;
    int cellType = 3;

    //compute vessel radius at this point
    float r = this->GetRadius( idx );
  
    if ( (*nodeItr)->parent == NULL )
      {
      parentID = -1;
      cellType = 1;
      }
    else
      {
      parentID = (*nodeItr)->parent->ID;
      }
    
    output << (*nodeItr)->ID << " " << cellType << " " << x << " " << y << " " << z << " " << "1" << " " << parentID << std::endl; 
    //output << (*nodeItr)->ID << " " << cellType << " " << x << " " << y << " " << z << " " << r << " " << parentID << std::endl; 
    }
  output.close();
}

////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::MaskAwaySomas()
{
  typedef itk::BinaryBallStructuringElement<
      CharImageType3D::PixelType,3>                  StructuringElementType;
  StructuringElementType structuringElement;
  structuringElement.SetRadius(5);
  structuringElement.CreateStructuringElement();

  typedef itk::BinaryDilateImageFilter <CharImageType3D, CharImageType3D, StructuringElementType>
          BinaryDilateImageFilterType;
  BinaryDilateImageFilterType::Pointer dilateFilter
          = BinaryDilateImageFilterType::New();
  dilateFilter->SetInput( this->SomaImage );
  dilateFilter->SetKernel(structuringElement);

  MaskFilterType::Pointer maskFilter = MaskFilterType::New();
  maskFilter->SetInput( this->InputImage );
  maskFilter->SetMaskImage( dilateFilter->GetOutput() );
  maskFilter->SetOutsideValue(0);
  maskFilter->Update();

  this->InputImage = maskFilter->GetOutput();
}

////////////////////////////////////////////////////////////////////////////////
void MicrogliaProcessTracer::ComputeAdjacencies( std::vector< Node * > nodes )
{
  std::vector< Node * >::iterator nodeItr, nbrItr;
  for(nodeItr = nodes.begin(); nodeItr != nodes.end(); ++nodeItr)
    {
    itk::Index<3> start = (*nodeItr)->index;
    std::list< std::pair< double, Node *> > neighbors;

    for(nbrItr = this->Open.begin(); nbrItr != this->Open.end(); ++nbrItr)
      {
      if( (*nodeItr)->ID == (*nbrItr)->ID )
        {
        continue;
        }

      itk::Index<3> end = (*nbrItr)->index;
      
      itk::Point<double,3> startPoint;
      itk::Point<double,3> endPoint;
      this->InputImage->TransformIndexToPhysicalPoint( start, startPoint );
      this->InputImage->TransformIndexToPhysicalPoint( end, endPoint );
      double euclideanDistance = startPoint.EuclideanDistanceTo(endPoint);
      double scaledDistance = this->GetDistanceBetweenPoints(start, end);

      double thresh = this->MaxDistance * 2; //doubled for intensity scaling
      if( (*nodeItr)->IsOpen == false )
        {
        //thresh += 10;
        }
      if( scaledDistance > this->MaxDistance )
      //if( euclideanDistance > thresh || scaledDistance > this->MaxDistance )
        {
        continue;
        }
      std::pair< double, Node * > nbr;
      nbr.first = euclideanDistance;
      nbr.second = (*nbrItr);
      neighbors.push_front(nbr);
      }
    neighbors.sort(CompareNeighbors);
    this->AdjacencyMap[(*nodeItr)] = neighbors;
    }
}

double MicrogliaProcessTracer::GetDistanceBetweenPoints(itk::Index<3> start, itk::Index<3> end)
{
  ImageType3D::SizeType imageSize = this->ThresholdedImage->GetLargestPossibleRegion().GetSize();

  //compute Euclidean distance between the two points
  itk::Point<double,3> startPoint;
  itk::Point<double,3> endPoint;
  this->InputImage->TransformIndexToPhysicalPoint( start, startPoint );
  this->InputImage->TransformIndexToPhysicalPoint( end, endPoint );
  double euclideanDistance = startPoint.EuclideanDistanceTo(endPoint);
 
  //analyze the pixels along the line between the start & end points.
  //the distance returned by this function is typically less than the actual Euclidean distance.
  //How much less is based on how many "foreground" pixels lie along this line.
  std::vector< itk::Index<3> > indices = this->Line.BuildLine(start, end);
  float scaledDistance = 0;
  float numPixels = 0;
  float numBackground = 0;
  for(unsigned int i = 0; i < indices.size(); i++)
    {
    float distanceWeight = 0;
    if(indices[i] == start)
      {
      continue;
      }
    if(indices[i] == end)
      {
      break;
      }

    //search for foreground pixels in each 3x3x3 neighborhood
    bool foregroundPixelFound = false;
    for(int x = indices[i][0] - 1; x < indices[i][0] + 2; x++)
      {
      if(x < 0 || x > imageSize[0] - 1) 
        {
        continue;
        }
      for(int y = indices[i][1] - 1; y < indices[i][1] + 2; y++)
        {
        if(y < 0 || y > imageSize[0] - 1) 
          {
          continue;
          }
        for(int z = indices[i][2] - 1; z < indices[i][2] + 2; z++)
          {
          if(z < 0 || z > imageSize[0] - 1)
            {
            continue;
            }
          //here's where we check for foreground pixels
          itk::Index<3> nbrIdx = {x, y, z};
          if(this->ThresholdedImage->GetPixel(nbrIdx) != 0)
            {
            foregroundPixelFound = true;
            break;
            }
          }
        if(foregroundPixelFound)
          {
          break;
          }
        }
      if(foregroundPixelFound)
        {
        break;
        }
      }
    if(!foregroundPixelFound)
      {
      distanceWeight = 1.0;
      numBackground++;
      }

    //take image spacing into account
    itk::Offset<3> diff = indices[i] - indices[i-1];
    if(diff[2] != 0)
      {
      scaledDistance += ( this->InputImage->GetSpacing()[2] * distanceWeight ); 
      }
    else if(diff[1] != 0)
      {
      scaledDistance += ( this->InputImage->GetSpacing()[1] * distanceWeight ); 
      }
    else
      {
      scaledDistance += ( this->InputImage->GetSpacing()[0] * distanceWeight ); 
      }
    ++numPixels;
    }

  return scaledDistance;
}
