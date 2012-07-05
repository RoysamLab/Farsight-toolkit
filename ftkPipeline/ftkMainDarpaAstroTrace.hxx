template< typename TINPUT >
void ftkMainDarpaAstroTrace::splitStore( typename TINPUT::Pointer ImageMontage, std::string nameInput )
{
	int found = nameInput.find_last_of("/\\");
	std::string nameInputNoPath = nameInput.substr(found+1);
	
	computeSplitConst( ImageMontage );
	itk::Size<3> ImageMontageSize = ImageMontage->GetLargestPossibleRegion().GetSize();
	
	int contadorSplit = 0;
#pragma omp parallel for collapse(3) //TEST
	for(int xco = 0; xco < _kx; xco++)
	{
		for(int yco = 0; yco < _ky; yco++)
		{
			for(int zco = 0; zco < _kz; zco++)
			{
// 				#pragma omp critical
// 				{

// 				}
				
				typename TINPUT::RegionType regionLocal_all = ComputeLocalRegionSplit( ImageMontageSize, xco, yco, zco );
				typename TINPUT::RegionType regionMontage_all = ComputeGlobalRegionSplit( ImageMontageSize, xco, yco, zco );
				
				typename TINPUT::Pointer imageLocal = TINPUT::New();
				imageLocal->SetRegions(regionLocal_all);
				imageLocal->Allocate();
				if(imageLocal->GetBufferPointer()==NULL)
					printf("Couldn't allocate memory - 4 .. going to crash now\n");

				typedef itk::ImageRegionIterator<TINPUT> IteratorType;
				IteratorType iterMontage(ImageMontage,regionMontage_all);
				IteratorType iterLocal(imageLocal,regionLocal_all);

				iterMontage.GoToBegin();
				iterLocal.GoToBegin();
				for(;!iterMontage.IsAtEnd();++iterMontage,++iterLocal)
				{
					iterLocal.Set(iterMontage.Get());
				}
				
				// Store local image
				std::stringstream out_x;
				std::stringstream out_y;
				std::stringstream out_z;
				out_x<<xco;
				out_y<<yco;
				out_z<<zco;
				std::string xStr = out_x.str();
				std::string yStr = out_y.str();
				std::string zStr = out_z.str();
				// !!! this has to be nrrd is just to test
				#pragma omp critical
				{
					++contadorSplit;
					std::cout<<std::endl<< "ImageSplit " << contadorSplit << " of " << _kx*_ky*_kz;
					std::string temp9 = _outPathTemp+"/"+nameInputNoPath+"_"+xStr+"_"+yStr+"_"+zStr+".nrrd";
					writeImage< TINPUT >(imageLocal,temp9.c_str());
				}
			}
		}
	}
}

template< typename TINPUT >
void ftkMainDarpaAstroTrace::computeSplitConst( typename TINPUT::Pointer ImageMontage )
{
	itk::Size<3> ImageMontageSize = ImageMontage->GetLargestPossibleRegion().GetSize();

	_kx = ImageMontageSize[0] /(_xTile-_xTileBor);
	_ky = ImageMontageSize[1] /(_yTile-_yTileBor);
	_kz = ImageMontageSize[2] /(_zTile-_zTileBor);
// 	std::cout << std::endl << _kx << " " << _ky << " " << _kz;

	int remx = ImageMontageSize[0] % (_xTile-_xTileBor);
	int remy = ImageMontageSize[1] % (_yTile-_yTileBor);
	int remz = ImageMontageSize[2] % (_zTile-_zTileBor);
// 	std::cout << std::endl << remx << " " << remy << " " << remz;

	itk::Size<3> size_test;
	size_test[0] =  MINNIC((_kx)*(_xTile-_xTileBor)+_xTile-1,ImageMontageSize[0]-1) -  (_kx) * (_xTile-_xTileBor) +1;
	size_test[1] =  MINNIC((_ky)*(_yTile-_yTileBor)+_yTile-1,ImageMontageSize[1]-1) -  (_ky) * (_yTile-_yTileBor) +1;
	size_test[2] =  MINNIC((_kz)*(_zTile-_zTileBor)+_zTile-1,ImageMontageSize[2]-1) -  (_kz) * (_zTile-_zTileBor) +1;
// 	std::cout << std::endl << size_test[0] << " " << size_test[1] << " " << size_test[2];

	if( size_test[0] > _xTileBor  )
	{
		if ( remx > 0 )
			_kx ++;
	}
	if( size_test[1] > _yTileBor  )
	{
		if ( remy > 0 )
			_ky ++;
	}
	if( size_test[2] > _zTileBor  )
	{
		if ( remz > 0 )
			_kz ++;
	}
	
	std::cout << std::endl << size_test[0] << " " << size_test[1] << " " << size_test[2];
	std::cout << std::endl << _kx << " " << _ky << " " << _kz;
	std::cout << std::endl << remx << " " << remy << " " << remz;
}