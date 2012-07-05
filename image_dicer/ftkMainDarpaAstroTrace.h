// ############################################################################################################################################################################
#ifndef _ftkMainDarpaAstroTrace_h_
#define _ftkMainDarpaAstroTrace_h_
// ############################################################################################################################################################################

#ifdef _OPENMP
#include "omp.h"
#endif

// GLOBAL INCLUDES
#include "ftkMainDarpaGlobalInclude.h"
#include "ftkMainDarpa.h"

// TEMPLATES
#include "ftkMainDarpaTemplates.h"

// DEFINITIONS
#include "ftkMainDarpaDeclaration.h"

//MACROS
 #define MINNIC(a,b) (((a) > (b))? (b) : (a))
 #define MAXNIC(a,b) (((a) < (b))? (b) : (a))

// ############################################################################################################################################################################
// MAIN CLASS
class ftkMainDarpaAstroTrace
{
public:
	
	void readParameters( std::string );
	void runPreprocesing(); // Before preprocess always run redParams before
	template< typename TINPUT >
	void splitStore( typename TINPUT::Pointer ImageMontage, std::string );
	void runSplitting();
	void runInterestPoints();
	void runTracing();
	
protected:
	void computeSplitConst( rawImageType_8bit::Pointer ImageMontage );
	std::vector< itk::Index<3> > getCentroidList();
	std::vector< itk::Index<3> > getSomaTable( std::vector< itk::Index<3> > , int , int , int );
	template< typename TINPUT >
	void computeSplitConst( typename TINPUT::Pointer ImageMontage );
	template< typename TINPUT >
	typename TINPUT::Pointer cropImages( typename TINPUT::Pointer , int , int , int );
	//void WriteCenterTrace(vtkSmartPointer< vtkTable > , int , int , int , std::string );

	void RemoveLabelNearBorder(rawImageType_8bit::RegionType, std::vector< rawImageType_16bit::Pointer >&, std::vector< vtkSmartPointer< vtkTable > >&, std::vector< std::map< unsigned int, itk::Index<3> > >& );
	rawImageType_8bit::RegionType ComputeLocalRegionSegment( itk::Size<3>, int, int, int );
	rawImageType_uint::RegionType ComputeGlobalRegionSegment( itk::Size<3>, int, int, int );
	rawImageType_8bit::RegionType ComputeLocalRegionSplit( itk::Size<3>, int, int, int );
	rawImageType_8bit::RegionType ComputeGlobalRegionSplit( itk::Size<3>, int, int, int );
	
private:
	std::string _segmentParams;
	int _xSize;
	int _ySize;
	int _zSize;
	int _xTile;
	int _yTile;
	int _zTile;
 	int _xTileBor;
 	int _yTileBor;
 	int _zTileBor;
	int _num_threads;
	std::string _Cy5_Image;
	std::string _TRI_Image;
	std::string _GFP_Image;
	std::string _DAP_Image;
	std::string _Dist_Map_Image;
	std::string _Soma_Centroids;
	std::string _Soma_Montage;
	int _isSmall;
	std::string _astroTraceParams;
// 	std::string _projectDefinition;
//	std::string _optionsASTR;
	std::string _outPath;
	std::string _outPathDebug;
	std::string _outPathDebugLevel2;
	std::string _outPathTemp;
	
	std::string _Cy5_ImageNRRD;
	std::string _TRI_ImageNRRD;
	std::string _GFP_ImageNRRD;
	std::string _DAP_ImageNRRD;
	std::string _Dist_Map_ImageNRRD;
	std::string _Soma_MontageNRRD;
// 	
// // 	std::vector< rawImageType_8bit::Pointer > _ImageMontages;
 	rawImageType_8bit::Pointer _ImageMontage_Cy5;
 	rawImageType_8bit::Pointer _ImageMontage_TRI;
 	rawImageType_8bit::Pointer _ImageMontage_GFP;
 	rawImageType_8bit::Pointer _ImageMontage_DAP;
	rawImageType_16bit::Pointer _ImageMontage_Dist_Map;
// 	
 	itk::Size<3> _ImageMontage_Cy5Size;
 	itk::Size<3> _ImageMontage_TRISize;
 	itk::Size<3> _ImageMontage_GFPSize;
 	itk::Size<3> _ImageMontage_DAPSize;
	itk::Size<3> _ImageMontage_Dist_MapSize;
 	
 	// Image Split
 	int _kx;
 	int _ky;
 	int _kz;

	// MNT
	std::string _TRI_ImagePREPMNT;
	
	// Division
	int _numDivisionsInRowCEN;
	std::vector< itk::Index<3> > _initialBigTileLOG;
	std::vector< itk::Size<3> > _sizeOfBigTilesLOG;
	std::vector< itk::Index<3> > _initialBigTileCEN;
	std::vector< itk::Size<3> > _sizeOfBigTilesCEN;
	
	// BIG IMAGES
	rawImageType_flo::Pointer _img_astroTraceDesiredRegion;
	rawImageType_uint::Pointer _somaMontageDesiredRegion;
};

#include "ftkMainDarpaAstroTrace.hxx"

#endif