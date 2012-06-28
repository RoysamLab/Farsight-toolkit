// ############################################################################################################################################################################
#ifndef _ftkMainDarpaSegment_h_
#define _ftkMainDarpaSegment_h_
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
class ftkMainDarpaSegment
{
public:
	
	void readParameters( std::string );
	void runSegment();
	void runSpliting( );
	void runStich( );
	void runStichOneAtTheTime(  ); // HAVE NOT FINISH THE TABLE PART, THE IDEA IS TO RUN STICHIN WITHOUT LOADING THE TWO IMAGES AT THE SAME TIME, ALSO ONLY LOCAL MAPS
	void splitStore( rawImageType_8bit::Pointer ImageMontage, std::string );
	
protected:
	void computeSplitConst( rawImageType_8bit::Pointer ImageMontage );
	void RunSegmentation(rawImageType_8bit::RegionType, std::vector< rawImageType_8bit::Pointer >&, std::vector< rawImageType_16bit::Pointer >&, std::vector< vtkSmartPointer< vtkTable > >&, std::vector< std::map< unsigned int, itk::Index<3> > > &);
	
	rawImageType_16bit::Pointer RunNuclearSegmentation(rawImageType_8bit::Pointer );
	vtkSmartPointer< vtkTable > ComputeFeaturesAndAssociations( std::vector< rawImageType_8bit::Pointer >&, std::vector< rawImageType_16bit::Pointer >& );
	std::map< unsigned int, itk::Index<3> > GetLabelToCentroidMap( vtkSmartPointer< vtkTable > );
	
	void RemoveLabelNearBorder(rawImageType_8bit::RegionType, std::vector< rawImageType_16bit::Pointer >&, std::vector< vtkSmartPointer< vtkTable > >&, std::vector< std::map< unsigned int, itk::Index<3> > >& );
	
	rawImageType_8bit::RegionType ComputeLocalRegionSegment( itk::Size<3>, int, int, int );
	rawImageType_uint::RegionType ComputeGlobalRegionSegment( itk::Size<3>, int, int, int );
	rawImageType_8bit::RegionType ComputeLocalRegionSplit( itk::Size<3>, int, int, int );
	rawImageType_8bit::RegionType ComputeGlobalRegionSplit( itk::Size<3>, int, int, int );
	
	
private:
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
	int _isSmall;
	std::string _segParams;
	std::string _projectDefinition;
	std::string _optionsMNT;
	std::string _outPath;
	std::string _outPathDebug;
	std::string _outPathDebugLevel2;
	std::string _outPathTemp;
	std::string _outPathData;
	
	std::string _Cy5_ImageNRRD;
	std::string _TRI_ImageNRRD;
	std::string _GFP_ImageNRRD;
	std::string _DAP_ImageNRRD;
	
// 	std::vector< rawImageType_8bit::Pointer > _ImageMontages;
	rawImageType_8bit::Pointer _ImageMontage_Cy5;
	rawImageType_8bit::Pointer _ImageMontage_TRI;
	rawImageType_8bit::Pointer _ImageMontage_GFP;
	rawImageType_8bit::Pointer _ImageMontage_DAP;
	
	itk::Size<3> _ImageMontage_Cy5Size;
	itk::Size<3> _ImageMontage_TRISize;
	itk::Size<3> _ImageMontage_GFPSize;
	itk::Size<3> _ImageMontage_DAPSize;
	
	// Image Split
	int _kx;
	int _ky;
	int _kz;
};

// #include "ftkMainDarpaSegment.hxx"

#endif