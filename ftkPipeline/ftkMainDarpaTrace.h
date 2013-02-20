// ############################################################################################################################################################################
#ifndef _ftkMainDarpaTrace_h_
#define _ftkMainDarpaTrace_h_
// ############################################################################################################################################################################

#ifdef _OPENMP
#include "omp.h"
#endif

// GLOBAL INCLUDES
#include "ftkMainDarpaGlobalInclude.h"
#include "ftkMainDarpa.h"
#include <boost/lexical_cast.hpp>

// TEMPLATES
#include "ftkMainDarpaTemplates.h"

// DEFINITIONS
#include "ftkMainDarpaDeclaration.h"

//MACROS
// #define MINNIC(a,b) (((a) > (b))? (b) : (a))
// #define MAXNIC(a,b) (((a) < (b))? (b) : (a))

// ############################################################################################################################################################################
// MAIN CLASS
class ftkMainDarpaTrace
{
  public:

    void readParameters( std::string );
    void runPreprocesing(); // Before preprocess always run redParams before
    void runTracing();
    void calcLMeasures(int argc, char *argv[]);
    float getCalcThreshold(std::vector<float> &features, std::string type);
    std::vector<float> computeFeatures(rawImageType_flo::Pointer &image);
    void computeTileGVFAndVesselness();



  protected:
    void computeSplitConst( );
    std::vector< itk::Index<3> > getCentroidList();
    std::vector< itk::Index<3> > getSomaTable( std::vector< itk::Index<3> > , int , int , int );
    template<typename TINPUT >
      typename TINPUT::Pointer cropImages( typename TINPUT::Pointer , int , int , int );
    void WriteCenterTrace(vtkSmartPointer< vtkTable > , int , int , int , std::string );

  private:
    std::string _segmentParams;
    int _xSize;
    int _ySize;
    int _zSize;
    int _xTile;
    int _yTile;
    int _zTile;
    // 	int _xTileBor;
    // 	int _yTileBor;
    // 	int _zTileBor;
    int _num_threads;
    std::string _Cy5_Image;
    std::string _TRI_Image;
    std::string _GFP_Image;
    std::string _DAP_Image;
    std::string _Soma_Centroids;
    std::string _Soma_Montage;
    int _isSmall;
    std::string _traceParams;
    // 	std::string _projectDefinition;
    // 	std::string _optionsMNT;
    std::string _outPath;
    std::string _outPathDebug;
    std::string _outPathDebugLevel2;
    std::string _outPathTemp;

    std::string _Cy5_ImageNRRD;
    std::string _TRI_ImageNRRD;
    std::string _GFP_ImageNRRD;
    std::string _DAP_ImageNRRD;
    std::string _Soma_MontageNRRD;
    std::string _overridedefaultsTraceParams;

    bool _optimizeCoverage;
    //
    // // 	std::vector< rawImageType_8bit::Pointer > _ImageMontages;
    // 	rawImageType_8bit::Pointer _ImageMontage_Cy5;
    // 	rawImageType_8bit::Pointer _ImageMontage_TRI;
    // 	rawImageType_8bit::Pointer _ImageMontage_GFP;
    // 	rawImageType_8bit::Pointer _ImageMontage_DAP;
    //
    // 	itk::Size<3> _ImageMontage_Cy5Size;
    // 	itk::Size<3> _ImageMontage_TRISize;
    // 	itk::Size<3> _ImageMontage_GFPSize;
    // 	itk::Size<3> _ImageMontage_DAPSize;
    //
    // 	// Image Split
    // 	int _kx;
    // 	int _ky;
    // 	int _kz;

    // MNT
    std::string _GFP_ImagePREPMNT;
    std::string _GVF_ImagePREMNT;
    std::string _Vesselness_ImagePREMNT;

    // Division
    int _numDivisionsInRowCEN;
    std::vector< itk::Index<3> > _initialBigTileLOG;
    std::vector< itk::Size<3> > _sizeOfBigTilesLOG;
    std::vector< itk::Index<3> > _initialBigTileCEN;
    std::vector< itk::Size<3> > _sizeOfBigTilesCEN;

    // BIG IMAGES
    rawImageType_flo::Pointer _img_traceDesiredRegion;
    rawImageType_uint::Pointer _somaMontageDesiredRegion;
    rawImageType_flo::Pointer _img_VesselDesiredRegion;
    GradientImageType::Pointer _img_GVFDesiredRegion;

};

#include "ftkMainDarpaTrace.hxx"

#endif
