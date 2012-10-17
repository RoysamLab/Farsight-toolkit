// ############################################################################################################################################################################
#ifndef _ftkMainDarpa_h_
#define _ftkMainDarpa_h_
// ############################################################################################################################################################################

#ifdef _OPENMP
#include "omp.h"
#endif

// GLOBAL INCLUDES
#include "ftkMainDarpaGlobalInclude.h"

// TEMPLATES
#include "ftkMainDarpaTemplates.h"

// DEFINITIONS
#include "ftkMainDarpaDeclaration.h"


// ############################################################################################################################################################################
// MAIN CLASS
class ftkMainDarpa
{
public:
	template<typename TINPUT, typename TOUTPUT >
	void projectImage( std::string, std::string, std::string, std::string );
	
	// To merge with the other one
	template<typename TINPUT, typename TOUTPUT >
	void projectImage( typename TINPUT::Pointer, std::string, std::string, std::string, std::string );

	// Assumed to be 8bit in general, for 16 bits some changes has to be made I thinkg
	template<typename TGFP, typename TSOMA >
	void projectImageRGB( std::string inputImageName1 /*GFP*/, std::string inputImageName2 /*Soma*/, std::string inputImageName3 /*Other */, std::string outputPath, std::string NAME );
	
	template<typename TINPUT, typename TOUTPUT >
	std::vector< typename TOUTPUT::Pointer > getProjectImage( typename TINPUT::Pointer inputImage, std::string projectOptions );

	template<typename TINPUT, typename TOUTPUT >
	void rescaleImage( std::string, std::string );

	template<typename TINPUT, typename TOUTPUT >
	void computeDistMap( std::string inputImageName, std::string outputPath, std::string imageType );
	
	template<typename TINPUT, typename TOUTPUT >
	void computeMedianFilter( std::string inputImageName, std::string outputImageName, std::string imageType );
	
	template<typename TINPUT >
	void saveNRRD( std::string );

	template<typename TINPUT >
	void cropImageDarpa( std::string, std::string, std::string );

	
private:
	
};

#include "ftkMainDarpa.hxx"

#endif