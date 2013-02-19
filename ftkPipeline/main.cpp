

#include "ftkMainDarpaGlobalInclude.h"

#include "ftkMainDarpa.h"
#include "ftkMainDarpaSegment.h"
#include "ftkMainDarpaTrace.h"
#include "ftkMainDarpaAstroTrace.h"

enum STEPS { CROPIMAGE, SAVENRRD, PROJECTION, PROJECTION_8BIT, PROJECTIONRGB, PROJECTIONFLO, MEDIAN, RESCALE, RESCALE_8BIT, DISTANCE_MAP, SEGMENT, TRACE, ASTRO_TRACE};
std::map< std::string, STEPS> stepsmap;

void register_stepsmap()
{
  stepsmap["CROPIMAGE"] = CROPIMAGE;
  stepsmap["SAVENRRD"] = SAVENRRD;
  stepsmap["PROJECTION"] = PROJECTION;
  stepsmap["PROJECTION_8BIT"] = PROJECTION_8BIT;
  stepsmap["PROJECTIONRGB"] = PROJECTIONRGB;
  stepsmap["PROJECTIONFLO"] = PROJECTIONFLO;
  stepsmap["MEDIAN"] = MEDIAN;
  stepsmap["RESCALE"] = RESCALE;
  stepsmap["RESCALE_8BIT"] = RESCALE_8BIT;
  stepsmap["DISTANCE_MAP"] = DISTANCE_MAP;
  stepsmap["SEGMENT"] = SEGMENT;
  stepsmap["TRACE"] = TRACE;
  stepsmap["ASTRO_TRACE"] = ASTRO_TRACE;
}

int main(int argc, char *argv[])
{
  ftkMainDarpa *objftkMainDarpa;
  objftkMainDarpa = new ftkMainDarpa;
  // 	ftkMainDarpaSegment *objftkMainDarpaSegment_1;
  // 	objftkMainDarpaSegment_1 = new ftkMainDarpaSegment;

  register_stepsmap();
  switch( stepsmap[argv[1]] )
  {
    case CROPIMAGE:
      {
        std::cout << "CROPIMAGE!";
        std::string imageInput = argv[2];
        std::string tableInput = argv[3];
        std::string coordinate = argv[4];

        objftkMainDarpa->cropImageDarpa< rawImageType_uint >( imageInput, tableInput, coordinate );
        break;
      }
    case SAVENRRD:
      {
        std::cout << "SAVENRRD!";
        std::string imageInputName = argv[2];

        objftkMainDarpa->saveNRRD< rawImageType_16bit >( imageInputName );
        break;
      }
    case PROJECTION:
      {
        std::cout << "PROJECTION!";
        std::string imageInputName = argv[2];
        std::string projectImagePath = argv[3];
        std::string projectOptions = argv[4];
        std::string imageType = argv[5];

        objftkMainDarpa->projectImage<rawImageType_16bit, rawImageType_16bit>( imageInputName, projectImagePath, projectOptions, imageType );
        break;
      }

    case PROJECTION_8BIT:
      {
        std::cout << "PROJECTION_8BIT!";
        std::string imageInputName = argv[2];
        std::string projectImagePath = argv[3];
        std::string projectOptions = argv[4];
        std::string imageType = argv[5];

        objftkMainDarpa->projectImage<rawImageType_8bit, rawImageType_8bit>( imageInputName, projectImagePath, projectOptions, imageType );
        break;
      }

    case PROJECTIONRGB:
      {
        std::cout << "PROJECTIONRGB!";
        std::string imageInputGFPName = argv[2];
        std::string imageInputSOMAName = argv[3];
        std::string projectImagePath = argv[4];
        std::string name = argv[5];
        std::string thirdChannel = "no3chanel";


        objftkMainDarpa->projectImageRGB<rawImageType_8bit, rawImageType_uint>( imageInputGFPName, imageInputSOMAName, thirdChannel, projectImagePath, name );
        break;
      }

    case PROJECTIONFLO:
      {
        std::cout << "PROJECTIONFLO!";
        std::string imageInputName = argv[2];
        std::string projectImagePath = argv[3];
        std::string projectOptions = argv[4];
        std::string imageType = argv[5];


        objftkMainDarpa->projectImage<rawImageType_flo, rawImageType_flo>( imageInputName, projectImagePath, projectOptions, imageType );
        break;
      }

    case MEDIAN:
      {
        std::cout << "MEDIAN!";
        std::string imageInputName = argv[2];
        std::string imageOutputName = argv[3];
        std::string imageType = argv[4];

        objftkMainDarpa->computeMedianFilter<rawImageType_16bit, rawImageType_16bit>( imageInputName, imageOutputName, imageType );
        break;
      }

    case RESCALE:
      {
        std::cout << "RESCALE!";
        std::string imageInputName = argv[2];
        std::string imageOutputName = argv[3];

        objftkMainDarpa->rescaleImage<rawImageType_16bit, rawImageType_16bit>( imageInputName, imageOutputName );
        break;
      }

    case RESCALE_8BIT:
      {
        std::cout << "RESCALE_8BIT!";
        std::string imageInputName = argv[2];
        std::string imageOutputName = argv[3];

        objftkMainDarpa->rescaleImage<rawImageType_16bit, rawImageType_8bit>( imageInputName, imageOutputName );
        break;
      }

    case DISTANCE_MAP:
      {
        std::cout << "DISTANCE_MAP!";
        std::string imageInputName = argv[2];
        std::string imageOutputName = argv[3];
        std::string imageType = argv[4];

        objftkMainDarpa->computeDistMap<rawImageType_uint, rawImageType_flo>( imageInputName, imageOutputName, imageType );
        break;
      }

    case SEGMENT:
      {
        ftkMainDarpaSegment *objftkMainDarpaSegment_1;
        objftkMainDarpaSegment_1 = new ftkMainDarpaSegment;
        std::cout << "SEGMENT!";
        std::string segmentParams = argv[2];

        objftkMainDarpaSegment_1->readParameters( segmentParams );
        objftkMainDarpaSegment_1->runSpliting();
        objftkMainDarpaSegment_1->runSegment();
        objftkMainDarpaSegment_1->runStich();
        delete objftkMainDarpaSegment_1;
        break;
      }

    case TRACE:
      {
        ftkMainDarpaTrace *objftkMainDarpaTrace_1;
        objftkMainDarpaTrace_1 = new ftkMainDarpaTrace;
        std::cout << "TRACE!";
        std::string segmentParams = argv[2];

        objftkMainDarpaTrace_1->readParameters( segmentParams );
        objftkMainDarpaTrace_1->runPreprocesing();
        //objftkMainDarpaTrace_1->computeTileGVFAndVesselness();
        objftkMainDarpaTrace_1->runTracing();
        delete objftkMainDarpaTrace_1;
        break;
      }

    case ASTRO_TRACE:
      {
        ftkMainDarpaAstroTrace *objftkMainDarpaAstroTrace_1;
        objftkMainDarpaAstroTrace_1 = new ftkMainDarpaAstroTrace;
        std::cout << "ASTRO_TRACE!";
        std::string astroTraceParams = argv[2];

        objftkMainDarpaAstroTrace_1->readParameters( astroTraceParams );
        objftkMainDarpaAstroTrace_1->runPreprocesing();
        objftkMainDarpaAstroTrace_1->runSplitting();
        objftkMainDarpaAstroTrace_1->runInterestPoints();
        objftkMainDarpaAstroTrace_1->runStitchRoots();
        objftkMainDarpaAstroTrace_1->computeRootFeaturesForNuclei();
        delete objftkMainDarpaAstroTrace_1;
        break;
      }

  }
  return 0;
}








