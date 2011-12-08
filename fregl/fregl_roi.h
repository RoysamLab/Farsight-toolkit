#include "itkImage.h"
//#include "fregl_image_manager.h"
class fregl_roi
{
public:
	typedef itk::Image< unsigned char, 3 >  ImageType;

	fregl_roi(std::string xml_file, std::string img_path, std::string anchor_image);
	~fregl_roi();

	void setROI(ImageType::IndexType roi_origin, ImageType::SizeType roi_size);
	void getMontage();
	
private:
	ImageType::IndexType roi_origin;
	ImageType::SizeType roi_size;

	std::string xml_file;
	std::string img_path;
	std::string anchor_image;
};

