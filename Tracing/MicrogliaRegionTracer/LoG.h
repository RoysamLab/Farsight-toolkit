#ifndef LoG_H
#define LoG_H

#include "Cell.h"

class LoG
{
private:
	typedef Cell::ImageType ImageType;
	typedef Cell::LoGImageType LoGImageType;

private:
    LoG();
	~LoG();
    
    static LoGImageType::Pointer RunLoG(const ImageType::Pointer & image, float scale);
    
public:
	static void WriteLoGImage(const std::string & filename, const LoGImageType::Pointer & image);
    static LoGImageType::Pointer RunMultiScaleLoG(const Cell & cell);
};

#endif
