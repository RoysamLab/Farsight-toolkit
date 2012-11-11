#ifndef LoG_H
#define LoG_H

#include "Cell.h"

class LoG
{
private:
	typedef Cell::ImageType ImageType;
	typedef Cell::LoGImageType LoGImageType;

public:
	explicit LoG();
	~LoG();

    LoGImageType::Pointer RunLoG(const ImageType::Pointer & image, float scale) const;
	void WriteLoGImage(const std::string & filename, const LoGImageType::Pointer & image) const;

    LoGImageType::Pointer RunMultiScaleLoG(const Cell & cell) const;
};

#endif
