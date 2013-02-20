#include <string>
#include <vector>

class ImageTransformData
{
public:
	ImageTransformData(std::string fileName, double tX, double tY, double tZ)
	{
		this->fileName = fileName;
		this->tX = tX;
		this->tY = tY;
		this->tZ = tZ;
	}
	~ImageTransformData();	

	std::string getFileName()
	{
		return fileName;
	}
	
	double gettX()
	{
		return tX;
	}
	
	double gettY()
	{
		return tY;
	}
	
	double gettZ()
	{	
		return tZ;
	}

private:
	std::string fileName;
	double tX;
	double tY;
	double tZ;
};
