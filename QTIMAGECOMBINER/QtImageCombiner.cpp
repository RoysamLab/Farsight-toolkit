
#include "QtImageCombiner.h"





QtImageCombiner::QtImageCombiner(QImage* image1 , QImage* image2)
{
	this->image1 = image1;
	this->image2 = image2;
	

}

QtImageCombiner::QtImageCombiner()
{

}

QtImageCombiner::~QtImageCombiner()
{

}

int QtImageCombiner::getpixelvalue(QRgb pixel)
{
	
	QColor c = QColor:: fromRgb(pixel);

	if(c == Qt::white)
		return 255;

	else if(c == Qt::black)
		return 0;
}
	

QImage* QtImageCombiner::CombineImages()
{
	int width = image1->width();
	int height = image1->height();


	QImage* image = new QImage(width,height,QImage::Format_Indexed8);
	QRgb value;

	for(int i = 0;i<width;i++)
	{
		for(int j = 0;j<height ; j++)
		{
			int value1 = getpixelvalue(image1->pixel(i,j));
			int value2 = getpixelvalue(image2->pixel(i,j));

			value = qRgb(0,0,0);
			image->setColor(0, value);			
			value = qRgb(255,255,255);
			image->setColor(1, value);

			
			if((value1 !=0) || (value2 !=0))
				image->setPixel(i, j, 1);

			else 
				image->setPixel(i, j,0);

		}
	}
	return image;
}




