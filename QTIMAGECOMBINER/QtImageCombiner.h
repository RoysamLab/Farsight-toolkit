#include <QtCore/QFileInfo>
#include <QtCore/QThread>
#include <QtCore/QSettings>
#include <QtCore/QSignalMapper>
# include <QImage>
#include <QColor>
#include <QImageIOHandler>
#include <QImageWriter>
#include <vector>
#include <iostream>

class QtImageCombiner
{
	public:

		QImage* image1 ;
		QImage* image2 ;

		QtImageCombiner();
		QtImageCombiner(QImage* image1 , QImage* image2 );
		~QtImageCombiner();
		QImage*  CombineImages();
		
	private:

		
		int getpixelvalue(QRgb pixel);
};



