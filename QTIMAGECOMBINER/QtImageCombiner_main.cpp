
#include "QtImageCombiner.h"



int main(int argc, char * argv[])
{

	argv[1]= "C:\\SACHIN\\Farsight_bin\\exe\\Release\\4891.tif";
	argv[2]= "C:\\SACHIN\\Farsight_bin\\exe\\Release\\ho.tif";

	//const char *fileName, const char *format = 0

	QImage image1 = QImage ( argv[1],  0 );
	QImage image2 = QImage ( argv[2],  0 );
	QtImageCombiner* mycombiner = new QtImageCombiner( &image1,&image2);
	QImage* image = mycombiner->CombineImages();
	
	QImageWriter writer;	
	writer.setFileName ( "sach.tif" );
	writer.write (*image );
	writer.error () ;
	return 0;

}
