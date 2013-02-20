/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/
#include "GalleryDialog.h"

//Constructors:
GalleryDialog::GalleryDialog(std::vector<std::pair<QImage,int> > gallery,QVector<QColor> colorVec,QWidget *parent)
: QDialog(parent)
{
	
	int counter =0;
	int row =0;

	this->setWindowTitle(tr("Gallery"));
	this->setModal(false);

	//Master Layout
	QGridLayout * layout = new QGridLayout;

	//Top-row of the window 

	while(counter<gallery.size())
	{
		QHBoxLayout *topRow = new QHBoxLayout;
		for(int i=0; i<MIN(5,gallery.size());++i)
		{	
			if(counter==gallery.size())
				break;
			QImage cellImg = gallery[counter].first;
			QImage image(80, 80, QImage::Format_RGB32);

			image.fill(colorVec.at(gallery[counter].second - 1).rgb());

			//image.fill(colorVec.at(0).rgb());
			QPainter painter(&image);
			painter.setCompositionMode(QPainter::CompositionMode_Source);
			painter.drawImage(10, 10, cellImg.copy(20,20,60,60));

			QLabel *imageLabel = new QLabel(this);
			imageLabel->setPixmap(QPixmap::fromImage(image));
			topRow->addWidget(imageLabel,1,0);
			counter++;
		}
		layout->addLayout(topRow,row,0,0);
		row++;
	}
	
	// Info Label 
	QLabel *infoLabel;
	if(gallery.size() == 0)
		infoLabel = new QLabel("No actively labeled queries ", this);
	else
		infoLabel = new QLabel("The colors correspond to the classes (centroids) ", this);

		QHBoxLayout *topRow = new QHBoxLayout;
		topRow->addWidget(infoLabel,1,0);
		layout->addLayout(topRow,row,0,0);


	this->setLayout(layout);
}