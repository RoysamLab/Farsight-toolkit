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
#include "ExclusionDialog.h"

//Constructors:
ExclusionDialog::ExclusionDialog(QImage *pImage, QWidget *parent)
: QDialog(parent)
{
	masterLayout = new QVBoxLayout();
	QLabel * header = new QLabel(tr("Please set parameters for the Region of Interest:"));
	masterLayout->addWidget(header);

	spinLayout = new QVBoxLayout();
	this->addSpin(tr("Left Margin: "),0,10,1000,tr("pixels"));
	this->addSpin(tr("Right Margin: "),0,10,1000,tr("pixels"));
	this->addSpin(tr("Top Margin: "),0,10,1000,tr("pixels"));
	this->addSpin(tr("Bottom Margin: "),0,10,1000,tr("pixels"));
	this->addSpin(tr("Lower Z Margin: "),0,4,100,tr("slices"));
	this->addSpin(tr("Upper Z Margin: "),0,4,100,tr("slices"));
	masterLayout->addLayout(spinLayout);

	baseImage = pImage;
	preview = new QLabel();
	preview->setMaximumSize(256,256);
	preview->setScaledContents(true);	//Image will be stretched to fill label
	masterLayout->addWidget(preview);

	QHBoxLayout * okLayout = new QHBoxLayout();
	okLayout->addStretch(50);
	okButton = new QPushButton(tr("OK"),this);
	connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
	okLayout->addWidget(okButton);
	masterLayout->addLayout(okLayout);

	this->setLayout(masterLayout);
	this->setWindowTitle(tr("Apply Exclusion Margin"));

	this->updatePreview();
}

void ExclusionDialog::updatePreview()
{
	if(!baseImage)
	{
		preview->setText(tr("ROI Preview will be displayed here"));
		return;
	}

	QImage prevImg = baseImage->copy();

	int x1 = this->getMargin(0);
	int x2 = prevImg.rect().width() - this->getMargin(1) - x1;
	int y1 = this->getMargin(2);
	int y2 = prevImg.rect().height() - this->getMargin(3) - y1;

	QPainter painter(&prevImg);
	painter.setPen(Qt::gray);
	painter.drawRect(x1,y1,x2,y2);

	preview->setPixmap(QPixmap::fromImage(prevImg));
	preview->adjustSize();

}

int ExclusionDialog::getMargin(int x)
{
	if( x < (int)spins.size() )
		return spins.at(x)->value();
	else
		return 0;
}

void ExclusionDialog::addSpin(QString label, int min, int deflt, int max, QString units)
{
	QHBoxLayout * hLayout = new QHBoxLayout();
	QLabel * mLabel = new QLabel(label);
	QSpinBox * mSpin = new QSpinBox();
	mSpin->setMinimum(min);
	mSpin->setMaximum(max);
	mSpin->setValue(deflt);
	QLabel * unitLabel = new QLabel(units);

	spins.push_back(mSpin);
	//int row = spins.size();

	hLayout->addWidget(mLabel);
	hLayout->addWidget(mSpin);
	hLayout->addWidget(unitLabel);
	hLayout->addStretch(50);
	spinLayout->addLayout(hLayout);

	connect(mSpin, SIGNAL(valueChanged(int)), this, SLOT(updatePreview()));
}