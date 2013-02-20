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
#include "PreferencesDialog.h"

//Constructors:
PreferencesDialog::PreferencesDialog(QMap<QString, QColor> colorItemsMap, QWidget *parent)
: QDialog(parent)
{
	// LEFT THESE AS EXAMPLE
	//colorItemsMap["Selected Objects"] = Qt::yellow;
	//colorItemsMap["Object Boundaries"] = Qt::cyan;
	//colorItemsMap["Object IDs"] = Qt::green;
	//colorItemsMap["ROI Boundary"] = Qt::gray;

	m_colorItemsMap = colorItemsMap;

	signalMapper = new QSignalMapper(this);

	QGridLayout * layout = new QGridLayout;

	QMap<QString, QColor>::const_iterator it = m_colorItemsMap.constBegin();
	int i = 0;
	for ( ; it != m_colorItemsMap.constEnd(); it++ )
	{
		QLabel * label = new QLabel( it.key() );
		label->setAutoFillBackground(true);
		label->setAlignment(Qt::AlignCenter);
		label->setPalette(QPalette( it.value() ));
		label->setFrameStyle( QFrame::StyledPanel | QFrame::Plain );
		label->setMargin(1);
		colorLabelsMap[ it.key() ] = label;
		layout->addWidget(label, i, 0, 1, 5);

		QPushButton * button = new QPushButton(tr("Change"));
		connect(button, SIGNAL(clicked()), signalMapper, SLOT(map()));
		signalMapper->setMapping( button, it.key() );
		layout->addWidget(button, i, 5, 1, 1);
		++i;
	}
	connect(signalMapper, SIGNAL(mapped(const QString &)), this, SLOT(chooseColor(const QString &)));

	QHBoxLayout * okLayout = new QHBoxLayout();
	okLayout->addStretch(50);
	QPushButton * okButton = new QPushButton(tr("OK"),this);
	connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
	okLayout->addWidget(okButton);

	QVBoxLayout * masterLayout = new QVBoxLayout;
	masterLayout->addLayout(layout);
	masterLayout->addStretch(50);
	masterLayout->addLayout(okLayout);
	this->setLayout(masterLayout);
	this->setWindowTitle(tr("Choose Colors"));
}

void PreferencesDialog::chooseColor(const QString & colorItem)
{
	QLabel * label = colorLabelsMap.value(colorItem);
	QColor color = QColorDialog::getColor(label->palette().window().color(), this, colorItem, QColorDialog::ShowAlphaChannel);
    if (color.isValid())
	{
        label->setPalette(QPalette(color));
		m_colorItemsMap.insert(colorItem, color);
    }
}
