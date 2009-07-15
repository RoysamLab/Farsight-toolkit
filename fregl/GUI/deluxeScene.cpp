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

#include "deluxescene.h"
#include <QRectF>
#include <iostream>

DeluxeScene::DeluxeScene(QObject *parent)
:QGraphicsScene(parent)
{
  lastBox = NULL;
  scaling = 1;
}

void DeluxeScene::mousePressEvent(QGraphicsSceneMouseEvent *mEvent)
{
  // Update the color
  for (int i = 0; i<myBoxes.size(); i++) {
    QPen pen = myBoxes[i]->pen();
    pen.setBrush(Qt::yellow);
    myBoxes[i]->setPen(pen);
  }

  // Add the new item
  lastMousePressEvent = mEvent->scenePos();
  qreal x = lastMousePressEvent.x();
  qreal y = lastMousePressEvent.y();
  lastBox = new QGraphicsRectItem(x,y,0,0);
  lastBox->setPen(QPen(Qt::magenta,2,Qt::SolidLine,Qt::RoundCap,Qt::RoundJoin));
  lastBox->setZValue(1);
  this->addItem(lastBox);
}

void DeluxeScene::mouseReleaseEvent(QGraphicsSceneMouseEvent *mEvent)
{
  lastMouseReleaseEvent = mEvent->scenePos();
  
  qreal x1 = lastMousePressEvent.x();
  qreal y1 = lastMousePressEvent.y();
  qreal x2 = lastMouseReleaseEvent.x();
  qreal y2 = lastMouseReleaseEvent.y();
  
  if ( x1 == x2 && y1 == y2 )
    {
      emit mouseClicked();
    }
  else
    {
      emit mouseDragged();
      /*
      QRectF rect;
      // scale to the original image space
      rect.setLeft(lastBox->rect().left()/scaling);
      rect.setRight(lastBox->rect().right()/scaling);
      rect.setTop(lastBox->rect().top()/scaling);
      rect.setBottom(lastBox->rect().bottom()/scaling);
      lastBox->setRect( rect );
      */
      myBoxes.push_back( lastBox );
      zooms.push_back(scaling);
    }
}

void DeluxeScene::mouseMoveEvent(QGraphicsSceneMouseEvent *mEvent)
{  
  qreal x1 = lastMousePressEvent.x();
  qreal y1 = lastMousePressEvent.y();
  qreal x2 = mEvent->scenePos().x();
  qreal y2 = mEvent->scenePos().y();
  lastBox->setRect(x1,y1,x2-x1,y2-y1);
  lastBox->setVisible(true);
  lastBox->show();
  this->update();
}


