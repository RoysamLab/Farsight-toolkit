#ifndef DELUXESCENE_H
#define DELUXESCENE_H

#include <QWidget>
#include <QGraphicsScene>
#include <QGraphicsPixmapItem>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QWidget>
#include <QObject>
#include <QString>
#include <QLabel>
#include <QFrame>
#include <QtGui>
#include <QWidget>
#include <QVector>

class DeluxeScene : public QGraphicsScene
{
  Q_OBJECT;
  
public:
  DeluxeScene( QObject *parent = 0);
  
  ~DeluxeScene(){};

protected:
  void mousePressEvent(QGraphicsSceneMouseEvent *mEvent);
  void mouseReleaseEvent(QGraphicsSceneMouseEvent *mEvent);
  void mouseMoveEvent(QGraphicsSceneMouseEvent *mEvent);
  
signals:
  void mouseClicked(void);
  void mouseDragged(void);
  void boxSelected(void);
  
public:
  QPointF lastMousePressEvent;
  QPointF lastMouseReleaseEvent;
  QVector<QGraphicsRectItem*> myBoxes;
  QVector<float> zooms;
  QGraphicsRectItem* lastBox; //the box last clicked
  float scaling;
};

#endif
