#ifndef IMAGEWINDOW_H
#define IMAGEWINDOW_H

#include <QMainWindow>

#include "gnrcimagedialog.h"

class QAction;
class QActionGroup;
class QLabel;
class QMenu;
class QToolBar;
class QWorkspace;
class TransactionThread;


class MainWindow : public QMainWindow
{
  Q_OBJECT;
public:
  MainWindow();
  
protected:
  void closeEvent(QCloseEvent *event);
                                     
private slots:
  void    loadImage();
  void    saveImageAs();
  void    saveROIAs();
  void    deleteROI();
  void    defineROI();

private:
  void    createActions();
  void    createMenus();
  void    createStatusBar();
  
  QMenu   *fileMenu;
  QAction *openAction;
  QAction *saveAsAction;
  QAction *deleteAction;
  QAction *defineAction;
  QAction *exitAction;
  
  QAction *aboutAction;
  QLabel  *readyLabel;
  
  QString  fileFilters;
  
  QWorkspace   *workspace;
  QWidgetList   windows;
  QToolBar     *fileToolBar;
  QActionGroup *windowActionGroup;
  GNRCImageDialog *inDlg;

private: 
  GNRCImageDialog *createGNRCImage();

};

#endif
