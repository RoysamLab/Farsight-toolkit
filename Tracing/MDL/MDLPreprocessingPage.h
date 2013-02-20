#ifndef MDLPREPROCESSINGPAGE_H
#define MDLPREPROCESSINGPAGE_H
#include <QtGui>

class MDLPreprocessingPage : public QWizardPage
{
  Q_OBJECT

  public:
  MDLPreprocessingPage(QWidget *parent = 0);
  bool isComplete() const;
  void CheckIfComplete();

  private:
};
#endif

