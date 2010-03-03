#ifndef MDL2PREPROCESSINGPAGE_H
#define MDL2PREPROCESSINGPAGE_H
#include <QtGui>

class MDL2PreprocessingPage : public QWizardPage
{
  Q_OBJECT

  public:
  MDL2PreprocessingPage(QWidget *parent = 0);
  bool isComplete() const;
  void CheckIfComplete();

  private:
};
#endif

