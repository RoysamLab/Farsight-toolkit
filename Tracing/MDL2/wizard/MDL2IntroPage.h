#ifndef MDL2INTROPAGE_H
#define MDL2INTROPAGE_H
#include <QtGui>

class MDL2IntroPage : public QWizardPage
{
  Q_OBJECT

  public:
  MDL2IntroPage(QWidget *parent = 0);
  bool isComplete() const;
  void CheckIfComplete();

  private:
};
#endif

