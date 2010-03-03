#ifndef MDL2PHASEONEPAGE_H
#define MDL2PHASEONEPAGE_H
#include <QtGui>

class MDL2PhaseOnePage : public QWizardPage
{
  Q_OBJECT

  public:
  MDL2PhaseOnePage(QWidget *parent = 0);
  bool isComplete() const;
  void CheckIfComplete();

  private:
};
#endif

