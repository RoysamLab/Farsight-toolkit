#ifndef MDL2PHASETWOPAGE_H
#define MDL2PHASETWOPAGE_H
#include <QtGui>

class MDL2PhaseTwoPage : public QWizardPage
{
  Q_OBJECT

  public:
  MDL2PhaseTwoPage(QWidget *parent = 0);
  bool isComplete() const;
  void CheckIfComplete();

  private:
};
#endif

