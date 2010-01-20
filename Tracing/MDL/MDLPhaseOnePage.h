#ifndef MDLPHASEONEPAGE_H
#define MDLPHASEONEPAGE_H
#include <QtGui>

class MDLPhaseOnePage : public QWizardPage
{
  Q_OBJECT

  public:
  MDLPhaseOnePage(QWidget *parent = 0);
  bool isComplete() const;
  void CheckIfComplete();

  private:
};
#endif

