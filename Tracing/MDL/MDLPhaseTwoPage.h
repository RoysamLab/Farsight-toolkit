#ifndef MDLPHASETWOPAGE_H
#define MDLPHASETWOPAGE_H
#include <QtGui>

class MDLPhaseTwoPage : public QWizardPage
{
  Q_OBJECT

  public:
  MDLPhaseTwoPage(QWidget *parent = 0);
  bool isComplete() const;
  void CheckIfComplete();

  private:
};
#endif

