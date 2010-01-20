#ifndef MDLINTROPAGE_H
#define MDLINTROPAGE_H
#include <QtGui>

class MDLIntroPage : public QWizardPage
{
  Q_OBJECT

  public:
  MDLIntroPage(QWidget *parent = 0);
  bool isComplete() const;
  void CheckIfComplete();

  private:
};
#endif

