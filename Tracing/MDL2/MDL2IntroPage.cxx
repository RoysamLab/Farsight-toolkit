#include "MDL2IntroPage.h"
#include "MDL2Wizard.h"

#include <iostream>
using std::cout;
using std::endl;

MDL2IntroPage::MDL2IntroPage(QWidget *parent)
{
}

bool MDL2IntroPage::isComplete() const
{
  MDL2Wizard *wiz = static_cast<MDL2Wizard*>(this->wizard());
  if(wiz->InputImageLabel->text() != "" &&
     wiz->BackboneOutputLabel->text() != "" &&
     wiz->SkeletonOutputLabel->text() != "")
    {
    return true;
    }
  return false;
}

void MDL2IntroPage::CheckIfComplete()
{
  emit this->completeChanged();
}

