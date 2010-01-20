#include "MDLIntroPage.h"
#include "MDLWizard.h"

#include <iostream>
using std::cout;
using std::endl;

MDLIntroPage::MDLIntroPage(QWidget *parent)
{
}

bool MDLIntroPage::isComplete() const
{
  MDLWizard *wiz = static_cast<MDLWizard*>(this->wizard());
  if(wiz->InputImageLabel->text() != "" &&
     wiz->BackboneOutputLabel->text() != "" &&
     wiz->SpinesOutputLabel->text() != "")
    {
    return true;
    }
  return false;
}

void MDLIntroPage::CheckIfComplete()
{
  emit this->completeChanged();
}

