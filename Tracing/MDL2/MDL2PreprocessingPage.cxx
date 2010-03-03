#include "MDL2PreprocessingPage.h"
#include "MDL2Wizard.h"

MDL2PreprocessingPage::MDL2PreprocessingPage(QWidget *parent)
{
}

bool MDL2PreprocessingPage::isComplete() const
{
  MDL2Wizard *wiz = static_cast<MDL2Wizard*>(this->wizard());
  if(wiz->PreprocessingDone)
    {
    return true;
    }
  return false;
}

void MDL2PreprocessingPage::CheckIfComplete()
{
  emit this->completeChanged();
}

