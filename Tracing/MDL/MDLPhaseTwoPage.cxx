#include "MDLPhaseTwoPage.h"
#include "MDLWizard.h"

MDLPhaseTwoPage::MDLPhaseTwoPage(QWidget *parent)
{
}

bool MDLPhaseTwoPage::isComplete() const
{
  MDLWizard *wiz = static_cast<MDLWizard*>(this->wizard());
  if(wiz->PhaseTwoDone)
    {
    return true;
    }
  return false;
}

void MDLPhaseTwoPage::CheckIfComplete()
{
  emit this->completeChanged();
}

