#include "MDL2PhaseTwoPage.h"
#include "MDL2Wizard.h"

MDL2PhaseTwoPage::MDL2PhaseTwoPage(QWidget *parent)
{
}

bool MDL2PhaseTwoPage::isComplete() const
{
  MDL2Wizard *wiz = static_cast<MDL2Wizard*>(this->wizard());
  if(wiz->PhaseTwoDone)
    {
    return true;
    }
  return false;
}

void MDL2PhaseTwoPage::CheckIfComplete()
{
  emit this->completeChanged();
}

