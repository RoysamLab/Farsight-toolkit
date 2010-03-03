#include "MDL2PhaseOnePage.h"
#include "MDL2Wizard.h"

MDL2PhaseOnePage::MDL2PhaseOnePage(QWidget *parent)
{
}

bool MDL2PhaseOnePage::isComplete() const
{
  MDL2Wizard *wiz = static_cast<MDL2Wizard*>(this->wizard());
  if(wiz->PhaseOneDone)
    {
    return true;
    }
  return false;
}

void MDL2PhaseOnePage::CheckIfComplete()
{
  emit this->completeChanged();
}

