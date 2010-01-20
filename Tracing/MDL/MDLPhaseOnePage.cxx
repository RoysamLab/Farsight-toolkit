#include "MDLPhaseOnePage.h"
#include "MDLWizard.h"

MDLPhaseOnePage::MDLPhaseOnePage(QWidget *parent)
{
}

bool MDLPhaseOnePage::isComplete() const
{
  MDLWizard *wiz = static_cast<MDLWizard*>(this->wizard());
  if(wiz->PhaseOneDone)
    {
    return true;
    }
  return false;
}

void MDLPhaseOnePage::CheckIfComplete()
{
  emit this->completeChanged();
}

