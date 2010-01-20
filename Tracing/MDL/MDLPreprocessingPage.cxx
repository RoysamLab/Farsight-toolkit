#include "MDLPreprocessingPage.h"
#include "MDLWizard.h"

MDLPreprocessingPage::MDLPreprocessingPage(QWidget *parent)
{
}

bool MDLPreprocessingPage::isComplete() const
{
  MDLWizard *wiz = static_cast<MDLWizard*>(this->wizard());
  if(wiz->PreprocessingDone)
    {
    return true;
    }
  return false;
}

void MDLPreprocessingPage::CheckIfComplete()
{
  emit this->completeChanged();
}

