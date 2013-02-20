#include <QtGui/QApplication>
#include "ActiveValidationWizard.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    Active_Validation_Wizard avw;
    avw.show();

    return app.exec();
}