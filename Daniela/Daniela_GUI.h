#include <QtGui/QMainWindow>
#include <QtGui/QMenubar>
#include "ftkGUI/LabelImageViewQT.h"

class Daniela_GUI : public QMainWindow
{
    Q_OBJECT;


private slots:
	void loadImage(void);
	void loadResult(void);


public:
	LabelImageViewQT *segView;
private:
	QSettings settings;
	ObjectSelection *selection;

public:
	Daniela_GUI(QWidget * parent = 0, Qt::WindowFlags flags = 0);

private:
	
	void createMenus(void);
   
};