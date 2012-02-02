#include "GUITester.h"
#include "GUITesterHelper.h"

#include "pqTestUtility.h"
#include "pqEventObserver.h"
#include "pqEventSource.h"

#include <QFileDialog>
#include <QTextStream>
#include <QXmlStreamAttributes>
#include <QXmlStreamReader>
#include <QXmlStreamWriter>
#include <QtDebug>

#include "vtkRenderWindow.h"
#include "vtkTesting.h"

#include "itkTestDriverInclude.h"

class XMLEventObserver : public pqEventObserver
{
  QXmlStreamWriter* XMLStream;
  QString XMLString;

public:
  XMLEventObserver(QObject* p) : pqEventObserver(p)
  {
  this->XMLStream = NULL;
  }

  ~XMLEventObserver()
    {
    delete this->XMLStream;
    }

protected:
  virtual void setStream(QTextStream* stream)
    {
    if (this->XMLStream)
      {
      this->XMLStream->writeEndElement();
      this->XMLStream->writeEndDocument();
      delete this->XMLStream;
      this->XMLStream = NULL;
      }
    if (this->Stream)
      {
      *this->Stream << this->XMLString;
      }
    this->XMLString = QString();
    pqEventObserver::setStream(stream);
    if (this->Stream)
      {
      this->XMLStream = new QXmlStreamWriter(&this->XMLString);
      this->XMLStream->setAutoFormatting(true);
      this->XMLStream->writeStartDocument();
      this->XMLStream->writeStartElement("events");
      }
    }

  virtual void onRecordEvent(const QString& widget, const QString& command,
    const QString& arguments)
    {
    if(this->XMLStream)
      {
      this->XMLStream->writeStartElement("event");
      this->XMLStream->writeAttribute("widget", widget);
      this->XMLStream->writeAttribute("command", command);
      this->XMLStream->writeAttribute("arguments", arguments);
      this->XMLStream->writeEndElement();
      }
    }
};

class XMLEventSource : public pqEventSource
{
  typedef pqEventSource Superclass;
  QXmlStreamReader *XMLStream;

public:
  XMLEventSource(QObject* p): Superclass(p) { this->XMLStream = NULL;}
  ~XMLEventSource() { delete this->XMLStream; }

protected:
  virtual void setContent(const QString& xmlfilename)
    {
    delete this->XMLStream;
    this->XMLStream = NULL;

    QFile xml(xmlfilename);
    if (!xml.open(QIODevice::ReadOnly))
      {
      qDebug() << "Failed to load " << xmlfilename;
      return;
      }
    QByteArray data = xml.readAll();
    this->XMLStream = new QXmlStreamReader(data);
    /* This checked for valid event objects, but also caused the first event
     * to get dropped. Commenting this out in the example. If you wish to report
     * empty XML test files a flag indicating whether valid events were found is
     * probably the best way to go.
    while (!this->XMLStream->atEnd())
      {
      QXmlStreamReader::TokenType token = this->XMLStream->readNext();
      if (token == QXmlStreamReader::StartElement)
        {
        if (this->XMLStream->name() == "event")
          {
          break;
          }
        }
      } */
    if (this->XMLStream->atEnd())
      {
      qDebug() << "Invalid xml" << endl;
      }
    }

  int getNextEvent(QString& widget, QString& command, QString&
    arguments)
    {
    if (this->XMLStream->atEnd())
      {
      return DONE;
      }
    while (!this->XMLStream->atEnd())
      {
      QXmlStreamReader::TokenType token = this->XMLStream->readNext();
      if (token == QXmlStreamReader::StartElement)
        {
        if (this->XMLStream->name() == "event")
          {
          break;
          }
        }
      }
    if (this->XMLStream->atEnd())
      {
      return DONE;
      }
    widget = this->XMLStream->attributes().value("widget").toString();
    command = this->XMLStream->attributes().value("command").toString();
    arguments = this->XMLStream->attributes().value("arguments").toString();
    return SUCCESS;
    }
};

//-----------------------------------------------------------------------------
GUITester::GUITester(QWidget *parent) : QWidget(parent)
{
  this->TestUtility = new pqTestUtility(this);
  this->TestUtility->addEventObserver("xml", new XMLEventObserver(this));
  this->TestUtility->addEventSource("xml", new XMLEventSource(this));
  
  this->Testing = vtkSmartPointer< vtkTesting >::New();
  this->Testing->SetValidImageFileName(NULL);
  this->Testing->SetRenderWindow(NULL);
  this->Testing->AddArgument("-T");
  this->Testing->AddArgument(TEST_OUTPUT_DIR);

  this->Threshold = 100.0;
  this->BaselineSet = false;
}

//-----------------------------------------------------------------------------
GUITester::~GUITester()
{
  delete this->TestUtility;
}

//-----------------------------------------------------------------------------
void GUITester::record()
{
  QString filename = QFileDialog::getSaveFileName (this, "Test File Name",
    QString(), "XML Files (*.xml)");
  if (!filename.isEmpty())
    {
    this->TestUtility->recordTests(filename);
    }
}

//-----------------------------------------------------------------------------
void GUITester::play()
{
  QString filename = QFileDialog::getOpenFileName (this, "Test File Name",
    QString(), "XML Files (*.xml)");
  if (!filename.isEmpty())
    {
    this->playTestFile(filename);
    }
}

//-----------------------------------------------------------------------------
bool GUITester::playTestAndCompareResults( QString filename )
{
  this->playTestFile(filename);
  return this->compareResults();
}

//-----------------------------------------------------------------------------
void GUITester::playTestFile( QString filename )
{
  this->TestUtility->playTests(filename);
}

//-----------------------------------------------------------------------------
bool GUITester::compareResults()
{
  if(this->BaselineSet && this->Testing->GetRenderWindow() != NULL)
    {
    int res = this->Testing->RegressionTest( this->Threshold );
    if(res == vtkTesting::PASSED )
      {
      return true;
      }
    }
  return false;
}

//-----------------------------------------------------------------------------
bool GUITester::compareResults( QString testImgFileName )
{
  if(this->BaselineSet)
    {
    int res = RegressionTestImage(
      testImgFileName.toStdString().c_str(),  
      this->BaselineImage.toStdString().c_str(),
      1, 0.0, (unsigned int)this->Threshold, 0);
    if(res == 0 )
      {
      return true;
      }
    }
  return false;
}

//-----------------------------------------------------------------------------
void GUITester::SetBaselineImage(const char *fn)
{
  this->Testing->AddArgument("-V");
  this->Testing->AddArgument(fn);
  this->BaselineImage = fn;
  this->BaselineSet = true;
}

//-----------------------------------------------------------------------------
void GUITester::SetRenderWindow(vtkRenderWindow *rw)
{
  this->Testing->SetRenderWindow(rw);
}

//-----------------------------------------------------------------------------
void GUITester::SetThreshold(double t)
{
  this->Threshold = t;
}
