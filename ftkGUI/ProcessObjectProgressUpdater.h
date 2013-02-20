#ifndef PROCESSOBJECTPROGRESSUPDATER_H
#define PROCESSOBJECTPROGRESSUPDATER_H

#include <QProgressBar>
#include <QLabel>
#include "itkCommand.h"

/** Update a QProgressBar when an itk::ProcessObject is running GenerateData().
 * */
class ProcessObjectProgressUpdater: public itk::Command
{
public:
	typedef ProcessObjectProgressUpdater    Self;
	typedef itk::Command                    Superclass;
	typedef itk::SmartPointer< Self >       Pointer;
	typedef itk::SmartPointer< const Self > ConstPointer;

	/** Creates ::New() */
	itkNewMacro( Self );

	/** Set/Get the QProgressBar used to display the progress. */
	QProgressBar * GetProgressBar();
	void SetProgressBar( QProgressBar * );

	/** Set/Get the Qt widget used to display a textual description of the
	 * process. */
	QLabel * GetTextWidget();
	void SetTextWidget( QLabel * widget );

	/** Set/Get the process description. */
	std::string GetDescription() const;
	void SetDescription( const std::string & description );
	
	/** Set/Get whether we should display the percent completion text in
	 * the progress bar.  Default to true.  May want to set this to false if
	 * no itk::ProgressEvents will be emitted. */
	void SetDisplayProgressBarText( bool display );
	bool GetDisplayProgressBarText() const;
	
	void Execute( itk::Object * caller, const itk::EventObject & event );
	void Execute( const itk::Object * caller, const itk::EventObject & event );

private:
	ProcessObjectProgressUpdater();

	QProgressBar * m_ProgressBar;
	QLabel *       m_TextWidget;
	std::string    m_Description;
	float          m_LastProgress;
	bool           m_DisplayProgressBarText;
};

#endif
