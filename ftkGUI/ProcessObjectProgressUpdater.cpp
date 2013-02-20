#include "ProcessObjectProgressUpdater.h"

#include "itkProcessObject.h"
#include "itkNumericTraits.h"

#include <QApplication>

ProcessObjectProgressUpdater::
ProcessObjectProgressUpdater():
  m_ProgressBar( NULL ),
  m_TextWidget( NULL ),
  m_DisplayProgressBarText( true )
{
}

QProgressBar *
ProcessObjectProgressUpdater::
GetProgressBar()
{
	return this->m_ProgressBar;
}

void
ProcessObjectProgressUpdater::
SetProgressBar( QProgressBar * progressBar )
{
	this->m_ProgressBar = progressBar;
}

QLabel *
ProcessObjectProgressUpdater::
GetTextWidget()
{
	return this->m_TextWidget;
}

void
ProcessObjectProgressUpdater::
SetTextWidget( QLabel * progressBar )
{
	this->m_TextWidget = progressBar;
}

std::string
ProcessObjectProgressUpdater::
GetDescription() const
{
	return this->m_Description;
}

void
ProcessObjectProgressUpdater::
SetDescription( const std::string & description )
{
	this->m_Description = description;
}

bool
ProcessObjectProgressUpdater::
GetDisplayProgressBarText() const
{
	return this->m_DisplayProgressBarText;
}

void
ProcessObjectProgressUpdater::
SetDisplayProgressBarText( bool display )
{
	this->m_DisplayProgressBarText = display;
}

void
ProcessObjectProgressUpdater::
Execute( itk::Object * caller, const itk::EventObject & event )
{
	this->Execute( const_cast< const itk::Object * >( caller ), event );
}

void
ProcessObjectProgressUpdater::
Execute( const itk::Object * caller, const itk::EventObject & event )
{
	if( ! this->m_ProgressBar )
		{
		itkExceptionMacro( << "ProcessObjectProgressUpdater: ProgressBar has not been set." );
		}
	if( ! this->m_TextWidget )
		{
		itkExceptionMacro( << "ProcessObjectProgressUpdater: TextWidget has not been set." );
		}

	if( itk::StartEvent().CheckEvent( &event ) )
		{
		const QString description( this->m_Description.c_str() );
		this->m_TextWidget->setText( description );
		this->m_ProgressBar->setRange( 0, itk::NumericTraits< int >::max() );
		this->m_LastProgress = 0.0f;
		this->m_ProgressBar->setValue( 0 );
		this->m_ProgressBar->setTextVisible( this->m_DisplayProgressBarText );
		QApplication::processEvents();
		}
	else if( itk::ProgressEvent().CheckEvent( &event ) )
		{
		const itk::ProcessObject * filter = dynamic_cast< const itk::ProcessObject * >( caller );
		if( ! filter )
		  {
	          itkExceptionMacro( << "ProcessObjectProgressUpdater only works with itk::ProcessObject's." );
		  }
		const float progress = filter->GetProgress();
		if( progress - m_LastProgress > 0.2f )
		  {
		  this->m_ProgressBar->setValue( static_cast< int >( progress * itk::NumericTraits< int >::max() ));
		  QApplication::processEvents();
		  this->m_LastProgress = progress;
		  }
		}
	else if( itk::EndEvent().CheckEvent( &event ) )
		{
		this->m_TextWidget->setText( "" );
		this->m_ProgressBar->setTextVisible( false );
		QApplication::processEvents();
		}
}
