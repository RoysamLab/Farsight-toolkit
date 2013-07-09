SET( TracingSysLibraries
    ${ITK_LIBRARIES}
    rtvl
    rgtl
    ${VTK_LIBRARIES}
)

SET( GUI_RCS
  icons.qrc
)

QT4_ADD_RESOURCES( GUI_RC_SRCS ${GUI_RCS} )
QT4_WRAP_CPP(GUI_MOC_SRCS TracingGUI.h TracingView.h ObjectSelection.h ParametersGroup.h SlidersGroup.h OpenSnakeTracer.h dialogs.h dialogs_montage.h ScatterView.h dialogs_scatter.h)

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} )

SET(TracingCore
	TracingCore/PointOperation.h
    TracingCore/PointOperation.cpp
    TracingCore/ImageOperation.h
    TracingCore/ImageOperation.cpp
    TracingCore/SnakeOperation.h
    TracingCore/SnakeOperation.cpp
    TracingCore/Montage/Register.h
    TracingCore/Montage/Register.cpp
	OpenSnakeTracer.h
	OpenSnakeTracer.cpp
    )

ADD_LIBRARY( TraceCoreLibrary ${TracingCore} )

SET(TracingSystem_SOURCES
  main.cpp
  SlidersGroup.cpp
  ParametersGroup.cpp
  dialogs.cpp
  dialogs_montage.cpp
  dialogs_scatter.cpp
  TracingGUI.cpp
  TracingView.cpp
  ObjectSelection.cpp
  OpenSnakeTracer.cpp
  ScatterView.cpp
  Rendering/ImageActors.h
  Rendering/ImageActors.cxx
  ${GUI_MOC_SRCS}
  ${GUI_RC_SRCS}
  )

SET(TracingSystem_LIBRARIES
  ${QT_LIBRARIES} TraceCoreLibrary ftkCommon ${TracingSysLibraries} TinyXML)

add_farsight_executable(TracingSystem "${TracingSystem_SOURCES}" "${TracingSystem_LIBRARIES}")
