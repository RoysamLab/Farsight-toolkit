if(NOT DEFINED qt_dir)
  message(FATAL_ERROR "error: required variable qt_dir not defined...")
endif()

#
# Windows builds with nmake are "in-source/in-place" builds.
#   i.e. source == build == install dir
#   and prefix is irrelevant...
#
# Linux/Mac builds with make honor -prefix. They do not necessarily have to
# be "in-source/in-place" builds, although they can be... So we do in-source
# everywhere for consistency's sake.
#
if(WIN32)
  set(cmd configure)
  set(args "")
else()
  set(cmd ./configure)
  set(args -prefix ${qt_dir})
endif()

if(APPLE)
  set(args ${args} -cocoa -sdk /Developer/SDKs/MacOSX10.5.sdk -arch x86_64)
  set(ENV{MACOSX_DEPLOYMENT_TARGET} "10.5")
endif()

set(args ${args} -opensource -confirm-license -nomake demos -nomake examples)

message(STATUS "cmd='${cmd}'")
message(STATUS "args='${args}'")
message(STATUS "qt_dir='${qt_dir}'")

execute_process(COMMAND ${cmd} ${args}
  WORKING_DIRECTORY ${qt_dir}
  RESULT_VARIABLE rv
)

if(NOT "${rv}" STREQUAL "0")
  message(FATAL_ERROR "error: problem configuring Qt: rv='${rv}'")
endif()
