INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )
# A program to run regression tests on Comma Separated Value (CSV) files.
# For usage information, do
#
#   CSVRegressionTest -h
ADD_EXECUTABLE( CSVRegressionTest CSVRegressionTest.cpp
  CSVRegressionTestCompareTestBaseline.cpp
  CSVRegressionTestparseArgs.cpp
  )
TARGET_LINK_LIBRARIES( CSVRegressionTest ${ITK_LIBRARIES} )

IF( NOT FARSIGHT_DATA_ROOT )
  MESSAGE( SEND_ERROR "The CSVRegressionTest requires FARSIGHT_DATA_ROOT
  to run its own tests." )
ENDIF()
SET( TEST_INPUT_DIR "${FARSIGHT_DATA_ROOT}/Baseline" )
SET( TEST_OUTPUT_DIR "${Farsight_BINARY_DIR}/Testing/Temporary" )

ADD_TEST( CSVRegressionTestTest ${EXE_DIR}/CSVRegressionTest
  -o "${TEST_OUTPUT_DIR}/CSVRegressionTestTest"
  ${TEST_INPUT_DIR}/CSVRegressionTestInput1Modified.csv
  ${TEST_INPUT_DIR}/CSVRegressionTestInput1.csv
  )
ADD_TEST( CSVRegressionTestAbsoluteToleranceTest ${EXE_DIR}/CSVRegressionTest
  -a 1.0
  -o "${TEST_OUTPUT_DIR}/CSVRegressionTestAbsoluteToleranceTest"
  ${TEST_INPUT_DIR}/CSVRegressionTestInput1Modified.csv
  ${TEST_INPUT_DIR}/CSVRegressionTestInput1.csv
  )
ADD_TEST( CSVRegressionTestFractionalToleranceTest ${EXE_DIR}/CSVRegressionTest
  -f 0.01
  -o "${TEST_OUTPUT_DIR}/CSVRegressionTestFractionalToleranceTest"
  ${TEST_INPUT_DIR}/CSVRegressionTestInput1Modified.csv
  ${TEST_INPUT_DIR}/CSVRegressionTestInput1.csv
  )
ADD_TEST( CSVRegressionTestStringErrorTest ${EXE_DIR}/CSVRegressionTest
  -o "${TEST_OUTPUT_DIR}/CSVRegressionTestStringErrorTest"
  ${TEST_INPUT_DIR}/CSVRegressionTestInput1StringError.csv
  ${TEST_INPUT_DIR}/CSVRegressionTestInput1.csv
  )
SET_TESTS_PROPERTIES( CSVRegressionTestAbsoluteToleranceTest
  CSVRegressionTestFractionalToleranceTest
  CSVRegressionTestStringErrorTest
  PROPERTIES WILL_FAIL TRUE )

ADD_TEST( CSVRegressionTestShowHelpTest ${EXE_DIR}/CSVRegressionTest -h )
ADD_TEST( CSVRegressionTestBadArgsTest  ${EXE_DIR}/CSVRegressionTest -c )
SET_TESTS_PROPERTIES( CSVRegressionTestBadArgsTest PROPERTIES WILL_FAIL TRUE )
ADD_TEST( CSVRegressionTestBadFilenameTest ${EXE_DIR}/CSVRegressionTest
  fileThatDoesntExist deliciousBritishFood
  )
SET_TESTS_PROPERTIES( CSVRegressionTestBadFilenameTest PROPERTIES WILL_FAIL TRUE )

ADD_TEST( CSVRegressionTestMissingRowTest ${EXE_DIR}/CSVRegressionTest
  -o "${TEST_OUTPUT_DIR}/CSVRegressionTestMissingRowTest"
  ${TEST_INPUT_DIR}/CSVRegressionTestInput1MissingRow.csv
  ${TEST_INPUT_DIR}/CSVRegressionTestInput1.csv
  )
ADD_TEST( CSVRegressionTestMissingColumnTest ${EXE_DIR}/CSVRegressionTest
  -o "${TEST_OUTPUT_DIR}/CSVRegressionTestMissingColumnTest"
  ${TEST_INPUT_DIR}/CSVRegressionTestInput1MissingColumn.csv
  ${TEST_INPUT_DIR}/CSVRegressionTestInput1.csv
  )
SET_TESTS_PROPERTIES( CSVRegressionTestMissingRowTest
  CSVRegressionTestMissingColumnTest
  PROPERTIES WILL_FAIL TRUE )

ADD_TEST( CSVRegressionTestNanTest ${EXE_DIR}/CSVRegressionTest
  -o "${TEST_OUTPUT_DIR}/CSVRegressionTestNanTest"
  ${TEST_INPUT_DIR}/CSVRegressionTestInput1Nan.csv
  ${TEST_INPUT_DIR}/CSVRegressionTestInput1.csv
  )
SET_TESTS_PROPERTIES( CSVRegressionTestNanTest
  PROPERTIES WILL_FAIL TRUE )

ADD_TEST( CSVRegressionTestTabDelimiterTest ${EXE_DIR}/CSVRegressionTest -d '\t'
  -o "${TEST_OUTPUT_DIR}/CSVRegressionTestTabDelimiterTest"
  ${TEST_INPUT_DIR}/CSVRegressionTestInput2Modified.csv
  ${TEST_INPUT_DIR}/CSVRegressionTestInput2.csv
  )
ADD_TEST( CSVRegressionTestExtraWhiteSpaceTest ${EXE_DIR}/CSVRegressionTest
  -o "${TEST_OUTPUT_DIR}/CSVRegressionTestExtraWhiteSpaceTest"
  ${TEST_INPUT_DIR}/CSVRegressionTestInput3Modified.csv
  ${TEST_INPUT_DIR}/CSVRegressionTestInput3.csv
  )
