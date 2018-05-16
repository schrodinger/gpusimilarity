function(add_boost_test SOURCE_FILE_NAME LINK_LIB BOOST_RUNTIME_PATH)
  get_filename_component(TEST_EXE ${SOURCE_FILE_NAME} NAME_WE)
  add_executable(${TEST_EXE} ${SOURCE_FILE_NAME})
  target_link_libraries(${TEST_EXE}
                        ${LINK_LIB} Boost::unit_test_framework)

  file(STRINGS "${SOURCE_FILE_NAME}" TEST_LIST REGEX "BOOST_AUTO_TEST_CASE\\( *([A-Za-z_0-9]+) *\\)")
  foreach(TEST_STRING ${TEST_LIST})
    string(REGEX REPLACE ".*\\( *([A-Za-z_0-9]+) *\\).*" "\\1" SUBTEST ${TEST_STRING})
    set(SUBTEST_NAME "${TEST_EXE}.${SUBTEST}")
    add_test(NAME "${SUBTEST_NAME}"
             COMMAND ${TEST_EXE}
             --run_test=${SUBTEST} --catch_system_error=yes)
    set_tests_properties("${SUBTEST_NAME}" PROPERTIES ENVIRONMENT "LD_LIBRARY_PATH=${BOOST_RUNTIME_PATH}:$ENV{LD_LIBRARY_PATH}")
  endforeach()
endfunction()
