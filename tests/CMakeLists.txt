#find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)

# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
include_directories (${TEST_SOURCE_DIR}/../src) 


#ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK) 

add_executable (tests tests.cpp)
add_executable (systests sysTests.cpp) 
add_executable (overlaptests overlapTests.cpp)
add_executable (importHicupTests importHicupTests.cpp)
#add_executable (tests tests.cpp ${TEST_SOURCE_DIR}/../src/hicupData.cpp) 

#target_link_libraries(test
#                      ${Boost_FILESYSTEM_LIBRARY}
#                      ${Boost_SYSTEM_LIBRARY}
#                      ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
#                      )

target_link_libraries(tests hicupData BinomData Utils gtest gtest_main)
target_link_libraries(overlaptests hicupData BinomData Utils gtest gtest_main)
target_link_libraries(importHicupTests hicupData BinomData Utils gtest gtest_main)
target_link_libraries(systests gtest gtest_main)

add_test(Sys_Tests systests)
add_test(HiCUP_Tests tests)
add_test(Overlap_Tests overlaptests)
add_test(importHiCUP_Tests importHicupTests)

