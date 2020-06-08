#find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)

# Add library locations to help find gtest
list(APPEND CMAKE_PREFIX_PATH "/usr/local/lib64/")
list(APPEND CMAKE_PREFIX_PATH "/usr/local/include/")

# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
include_directories (${TEST_SOURCE_DIR}../src) 

#ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK) 

set(LIBS hicupTest sysTests overlapTest importHicupTest binInterTest calcFreqTest pBhAdjustTest)
set(LIBS ${LIBS} sortPositionsTest countDupsTest removeDupsTest findCisTest binomialTest binomTestTest)
set(LIBS ${LIBS} findTransTest fixChromoTest halfInteractionTest hicupCompInterTest)
set(LIBS ${LIBS} sumSquareTest) # sumSquareCompTest)

foreach(lib ${LIBS})
  add_executable(${lib} ${lib}.cpp)
  target_link_libraries (${lib} hicupData BinomData Utils Interactions gtest gtest_main)
  string(REPLACE "Test" "_Test" name ${lib})
  message(STATUS ${name})
  add_test(${name} ${lib})
endforeach(lib)

#target_link_libraries(test
#                      ${Boost_FILESYSTEM_LIBRARY}
#                      ${Boost_SYSTEM_LIBRARY}
#                      ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
#                      )

#target_link_libraries(hicupTests hicupData BinomData Utils gtest gtest_main)
#target_link_libraries(overlaptests hicupData BinomData Utils gtest gtest_main)
#target_link_libraries(importHicupTests hicupData BinomData Utils gtest gtest_main)
#target_link_libraries(binInterTests hicupData BinomData Utils gtest gtest_main)
#target_link_libraries(sortPositionsTest hicupData BinomData Utils gtest gtest_main)
#target_link_libraries(countDupsTests hicupData BinomData Utils gtest gtest_main)
#target_link_libraries(systests gtest gtest_main)


#add_test(Sys_Tests systests)
#add_test(HiCUP_Tests hicupTests)
#add_test(Overlap_Tests overlaptests)
#add_test(importHiCUP_Tests importHicupTests)
#add_test(binInteractions_Tests binInterTests)
#add_test(sortPositions_Tests sortPositionsTest)
#add_test(countDuplicates_Tests countDupsTests)
#add_test(removeDuplicates_Tests removeDupsTest)
#add_test(findCis_Test findCisTest)
#add_test(findTrans_Test findTransTest)


