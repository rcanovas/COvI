cmake_minimum_required(VERSION 3.5.1)
project(COvI)

set(CMAKE_BUILD_TYPE "Release")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#Add sdsl library
#change to the location of yours sdsl library
include_directories(/home/rcanovas/include)
link_directories(/home/rcanovas/lib)
link_libraries(sdsl divsufsort divsufsort64)

set(SOURCE_FILES2 tests/CreateCOvI.cpp include/covi/covi.h)
set(SOURCE_FILES3 tests/TestOperations.cpp include/covi/covi.h)


add_executable(CreateCOvI ${SOURCE_FILES2})
add_executable(TestOperations ${SOURCE_FILES3})
