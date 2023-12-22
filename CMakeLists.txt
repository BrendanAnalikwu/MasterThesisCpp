cmake_minimum_required(VERSION 3.5)

PROJECT(benchmark CXX C)


#### Set Standard options
SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Configuration type" )
SET(CMAKE_CXX_FLAGS_RELEASE "-O3" CACHE STRING "Configuration type" )

#### Include setup for all Documented Examples
INCLUDE(../CMakeLists_Examples.txt)


LINK_LIBRARIES(${GASCOIGNE_LIBRARIES})

#### INCLUDE DIRS
INCLUDE_DIRECTORIES(${GASCOIGNE_INCLUDE_DIRS})


#### Executables to be generated
ADD_EXECUTABLE("bench" loop.cc main.cc myequation.cc solver.cc)


