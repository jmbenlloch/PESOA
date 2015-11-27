# -*- cmake -*-
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(GATE)

MESSAGE("CMAKE_SYSTEM_NAME: ${CMAKE_SYSTEM_NAME}")
#Set compilier flags
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
   # Mac OS X specific code
   SET(CMAKE_CXX_FLAGS "-Wall -g -O2 -std=c++11" )
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    # Linux specific code
    SET(CMAKE_CXX_FLAGS "-Wall -g -O2 -std=c++0x" )
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")


SET(CMAKE_VERBOSE_MAKEFILE ON)
#SET(CMAKE_BUILD_TYPE ${CMAKE_CXX_FLAGS})

MESSAGE("CMAKE source dir: ${CMAKE_SOURCE_DIR}")
MESSAGE("CMAKE binary dir: ${CMAKE_BINARY_DIR}")

MESSAGE("CMAKE compiler: ${CMAKE_CXX_COMPILER}")
MESSAGE("CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
MESSAGE("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# Set where to look for external packages 
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/") 

MESSAGE("CMAKE_MODULE_PATH source dir: ${CMAKE_MODULE_PATH}")


# Look for the needed libraries

find_package(ROOT REQUIRED)

MESSAGE("ROOT include directories: ${ROOT_INCLUDE_DIR}")


MESSAGE("ROOT libraries: ${ROOT_LIBRARIES}")

SET(EMODEL ${CMAKE_SOURCE_DIR}/emodel) # Define the directories
SET(IO ${CMAKE_SOURCE_DIR}/IO) # Define the directories
SET(UTILS ${CMAKE_SOURCE_DIR}/UTILS) # Define the directories

#Make symbolic links in GATE/ to all header files
EXECUTE_PROCESS(
  COMMAND mkdir -p GATE)

FILE(GLOB gate_headers . ${EMODEL}/*.h ${IO}/*.h ${UTILS}/*.h)

FOREACH(file ${gate_headers})
  EXECUTE_PROCESS(
    COMMAND ln -s -f ${file} GATE)
ENDFOREACH(file)

# Directories that will be searched for headers during compilation
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR} ${EMODEL}
					${IO} ${UTILS} ${ROOT_INCLUDE_DIR})


set(EMODEL_DICTIONARY ${CMAKE_BINARY_DIR}/emodelDict.cxx)  # set dictionary path

MESSAGE("Generate emodel dictionary")
FILE(GLOB EMODEL_HEADERS . ${EMODEL}/*.h)

set(IO_DICTIONARY ${CMAKE_BINARY_DIR}/ioDict.cxx)  # set dictionary path

MESSAGE("Generate io dictionary")
FILE(GLOB IO_HEADERS . ${IO}/*.h)

set(UTILS_DICTIONARY ${CMAKE_BINARY_DIR}/utilsDict.cxx)  # set dictionary path

MESSAGE("Generate utils dictionary")
FILE(GLOB UTILS_HEADERS . ${UTILS}/*.h)

#SET(EMODEL_HEADERS ${EMODEL}/*.h)

ROOT_GENERATE_DICTIONARY(
           "${EMODEL_DICTIONARY}"            # path to dictionary to generate
           "${EMODEL}"  # list of includes
           "${EMODEL_HEADERS}" 
           "${EMODEL}/LinkDef.hh"   # ROOT linkDef file
         )

ROOT_GENERATE_DICTIONARY(
           "${IO_DICTIONARY}"            # path to dictionary to generate
           "${IO}"  # list of includes
           "${IO_HEADERS}" 
           "${IO}/LinkDef.hh"   # ROOT linkDef file
         )

ROOT_GENERATE_DICTIONARY(
           "${UTILS_DICTIONARY}"            # path to dictionary to generate
           "${UTILS}"  # list of includes
           "${UTILS_HEADERS}" 
           "${UTILS}/LinkDef.hh"   # ROOT linkDef file
         )

FILE(GLOB gate_src . ${EMODEL}/*.cc ${IO}/*.cc ${UTILS}/*.cc)


# # # # Finally, create it!

MESSAGE("Creating GATE library")

ADD_LIBRARY(GATE SHARED ${gate_src} ${EMODEL_DICTIONARY} ${IO_DICTIONARY} 
  ${UTILS_DICTIONARY})

SET_TARGET_PROPERTIES(GATE PROPERTIES VERSION 1.1.0 SOVERSION 1)
TARGET_LINK_LIBRARIES(GATE ${ROOT_LIBRARIES})

INSTALL(TARGETS GATE DESTINATION lib)
INSTALL(FILES ${gate_headers} DESTINATION include/GATE)
INSTALL(FILES ${CMAKE_BINARY_DIR}/emodelDict_rdict.pcm DESTINATION lib)
INSTALL(FILES ${CMAKE_BINARY_DIR}/ioDict_rdict.pcm DESTINATION lib)
INSTALL(FILES ${CMAKE_BINARY_DIR}/utilsDict_rdict.pcm DESTINATION lib)


# # Installation destinations (just for the lib, not the python module)
#INSTALL(TARGETS GATE LIBRARY DESTINATION lib)
#INSTALL(FILES ${GATE} DESTINATION include/GATE)

# # Get ready for creating a package
 #INCLUDE(CPack)
