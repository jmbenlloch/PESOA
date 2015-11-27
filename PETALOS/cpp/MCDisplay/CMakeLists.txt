#
SET(SRC ${CMAKE_SOURCE_DIR}/MCDisplay)


MESSAGE("Entering MCDisplay directory, SRC=${SRC}")

INCLUDE_DIRECTORIES(${SRC})


MESSAGE("Creating MCDisplay library")
FILE(GLOB algo_src . ${SRC}/*.cc)


#SET_TARGET_PROPERTIES(MCDisplay PROPERTIES VERSION 1.1.0 SOVERSION 1)

set(MCDisplay_DICTIONARY ${CMAKE_BINARY_DIR}/MCDisplay/MCDisplayDict.cxx)  # set dictionary

MESSAGE("Generate MCDisplay dictionary = ${MCDisplay_DICTIONARY}")

FILE(GLOB MCDisplay_HEADERS . ${SRC}/*.h ${GATE_INCLUDE_DIR}/GATE/Centella.h)
FILE(GLOB MCDisplay_SRC . ${SRC} ${GATE_INCLUDE_DIR}/GATE)

MESSAGE("MCDisplay_HEADERS = ${MCDisplay_HEADERS}")
MESSAGE("MCDisplay_SRC = ${MCDisplay_SRC}")

ROOT_GENERATE_DICTIONARY(
           "${MCDisplay_DICTIONARY}"            # path to dictionary to generate
           "${MCDisplay_SRC}"  # list of includes
           "${MCDisplay_HEADERS}" 
           "${SRC}/LinkDef.hh"   # ROOT linkDef file
         )

ADD_LIBRARY(MCDisplay SHARED ${algo_src} ${MCDisplay_DICTIONARY})
TARGET_LINK_LIBRARIES(MCDisplay ${ROOT_LIBRARIES} ${GATE_LIBRARIES})


# MESSAGE("Adding RunMCDisplay executable")
# add_executable(RunMCDisplay ${BIN}/RunMCDisplay.cc)

# target_link_libraries (RunMCDisplay MCDisplay ${ROOT_LIBRARIES} ${GATE_LIBRARIES})


INSTALL(TARGETS MCDisplay DESTINATION lib)
INSTALL(FILES ${CMAKE_BINARY_DIR}/MCDisplay/MCDisplayDict_rdict.pcm DESTINATION lib)
#install (TARGETS RunMCDisplay DESTINATION bin)


# add_subdirectory(bin)

# set (CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})
