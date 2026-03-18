# Minimum CMake version required
cmake_minimum_required(VERSION 3.10)

# Project name and C++ standard
project(gen_tcs_bh CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fPIC -ansi -Wall -Wno-overloaded-virtual -Wno-long-long -fno-common -O2")

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS RIO Net Hist Tree)
include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIRS})
add_definitions(${ROOT_DEFINITIONS})

# Source files
set(ANA_C
    FillTables.cc
    RadProcess.cc
    TableProvider.cc
    TabUtils.cc
    interpol.cc
    Utils.cc
    Options_tcs.cc
    FormFactors.cc
    TreatOptions.cc
    PDFparam.cc
    Constants.cc
    Polar.cc
    PDGinfo.cc
    Kinematics.cc
    LeptonPair.cc
    PartUtils.cc
    ReactionKinTwo.cc
)

# Main executable
add_executable(gen_tcs_bh main.cc ${ANA_C})

# Link ROOT libraries
target_link_libraries(gen_tcs_bh ${ROOT_LIBRARIES})


# Custom target for cleaning files, copying scripts, and executing the script
add_custom_target(custom_clean_and_setup
    # Cleaning up generated files
    COMMAND ${CMAKE_COMMAND} -E remove -f *.o *~ main includes/*.o *.o Data/*
    COMMENT "Cleaning up generated files"

    # Copying files
    COMMAND ${CMAKE_COMMAND} -E copy ../set.csh ${CMAKE_BINARY_DIR}/set.csh
#    COMMAND ${CMAKE_COMMAND} -E copy ../TCS.input ${CMAKE_BINARY_DIR}/TCS.input
#    COMMAND ${CMAKE_COMMAND} -E copy ../PS_EEPHOTO_FIX.input ${CMAKE_BINARY_DIR}/PS_EEPHOTO_FIX.input
    # Executing the script
    COMMAND ./set.csh
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Copying set.csh, then executing set.csh"
)