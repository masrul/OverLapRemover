cmake_minimum_required (VERSION 3.0)
project (overlap_remover)
enable_language (CXX)


file(GLOB_RECURSE CXXFiles 
    src/pbc.cpp
    src/remover.cpp
    src/string_manipulate.cpp
    src/main.cpp
    )

file(GLOB_RECURSE HeaderFiles 
    src/remover.h
    src/util.h
    src/common.h
    )

add_executable(overlap_remover ${CXXFiles} ${HeaderFiles})
install(TARGETS overlap_remover DESTINATION ${CMAKE_INSTALL_PREFIX})
