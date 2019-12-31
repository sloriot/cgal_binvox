# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.14)
project( cgal_binvox.git )

find_package(CGAL REQUIRED QUIET)

create_single_source_cgal_program( "binvox.cpp" )

