#!/bin/bash

make clean
rm CMakeCache.txt
rm -r CMakeFiles/
rm cmake_install.cmake
rm Makefile
cd sources
rm CMakeCache.txt
rm -r CMakeFiles/
rm cmake_install.cmake
rm Makefile
cd core_library
rm CMakeCache.txt
rm -r CMakeFiles/
rm cmake_install.cmake
rm Makefile
cd sources
rm -r CMakeFiles/
rm cmake_install.cmake
rm Makefile
cd ..
cd ..