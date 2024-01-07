#!/bin/bash

# Remove the "build" directory and its contents
rm -rf build

# Create the "build" directory
mkdir build

# Change to the "build" directory
cd build

# Run CMake to generate Makefiles
if [ "$1" == "debug" ]; then
    # Run cmake in the parent directory
    cmake -DCMAKE_BUILD_TYPE=Debug ..
else
    cmake ..
fi

make

# Pause to keep the terminal open (optional)
read -p "Press Enter to exit..."
