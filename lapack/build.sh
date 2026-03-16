#!/bin/bash

# Check if a folder name is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <project_folder>"
  exit 1
fi

# Input parameter
PROJECT_FOLDER="$1"

# Check if the folder exists
if [ ! -d "$PROJECT_FOLDER" ]; then
  echo "Error: Folder '$PROJECT_FOLDER' does not exist."
  exit 1
fi

# Check if the folder contains a CMakeLists.txt
if [ ! -f "$PROJECT_FOLDER/CMakeLists.txt" ]; then
  echo "Error: No CMakeLists.txt found in '$PROJECT_FOLDER'."
  exit 1
fi

# Build directory (nested inside the project folder)
BUILD_DIR="$PROJECT_FOLDER/build"

# rm -rf "$BUILD_DIR"
# mkdir -p "$BUILD_DIR"

# Step 1: Generate build files
echo "Generating build files for project in '$PROJECT_FOLDER'..."
# cmake --preset=default -DCMAKE_TOOLCHAIN_FILE=/opt/programs/vcpkg/scripts/buildsystems/vcpkg.cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug -B "$BUILD_DIR" -S "$PROJECT_FOLDER"
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug -B "$BUILD_DIR" -S "$PROJECT_FOLDER"

# Step 2: Build the project
echo "Building the project..."
cmake --build "$BUILD_DIR" --config Debug

# Step 3: Find the executable and run it
EXECUTABLE=$(find "$BUILD_DIR" -maxdepth 1 -type f -executable)

if [ -n "$EXECUTABLE" ]; then
  echo "Running the executable: $EXECUTABLE"
  "$EXECUTABLE"
else
  echo "Error: No executable found in '$BUILD_DIR'."
  exit 1
fi
