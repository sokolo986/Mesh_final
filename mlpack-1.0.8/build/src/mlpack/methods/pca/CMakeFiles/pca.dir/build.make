# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build

# Include any dependencies generated for this target.
include src/mlpack/methods/pca/CMakeFiles/pca.dir/depend.make

# Include the progress variables for this target.
include src/mlpack/methods/pca/CMakeFiles/pca.dir/progress.make

# Include the compile flags for this target's objects.
include src/mlpack/methods/pca/CMakeFiles/pca.dir/flags.make

src/mlpack/methods/pca/CMakeFiles/pca.dir/pca_main.cpp.o: src/mlpack/methods/pca/CMakeFiles/pca.dir/flags.make
src/mlpack/methods/pca/CMakeFiles/pca.dir/pca_main.cpp.o: ../src/mlpack/methods/pca/pca_main.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/mlpack/methods/pca/CMakeFiles/pca.dir/pca_main.cpp.o"
	cd /home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/pca && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/pca.dir/pca_main.cpp.o -c /home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/src/mlpack/methods/pca/pca_main.cpp

src/mlpack/methods/pca/CMakeFiles/pca.dir/pca_main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pca.dir/pca_main.cpp.i"
	cd /home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/pca && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/src/mlpack/methods/pca/pca_main.cpp > CMakeFiles/pca.dir/pca_main.cpp.i

src/mlpack/methods/pca/CMakeFiles/pca.dir/pca_main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pca.dir/pca_main.cpp.s"
	cd /home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/pca && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/src/mlpack/methods/pca/pca_main.cpp -o CMakeFiles/pca.dir/pca_main.cpp.s

src/mlpack/methods/pca/CMakeFiles/pca.dir/pca_main.cpp.o.requires:
.PHONY : src/mlpack/methods/pca/CMakeFiles/pca.dir/pca_main.cpp.o.requires

src/mlpack/methods/pca/CMakeFiles/pca.dir/pca_main.cpp.o.provides: src/mlpack/methods/pca/CMakeFiles/pca.dir/pca_main.cpp.o.requires
	$(MAKE) -f src/mlpack/methods/pca/CMakeFiles/pca.dir/build.make src/mlpack/methods/pca/CMakeFiles/pca.dir/pca_main.cpp.o.provides.build
.PHONY : src/mlpack/methods/pca/CMakeFiles/pca.dir/pca_main.cpp.o.provides

src/mlpack/methods/pca/CMakeFiles/pca.dir/pca_main.cpp.o.provides.build: src/mlpack/methods/pca/CMakeFiles/pca.dir/pca_main.cpp.o

# Object files for target pca
pca_OBJECTS = \
"CMakeFiles/pca.dir/pca_main.cpp.o"

# External object files for target pca
pca_EXTERNAL_OBJECTS =

bin/pca: src/mlpack/methods/pca/CMakeFiles/pca.dir/pca_main.cpp.o
bin/pca: lib/libmlpack.so.1.0
bin/pca: /usr/lib/libarmadillo.so
bin/pca: /usr/lib/libboost_program_options-mt.so
bin/pca: /usr/lib/libboost_unit_test_framework-mt.so
bin/pca: /usr/lib/libboost_random-mt.so
bin/pca: /usr/lib/libboost_program_options-mt.so
bin/pca: /usr/lib/libboost_unit_test_framework-mt.so
bin/pca: /usr/lib/libboost_random-mt.so
bin/pca: /usr/lib/i386-linux-gnu/libxml2.so
bin/pca: src/mlpack/methods/pca/CMakeFiles/pca.dir/build.make
bin/pca: src/mlpack/methods/pca/CMakeFiles/pca.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../../../bin/pca"
	cd /home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/pca && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pca.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/mlpack/methods/pca/CMakeFiles/pca.dir/build: bin/pca
.PHONY : src/mlpack/methods/pca/CMakeFiles/pca.dir/build

src/mlpack/methods/pca/CMakeFiles/pca.dir/requires: src/mlpack/methods/pca/CMakeFiles/pca.dir/pca_main.cpp.o.requires
.PHONY : src/mlpack/methods/pca/CMakeFiles/pca.dir/requires

src/mlpack/methods/pca/CMakeFiles/pca.dir/clean:
	cd /home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/pca && $(CMAKE_COMMAND) -P CMakeFiles/pca.dir/cmake_clean.cmake
.PHONY : src/mlpack/methods/pca/CMakeFiles/pca.dir/clean

src/mlpack/methods/pca/CMakeFiles/pca.dir/depend:
	cd /home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8 /home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/src/mlpack/methods/pca /home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build /home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/pca /home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/pca/CMakeFiles/pca.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/mlpack/methods/pca/CMakeFiles/pca.dir/depend

