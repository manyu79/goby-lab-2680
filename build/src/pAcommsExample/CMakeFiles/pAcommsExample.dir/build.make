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
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /opt/local/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/abelani/moos-ivp-abelani-acomms

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/abelani/moos-ivp-abelani-acomms/build

# Include any dependencies generated for this target.
include src/pAcommsExample/CMakeFiles/pAcommsExample.dir/depend.make

# Include the progress variables for this target.
include src/pAcommsExample/CMakeFiles/pAcommsExample.dir/progress.make

# Include the compile flags for this target's objects.
include src/pAcommsExample/CMakeFiles/pAcommsExample.dir/flags.make

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.o: src/pAcommsExample/CMakeFiles/pAcommsExample.dir/flags.make
src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.o: ../src/pAcommsExample/AcommsExample.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/abelani/moos-ivp-abelani-acomms/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.o"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/pAcommsExample && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.o -c /Users/abelani/moos-ivp-abelani-acomms/src/pAcommsExample/AcommsExample.cpp

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.i"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/pAcommsExample && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/abelani/moos-ivp-abelani-acomms/src/pAcommsExample/AcommsExample.cpp > CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.i

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.s"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/pAcommsExample && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/abelani/moos-ivp-abelani-acomms/src/pAcommsExample/AcommsExample.cpp -o CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.s

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.o.requires:
.PHONY : src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.o.requires

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.o.provides: src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.o.requires
	$(MAKE) -f src/pAcommsExample/CMakeFiles/pAcommsExample.dir/build.make src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.o.provides.build
.PHONY : src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.o.provides

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.o.provides.build: src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.o

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.o: src/pAcommsExample/CMakeFiles/pAcommsExample.dir/flags.make
src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.o: ../src/pAcommsExample/AcommsExample_ExampleConfig.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/abelani/moos-ivp-abelani-acomms/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.o"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/pAcommsExample && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.o -c /Users/abelani/moos-ivp-abelani-acomms/src/pAcommsExample/AcommsExample_ExampleConfig.cpp

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.i"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/pAcommsExample && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/abelani/moos-ivp-abelani-acomms/src/pAcommsExample/AcommsExample_ExampleConfig.cpp > CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.i

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.s"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/pAcommsExample && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/abelani/moos-ivp-abelani-acomms/src/pAcommsExample/AcommsExample_ExampleConfig.cpp -o CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.s

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.o.requires:
.PHONY : src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.o.requires

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.o.provides: src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.o.requires
	$(MAKE) -f src/pAcommsExample/CMakeFiles/pAcommsExample.dir/build.make src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.o.provides.build
.PHONY : src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.o.provides

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.o.provides.build: src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.o

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/main.cpp.o: src/pAcommsExample/CMakeFiles/pAcommsExample.dir/flags.make
src/pAcommsExample/CMakeFiles/pAcommsExample.dir/main.cpp.o: ../src/pAcommsExample/main.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/abelani/moos-ivp-abelani-acomms/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/pAcommsExample/CMakeFiles/pAcommsExample.dir/main.cpp.o"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/pAcommsExample && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/pAcommsExample.dir/main.cpp.o -c /Users/abelani/moos-ivp-abelani-acomms/src/pAcommsExample/main.cpp

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pAcommsExample.dir/main.cpp.i"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/pAcommsExample && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/abelani/moos-ivp-abelani-acomms/src/pAcommsExample/main.cpp > CMakeFiles/pAcommsExample.dir/main.cpp.i

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pAcommsExample.dir/main.cpp.s"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/pAcommsExample && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/abelani/moos-ivp-abelani-acomms/src/pAcommsExample/main.cpp -o CMakeFiles/pAcommsExample.dir/main.cpp.s

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/main.cpp.o.requires:
.PHONY : src/pAcommsExample/CMakeFiles/pAcommsExample.dir/main.cpp.o.requires

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/main.cpp.o.provides: src/pAcommsExample/CMakeFiles/pAcommsExample.dir/main.cpp.o.requires
	$(MAKE) -f src/pAcommsExample/CMakeFiles/pAcommsExample.dir/build.make src/pAcommsExample/CMakeFiles/pAcommsExample.dir/main.cpp.o.provides.build
.PHONY : src/pAcommsExample/CMakeFiles/pAcommsExample.dir/main.cpp.o.provides

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/main.cpp.o.provides.build: src/pAcommsExample/CMakeFiles/pAcommsExample.dir/main.cpp.o

# Object files for target pAcommsExample
pAcommsExample_OBJECTS = \
"CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.o" \
"CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.o" \
"CMakeFiles/pAcommsExample.dir/main.cpp.o"

# External object files for target pAcommsExample
pAcommsExample_EXTERNAL_OBJECTS =

../bin/pAcommsExample: src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.o
../bin/pAcommsExample: src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.o
../bin/pAcommsExample: src/pAcommsExample/CMakeFiles/pAcommsExample.dir/main.cpp.o
../bin/pAcommsExample: src/pAcommsExample/CMakeFiles/pAcommsExample.dir/build.make
../bin/pAcommsExample: /Users/abelani/moos-ivp/MOOS/MOOSCore/lib/libMOOS.a
../bin/pAcommsExample: ../lib/liblab12codecs.dylib
../bin/pAcommsExample: ../lib/liblab12messages.dylib
../bin/pAcommsExample: /Users/abelani/goby/lib/libgoby_acomms.dylib
../bin/pAcommsExample: /Users/abelani/goby/lib/libgoby_moos.dylib
../bin/pAcommsExample: /Users/abelani/goby/lib/libgoby_common.dylib
../bin/pAcommsExample: /Users/abelani/goby/lib/libgoby_util.dylib
../bin/pAcommsExample: /opt/local/lib/libprotobuf.dylib
../bin/pAcommsExample: /opt/local/lib/libboost_signals-mt.dylib
../bin/pAcommsExample: /opt/local/lib/libboost_system-mt.dylib
../bin/pAcommsExample: /Users/abelani/goby/lib/libgoby_acomms.dylib
../bin/pAcommsExample: /Users/abelani/goby/lib/libgoby_moos.dylib
../bin/pAcommsExample: /Users/abelani/goby/lib/libgoby_common.dylib
../bin/pAcommsExample: /Users/abelani/goby/lib/libgoby_util.dylib
../bin/pAcommsExample: /opt/local/lib/libprotobuf.dylib
../bin/pAcommsExample: src/pAcommsExample/CMakeFiles/pAcommsExample.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../../bin/pAcommsExample"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/pAcommsExample && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pAcommsExample.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/pAcommsExample/CMakeFiles/pAcommsExample.dir/build: ../bin/pAcommsExample
.PHONY : src/pAcommsExample/CMakeFiles/pAcommsExample.dir/build

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/requires: src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample.cpp.o.requires
src/pAcommsExample/CMakeFiles/pAcommsExample.dir/requires: src/pAcommsExample/CMakeFiles/pAcommsExample.dir/AcommsExample_ExampleConfig.cpp.o.requires
src/pAcommsExample/CMakeFiles/pAcommsExample.dir/requires: src/pAcommsExample/CMakeFiles/pAcommsExample.dir/main.cpp.o.requires
.PHONY : src/pAcommsExample/CMakeFiles/pAcommsExample.dir/requires

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/clean:
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/pAcommsExample && $(CMAKE_COMMAND) -P CMakeFiles/pAcommsExample.dir/cmake_clean.cmake
.PHONY : src/pAcommsExample/CMakeFiles/pAcommsExample.dir/clean

src/pAcommsExample/CMakeFiles/pAcommsExample.dir/depend:
	cd /Users/abelani/moos-ivp-abelani-acomms/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/abelani/moos-ivp-abelani-acomms /Users/abelani/moos-ivp-abelani-acomms/src/pAcommsExample /Users/abelani/moos-ivp-abelani-acomms/build /Users/abelani/moos-ivp-abelani-acomms/build/src/pAcommsExample /Users/abelani/moos-ivp-abelani-acomms/build/src/pAcommsExample/CMakeFiles/pAcommsExample.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/pAcommsExample/CMakeFiles/pAcommsExample.dir/depend
