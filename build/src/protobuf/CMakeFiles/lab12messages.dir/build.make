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
include src/protobuf/CMakeFiles/lab12messages.dir/depend.make

# Include the progress variables for this target.
include src/protobuf/CMakeFiles/lab12messages.dir/progress.make

# Include the compile flags for this target's objects.
include src/protobuf/CMakeFiles/lab12messages.dir/flags.make

src/protobuf/acomms_example.pb.cc: ../src/protobuf/acomms_example.proto
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/abelani/moos-ivp-abelani-acomms/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Running C++ protocol buffer compiler on acomms_example.proto"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf && /opt/local/bin/protoc --cpp_out /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf /Users/abelani/moos-ivp-abelani-acomms/src/protobuf/acomms_example.proto -I/opt/local/include -I/Users/abelani/moos-ivp-abelani-acomms/src/protobuf -I/Users/abelani/goby/include -I/Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf -I/Users/abelani/moos-ivp-abelani-acomms/src/protobuf

src/protobuf/acomms_example.pb.h: src/protobuf/acomms_example.pb.cc

src/protobuf/mini_command.pb.cc: ../src/protobuf/mini_command.proto
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/abelani/moos-ivp-abelani-acomms/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Running C++ protocol buffer compiler on mini_command.proto"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf && /opt/local/bin/protoc --cpp_out /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf /Users/abelani/moos-ivp-abelani-acomms/src/protobuf/mini_command.proto -I/opt/local/include -I/Users/abelani/moos-ivp-abelani-acomms/src/protobuf -I/Users/abelani/goby/include -I/Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf -I/Users/abelani/moos-ivp-abelani-acomms/src/protobuf

src/protobuf/mini_command.pb.h: src/protobuf/mini_command.pb.cc

src/protobuf/ctd_default.pb.cc: ../src/protobuf/ctd_default.proto
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/abelani/moos-ivp-abelani-acomms/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Running C++ protocol buffer compiler on ctd_default.proto"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf && /opt/local/bin/protoc --cpp_out /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf /Users/abelani/moos-ivp-abelani-acomms/src/protobuf/ctd_default.proto -I/opt/local/include -I/Users/abelani/moos-ivp-abelani-acomms/src/protobuf -I/Users/abelani/goby/include -I/Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf -I/Users/abelani/moos-ivp-abelani-acomms/src/protobuf

src/protobuf/ctd_default.pb.h: src/protobuf/ctd_default.pb.cc

src/protobuf/ctd.pb.cc: ../src/protobuf/ctd.proto
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/abelani/moos-ivp-abelani-acomms/build/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Running C++ protocol buffer compiler on ctd.proto"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf && /opt/local/bin/protoc --cpp_out /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf /Users/abelani/moos-ivp-abelani-acomms/src/protobuf/ctd.proto -I/opt/local/include -I/Users/abelani/moos-ivp-abelani-acomms/src/protobuf -I/Users/abelani/goby/include -I/Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf -I/Users/abelani/moos-ivp-abelani-acomms/src/protobuf

src/protobuf/ctd.pb.h: src/protobuf/ctd.pb.cc

src/protobuf/CMakeFiles/lab12messages.dir/acomms_example.pb.cc.o: src/protobuf/CMakeFiles/lab12messages.dir/flags.make
src/protobuf/CMakeFiles/lab12messages.dir/acomms_example.pb.cc.o: src/protobuf/acomms_example.pb.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/abelani/moos-ivp-abelani-acomms/build/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/protobuf/CMakeFiles/lab12messages.dir/acomms_example.pb.cc.o"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/lab12messages.dir/acomms_example.pb.cc.o -c /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf/acomms_example.pb.cc

src/protobuf/CMakeFiles/lab12messages.dir/acomms_example.pb.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lab12messages.dir/acomms_example.pb.cc.i"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf/acomms_example.pb.cc > CMakeFiles/lab12messages.dir/acomms_example.pb.cc.i

src/protobuf/CMakeFiles/lab12messages.dir/acomms_example.pb.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lab12messages.dir/acomms_example.pb.cc.s"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf/acomms_example.pb.cc -o CMakeFiles/lab12messages.dir/acomms_example.pb.cc.s

src/protobuf/CMakeFiles/lab12messages.dir/acomms_example.pb.cc.o.requires:
.PHONY : src/protobuf/CMakeFiles/lab12messages.dir/acomms_example.pb.cc.o.requires

src/protobuf/CMakeFiles/lab12messages.dir/acomms_example.pb.cc.o.provides: src/protobuf/CMakeFiles/lab12messages.dir/acomms_example.pb.cc.o.requires
	$(MAKE) -f src/protobuf/CMakeFiles/lab12messages.dir/build.make src/protobuf/CMakeFiles/lab12messages.dir/acomms_example.pb.cc.o.provides.build
.PHONY : src/protobuf/CMakeFiles/lab12messages.dir/acomms_example.pb.cc.o.provides

src/protobuf/CMakeFiles/lab12messages.dir/acomms_example.pb.cc.o.provides.build: src/protobuf/CMakeFiles/lab12messages.dir/acomms_example.pb.cc.o

src/protobuf/CMakeFiles/lab12messages.dir/mini_command.pb.cc.o: src/protobuf/CMakeFiles/lab12messages.dir/flags.make
src/protobuf/CMakeFiles/lab12messages.dir/mini_command.pb.cc.o: src/protobuf/mini_command.pb.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/abelani/moos-ivp-abelani-acomms/build/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/protobuf/CMakeFiles/lab12messages.dir/mini_command.pb.cc.o"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/lab12messages.dir/mini_command.pb.cc.o -c /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf/mini_command.pb.cc

src/protobuf/CMakeFiles/lab12messages.dir/mini_command.pb.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lab12messages.dir/mini_command.pb.cc.i"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf/mini_command.pb.cc > CMakeFiles/lab12messages.dir/mini_command.pb.cc.i

src/protobuf/CMakeFiles/lab12messages.dir/mini_command.pb.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lab12messages.dir/mini_command.pb.cc.s"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf/mini_command.pb.cc -o CMakeFiles/lab12messages.dir/mini_command.pb.cc.s

src/protobuf/CMakeFiles/lab12messages.dir/mini_command.pb.cc.o.requires:
.PHONY : src/protobuf/CMakeFiles/lab12messages.dir/mini_command.pb.cc.o.requires

src/protobuf/CMakeFiles/lab12messages.dir/mini_command.pb.cc.o.provides: src/protobuf/CMakeFiles/lab12messages.dir/mini_command.pb.cc.o.requires
	$(MAKE) -f src/protobuf/CMakeFiles/lab12messages.dir/build.make src/protobuf/CMakeFiles/lab12messages.dir/mini_command.pb.cc.o.provides.build
.PHONY : src/protobuf/CMakeFiles/lab12messages.dir/mini_command.pb.cc.o.provides

src/protobuf/CMakeFiles/lab12messages.dir/mini_command.pb.cc.o.provides.build: src/protobuf/CMakeFiles/lab12messages.dir/mini_command.pb.cc.o

src/protobuf/CMakeFiles/lab12messages.dir/ctd_default.pb.cc.o: src/protobuf/CMakeFiles/lab12messages.dir/flags.make
src/protobuf/CMakeFiles/lab12messages.dir/ctd_default.pb.cc.o: src/protobuf/ctd_default.pb.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/abelani/moos-ivp-abelani-acomms/build/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/protobuf/CMakeFiles/lab12messages.dir/ctd_default.pb.cc.o"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/lab12messages.dir/ctd_default.pb.cc.o -c /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf/ctd_default.pb.cc

src/protobuf/CMakeFiles/lab12messages.dir/ctd_default.pb.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lab12messages.dir/ctd_default.pb.cc.i"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf/ctd_default.pb.cc > CMakeFiles/lab12messages.dir/ctd_default.pb.cc.i

src/protobuf/CMakeFiles/lab12messages.dir/ctd_default.pb.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lab12messages.dir/ctd_default.pb.cc.s"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf/ctd_default.pb.cc -o CMakeFiles/lab12messages.dir/ctd_default.pb.cc.s

src/protobuf/CMakeFiles/lab12messages.dir/ctd_default.pb.cc.o.requires:
.PHONY : src/protobuf/CMakeFiles/lab12messages.dir/ctd_default.pb.cc.o.requires

src/protobuf/CMakeFiles/lab12messages.dir/ctd_default.pb.cc.o.provides: src/protobuf/CMakeFiles/lab12messages.dir/ctd_default.pb.cc.o.requires
	$(MAKE) -f src/protobuf/CMakeFiles/lab12messages.dir/build.make src/protobuf/CMakeFiles/lab12messages.dir/ctd_default.pb.cc.o.provides.build
.PHONY : src/protobuf/CMakeFiles/lab12messages.dir/ctd_default.pb.cc.o.provides

src/protobuf/CMakeFiles/lab12messages.dir/ctd_default.pb.cc.o.provides.build: src/protobuf/CMakeFiles/lab12messages.dir/ctd_default.pb.cc.o

src/protobuf/CMakeFiles/lab12messages.dir/ctd.pb.cc.o: src/protobuf/CMakeFiles/lab12messages.dir/flags.make
src/protobuf/CMakeFiles/lab12messages.dir/ctd.pb.cc.o: src/protobuf/ctd.pb.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/abelani/moos-ivp-abelani-acomms/build/CMakeFiles $(CMAKE_PROGRESS_8)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/protobuf/CMakeFiles/lab12messages.dir/ctd.pb.cc.o"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/lab12messages.dir/ctd.pb.cc.o -c /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf/ctd.pb.cc

src/protobuf/CMakeFiles/lab12messages.dir/ctd.pb.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lab12messages.dir/ctd.pb.cc.i"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf/ctd.pb.cc > CMakeFiles/lab12messages.dir/ctd.pb.cc.i

src/protobuf/CMakeFiles/lab12messages.dir/ctd.pb.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lab12messages.dir/ctd.pb.cc.s"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf/ctd.pb.cc -o CMakeFiles/lab12messages.dir/ctd.pb.cc.s

src/protobuf/CMakeFiles/lab12messages.dir/ctd.pb.cc.o.requires:
.PHONY : src/protobuf/CMakeFiles/lab12messages.dir/ctd.pb.cc.o.requires

src/protobuf/CMakeFiles/lab12messages.dir/ctd.pb.cc.o.provides: src/protobuf/CMakeFiles/lab12messages.dir/ctd.pb.cc.o.requires
	$(MAKE) -f src/protobuf/CMakeFiles/lab12messages.dir/build.make src/protobuf/CMakeFiles/lab12messages.dir/ctd.pb.cc.o.provides.build
.PHONY : src/protobuf/CMakeFiles/lab12messages.dir/ctd.pb.cc.o.provides

src/protobuf/CMakeFiles/lab12messages.dir/ctd.pb.cc.o.provides.build: src/protobuf/CMakeFiles/lab12messages.dir/ctd.pb.cc.o

# Object files for target lab12messages
lab12messages_OBJECTS = \
"CMakeFiles/lab12messages.dir/acomms_example.pb.cc.o" \
"CMakeFiles/lab12messages.dir/mini_command.pb.cc.o" \
"CMakeFiles/lab12messages.dir/ctd_default.pb.cc.o" \
"CMakeFiles/lab12messages.dir/ctd.pb.cc.o"

# External object files for target lab12messages
lab12messages_EXTERNAL_OBJECTS =

../lib/liblab12messages.dylib: src/protobuf/CMakeFiles/lab12messages.dir/acomms_example.pb.cc.o
../lib/liblab12messages.dylib: src/protobuf/CMakeFiles/lab12messages.dir/mini_command.pb.cc.o
../lib/liblab12messages.dylib: src/protobuf/CMakeFiles/lab12messages.dir/ctd_default.pb.cc.o
../lib/liblab12messages.dylib: src/protobuf/CMakeFiles/lab12messages.dir/ctd.pb.cc.o
../lib/liblab12messages.dylib: src/protobuf/CMakeFiles/lab12messages.dir/build.make
../lib/liblab12messages.dylib: /opt/local/lib/libprotobuf.dylib
../lib/liblab12messages.dylib: /Users/abelani/goby/lib/libgoby_acomms.dylib
../lib/liblab12messages.dylib: /Users/abelani/goby/lib/libgoby_moos.dylib
../lib/liblab12messages.dylib: /Users/abelani/goby/lib/libgoby_common.dylib
../lib/liblab12messages.dylib: /Users/abelani/goby/lib/libgoby_util.dylib
../lib/liblab12messages.dylib: src/protobuf/CMakeFiles/lab12messages.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared library ../../../lib/liblab12messages.dylib"
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lab12messages.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/protobuf/CMakeFiles/lab12messages.dir/build: ../lib/liblab12messages.dylib
.PHONY : src/protobuf/CMakeFiles/lab12messages.dir/build

src/protobuf/CMakeFiles/lab12messages.dir/requires: src/protobuf/CMakeFiles/lab12messages.dir/acomms_example.pb.cc.o.requires
src/protobuf/CMakeFiles/lab12messages.dir/requires: src/protobuf/CMakeFiles/lab12messages.dir/mini_command.pb.cc.o.requires
src/protobuf/CMakeFiles/lab12messages.dir/requires: src/protobuf/CMakeFiles/lab12messages.dir/ctd_default.pb.cc.o.requires
src/protobuf/CMakeFiles/lab12messages.dir/requires: src/protobuf/CMakeFiles/lab12messages.dir/ctd.pb.cc.o.requires
.PHONY : src/protobuf/CMakeFiles/lab12messages.dir/requires

src/protobuf/CMakeFiles/lab12messages.dir/clean:
	cd /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf && $(CMAKE_COMMAND) -P CMakeFiles/lab12messages.dir/cmake_clean.cmake
.PHONY : src/protobuf/CMakeFiles/lab12messages.dir/clean

src/protobuf/CMakeFiles/lab12messages.dir/depend: src/protobuf/acomms_example.pb.cc
src/protobuf/CMakeFiles/lab12messages.dir/depend: src/protobuf/acomms_example.pb.h
src/protobuf/CMakeFiles/lab12messages.dir/depend: src/protobuf/mini_command.pb.cc
src/protobuf/CMakeFiles/lab12messages.dir/depend: src/protobuf/mini_command.pb.h
src/protobuf/CMakeFiles/lab12messages.dir/depend: src/protobuf/ctd_default.pb.cc
src/protobuf/CMakeFiles/lab12messages.dir/depend: src/protobuf/ctd_default.pb.h
src/protobuf/CMakeFiles/lab12messages.dir/depend: src/protobuf/ctd.pb.cc
src/protobuf/CMakeFiles/lab12messages.dir/depend: src/protobuf/ctd.pb.h
	cd /Users/abelani/moos-ivp-abelani-acomms/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/abelani/moos-ivp-abelani-acomms /Users/abelani/moos-ivp-abelani-acomms/src/protobuf /Users/abelani/moos-ivp-abelani-acomms/build /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf /Users/abelani/moos-ivp-abelani-acomms/build/src/protobuf/CMakeFiles/lab12messages.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/protobuf/CMakeFiles/lab12messages.dir/depend
