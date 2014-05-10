#
# 'make'        build executable file
# 'make clean'  removes all .o and executable files
#

# Executables to build
#SDLEXEC += viewer
#SDLEXEC += subgraph
#SDLEXEC += mass_spring
SDLEXEC += shallow_water_ext
#SDLEXEC += test_nodes
#SDLEXEC += test_edges
#SDLEXEC += test

#SDLEXEC += viewer
#SDLEXEC += subgraph
#SDLEXEC += shortest_path
SDLEXEC += mass_spring
#SDLEXEC += poisson
#SDLEXEC += shallow_water
#SDLEXEC += cp_shallow_water
SDLEXEC += OPENMPshallow
SDLEXEC += OPENMPshallow2
SDLEXEC += openmp_mass_spring

# Get the shell name to determine the OS
UNAME := $(shell uname)

# Define the C++ compiler to use
CXX := $(shell which g++) -std=gnu++0x
ifeq ($(UNAME), Darwin)
CC := $(shell which gcc) -O3
SDLOBJS := CS207/SDLMain.o
endif

# Dependency directory and flags
DEPSDIR := $(shell mkdir -p .deps; echo .deps)
# MD: Dependency as side-effect of compilation
# MF: File for output
# MP: Include phony targets
DEPSFILE = $(DEPSDIR)/$(notdir $*.d)
DEPSFLAGS = -MD -MF $(DEPSFILE) #-MP

# Define any directories containing header files
#   To include directories use -Ipath/to/files

#INCLUDES += -I.
INCLUDES += -I. -I/usr/include -I/usr/include/libxml2 #added
INCLUDES += -I./mlpack-1.0.8/build/include #added

# Define CXX compile flags
CXXFLAGS += -fopenmp -funroll-loops -O3 -W -Wall -Wextra #-Wfatal-errors

# Define any directories containing libraries
#   To include directories use -Lpath/to/files
LDFLAGS += -L./mlpack-1.0.8/build/lib

# Define any libraries to link into executable
#   To link in libraries (libXXX.so or libXXX.a) use -lXXX
ifeq ($(UNAME), Linux)
LDLIBS += -lSDL -lGL -lGLU -lmlpack
endif
ifeq ($(UNAME), Darwin)
LDLIBS += -framework SDL -framework OpenGL -framework Cocoa -lmlpack
endif

##################
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
##################

# 'make' - default rule
all: $(EXEC) $(SDLEXEC)

# Default rule for creating an exec of $(EXEC) from a .o file
$(EXEC): % : %.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# Default rule for creating an exec of $(EXEC) from a .o file
$(SDLEXEC): % : %.o $(SDLOBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# Default rule for creating a .o file from a .cpp file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DEPSFLAGS) -c -o $@ $<

# Extra dependencies for executables
#   Nothing here

# 'make clean' - deletes all .o files, exec, and dependency files
clean:
	-$(RM) *.o $(EXEC) $(SDLEXEC) $(SDLOBJS)
	$(RM) -r $(DEPSDIR)

# Define rules that do not actually generate the corresponding file
.PHONY: clean all

# Include the dependency files
-include $(wildcard $(DEPSDIR)/*.d)
