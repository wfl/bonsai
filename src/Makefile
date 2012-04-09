################################################################################
# Makefile for BONSAI - INDEL Haplotype Local Realigner
# Wen Fung Leong
# Boston College, 2012
################################################################################

# -------------------
# BAMTOOLs
# -------------------
BAMTOOLS_DIR=../bamtools
BAMTOOLS_INCLUDE_DIR = $(BAMTOOLS_DIR)/include/
BAMTOOLS_LIB_DIR     = $(BAMTOOLS_DIR)/lib

# --------------------------------
# # Defining the compiler:
# --------------------------------
CXX=g++

# --------------------------------------
# Defining compiler flags
# --------------------------------------
CXXFLAGS = -Wall -O2 -g -D_FILE_OFFSET_BITS=64 -Wl,-rpath,$(BAMTOOLS_LIB_DIR)

OBJ_DIR = .
BIN_DIR = ../bin

PROGRAM=bonsai
all: $(BIN_DIR)/$(PROGRAM)
		
# --------------------------------------
# Defining includes
# --------------------------------------
LIBS     = -lbamtools -lz
INCLUDES = -I$(BAMTOOLS_INCLUDE_DIR)

# ----------------------------------
# Defining source and object files
# ----------------------------------
SOURCES=Fasta.cpp \
		ModifiedSmithWatermanGotoh.cpp \
		split.cpp \
		GenFunc.cpp \
		Realign.cpp \
		realigner.cpp

OBJECTS= $(SOURCES:.cpp=.o)
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

# -----------------------------------------------------
# Compile Executable file
# -----------------------------------------------------
$(PROGRAM) ../bin/$(PROGRAM): $(BUILT_OBJECTS)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $(BIN_DIR)/$@ $^ -L$(BAMTOOLS_LIB_DIR) $(LIBS)

# -----------------------------------------------------
# Compile Objects files
# -----------------------------------------------------
$(BUILT_OBJECTS): $(SOURCES)
	$(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(INCLUDES)


# -----------------------------------------------------
# Removing the executable and the object files
# -----------------------------------------------------
clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/*.o $(BIN_DIR)/$(PROGRAM)
	cd $(BAMTOOLS_DIR)/build && make clean

.PHONY: clean all


