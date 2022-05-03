
##########################################################

## CC COMPILER OPTIONS ##

# CC compiler options:
CC=nvcc
CC_FLAGS= -g
CC_LIBS=

# CUDA linking libraries:
CUDA_LINK_LIBS= -lcudart

##########################################################

## Project file structure ##

# Source file directory:
SRC_DIR = src

# Object file directory:
OBJ_DIR = bin

# Include header file diretory:
INC_DIR = include

##########################################################

## Make variables ##

# Target executable name:
EXE = ssp

SOURCES_C := $(wildcard $(SRC_DIR)/*.cpp)
SOURCES_CUDA := $(wildcard $(SRC_DIR)/*.cu)

# Object files:
OBJS := $(OBJ_DIR)/main.o $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SOURCES_C)) $(patsubst $(SRC_DIR)/%.cu, $(OBJ_DIR)/%.o, $(SOURCES_CUDA))

##########################################################

## Compile ##

# Link C and CUDA compiled object files to target executable:
$(EXE) : $(OBJS)
	$(CC) $(CC_FLAGS) $(OBJS) -o $@ $(CUDA_LINK_LIBS)

# Compile main .c file to object files:
$(OBJ_DIR)/%.o : %.cpp
	$(CC) $(CC_FLAGS) -c $< -o $@

# Compile C++ source files to object files:
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp include/%.h 
	$(CC) $(CC_FLAGS) -c $< -o $@

# Compile CUDA source files to object files:
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cu $(INC_DIR)/%.cuh
	$(CC) $(CC_FLAGS) -c $< -o $@

# Clean objects in object directory.
clean:
	$(RM) bin/* *.o $(EXE)
