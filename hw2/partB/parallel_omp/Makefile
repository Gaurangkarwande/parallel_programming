.PHONY: all clean

CC = gcc
INCLUDES = -I .
TARGET = strassen_parallel
MAP_FILE = strassen_parallel.map

CC_OPTIONS = -fopenmp -g -Wall
LD_FLAGS = -Xlinker -Map=$(MAP_FILE)
LD_LIBS = -lpthread

OBJ_FILES = strassen_parallel.o 

.SILENT: clean all %.o

all: ${OBJ_FILES}
	@echo "Building ${TARGET} ..."
	${CC} ${CC_OPTIONS} ${LD_FLAGS} ${OBJ_FILES} -o ${TARGET} ${LD_LIBS}

%.o: %.c
	@echo "Compiling $*.o ..."
	${CC} ${CC_OPTIONS} $(LD_FLAGS) -c $*.c ${INCLUDES}

clean:
	@echo "Cleaning $(OBJ_FILES) $(TARGET) $(MAP_FILE)..."
	-rm -rf	 *.o
	-rm -rf	 strassen_parallel.map
	-rm -rf	 strassen_parallel