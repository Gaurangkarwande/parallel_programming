.PHONY: all clean

CC = mpicc
INCLUDES = -I .
TARGET = bfs
MAP_FILE = bfs.map

CC_OPTIONS = -Wall -O2 -g -fopenmp
LD_FLAGS = -Xlinker -Map=$(MAP_FILE)
LD_LIBS = -lblas -lm

OBJ_FILES = bfs.o #utils.o quicksort.o

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
	-rm -rf	 bfs.map
	-rm -rf	 bfs