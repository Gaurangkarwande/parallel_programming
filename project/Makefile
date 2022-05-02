.PHONY: all clean

CC = nvcc
INCLUDES = -I .
TARGET = ssp
MAP_FILE = ssp.map

CC_OPTIONS = -g
LD_FLAGS = -Xlinker -Map=$(MAP_FILE)
LD_LIBS = -lblas -lm

SRC = $(wildcard *.c)
OBJ_FILES = $(SRC:.c=.o) #utils.o quicksort.o

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
	-rm -rf	 ssp.map
	-rm -rf	 ssp