MODE = INTEL
EXE = bin/amp.x
#-Wl,-stack_size,0x10000000 /256MB
#-Wl,-stack_size,0x40000000,-stack_addr,0xf0000000 /1GB
ifeq "$(MODE)" "INTEL"
	CC = icc
	#CCFLAGS = -O3 -openmp -vec-report=4 -g -traceback -check-pointers=rw
	#CCFLAGS = -O3 -openmp -vec-report=1
	CCFLAGS = -O3 -qopenmp -g -Wl,-stack_size,0x40000000,-stack_addr,0xf0000000
	#CCFLAGS = -O3 -qopenmp -restrict -g -profile-functions -profile-loops=all
	#CCFLAGS = -debug all -traceback -g 
endif

INCLUDES = -I.

SOURCE = $(wildcard *.cpp)
OBJS = $(SOURCE:.cpp=.o)

all:	$(EXE)
	@echo Compiling

$(EXE): $(OBJS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o $(EXE) $(OBJS) $(LFLAGS) $(LIBS)

.cpp.o:
	$(CC) $(CCFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	rm -f *.o *.optrpt *.vtk $(EXE)
