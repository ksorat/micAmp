ARCH = MAC
#MODE = PHI

#MODE = DEBUG
#MODE = PRODUCTION

EXE = bin/amp.x
CC = icc

ifeq "$(ARCH)" "MAC"
	CCFLAGS = -O3 -qopenmp -g -Wl,-stack_size,0x40000000,-stack_addr,0xf0000000

endif
ifeq "$(ARCH)" "PHI"
	CCFLAGS = -O3 -qopenmp -g -D DOPHI

endif
#-Wl,-stack_size,0x10000000 /256MB
#-Wl,-stack_size,0x40000000,-stack_addr,0xf0000000 /1GB
#CCFLAGS = -debug all -traceback -g 
#CCFLAGS = -restrict -vec-report=4 -profile-functions -profile-loops=all

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
