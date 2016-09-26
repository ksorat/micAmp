HOST=$(shell hostname)


ARCH = EMPTY
EXE = bin/amp.x
CC = icc

ifneq (, $(findstring soratka1,$(HOST)))

	ARCH = MAC
endif

ifneq (, $(findstring smic,$(HOST)))
	ARCH = PHI
endif

ifeq "$(ARCH)" "MAC"
	CCFLAGS = -O3 -qopenmp -g -Wl,-stack_size,0x40000000,-stack_addr,0xf0000000

endif
ifeq "$(ARCH)" "PHI"
	CCFLAGS = -O3 -openmp -g -opt-report-phase=offload -D DOPHI

endif
#-Wl,-stack_size,0x10000000 /256MB
#-Wl,-stack_size,0x40000000,-stack_addr,0xf0000000 /1GB
#CCFLAGS = -debug all -traceback -g 
#CCFLAGS = -restrict -vec-report=4 -profile-functions -profile-loops=all

INCLUDES = -I.

SOURCE = $(wildcard *.cpp)
OBJS = $(SOURCE:.cpp=.o)

all:	$(EXE)
	@echo "Compiled for $(ARCH) on $(HOST)"

$(EXE): $(OBJS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o $(EXE) $(OBJS) $(LFLAGS) $(LIBS)

.cpp.o:
	$(CC) $(CCFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	rm -f *.o *.optrpt *.vtk $(EXE)
