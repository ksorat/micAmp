HOST=$(shell hostname)


ARCH = EMPTY
EXE = bin/amp.x
CC = icc

#Values of thread decomposition for NON-Phi runs
#Number of simultaneous sub-blocks
#Number of threads per sub-block
NUMSBS = 2
TPSB = 4

#Thread decomposition for Phi runs
#Threads/core on device (1-4)
#Sub-blocks per device
TPC = 3
SBPDEV = 4

ifneq (, $(findstring soratka1,$(HOST)))
	ARCH = MAC
endif

ifneq (, $(findstring smic,$(HOST)))
	ARCH = PHI
endif

ifeq "$(ARCH)" "MAC"
	CCFLAGS = -O3 -qopenmp -g -Wl,-stack_size,0x40000000,-stack_addr,0xf0000000
	CCFLAGS += -DNUMSBS=$(NUMSBS) -DTPSB=$(TPSB)
endif
ifeq "$(ARCH)" "PHI"
	CCFLAGS = -O3 -openmp -opt-report-phase=offload -D DOPHI
	#CCFLAGS = -O3 -opt-report3 -openmp -g -opt-report-phase=offload -D DOPHI
	CCFLAGS += -DTPC=$(TPC) -DSBPDEV=$(SBPDEV)
endif

#-Wl,-stack_size,0x10000000 /256MB
#-Wl,-stack_size,0x40000000,-stack_addr,0xf0000000 /1GB
#CCFLAGS = -debug all -traceback -g 
#CCFLAGS = -restrict -vec-report=4 -profile-functions -profile-loops=all

INCLUDES = -I.

SOURCE = $(wildcard *.cpp)
OBJS = $(SOURCE:.cpp=.o)

all:	$(EXE)
	@echo ""
	@echo "Compiled for $(ARCH) on $(HOST)"

$(EXE): $(OBJS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o $(EXE) $(OBJS) $(LFLAGS) $(LIBS)

.cpp.o:
	$(CC) $(CCFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	rm -f *.o *.optrpt *.vtk $(EXE)
