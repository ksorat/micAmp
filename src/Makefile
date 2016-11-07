HOST=$(shell hostname)
#DOHOST = YES
EXE = bin/amp.x
CC = icc

#Values of thread decomposition for NON-Phi runs
#Number of simultaneous sub-blocks
#Number of threads per sub-block
NUMSBS = 8
TPSB = 4

#Thread decomposition for Phi runs
#Threads/core on device (1-4)
#Sub-blocks per device
TPC = 4
SBPDEV = 4

ARCH = EMPTY
ifneq (, $(findstring soratka1,$(HOST)))
	ARCH = MAC
endif

ifneq (, $(findstring smic,$(HOST)))
	ARCH = PHI
endif

ifneq (, $(findstring ys,$(HOST)))
	ARCH = YELLOWSTONE
endif

ifeq "$(DOHOST)" "YES"
	ARCH = SMICHOST
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
ifeq "$(ARCH)" "SMICHOST"
	CCFLAGS = -O3 -openmp
	#CCFLAGS = -O3 -opt-report3 -openmp -g -opt-report-phase=offload -D DOPHI
	CCFLAGS += -DNUMSBS=$(NUMSBS) -DTPSB=$(TPSB)
endif

ifeq "$(ARCH)" "YELLOWSTONE"
	NUMSBS = 4
	TPSB = 4
	CCFLAGS = -O3 -qopenmp 
	CCFLAGS += -DNUMSBS=$(NUMSBS) -DTPSB=$(TPSB)
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
