MODE = INTEL
SCHEME = MUSCL
EXE = ../bin/amp.x

ifeq "$(MODE)" "INTEL"
	CC = icc
	#CCFLAGS = -O3 -openmp -vec-report=4 -g -traceback -check-pointers=rw
	#CCFLAGS = -O3 -openmp -vec-report=1
	CCFLAGS = -O3 -restrict -qopenmp 
	#CCFLAGS = -O3 -qopenmp -restrict -g -profile-functions -profile-loops=all
	#CCFLAGS = -debug all -traceback -g -qopt-report=4
endif

INCLUDES = -I.

#CCFLAGS += -D$(SCHEME) -DVISCOSITY
SOURCE = $(wildcard *.cpp)
OBJS = $(SOURCE:.cpp=.o)

all:	$(EXE)
	@echo Compiling

$(EXE): $(OBJS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o $(EXE) $(OBJS) $(LFLAGS) $(LIBS)

.cpp.o:
	$(CC) $(CCFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	rm -f *.o *.optrpt $(EXE)
