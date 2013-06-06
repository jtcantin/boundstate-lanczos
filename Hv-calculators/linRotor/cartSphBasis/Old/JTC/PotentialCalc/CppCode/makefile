CC=g++
SUP_PATH1= ./vectorClasses/
SUP_PATH2= ./supportingRoutines/
INCLUDES= -I $(SUP_PATH1) -I $(SUP_PATH2)
CFLAGS=-c -Wall -fopenmp -DLINUX -ggdb -O0 $(INCLUDES)
LIB_SOURCES= $(SUP_PATH1)matvec2.cpp $(SUP_PATH1)Cvectors.cpp $(SUP_PATH2)gridFcns.cpp $(SUP_PATH2)energyRoutines.cpp $(SUP_PATH2)TIP4P_AH2_EngRoutines_Calc.cpp $(SUP_PATH2)rotmat.cpp $(SUP_PATH2)Alavi_H2_Routines.cpp
LDFLAGS=-fopenmp -lgfortran -llapack
SOURCES=PointPotentials_TIP4P_AH2_sep9.cpp $(LIB_SOURCES)
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=PointPotentials_TIP4P_AH2_sep9

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)


