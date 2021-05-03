# Declare compiler, required flags, and dependencies
CC = mpicxx
CXXFLAGS = -std=c++11 -Wall -pedantic -O2
LDLIBS = -lboost_program_options -llapack -lblas -lscalapack-openmpi
HDRS = LidDrivenCavity.h PoissonSolver.h
TARGET = Solve
OBJS = LidDrivenCavitySolver.o LidDrivenCavity.o PoissonSolver.o

# Define make target
default: $(TARGET)
all: $(TARGET)

# Generate object files
%.o : %.cpp $(HDRS)
	$(CC) $(CXXFLAGS) -o $@ -c $< $(LDLIBS)

# Link object file and produce executable
$(TARGET): $(OBJS)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

.PHONY: testParallel testSerial clean cleanPlot cleanAll

# Parallel test case for the code
test:$(TARGET)
	mpiexec -np 4 ./$(TARGET) --Nx=5 --Ny=5 --Px=2 --Py=2

# Serial test case for the code
testSerial:$(TARGET)
	mpiexec -np 1 ./$(TARGET) --Nx=5 --Ny=5

# Clean commands
clean:
	rm -f $(TARGET) *.o

cleanPlot:
	rm -f Plot/*

cleanAll:
	rm -f Plot/* *.o $(TARGET)
