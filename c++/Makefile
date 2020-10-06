# Lattice Boltzmann Makefile

# Compiler
COMPILE= g++

# Aditional user compiler flags
USER_FLAGS= -O3
OBJFILES= main.o \
          LatticeBoltzmann.o 

# ===================================================================

CFLAGS= $(USER_FLAGS)
PROGRAM= main

$(PROGRAM): $(OBJFILES)
	@echo Linking object files ...
	$(COMPILE) $(CFLAGS) main.o LatticeBoltzmann.o -o $(PROGRAM)

main.o: LatticeBoltzmann.cpp LatticeBoltzmann.h main.cpp
	@echo Compile main program ...
	$(COMPILE) $(CFLAGS) -c main.cpp

LatticeBoltzmann.o: LatticeBoltzmann.cpp LatticeBoltzmann.h 
	@echo Compile Lattice Boltzmann Class ...
	$(COMPILE) $(CFLAGS) -c LatticeBoltzmann.cpp

clean :
	rm -f DATA/*.o 
	rm -f DATA/$(PROGRAM)
	rm -f DATA/*.dat
