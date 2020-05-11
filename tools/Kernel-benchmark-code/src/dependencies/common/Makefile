all: mathconstants.mod mathlib.mod

mathconstants.mod: mathconstants.f90
	$(FC) $(FCFLAGS) -I ../base -c mathconstants.f90

mathlib.mod: mathlib.f90
	$(FC) $(FCFLAGS) -I ../base -c mathlib.f90

clean:
	rm -f *.o *.mod

