all: cell_types.mod

cell_types.mod: cell_types.f90
	$(FC) $(FCFLAGS) -I ../base -I ../common -c cell_types.f90

clean:
	rm -f *.o *.mod

