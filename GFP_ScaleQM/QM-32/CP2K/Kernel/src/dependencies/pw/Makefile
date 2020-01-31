all: pw_grid_types.mod pw_types.mod

pw_grid_types.mod: pw_grid_types.f90
	$(FC) $(FCFLAGS) -I ../base -c pw_grid_types.f90

pw_types.mod: pw_types.f90
	$(FC) $(FCFLAGS) -I ../base -c pw_types.f90

clean:
	rm -f *.o *.mod

