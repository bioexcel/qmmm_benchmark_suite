all: build_base build_common build_pw build_subsys build_qmmm_util

build_base:
	( cd base; make )

build_common:
	( cd common; make )

build_pw:
	( cd pw; make )

build_subsys:
	( cd subsys; make )

build_qmmm_util: qmmm_util.f90
	$(FC) $(FCFLAGS) -I base -I common -c qmmm_util.f90

clean:
	rm -f *.o *.mod
	( cd base; make clean )
	( cd common; make clean )
	( cd pw; make clean )
	( cd subsys; make clean )

