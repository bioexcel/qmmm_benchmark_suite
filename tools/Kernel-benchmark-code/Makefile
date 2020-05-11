export FC = gfortran
export FCFLAGS = -O3 -fopenmp
export LDFLAGS = $(FCFLAGS)

all:
	( cd src; make )

clean:
	rm -f kernel_benchmark
	( cd src; make clean )
