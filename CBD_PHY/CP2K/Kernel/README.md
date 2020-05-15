# CBD_PHY system

## QM/MM Kernel Benchmark

### About

The files within this directory, along with the benchmark source code in /tools/kernel-benchmark-code can be used to benchamrk ther performance of the periodic qmmm_force_with_gaussian_LG kernel used within the GEEP QMMM routines. 



``/data`` - Contains the input data required for running the kernel. The size of the data to be read in is dependent on the number of processes it was generated on. In this case 24 proceses were used.

``/outputs`` - Contains the output of the benchmark Forces.out which can be used to check correctness.




### To Run: 

* Build the kernel_benchmark in ``~/tools/kernel-benchmark-code``
* Run execubtable in the current working directory
* Correctness can be checked by doing a diff of the generate ``Forces.out`` with ``/output/Forces.out``
* On 1 thread the benchmark should take ~60 seconds on ARCHER
