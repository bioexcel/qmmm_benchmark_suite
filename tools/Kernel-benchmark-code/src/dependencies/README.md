## Dependencies

This directory contains the minimal set of additional source files needed to build the kernel benchmark code. These were taken from the CP2K v6.1 codebase, with the original directory structure under the CP2K /src/ dir respected. Files have been stripped of any code not strictly needed by the kernel benchmark, and of any further dependencies. Further minor modifications have been made, in particular to the cell_create subroutine in /subsys/cell_types.f. CP2K's internal pointer reference counting calls and CPASSERT checks have been removed.

### Building

These dependencies will be built and module files created as part of the make procedure initiated in the parent directory 

