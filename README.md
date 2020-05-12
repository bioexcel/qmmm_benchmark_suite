# QM/MM Benchmark Suite

This benchmark suite contains a number of QM/MM systems designed for benchmarking CP2K and the Grm+CP2K interface.
Each system may contain benchmarks for three different test cases


## Systems

|Name      |QM atoms       |Total atoms |Functional  |Basis set       |MD run type  |Periodic?|
|----------|---------------|------------|------------|----------------|-------------|---------|
|MQAE      |34             |16,396      |BLYP	     |DZVP-MOLOPT-GTH |NVE          |Y        |
|ClC-19    |19             |150,925	    |BLYP	     |DZVP-MOLOPT-GTH |NVE	        |Y        |
|ClC-253   |253            |150,925	    |BLYP	     |DZVP-MOLOPT-GTH |NVE          |Y        |
|CBD_PHY   |68             |167,922	    |PBE         |DZVP-MOLOPT-GTH |NVE          |Y        |
|GFP_QM-77 |20, 32, 54, 77 |28,264      |BLYP        |DZVP-GTH-BLYP   |NVT          |N        |


### MQAE

### ClC-19

### ClC-253

### CBD_PHY

### QMScale

## Test cases

### CP2K MD

Contains the input files to run an MD simulation in native CP2K
