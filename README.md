# rfEDIPIC

This is a modified version of the EDIPIC code from Dmytro Sydorenko. This version of the code directly incorporates
into input files the functionality to control an RF biased left wall and includes diagnostics designed to analyze
RF plasma processing devices.

For a complete list of diagnostics, operating instructions, and documentation, please refer to the [wiki](https://github.com/budjensen/rfEDIPIC/wiki).

### Quick Start Instructions

To compile, provided that you have fortran binded to MPI, do the following:

```console
mpif90 -c parMT19937.f90 parModules.f90 par*f90
mpif90 -o rfedipic par*o
```

You may need to specify input flags e.g. `-ffree-line-length-none`, depending on your specific fortran and mpi build. For ease, use the included Makefile.

To run, copy the executable file and the input data files (ssc_\*dat, inluding the cross sections for the selected ion species) into some directory, edit the input data files as necessary (note that the input is formatted, so one has to follow patterns provided in corresponding description lines), then do for example:

```console
mpirun -np 8 ./rfedipic > output.txt &  
```
