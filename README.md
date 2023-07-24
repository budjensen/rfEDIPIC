# rfEDIPIC

This is a modified version of the EDIPIC code from Dmytro Sydorenko. This version of the code directly incorporates
into the inpu files the functionality to have an RF biased left wall.

### Other Modifications

- (_in progress_) Output IADF data at each snapshot

### Compiling Instructions (from the original github repo)

To compile, if one has gcc and intel fortran binded with an MPI, for example, do the following:

```console
gcc -c WELL19937a\_new.c
mpif90 -c parModules.f90
mpif90 -o edipic1d.out par\*f90 WELL19937a\_new.o
```

To run, copy the executable file and the input data files (ssc_\*dat, inluding the cross sections for the selected ion species) into some directory, edit the input data files as necessary (note that the input is formatted, so one has to follow patterns provided in corresponding description lines), then do for example:

```console
mpirun -np 8 ./edipic1d.out > output.txt &  
```
