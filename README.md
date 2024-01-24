# rfEDIPIC

This is a modified version of the EDIPIC code from Dmytro Sydorenko. This version of the code directly incorporates
into the input files the functionality to have an RF biased left wall.

### Other Modifications

- Output EDF data at each snapshot
- Second model of electron-neutral excitation collisions
- (_in progress_) Output IADF data at each snapshot

### Compiling Instructions (from the original github repo)

To compile, if one has fortran binded with an MPI, for example, do the following:

```console
mpif90 -c parModules.f90
mpif90 -o edipic1d.out par\*f90
```

To run, copy the executable file and the input data files (ssc_\*dat, inluding the cross sections for the selected ion species) into some directory, edit the input data files as necessary (note that the input is formatted, so one has to follow patterns provided in corresponding description lines), then do for example:

```console
mpirun -np 8 ./edipic1d.out > output.txt &  
```
