# ParRADMeth
ParRADMeth, Parallel Regression Analysis of Differential Methilation is a parallel algorithm for the identification of differentially methylated regions in genomic analyses. Based on [RADMeth](http://smithlabresearch.org/software/radmeth/), from The Smith Lab.

Prerequisites
-------------

Before starting the installation, please confirm that the following software is available
in your system. Particular versions using during development are shown.

 - **GCC** v8.3.0
 - **GSL** v2.6
 - **make** v3.82
 - **MPI** compiler with support for OpenMP:
     - **OpenMPI** v3.1.4
 - **Git** (Optional)

Different versions of the software may work but they have not been tested.

Compilation
-----------

The project was design so a system level instalation is not needed, but user level compilation and execution can be done. Compilation of the tool for your system architecture can be done by following these steps:

1. **Download**. First, obtain project files by cloning this git repository.
2. **Compilation**. In this step go to the root folder of the project and use ```make```.
3. **Instalation**. In this steps, still on the root folder, use ```make install```. This will place the executable on the ```bin/``` folder.
4. **Cleaning** (Optional). Optionally, still on the root folder, use ```make clean``` to delete files generated during tool's compilation. Only unnecessary files are deleted, so executable files generated still work.

Execution
---------

ParRADMeth must be executed with MPI execution commands (```mpiexec``` (standard) or ```mpirun```).Therefore, the tool can be executed using the following command from the root folder of the proyect:

``` sh
mpiexec -n <numProcs> ./bin/ParRADMeth regression <options> <design_matrix> <proportion_table>
```

when ```numProcs``` is the number of MPI processes to execute, and ```options``` is a list of the following arguments:

 - **-h** (Optional). Print usage and exit.
 - **-v** (Optional). Verbose, print more run information.
 - **-factor \<f\>** (Compulsory). The factor to test.
 - **-o** (Optional). Output option.
 - **design_matrix** (Compulsory). Path to the file containing the design matrix.
 - **proportion_table** (Compulsory). Path to the file containing the proportion table.

License
-------

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
