# VortexHelper.jl

This package is designed to simplify the usage and analysis of the vortex-p code (https://github.com/dvallesp/vortex-p).
In particular, it contains julia wrapper function to execute the vortex program and read in the particle data.

## Preparation

In the beginning, the `test_runs` variable has to be set using `set_test_runs`.

All simulation output is assumed to be within directories test_runs/out_cluster_method/.
All vortex output will be located within directories test_runs/vortex_analysis/(filtered_)cluster_method/.

## Running vortex

Use the function `run_vortex`. This creates appropriate vortex parameter files and executes the program.
It requires an executable vortex_mfm/sph_(un)filtered in ./vortex-p/src (or ./vortex-GADGET/src if `new=false`) which has to be compiled before.

## Reading the output

The particle output can be read using `read_particle_data` making use of the FortranFiles package.
