Inspired by the clandpol package and the paper (https://doi.org/10.1021/acs.jctc.3c00278), I modified the `polarizer` for adding symmetrized bonds and write a `symLJ` to symmetrize the LJ parameters

The working flow is
1) use `clandp` and `fftool` packages (please refer to https://github.com/paduagroup/clandp and https://github.com/paduagroup/fftool) to generate an unpolarized lammps data file and pair.lmp (including LJ parameters)
2) use `polarizer` in this folder to generate symmetrized polarized lammps data file, this script will generate another file named `atom.info` for following steps to symmetrize LJ parameters.
3) use `scaleLJ` in this folder (forked from `clandpol` package) to scale LJ parameters by SAPT calculated parameters
4) use `symLJ` in this folder to symmetrize the LJ parameters gained from step 3.

Please reach out to me (sijiachen@uchicago.edu) if you find any bugs in the scripts.
