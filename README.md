# chemical_potential

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## About

This repository contains a bash/json interface to calculate chemical potentials over subdomains, which can be used to calculate quantities like chemical potential gradients. This is particularly useful for an out-of-equilibrium configuration, where chemical potential gradients might be a driving force towards equilibrium.

Mathematical details are in [chemical_potential_profile.pdf](https://github.com/jwjeffr/chemical_potential/blob/main/chemical_potential_profile.pdf).

## Tutorial

Most customizability features are accessible through a [JSON](https://en.wikipedia.org/wiki/JSON) configuration file and the bash script [chemical_potentials.sh](https://github.com/jwjeffr/chemical_potential/blob/main/chemical_potentials.sh). An example configuration file is provided in [config.json](https://github.com/jwjeffr/chemical_potential/blob/main/config.json). Run this example with:

```bash
git clone https://github.com/jwjeffr/chemical_potential.git
chmod 755 chemical_potentials.sh
./chemical_potentials.sh config.json
```

# The config file

The input configuration file defines necessary variables to run both LAMMPS and the associated analysis script to post-process the data output from LAMMPS.

The first key (``mpi_args``) defines a list of arguments passed to MPI. This can probably be left blank unless MPI is throwing errors at you. A common error is that the number of slots does not match the number of requested. This can be ignored by providing the ``--oversubscribe`` argument, which can be written into the configuration file as:

```json
"mpi_args": [
    "--oversubscribe"
]
```

Multiple arguments can be passed, hence the list notation.

The second key (``env_vars``) defines a dictionary of environment variables to be passed to the run. An example one is ``OMP_NUM_THREADS``, which specifies the number of OpenMP threads to use. This can be specified with:

```json
"env_vars": {
    "OMP_NUM_THREADS": 1
}
```

Multiple environment variables can be created, hence the dictionary notation.

Other variables define:

- ``np``: The number of MPI processors
- ``exec``: The name of the LAMMPS executable
- ``lmp_options``: Extra options to pass to LAMMPS
    - For example, to run with the [KOKKOS](https://docs.lammps.org/Speed_kokkos.html) accelerator with a GPU backend, one can pass ``"-k on g 1 -sf kk"``.
- ``log_file``: The desired location of the log file
- ``pair_style``: The [pair_style](https://docs.lammps.org/pair_style.html) command that will be provided to LAMMPS
- ``pair_coeff``: The [pair_coeff](https://docs.lammps.org/pair_coeff.html) command that will be provided to LAMMPS
- ``potential_file``: The desired location to write the pair_style and pair_coeff lines
- ``input_data_file``: The input LAMMPS data file to start the simulation
- ``relaxed_data_file``: The desired location of the initial data file after relaxing
- ``occupying_energies_file``: The desired location of the occupying energies file
- ``extra_vars``: Extra variables passed into LAMMPS, mostly related to minimization settings.

Then, the analysis script [chemical_potentials.py](https://github.com/jwjeffr/chemical_potential/blob/main/chemical_potentials.py) uses the final three variables specified in the config file. They are defined by:

- ``axis``: The desired axis to create the chemical potential profile over
- ``num_atoms_per_chunk``: The number of atoms to include in each chunk when averaging
- ``chemical_potentials_file``: The desired file to store the chemical potentials in. This is formatted like:

```
# center mu_1 mu_2 ...
z mu_1(z) mu_2(z) ...
.
.
.
```

where ``mu_i(z)`` is the chemical potential of ``i`` at position ``z``. The output file created by the example input file is at [example/chemical_potentials.txt](https://github.com/jwjeffr/chemical_potential/blob/main/example/chemical_potentials.txt).

## Acknowledgements

The  work  was  supported  by  the  grant  DE-SC0022980 funded by the U.S. Department of Energy,  Office of Science.

This material is based on work supported by the National Science Foundation under Grant Nos. MRI# 2024205, MRI# 1725573, and CRI# 2010270.

## Disclaimer

Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.