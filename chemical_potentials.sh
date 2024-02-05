#!/usr/bin/bash

config_file_name=$1

cfg () {
	jq -r $1 ${config_file_name}
}

lmp=$(cfg '.exec')
np=$(cfg '.np')

pair_style=$(cfg '.pair_style')
pair_coeff=$(cfg '.pair_coeff')

potential_file=$(cfg '.potential_file')
echo "pair_style ${pair_style}" > ${potential_file}
echo "pair_coeff ${pair_coeff}" >> ${potential_file}

input_data_file=$(cfg '.input_data_file')
ntypes=$(grep -oP '\d+(?=\s*atom types)' ${input_data_file})

# using cfg() doesn't work here? why?
mpi_args=$(jq -r '.mpi_args | join(" ")' ${config_file_name})

envvars_dict=$(cfg '.env_vars')
envvars=""
for key in $(echo ${envvars_dict} | jq -r 'keys_unsorted[]')
do
  value=$(echo ${envvars_dict} | jq -r ".${key}")
  envvars+=" ${key}=${value}"
done

options=$(cfg '.lmp_options')

occupying_energies_file=$(cfg '.occupying_energies_file')

var_str="-var ntypes ${ntypes} \
  -var potential_file ${potential_file} \
  -var input_data_file ${input_data_file} \
  -var occupying_energies_file ${occupying_energies_file} \
  -var relaxed_data_file $(cfg '.relaxed_data_file')"

more_vars=("units" "dmax" "pressure" "vmax" "etol" "ftol" "maxsteps")
vars_dict=$(cfg '.extra_vars')
for var in "${more_vars[@]}"
do
  var_str+=" -var ${var} $(echo ${vars_dict} | jq --arg keyvar "$var" '.[$keyvar]')"
done

env ${envvars} mpirun ${mpi_args} -np ${np} ${lmp} ${options} -in in.insertions -log $(cfg '.log_file') ${var_str}

python chemical_potentials.py ${occupying_energies_file} \
  $(cfg '.axis') \
  $(cfg '.num_atoms_per_chunk') \
  $(cfg '.chemical_potentials_file')