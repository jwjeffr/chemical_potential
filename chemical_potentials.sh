config_file_name=$1

cfg () {
	jq -r $1 ${config_file_name}
}

lmp=$(cfg '.exec')
np=$(cfg '.np')
units=$(cfg '.units')
input_data_file=$(cfg '.input_data')

pair_style=$(cfg '.pair_style')
pair_coeff=$(cfg '.pair_coeff')

potential_file=$(cfg '.potential_file')
echo "pair_style ${pair_style}" > ${potential_file}
echo "pair_coeff ${pair_coeff}" >> ${potential_file}

dmax=$(cfg '.dmax')
pressure=$(cfg '.pressure')
vmax=$(cfg '.vmax')
etol=$(cfg '.etol')
ftol=$(cfg '.ftol')
maxsteps=$(cfg '.maxsteps')

ntypes=$(grep -oP '\d+(?=\s*atom types)' ${input_data_file})
occupying=$(cfg '.occupying')

log_file=$(cfg '.log_file')

mpirun -np ${np} ${lmp} -in in.insertions \
  -var units ${units} \
  -var ntypes ${ntypes} \
  -var input_data_file ${input_data_file} \
  -var potential_file ${potential_file} \
  -var dmax ${dmax} \
  -var pressure ${pressure} \
  -var vmax ${vmax} \
  -var etol ${etol} \
  -var ftol ${ftol} \
  -var maxsteps ${maxsteps} \
  -var occupying_energies_file ${occupying} \
  -log ${log_file}