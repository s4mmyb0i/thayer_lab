#!/bin/bash

#SBATCH --job-name=create_dccm
#SBATCH --output=create_dccm.out
#SBATCH --error=create_dccm.err
#SBATCH --partition=exx96
#SBATCH --nodelist=n82
#SBATCH -B 1:1:1
#SBATCH -N 1
#SBATCH --mail-user=spsinghal@wesleyan.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# A script that creates a DCCM from a trajectory file
# Command ./create_dccm.sh /zfshomes/spsinghal/2025/summer/bharat/mdcrd /zfshomes/spsinghal/2025/summer/bharat/prmtop /zfshomes/spsinghal/2025/summer/bharat/dccm


show_help() {
    echo "Usage: $0 <mdcrd_dir> <prmtop_dir> <output_dir>"
}
if [[ "$1" == "--help" ]]; then
    show_help
    exit 0
fi

mdcrd_dir=$1
prmtop_dir=$2
output_dir=$3

if [[ -z "$mdcrd_dir" || -z "$prmtop_dir" || -z "$output_dir" ]]; then
    echo "Error: Missing arguments"
    show_help
    exit 1
fi

echo "Starting..."
module load amber/16
mkdir -p "$output_dir"

for mdcrd_file in ${mdcrd_dir}/*.mdcrd; do
    echo "Processing $mdcrd_file"
    base_name=$(basename "$mdcrd_file" .mdcrd)
    prmtop_file="${prmtop_dir}/${base_name}.prmtop"
    
    if [[ ! -f "$prmtop_file" ]]; then
        echo "Warning: No matching prmtop file for $mdcrd_file"
        continue
    fi
    
    output_file="${output_dir}/${base_name}_dccm.dat"
    input_script="dccm_${base_name}.in"
    
    cat << EOF > "$input_script"

parm $prmtop_file
trajin $mdcrd_file
reference $mdcrd_file 1

box auto
center @CA
rms reference mass out rmsd_${base_name}.dat
matrix correl @* out $output_file

run
EOF

    cpptraj -i "$input_script"
    if [[ $? -ne 0 ]]; then
        echo "Error processing $mdcrd_file"
        continue
    fi

    rm "$input_script"
    echo "DCCM created for $mdcrd_file with $prmtop_file at $output_file"
    ls "$output_dir"
done

echo "All trajectories have been processed"
echo "Stopped"