#!/bin/bash

# directory for all output for cell simulations
projectdir=~/project/cells

# directory for simulations specific to jamming
simtypedir=$projectdir/hopper

# make directories, unless they already exist
mkdir -p slurm
mkdir -p out

# inputs
NCELLS=$1
NT=$2
sizeDisp=$3
g=$4
w0=$5
w=$6
partition=$7
time=$8

# name strings
basestr=hopper2D_SP_N"$NCELLS"_NT"$NT"_sig"$sizeDisp"_g"$g"_w0"$w0"_w"$w"
runstr="$basestr"_PP

# label directory specific for this simulation
simdatadir=$simtypedir/$basestr

# save file
savef=$simdatadir/"$basestr".mat

# code to run matlab
MCODE="hopper2D_SP_PostProcess('$simdatadir','$savef'); quit;"

# setup slurm files
slurmf=slurm/"$runstr".slurm
job_name="$runstr"
runout=out/"$runstr".out
rm -f $slurmf

# echo about time
echo -- running time = $time for $partition

echo -- PRINTING SLURM FILE...
echo \#\!/bin/bash >> $slurmf
echo \#SBATCH --cpus-per-task=1 >> $slurmf
echo \#SBATCH -n 1 >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo \#SBATCH -o $runout >> $slurmf
echo module load MATLAB >> $slurmf
echo matlab -nodisplay -r "\"$MCODE\"" >> $slurmf
cat $slurmf

# run sbatch file
echo -- running on slurm in partition $partition
sbatch -t $time $slurmf


# ====================
#       INPUTS
# ====================
# 1. NCELLS
# 2. NT
# 3. size dispersion of particles
# 4. forcing strength g
# 5. hopper reservoir width w0
# 6. hopper orifice width w
# 7. partition
# 8. time




