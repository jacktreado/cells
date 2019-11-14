#!/bin/bash

# directory for all output for cell simulations
projectdir=~/project/cells

# directory for simulations specific to quasi-static compression protocols
simtypedir=$projectdir/jamming

# matfile directory
matfiledir=$simtypedir/matfiles

# make directories, unless they already exist
mkdir -p slurm
mkdir -p out
mkdir -p "$matfiledir"

# inputs
NCELLS=$1
NV=$2
calA=$3
kl=$4
ka=$5
partition=$6
time=$7

# name strings
basestr=qscomp_N"$NCELLS"_NV"$NV"_calA"$calA"_kl"$kl"_ka"$ka"
runstr="$basestr"_PROCESS

# label directory specific for this simulation
simdatadir=$simtypedir/$basestr

# create output files
savef="$basestr"_pp.mat

# create matlab code string
MCODE="qscompPostProcess('$simdatadir','$savef'); quit;"

# setup slurm files
slurmf=slurm/"$runstr".slurm
job_name="$runstr"
runout=out/"$runstr"-%a.out
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
echo "matlab -nodisplay -r \"$MCODE\"" >> $slurmf
cat $slurmf

# run sbatch file
echo -- running on slurm in partition $partition
sbatch -t $time $slurmf


# ====================
#       INPUTS
# ====================
# 1. NCELLS
# 2. NV
# 3. calA
# 4. kl
# 5. ka
# 6. partition
# 7. time
