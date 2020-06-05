#!/bin/bash

# directories with code
cellsdir=~/cells
srcdir=$cellsdir/src
maindir=$cellsdir/main

# directory for all output for cell simulations
outputdir=/gpfs/project/fas/ohern/jdt45/cells

# directory for simulations specific to bidisperse cell jamming (add confluence simulations to jamming directory)
simtypedir=$outputdir/bidcells

# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p $simtypedir
mkdir -p bin
mkdir -p tasks
mkdir -p slurm
mkdir -p out

# inputs
NCELLS=$1
NV=$2
calA0=$3
kl=$4
kb=$5
phiTarget=$6
partition=$7
time=$8
seedStart=$9
seedNum="${10}"

# other inputs
NOUTPUTS=60
let seedMax=$seedStart+$seedNum-1

# name strings
inputstr=bidcells_N"$NCELLS"_NV"$NV"_calA"$calA0"_kl"$kl"_kb"$kb"
basestr=bidconf_N"$NCELLS"_NV"$NV"_calA"$calA0"_kl"$kl"_kb"$kb"_pt"$phiTarget"
runstr="$basestr"_seedStart"$seedStart"_seedMax"$seedMax"

# make directory specific for this simulation
simdatadir=$simtypedir/$inputstr

# compile into binary using packing.h
binf=bin/"$runstr".o
mainf=$maindir/jamming/cellConfluence.cpp
echo Running $numSeeds compression to confluence sims of $NCELLS cells with $NV verts, calA0 = $calA0 , perimeter energy kl = $kl, and bending energy kb = $kb

# run compiler
rm -f $binf
g++ --std=c++11 -I $srcdir $mainf $srcdir/*.cpp -o $binf 
echo compiling with : g++ --std=c++11 -I $srcdir $mainf $srcdir/*.cpp -o $binf 

# check compilation
if [[ ! -f $binf ]]
then
    echo -- binary file does not exist, compilation failed.
    exit 1
fi

# create task file
taskf=tasks/"$runstr".task
rm -f $taskf

# loop over files
let fcount=0

# list of files
flist="$simdatadir"/"$inputstr"_seed*.jam

# loop over files initially to count
for f in $flist; do
    let fcount=$fcount+1
done

if [[ $fcount -eq 0 ]]
then
    echo no files found in flist, ending.
    exit 1
else
    echo $fcount files found in flist, looping...
fi

# reset fcount to 0
let fcount=0

# LOOP OVER FILES. 
for f in $flist; do
    # parse file name
    file=${f##*/}
    baseid=${file%%.jam}
    seed=${baseid#*seed*}
    echo seed = $seed, file = $file

    # echo to console
    echo On base seed $seed

    # check if seed is in correct range
    if [[ $seed -lt $seedStart ]]
    then
        echo seed = $seed too small, skipping...
        continue
    elif [[ $seed -gt $seedMax ]]
    then    
        echo seed = $seed too large, skipping
        continue
    else    
        # increment file count
        let fcount=$fcount+1
        echo seed = $seed is ready for primetime, adding to task file.
    fi

    # echo string of numSeedPerRun commands to task file
    runString="cd `pwd`"

    # append executable to run string
    enf="$simdatadir"/$basestr.en
    posf="$simdatadir"/$basestr.pos
    vdosf="$simdatadir"/$basestr.vdos

    # append to runString
    runString="$runString ; ./$binf $f $calA0 $kl $kb $NOUTPUTS $seed $enf $posf $vdosf"

    # finish off run string
    runString="$runString ;"

    # echo to task file
    echo "$runString" >> $taskf
done

# test if task file was created
if [[ ! -f "$taskf" ]]
then
    echo task file not created, ending before job submission
    exit 1
fi

# get number of jobs to submit to each array
let arraynum=$fcount
echo -- total number of array runs = $arraynum

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
echo \#SBATCH --array=1-$arraynum >> $slurmf
echo \#SBATCH -n 1 >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo \#SBATCH -o $runout >> $slurmf
echo sed -n \"\$\{SLURM_ARRAY_TASK_ID\}p\" "$taskf" \| /bin/bash >> $slurmf
cat $slurmf

# run sbatch file
echo -- running on slurm in partition $partition
echo sbatch -t $time $slurmf


# ====================
#       INPUTS
# ====================
# 1. NCELLS
# 2. NV
# 3. calA0
# 4. perimeter force scale (kl)
# 5. bending energy scale (kb)
# 6. target phi
# 7. partition
# 8. time
# 9. start seed (from list of files)
# 10. number of seeds




