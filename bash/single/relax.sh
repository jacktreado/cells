#!/bin/bash

# directories with code
cellsdir=~/cells
srcdir=$cellsdir/src
maindir=$cellsdir/main

# directory for all output for cell simulations
outputdir=~/project/cells

# directory for simulations specific to jamming
simtypedir=$outputdir/relax

# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p $simtypedir
mkdir -p bin
mkdir -p tasks
mkdir -p slurm
mkdir -p out

# inputs
NV=$1
calA0=$2
kl=$3
kb=$4
th0=$5
partition=$6
time=$7
numSeedsPerRun=$8
numRuns=$9
startSeed="${10}"

let numSeeds=$numSeedsPerRun*$numRuns
let endSeed=$startSeed+$numSeeds-1

# name strings
basestr=relax_NV"$NV"_calA"$calA0"_kl"$kl"_kb"$kb"_th0"$th0"
runstr="$basestr"_startseed"$startSeed"_endseed"$endSeed"

# make directory specific for this simulation
simdatadir=$simtypedir/$basestr
mkdir -p $simdatadir

# compile into binary using packing.h
binf=bin/"$runstr".o
mainf=$maindir/singleCellRelax.cpp
echo Running $numseed relaxation simulations for single cells with the following parameters:
echo . . . NV = $NV
echo . . . calA0 = $calA0
echo . . . kl = $kl
echo . . . kb = $kb
echo . . . th0 = $th0

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

# LOOP OVER FILES. 
for seed in `seq $startSeed $numSeedsPerRun $endSeed`; do
    # count files
    let fcount=$fcount+1

    # echo to console
    echo On base seed $seed

    # echo string of numSeedPerRun commands to task file
    runString="cd `pwd`"

    # loop over seeds to go into runString
    let ssMax=$numSeedsPerRun-1

    # make output directory
    specificdir=$simdatadir/"$basestr"
    mkdir -p $specificdir

    for ss in `seq 0 $ssMax`; do
        # get seed for actual run
        let runseed=$seed+ss

        # get file str
        filestr="$basestr"_seed"$seed"

        # create output files
        posf=$specificdir/$filestr.pos

        # append to runString
        runString="$runString ; ./$binf $NCELLS $NV $sizeDisp $calA0 $kl $ka $kb $runseed $posf $enf $jamf"
    done

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
# 1. NV
# 2. calA0
# 3. kl
# 4. kb
# 5. th0
# 6. partition
# 7. time
# 8. num seeds per run (for each entry in array)
# 9. number of runs (number of array entries, i.e. arraynum)
# 10. start seed (end seed determined by number of runs)




