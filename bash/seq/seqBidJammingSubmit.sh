#!/bin/bash

# directories with code
cellsdir=~/cells
srcdir=$cellsdir/src
maindir=$cellsdir/sequential

# directory for all output for cell simulations
outputdir=/gpfs/loomis/project/fas/ohern/jdt45/cells

# directory for simulations specific to jamming
simtypedir=$outputdir/seqBidCells

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
kbb=$6
partition=$7
time=$8
numSeedsPerRun=$9
numRuns="${10}"
startSeed="${11}"

# other variables
dphi=5e-4
Ptol=1e-7
Ftol=1e-11

let numSeeds=$numSeedsPerRun*$numRuns
let endSeed=$startSeed+$numSeeds-1

# name strings
basestr=seqbidcells_N"$NCELLS"_NV"$NV"_calA"$calA0"_kl"$kl"_kb"$kb"_kbb"$kbb"
runstr="$basestr"_startseed"$startSeed"_endseed"$endSeed"

# make directory specific for this simulation
simdatadir=$simtypedir/$basestr
mkdir -p $simdatadir

# compile into binary using packing.h
binf=bin/"$runstr".o
mainf=$maindir/jamming/bidRepulsiveCellJamming.cpp
echo Running $numSeeds jamming sims of $NCELLS cells with $NV verts, bidisperse , calA0 = $calA0 , and bending energy kb = $kb, kbb = $kbb

# run compiler
rm -f $binf
g++ --std=c++11 -O3 -I $srcdir $mainf -o $binf 
echo compiling with : g++ --std=c++11 -O3 -I $srcdir $mainf -o $binf 

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

    for ss in `seq 0 $ssMax`; do
        # get seed for actual run
        let runseed=$seed+ss

        # get file str
        filestr="$basestr"_seed"$seed"

        # create output files
        jamf=$simdatadir/$filestr.jam
        vdosf=$simdatadir/$filestr.vdos

        # append to runString
        runString="$runString ; ./$binf $NCELLS $NV $calA0 $dphi $kl $kb $kbb $Ptol $Ftol $runseed $jamf $vdosf"
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
sbatch -t $time $slurmf


# ====================
#       INPUTS
# ====================
# 1. NCELLS
# 2. NV
# 3. calA0
# 4. perimeter energy (kl)
# 5. bending energy (kb)
# 6. belt energy (kbb)
# 7. partition
# 8. time
# 9. num seeds per run (for each entry in array)
# 10. number of runs (number of array entries, i.e. arraynum)
# 11. start seed (end seed determined by number of runs)




