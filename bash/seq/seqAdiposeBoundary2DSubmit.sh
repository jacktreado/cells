#!/bin/bash

# directories with code
cellsdir=~/cells
srcdir=$cellsdir/src
maindir=$cellsdir/sequential

# directory for all output for cell simulations
outputdir=/gpfs/loomis/project/fas/ohern/jdt45/cells

# directory for simulations specific to sims
simtypedir=$outputdir/adiposeBoundary2D

# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p $simtypedir
mkdir -p bin
mkdir -p tasks
mkdir -p slurm
mkdir -p out

# inputs
aN=$1
tcalA0=$2
l1=$3
l2=$4
v0=$5
dDr=$6
partition=$7
time=$8
numRuns=$9
startSeed="${10}"

# other variables
areaRatio=25
NV=24
NT=1e7
NASKIP=5e4
numSeedsPerRun=1

let numSeeds=$numSeedsPerRun*$numRuns
let endSeed=$startSeed+$numSeeds-1

# name strings
basestr=invasion_aN"$aN"_tc"$tcalA0"_l1"$l1"_l2"$l2"_v0"$v0"_dDr"$dDr"
runstr="$basestr"_startseed"$startSeed"_endseed"$endSeed"

# make directory specific for this simulation
simdatadir=$simtypedir/$basestr
mkdir -p $simdatadir

# compile into binary using packing.h
binf=bin/"$runstr".o
mainf=$maindir/active/adiposeBoundary2D.cpp
echo Running $numSeeds adipose boundary invasion sims of $aN adipocytes 
echo - - - Area ratio                   $areaRatio
echo - - - tumor cells shape param      $tcalA0
echo - - - l1                           $l1
echo - - - l2                           $l2
echo - - - v0                           $v0
echo - - - NT                           $NT
echo - - - dDr                          $dDr
echo - - - NT                           $NT
echo - - - NASKIP                       $NASKIP

# run compiler
rm -f $binf
g++ --std=c++11 -O3 $mainf -o $binf 
echo compiling with : g++ --std=c++11 -O3 $mainf -o $binf  

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
        posf=$simdatadir/$filestr.pos

        # append to runString
        runString="$runString ; ./$binf $aN $areaRatio $NV $tcalA0 $l1 $l2 $v0 $dDr $NT $NASKIP $seed $posf"
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
# 1.    aN
# 2.    tumor calA0
# 3.    l1
# 4.    l2
# 5.    v0 
# 6.    dDr
# 7.    partition
# 8.    time
# 9.    number of runs (number of array entries, i.e. arraynum)
# 10.   start seed (end seed determined by number of runs)










