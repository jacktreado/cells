#!/bin/bash

# directories with code
cellsdir=~/cells
srcdir=$cellsdir/src
maindir=$cellsdir/main

# directory for all output for cell simulations
outputdir=~/project/cells

# directory for simulations specific to jamming
simtypedir=$outputdir/gelQS

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
plThresh=$3
kb=$4
aGelation=$5
partition=$6
time=$7
numSeedsPerRun=$8
numRuns=$9
startSeed="${10}"

# other parameters
sizeDisp=0
phiTarget=1.03
dl0=1e-3
kl=1.0
ka=1.0
del=1.0

let numSeeds=$numSeedsPerRun*$numRuns
let endSeed=$startSeed+$numSeeds-1

# name strings
basestr=gelQS_N"$NCELLS"_NV"$NV"_pt"$plThresh"_kb"$kb"_a"$aGelation"
runstr="$basestr"_startseed"$startSeed"_endseed"$endSeed"

# make directory specific for this simulation
simdatadir=$simtypedir/$basestr
mkdir -p $simdatadir

# compile into binary using packing.h
binf=bin/"$runstr".o
mainf=$maindir/gelQsRatchet.cpp
echo Running $numSeeds QS gelation sims of $NCELLS cells with $NV verts, starting from phi = $phiTarget using ratchet springs on perimeters with threshold $plThresh and bending = $kb, adhesion = $aGelation

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

    for ss in `seq 0 $ssMax`; do
        # get seed for actual run
        let runseed=$seed+ss

        # get file str
        filestr="$basestr"_seed"$seed"

        # make output directory
        specificdir=$simdatadir/"$filestr"
        mkdir -p $specificdir

        # create output files
        posf=$specificdir/$filestr.pos
        enf=$specificdir/$filestr.en
        contactf=$specificdir/$filestr.cm

        # append to runString
        runString="$runString ; ./$binf $NCELLS $NV $sizeDisp $phiTarget $plThresh $dl0 $kl $ka $kb $del $aGelation $runseed $posf $enf $contactf"
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
# 3. plThresh (deformation threshold)
# 4. bending energy parameter (kb)
# 5. attraction parameter (a)
# 6. partition
# 7. time
# 8. num seeds per run (for each entry in array)
# 9. number of runs (number of array entries, i.e. arraynum)
# 10. start seed (end seed determined by number of runs)




