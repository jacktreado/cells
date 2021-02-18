#!/bin/bash

# directory for all output for cell simulations
outputdir=/gpfs/loomis/project/fas/ohern/jdt45/cells

# directory for simulations specific to jamming
simtypedir=$outputdir/activeBrownianCrawler

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
v0=$5
Dr=$6
NT=$7
dt=$8
NFRAMES=$9
partition="${10}"
time="${11}"
numRuns="${12}"
startSeed="${13}"

# other variables
numSeedsPerRun=1

let numSeeds=$numSeedsPerRun*$numRuns
let endSeed=$startSeed+$numSeeds-1

# name strings
basestr=abc_NV"$NV"_ca"$calA0"_kl"$kl"_kb"$kb"_v0"$v0"_Dr"$Dr"
runstr="$basestr"_startseed"$startSeed"_endseed"$endSeed"

# make directory specific for this simulation
simdatadir=$simtypedir/$basestr
mkdir -p $simdatadir

# compile into binary using packing.h
binf=bin/"$runstr".o
mainf=$maindir/shapeMobility.cpp
echo Running $numSeeds sims of ABC particles 
echo - - - NV                           $NV
echo - - - calA0 						$calA0
echo - - - kl 							$kl
echo - - - kb                           $kb
echo - - - v0                           $v0
echo - - - Dr                           $Dr
echo - - - NT                           $NT

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
        savef=$simdatadir/$filestr.mat

        # append to runString
        runString="$runString ; matlab -nodisplay -r \"activeBrownianCrawler($NV,$calA0,$kl,$kb,$v0,$Dr,$NT,$dt,$NFRAMES,$seed,'$savef')\"; quit"
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
echo module load MATLAB >> $slurmf
echo sed -n \"\$\{SLURM_ARRAY_TASK_ID\}p\" "$taskf" \| /bin/bash >> $slurmf
cat $slurmf

# run sbatch file
echo -- running on slurm in partition $partition
sbatch -t $time $slurmf




# ====================
#       INPUTS
# ====================
# 1. NV
# 2. calA0
# 3. kl
# 4. kb
# 5. v0
# 6. Dr
# 7. NT
# 8. dt
# 9. NFRAMES
# 10. partition
# 11. time
# 12. number of runs (number of array entries, i.e. arraynum)
# 13. start seed (end seed determined by number of runs)












