#!/bin/bash

# directory for all output for cell simulations
projectdir=~/project/cells

# directory for simulations specific to jamming
simtypedir=$outputdir/gelQS

# make directories, unless they already exist
mkdir -p tasks
mkdir -p slurm
mkdir -p out

# inputs
NCELLS=$1
NV=$2
asphericity=$3
a=$4
partition=$5
time=$6
numSeeds=$7
startSeed=$8
movie=$9

# ending seed, 1 seed per array run
let numSeedsPerRun=1
let endSeed=$startSeed+$numSeeds-1

# name strings
basestr=gelQS_N"$NCELLS"_NV"$NV"_calA"$asphericity"_a"$a"
runstr="$basestr"_PROCESS_startseed"$startSeed"_endseed"$endSeed"

# label directory specific for this simulation
simdatadir=$simtypedir/$basestr

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

        # label specific data directory
        specificdir=$simdatadir/"$filestr"

        # create output files
        posf=$specificdir/$filestr.pos
        savef=$specificdir/$filestr.mat
        if [[ $movie -eq 1 ]]
        then
        	# name of movie file
        	movief=$specificdir/$filestr.mp4

        	# code to run MATLAB
        	MCODE="gel2DSimPostprocess('$posf','$savef','$movief'); quit;"
        else
        	MCODE="gel2DSimPostprocess('$posf','$savef'); quit;"
        fi

        # append to runString
        runString="$runString ; matlab -nodisplay -r \"$MCODE\""
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
# 3. asphericity
# 4. a (attraction)
# 5. partition
# 6. time
# 7. num seeds (i.e. number of array entries total)
# 8. start seed (end seed determined by number of runs)
# 9. write movie (1 = YES, 0 = NO)




