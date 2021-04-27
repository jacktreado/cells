#!/bin/bash

# directory for all output for cell simulations
outputdir=/gpfs/loomis/project/fas/ohern/jdt45/cells

# directory for simulations specific to jamming
simtypedir=$outputdir/singleDP

# make directories, unless they already exist
mkdir -p $outputdir
mkdir -p $simtypedir
mkdir -p bin
mkdir -p tasks
mkdir -p slurm
mkdir -p out

# inputs
NV=$1
kl=$2
kb=$3
partition=$4
time=$5
startDCalA0=$6
endDCalA0=$7
dDCalA0=$8

# create numerical counters for calA0
catmp=$(echo "scale=8; $startDCalA0" | bc)
camax=$(echo "scale=8; $endDCalA0" | bc)
dca=$(echo "scale=8; $dDCalA0" | bc)

# name strings
basestr=singleDP_NV"$NV"_kl"$kl"_kb"$kb"
runstr="$basestr"_startDCa"$startDCalA0"_endDCa"$endDCalA0"

# make directory specific for this simulation
simdatadir=$simtypedir/$basestr
mkdir -p $simdatadir

echo Running $numSeeds sims of single particle relaxation
echo - - - NV                           $NV
echo - - - kl 							$kl
echo - - - kb                           $kb
echo - - - start DCalA0                 $startDCalA0
echo - - - end DCalA0                   $endDCalA0
echo - - - step DCalA0                  $dDCalA0


# create task file
taskf=tasks/"$runstr".task
rm -f $taskf

# file indices
let fcount=0
let k=1
let kmax=1000

# LOOP OVER FILES
while [[ $(echo "$catmp < $camax" | bc -l) -eq 1  && $k -lt $kmax ]]; do
    # count files
    let fcount=$fcount+1

    # get calA0
    calA0=$(echo "scale=8; 1.0 + $catmp" | bc)

    # echo string command to task file
    runString="cd `pwd`"

    # get file str
    filestr="$basestr"_ca"$calA0"

    # create output files
    savef=$simdatadir/$filestr.mat

    # append to runString
    runString="$runString ; matlab -nodisplay -r \"singleParticleRelaxation($NV,$kl,$kb,$calA0,'$savef'); quit;\" ; "
    echo "$runString" >> $taskf

    # print to console
    echo k = $k, catmp = $catmp, calA0 = $calA0, condition = $(echo "$catmp < $camax" | bc -l)

    # update DCalA0
    catmp=$(echo "scale=8; $catmp*$dca" | bc)

    # update k
    let k=$k+1
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
# 2. kl
# 3. kb
# 4. partition
# 5. time
# 6. start DCalA0 (= 1 - calA0)
# 7. end DCalA0
# 8. step DCalA0












