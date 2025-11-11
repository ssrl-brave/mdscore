#!/bin/bash

mkdir trial_1_boltz trial_2_boltz trial_3_boltz

script=../../script_25/min.sh

PID0=$(sbatch $script|awk '{print $NF}')
echo $PID0


#script=../../script/min.sh

#PID=$(sbatch --dependency=afterany:$PID0 $script|awk '{print $NF}')
#echo $PID

for i in trial_1_boltz trial_2_boltz trial_3_boltz 
do
cd $i

script=../../../script_25/heat.sh

PID_1=$(sbatch --dependency=afterany:$PID0  $script| awk '{print $NF}')
echo $PID_1

script=../../../script_25/npt_cpu.sh

PID_2=$(sbatch --dependency=afterany:$PID_1  $script| awk '{print $NF}')
echo $PID_2

script=../../../script_25/md.sh

PID_3=$(sbatch --dependency=afterany:$PID_2  $script| awk '{print $NF}')
echo $PID_3

cd ..
done

