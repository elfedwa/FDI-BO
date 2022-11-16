#!/bin/bash

date

JOBID=`sbatch ~/HPC/run.bash  > vasp_aw.log`
JOBID=`echo $JOBID | awk {'print $4'}`

for i in {1..100};do

    sleep 1

    squeue | grep $JOBID > /dev/null

    if [[ $? -ne 0]];then

	echo "job completed"
	break;

    fi

done

date
