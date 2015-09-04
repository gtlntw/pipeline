#!/bin/bash
#default job_name = 'test'

while getopts i:j: opt; do
  case $opt in
 	i)
		d=$OPTARG
		;;
	j)
		job_name=$OPTARG
		;;
  esac
done

#command to run
[ -z "$job_name" ] && echo need to give a job name ./runmake.sh -j && exit 1
jobNum=100; #number of parallel jobs
python ./makefile.py -l slurm -j ${job_name};
nohup make -j ${jobNum} --keep-going > "runmake_${job_name}.log" 2>&1 &
