#!/bin/bash
#default: id=NA19625
id=""
launchMethod="slurm"
jobNo=10 #number of jobs in parallel

while getopts i:l: opt; do
  case $opt in
  i)
      id=$OPTARG
      ;;
  l)
      launchMethod=$OPTARG
      ;;  esac
done

#command to run
cmd='job_name=${id}_HGDP;
	 ./makefile_topmed_HGDP.pl -l ${launchMethod} -j ${job_name} -i ${id} -m makefile_${id};
	 nohup make -f makefile_${id} -j${jobNo} --keep-going > log/runmake_${id}.log 2>&1 &'

#if id is not present from the command line then read from the file id.txt 
if [ "$id" == "" ]; then
	while read -r line; do
		#skip comment lines
	    [[ "$line" =~ ^#.*$ ]] && continue
	    id="${line}"
	    #call local ancestry pipeline
	    eval $cmd
	done < "./id.txt"
else
	eval $cmd
fi
