#!/usr/bin/env python
#!/usr/bin/python
#!python

import sys
import os
import re
import getopt
import numpy
from string import Template

def usage():
	print """
	create makefile for the pipeline of simulations

	 makefile [-h] [-i <input filename>] [-o <output filename>]

	 -h                  ask help

	 -l                  launch by local or slurm

	 -j                  job name

	 """

#option variables
opts = {}
opts['help'] = ''
opts['verbose'] = ''
opts['debug'] = ''
opts['outputDir'] = '/net/frodo/home/khlin/p3'
opts['makeFile'] = 'makefile'
opts['launchMethod'] = 'local'
opts['jobName'] = 'test'


o, a = getopt.getopt(sys.argv[1:], 'l:j:h')
for k,v in o:
	opts[k] = v
	if opts.has_key('-h'):
	    usage(); sys.exit(0)
	elif k == "-l":
		opts['launchMethod'] = v
		if opts['launchMethod'] != "local" and opts['launchMethod'] != "slurm":
			print "Launch method has to be local or slurm"; exit(1)
	elif k == "-j":
		opts['jobName'] = v

##############
#print options
##############
print "Options\n"
print "output directory : {outputDir}".format(**opts)
print "launch method    : {launchMethod}".format(**opts)
print "output jobName   : {jobName}".format(**opts)
print "makefile         : {makeFile}".format(**opts)
print ""

#arrays for storing targets, dependencies and commands
tgts = []
deps = []
cmds = []
 
#temporary variables
tgt = ""
dep = ""
cmd = []
inputFiles = []
inputFilesOK = []
inputFile = ''
outputFile = ''

#create directory needed for slurm script 
if not os.path.exists(opts['outputDir']): os.makedirs(opts['outputDir'])
opts['slurmScriptsDir'] = '{outputDir}/slurm_scripts'.format(**opts)
if not os.path.exists(opts['slurmScriptsDir']): os.makedirs(opts['slurmScriptsDir'])
opts['slurmScriptNo'] = 0
if not os.path.exists(opts['jobName']): os.makedirs(opts['jobName'])

##########
#functions
##########
 
#run a job either locally or by slurm
def makeJob(method, tgt, dep, cmd):
	if method == "local":
		makeLocalStep(tgt, dep, cmd)
	elif method == "slurm":
		makeSlurm(tgt, dep, cmd)
 
#run slurm jobs
def makeSlurm(tgt, dep, cmd):
	tgts.append(tgt)
	deps.append(dep)
	cmd_tmp = []
	for c in cmd:
		opts.update({'command': c})
		#contains pipe or cd
		if re.search('\||^cd', c):
			opts['slurmScriptNo'] += 1
			opts['slurmScriptFile'] = '{slurmScriptsDir}/{slurmScriptNo}_{jobName}.sh'.format(**opts)
			IN = open('{slurmScriptFile}'.format(**opts), 'w')
			IN.write('#!/bin/bash\n')
			IN.write('set pipefail; {command}'.format(**opts))
			IN.close()
			os.chmod('{slurmScriptFile}'.format(**opts), 0755);

			cmd_tmp.append('\tsrun -p nomosix,main -J {jobName} -D {outputDir} {slurmScriptFile} \n'.format(**opts))
		else:
			cmd_tmp.append('\tsrun -p nomosix,main -J {jobName} -D {outputDir} {command} \n'.format(**opts))
	cmd_tmp.append('\ttouch {tgt}\n'.format(tgt=tgt))
	cmds.append(cmd_tmp)
 
#run a local job
def makeLocalStep(tgt, dep, cmd):
	tgts.append(tgt)
	deps.append(dep)
	cmd_tmp = []
	for c in cmd:
		cmd_tmp.append('\t{command}\n'.format(command=c))
	cmd_tmp.append('\ttouch {tgt}\n'.format(tgt=tgt))
	cmds.append(cmd_tmp)

################################################################
#Start the code
################################################################
##initialize
opts['seed']=1000 #number of replications
opts['n_rep_total']=1000
opts['n_rep']=200 #number of replications
opts['n_ite']=opts['n_rep_total']/opts['n_rep'] #number of replications
opts['n_family']=1000 #number of family
opts['p_dis']=0.05 #prevalence

######################
#1.0. log the start time
######################
tgt = '{outputDir}/start.runmake.{jobName}.OK'.format(**opts)
dep = ''
cmd = ['[ ! -f {outputDir}/runmake_{jobName}_time.log ] && echo > {outputDir}/runmake_{jobName}_time.log; date | awk \'{{print "\\nend: "$$0}}\' >> {outputDir}/runmake_{jobName}_time.log'.format(**opts)]
makeJob('local', tgt, dep, cmd)

######################
#1.1. run simulations by calling mainSim.R
######################
inputFilesOK = []
# opts['riskVariant'] = '39' #super rare variant
# for i in numpy.arange(1,2.8,0.1):
# 	for j in range(opts['n_ite']):
# 		opts['r'] = i
# 		tgt = 'callmainSim_{seed}.OK'.format(**opts)
# 		inputFilesOK.append(tgt)
# 		dep = ''
# 		cmd = ['R --vanilla --args seed {seed} n_rep {n_rep} r {r} n_family {n_family} p_dis {p_dis} risk.variant {riskVariant} < mainSim.R > mainSim_{n_rep}_{r}_{n_family}_{p_dis}_{riskVariant}.Rout{seed} 2>&1'.format(**opts)]
# 		makeJob('slurm', tgt, dep, cmd)
# 		opts['seed'] += 1

opts['riskVariant'] = '7' #rare
for i in numpy.arange(1,2.2,0.1):
	for j in range(opts['n_ite']):
		opts['r'] = i
		tgt = 'callmainSim_{seed}.OK'.format(**opts)
		inputFilesOK.append(tgt)
		dep = ''
		cmd = ['R --vanilla --args seed {seed} n_rep {n_rep} r {r} n_family {n_family} p_dis {p_dis} risk.variant {riskVariant} < mainSim.R > mainSim_{n_rep}_{r}_{n_family}_{p_dis}_{riskVariant}.Rout{seed} 2>&1'.format(**opts)]
		makeJob('slurm', tgt, dep, cmd)
		opts['seed'] += 1	

opts['riskVariant'] = '2' #common
for i in numpy.arange(1,1.6,0.1):
	for j in range(opts['n_ite']):
		opts['r'] = i
		tgt = 'callmainSim_{seed}.OK'.format(**opts)
		inputFilesOK.append(tgt)
		dep = ''
		cmd = ['R --vanilla --args seed {seed} n_rep {n_rep} r {r} n_family {n_family} p_dis {p_dis} risk.variant {riskVariant} < mainSim.R > mainSim_{n_rep}_{r}_{n_family}_{p_dis}_{riskVariant}.Rout{seed} 2>&1'.format(**opts)]
		makeJob('slurm', tgt, dep, cmd)
		opts['seed'] += 1	


######################
#1.2. combine the result
######################
tgt = 'pasteResults.OK'
dep = ' '.join(inputFilesOK)
cmd = ['Rscript paste_mainSim_results.R']
makeJob('local', tgt, dep, cmd)

######################
#1.3. copy the results to a folder with name $jobName
######################
tgt = 'backupResult.OK'
dep = 'pasteResults.OK'
cmd = ['cp * {jobName}/'.format(**opts)]
makeJob('local', tgt, dep, cmd)

######################
#1.4. log the end time
######################
tgt = '{outputDir}/end.runmake.{jobName}.OK'.format(**opts)
dep = 'pasteResults.OK'
cmd = ['date | awk \'{{print "Simulations pipeline\\n\\nstart: "$$0}}\' >> {outputDir}/runmake_{jobName}_time.log'.format(**opts)]
makeJob('local', tgt, dep, cmd)

#*******************
#Write out make file
#*******************
MAK = open('{makeFile}'.format(**opts), 'w')
MAK.write('.DELETE_ON_ERROR:\n\n')
MAK.write('all: {tgts}\n\n'.format(tgts=' '.join(tgts)))
 
#clean
tgts.append('clean')
deps.append('')
cmds.append('\trm -f *.OK *.log *.csv *.Rout* *.out* result*.txt mendelian*.txt data*.ped\n')

#clean_job
tgts.append('clean_job')
deps.append('')
cmds.append('\tps xu | awk \'{{if($$11=="make") print $$2}}\' | xargs kill; squeue -ukhlin | awk \'{{if($$3 ~ /p3/) print $$1}}\' | xargs -0 scancel\n'.format(**opts))
 
for tgt,dep,cmd in zip(tgts, deps, cmds):
	MAK.write('{tgt} : {dep}\n'.format(tgt=tgt, dep=dep))
	MAK.writelines(cmd)
	MAK.write('\n')

MAK.close()

