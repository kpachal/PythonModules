import os
import subprocess

class CondorHandler(object) :

  def __init__(self,log_path,batch_path) :

    self.log_path = log_path
    self.batch_path = batch_path

    self.job_length = "workday"

    self.email = 'katherine.pachal@cern.ch'

  def send_job(self,command,tag) :

    # make files
    bashfile = self.make_bash_file(command, tag)
    jobfile = self.make_job_file(bashfile, tag)
    # do submit thing
    subprocess.call("condor_submit {0}".format(jobfile),shell=True)

  def make_bash_file(self,command, tag) :

    runFile = self.batch_path+"batch_{0}.sh".format(tag)

    queue = 'short.q'

    with open(runFile,"w") as fr :
      fr.write('#!/bin/sh\n')
      fr.write('# '+tag+' batch run script\n')
      fr.write('#$ -cwd\n')
      fr.write('#$ -j y\n')
      fr.write('#$ -l cvmfs\n')
      fr.write('#$ -M '+self.email+'\n')
      fr.write("path="+os.getcwd()+"\n")      
      fr.write('export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase\n')
      fr.write('source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n')
      fr.write('pwd; ls -l\n')
      fr.write(command+"\n")
      fr.write('ls -l\n')

    print ("Made run file",runFile)
    return runFile

  def make_job_file(self,runFile, tag) :

    batchFile = self.batch_path+"batch_{0}.job".format(tag)
    with open(batchFile, "w") as fsubcondor :
      fsubcondor.write('Universe        = vanilla\n')
      fsubcondor.write('Executable      = '+runFile+'\n')
      fsubcondor.write('+JobFlavour     = "{0}"\n'.format(self.job_length)) # 8 hours is default
      fsubcondor.write('Output          = {0}/stdout_{1}.txt\n'.format(self.log_path,tag))
      fsubcondor.write('Error           = {0}/stderr_{1}.txt\n'.format(self.log_path,tag))
      fsubcondor.write('log             = {0}/batch_{1}.log\n'.format(self.log_path,tag))
      fsubcondor.write('\nqueue 1\n')

    print ("Made job file",batchFile)
    return batchFile
