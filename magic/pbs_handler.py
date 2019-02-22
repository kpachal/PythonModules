import os
import subprocess

class PBSHandler(object) :

  def __init__(self,log_path,batch_path) :

    self.log_path = log_path
    self.batch_path = batch_path

  def createBatchScript(self, command, tag):
    scriptName = self.batch_path+"batch_{0}.sh".format(tag)
    with open(scriptName,"w") as f:
      f.write("#!/bin/bash"+"\n")
      f.write("path="+os.getcwd()+"\n")
      f.write("export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase"+"\n")
      f.write("source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh"+"\n")
      f.write("cd $path\n")
      f.write("source setup.sh"+"\n")
      f.write("cd "+os.getcwd()+"\n")
      f.write(command+"\n")
      f.close()
    return scriptName
            
  def send_job(self, command, tag, queue="2nw",ncores=4):
    script = createBatchScript(command,tag)
    batch_cmd = "bsub -J "+tag+" -q "+queue+" -n "+str(ncores)+" < " + script
    print("\t > Submission of {0} on the lxbatch on {1} w/ {2} cores ".format(script,queue,ncores))
    subprocess.call(batch_cmd,shell=True) 