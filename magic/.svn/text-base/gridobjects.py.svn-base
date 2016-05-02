# This class helps the user submit jobs to the grid,
# determine when they have finished, download the output,
# and validate it as quickly and reliably as possible.
#
# Much of the structure stolen from the legendary
# Thomas Gillam, with thanks.


# Standard imports
import sys
import subprocess
import time
import os
import signal
from dq2.clientapi.DQ2 import DQ2
from os import environ, path

class GridJobset(object) :

  def __init__(self,command) :

    # Set defaults for statuses
    self.JobsetID = None
    self.JobIDs = []
    self.runCommand = command
    self.jobIdFile = '%s/pjobid.dat' % os.environ['PANDA_CONFIG_ROOT']

  def retrieveJobIDs(self) :
    '''Retrieve JobsetID, JobIDs of new job(s)'''

    jobIDList = []
    readableJobIdFile = open(self.jobIdFile)
    for line in readableJobIdFile :
      line.strip()
      jobIDList.append(long(line))
    readableJobIdFile.close()
    jobIDList.sort()
    self.JobsetID = jobIDList[0]-1
    self.JobIDs = jobIDList

  def submit(self) :
    '''Calls prun. Retrieves ID of newly created job(s).'''

    # Submit jobs
    subprocess.call(self.runCommand, shell=True)
    self.retrieveJobIDs()


class GridJob(object) :

  ## ----------------------------------------------------
  ## Initialisers

  def __init__(self, jobinfo, outputdir, dq2SetupScript, downloadTimeout, downloadLimit) :

    # jobinfo object has attributes:
    # ['id','JobID','PandaID','jobStatus','site','cloud','jobType',
    #  'jobName','inDS','outDS','libDS','provenanceID','creationTime',
    #  'lastUpdate','jobParams','dbStatus','buildStatus','retryID',
    #  'commandToPilot']
    self.JobInfo = jobinfo
    self.JobID = jobinfo.JobID
    self.inDS = jobinfo.inDS
    self.outDS = jobinfo.outDS
    self.site = jobinfo.site

    self.outputdir = outputdir
    self.dq2SetupScript = dq2SetupScript
    self.dq2process = None
    self.dq2processStart = 0
    self.status = ['None']
    self.statusSince = time.time()

    self.prunAttemptCount = 0
    self.dq2AttemptCount = 0
    self.defineDownloadAsStuck = downloadTimeout
    self.dq2RetryLimit = downloadLimit

    self.isCurrentlyDownloading = False
    self.isDownloadFinished = False
    self.isDownloadFailed = False
    self.corruptedFiles = []

    self.createDQ2Client()

  def createDQ2Client(self) :
    self.dq2Client = DQ2()

  ## ----------------------------------------------------
  ## Constituent functions

  def retrieveOutput(self) :
    '''Retrieves finished job via dq2-get. Starts a subprocess.'''

    # If already tried the maximum number of times, don't do this again.
    if self.dq2AttemptCount >= self.dq2RetryLimit :
      self.isDownloadFinished = True
      self.isDownloadFailed = True
      return

    self.dq2AttemptCount+=1
    self.isCurrentlyDownloading = True

    # Make sure directory for download exists
    if not os.path.isdir(self.outputdir):
      os.makedirs(self.outputdir)

    # Format dq2 command.
    command = '. '+self.dq2SetupScript+'; '
    command += 'dq2-get '+self.outDS
    print "Executing",command

    # This formatting stolen from Tom.
    process = subprocess.Popen(command, cwd=self.outputdir, \
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, \
        preexec_fn=os.setsid)
    self.dq2process = process
    self.dq2processStart = time.time()

  def validateOutput(self) :
    '''Make sure all information is in finished job.'''

    # Clear every time we re-validate.
    self.corruptedFiles = []

    # We can find out from dq2 what the output is supposed to look like.
    theoreticalFiles = {}
    datasets = self.dq2Client.listDatasetsInContainer(self.outDS)
    for dataset in datasets :
      dictofoutfiles = {}
      rootFiles = {}
      for item in self.dq2Client.listFilesInDataset(dataset) : # (content, timestamp)
        if isinstance(item, dict) :
          dictofoutfiles = item
      for handle in dictofoutfiles.keys() :
        filedata = dictofoutfiles[handle]
        filename = filedata['lfn']
        filesize = filedata['filesize']
        if "root" in filename:
          rootFiles[filename] = filesize
      theoreticalFiles[dataset] = rootFiles

    # Now allRootFiles should contain a dict of every root file expected
    # in download, with its size, sorted by containing dataset name.
    # Question is: are they all there?
    actualFiles = {}
    for dataset in theoreticalFiles.keys() :
      # Did it get downloaded at all?
      if dataset in os.listdir(self.outputdir) :
        Files = {}
        for fileName in [x for x in os.listdir(os.path.join(self.outputdir,dataset)) \
              if os.path.isfile(os.path.join(os.path.join(self.outputdir,dataset),x))]:
          Files[fileName] = os.path.getsize(os.path.join(os.path.join(self.outputdir,dataset),fileName))
        actualFiles[dataset] = Files

    # Possibilities:
    # Entire dataset missing.
    if len(actualFiles.keys()) < len(theoreticalFiles.keys()) :
      return False

    # Files within dataset missing.
    for datasetName in theoreticalFiles.keys() :
      theoreticalDataset = theoreticalFiles[datasetName]
      actualDataset = actualFiles[datasetName]
      shouldBeFiles = theoreticalDataset.keys()
      areFiles = actualDataset.keys()
      if len(areFiles) < len(shouldBeFiles) :
        print "Files in dataset", datasetName,"missing!"
        return False

    # Files present but wrong size.
      else :
        for file in shouldBeFiles :
          shouldBeFileSize = theoreticalDataset[file]
          actualFileSize = actualDataset[file]
          if actualFileSize != shouldBeFileSize :
            self.corruptedFiles.append([datasetName,file])
            print "File",file,"corrupted!"
            print "Actual size:",actualFileSize
            print "Should be size:",shouldBeFileSize
            return False

    # If we made it here, everything is there!
    return True

  def killThisDownload(self) :
    '''Stop a running subprocess'''
    if self.dq2process != None:
      os.killpg(self.dq2process.pid, signal.SIGTERM)
    self.dq2process = None
    self.dq2processStart = 0

  def cleanAndRestart(self) :
    '''Remove corrupted files and restart dq2'''
    # Reset flags
    self.isDownloadFinished=False
    # Delete corrupted files dataset folder
    for file in self.corruptedFiles :
      datasetname = file[0]
      filename = file[1]
      os.remove(os.path.join(os.path.join(self.outputdir,datasetname), filename))
    # Now restart dq2-get.
    self.retrieveOutput()

  def checkDQ2(self) :
    '''Check progress of dq2. Request restart if download stuck.'''
    isComplete = self.dq2process.poll()
    if isComplete == None :
      if time.time()-self.dq2processStart > self.defineDownloadAsStuck :
        print "Killing and retrying download attempt."
        self.killThisDownload()
        self.retrieveOutput()
    if isComplete != None :
      self.isDownloadFinished = True

  def checkDownloadProgress(self) :
    '''Report if download complete.'''

    # Check how download is going
    self.checkDQ2()

    # Validate. Report complete when finished.
    if self.isDownloadFinished==True :
      if self.isDownloadFailed==False :
        allFilesOK = self.validateOutput()
        if (allFilesOK):
          return 'complete'
        else :
          self.cleanAndRestart()
          return 'incomplete'
      else :
        return 'complete'

    return 'incomplete'

