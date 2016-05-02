# Kate Pachal, Feb 2013
#
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
from os import environ, path
from pandatools import PdbUtils
from gridobjects import GridJobset, GridJob

# Import pbook
tmp, __name__ = __name__, 'pbook'
execfile(path.join(environ["ATLAS_LOCAL_ROOT_BASE"], "x86_64/PandaClient/current/bin/pbook") )
__name__ = tmp

class Juggernaut(object) :

  ## ----------------------------------------------------
  ## Initialisers
 
  def __init__(self) :
    self.commandList = []
    self.outputdir = "/data/atlas/atlasdata/pachal/dijets"
    self.useNewSite = False
    self.downtime = 30 # Rest time in seconds
    self.defineJobAsStuck = 7200 # Two hours
    self.defineDownloadAsStuck = 3600 # One hour
    self.pandaRetryLimit =3 
    self.dq2RetryLimit = 3
    self.maxDQ2Streams = 3
    self.nDQ2StreamsInProgress = 0
    self.dq2SetupScript = "/home/pachal/scripts/dq2_setup.sh"

    self.syncToRunningJobs = False
    self.additionalJobs = []

    self.createPbook()
    self.createPdbUtils()

  def createPbook(self) :
    enforceEnter = False
    verbose   = False
    restoreDB = False
    self.pbookCore = PBookCore(enforceEnter, verbose, restoreDB)

  def createPdbUtils(self) :
    verbose = False
    PdbUtils.initialzieDB(verbose)

  ## ----------------------------------------------------
  ## Setters

  # Specify all the prun commands you will use
  def setCommandList(self, commandList) :
    self.commandList = commandList

  # Add all currently running jobs to monitored list
  def useRunningJobs(self, sync) :
    self.syncToRunningJobs = sync

  # Specify a list of previously-defined job numbers 
  # to monitor and download
  def syncToSpecifiedJobs(self, getJobsList) :
    self.additionalJobs = getJobsList

  # Specify name of your local directory where output will be sent
  def setOutputDirectory(self, outputdir) :
    self.outputdir = outputdir

  # You don't want this if you're using slims 
  # stored at only two or three sites
  def setUseNewSite(self, doUseNew) :
    self.useNewSite = doUseNew

  # Maximum number of subjobs dq2-getting things at once
  def setMaxDQ2Streams(self, maxStreams) :
    self.maxDQ2Streams = maxStreams

  # Time at which a job not yet running is considered stuck
  def setTimeToDefineStuckJob(self, timeout) :
    self.defineJobAsStuck = timeout

  # Maximum amount of time to try to dq2-get before killing subprocess
  def setDownloadTimeout(self, timeout) :
    self.defineDownloadAsStuck = timeout

  # Setup script for dq2
  # NOTE: Don't change, if possible. Script also needs to set up new
  # version of python compatible with dq2.
  # Source: http://atlas-sw.cern.ch/cgi-bin/cvsweb.cgi/dq2.clientapi/ 
  # lib/dq2/clientapi/DQ2.py?rev=1.44;content-type=text%2Fplain;
  # cvsroot=atlas-dq2;only_with_tag=MAIN
  def setDQ2SetupScript(self, script) :
    self.dq2SetupScript = script

  ## ----------------------------------------------------
  ## Main function

  def execute(self) :

    # Hold output jobs
    self.failedJobs = []
    self.successfulJobs = []

    # Run all commands to submit new jobsets to the grid.
    # Retrieve the JobIDs and stats of all newly created jobs.
    print "About to submit requested jobs"
    self.currentJobs = {}
    for item in self.commandList :
      newjobset = GridJobset(item)
      newjobset.submit()
      createdJobs = newjobset.JobIDs
      # Register newly existing JobIDS
      self.pbookCore.sync()
      PdbUtils.getListOfJobIDs(True,False)
      for jobID in createdJobs :
        self.addJobToList(jobID)
    print "Submitted jobs",self.currentJobs.keys()
    
    # If specified jobs to pick up, get those:
    if len(self.additionalJobs)!=0 :
      print "Adding specified jobs to list."
      for item in self.additionalJobs :
        self.addJobToList(item)

    # If we are supposed to include all running jobs, add those
    if self.syncToRunningJobs==True :
      print "Adding currently running jobs to list."
      # Param 'True' means only jobs not 'frozen' are kept
      runningJobs = PdbUtils.getListOfJobIDs(True,False)
      for item in runningJobs :
        self.addJobToList(item)

    print "Total list of jobs to monitor is now:",self.currentJobs.keys(),"\n"

    # Run this until all jobs are complete.
    while len(self.currentJobs.keys())>0 :

      # Synchronise pbook.
      self.pbookCore.sync()

      # Check each job and act accordingly.
      self.currentJobIDs = sorted(self.currentJobs.keys())
      for jobID in self.currentJobIDs :
        job = self.currentJobs[jobID]
        currentStatus = self.checkCurrentStatus(job)

        if currentStatus == 'stillRunning' :
          continue

        elif currentStatus == 'stuck' :
          self.unstick(job)

        elif currentStatus == 'failed' :
          if job.prunAttemptCount < self.pandaRetryLimit :
            self.retryFailed(job)
          else :
            self.failedJobs.append[jobID]
            del self.currentJobs[jobID]

        elif currentStatus == 'finished' :
          ## If running a test code which does not produce an
          ## output dataset, the outDS will be blank.
          if job.outDS == "" :
            del self.currentJobs[jobID]
            self.successfulJobs.append(jobID)
          ## dq2-get output.
          else :
            self.getOutput(job)

        else :
          print "Error!"
          self.currentJobs = {}
          break

      # Wait required gap time, then go to next iteration.
      print "\n"
      time.sleep(self.downtime)

    print "All jobs finished."
    print "Successful jobs:",self.successfulJobs
    print "Failed jobs:",self.failedJobs

    sys.exit(0)

  ## ----------------------------------------------------
  ## Constituent functions

  def addJobToList(self, jobID) :
    '''Get job info for jobID and add GridJob to self.currentJobs'''
    self.pbookCore.sync()
    if not isinstance(jobID, str) :
      jobID = '%d' % jobID
    if jobID not in self.currentJobs.keys() :
      jobinfo = PdbUtils.readJobDB(jobID,False)
      if jobinfo==None :
        print "No job found with ID",jobID,"!"
        self.failedJobs.append(jobID)
        return
      newjob = GridJob(jobinfo,self.outputdir,self.dq2SetupScript,\
               self.defineDownloadAsStuck,self.dq2RetryLimit)
      self.currentJobs[jobID] = newjob

  def performRetry(self,job) :
    '''Retry a previous job. Update stored information.'''
    # Get current info for this job
    oldJobID = job.JobID
    jobInfo = PdbUtils.readJobDB(oldJobID,False)
    # Determine whether we need to redo build
    if jobInfo.buildStatus in ['','finished'] :
      print "Build OK."
      self.pbookCore.retry(oldJobID,newSite=self.useNewSite)
    else :
      print "Retrying build."
      self.pbookCore.retry(oldJobID,retryBuild=True,newSite=self.useNewSite)
    # Retrieve ID of new job
    jobInfo = PdbUtils.readJobDB(oldJobID,False)
    newJobID = jobInfo.retryID
    if newJobID == 0 :
      newJobID == oldJobID
    # Remove old information in currentJobs
    if not isinstance(oldJobID,str) :
      oldJobID = '%d' % oldJobID
    del self.currentJobs[oldJobID]
    # Put current information in currentJobs
    self.addJobToList(newJobID)
    print "Retrying job",oldJobID,". New JobID:",newJobID

  def unstick(self,job) :
    '''Kill and retry stuck job. Replace old jobID with new one.'''
    idToKill = job.JobID
    self.pbookCore.kill(idToKill)
    self.performRetry(job)

  def retryFailed(self, job) :
    '''Retry a failed job once db entry is frozen.'''
    print "retrying failed job"
    #Is anything in the job still running? If so, kill it.
    jobInfo = PdbUtils.readJobDB(job.JobID,False)
    if jobInfo.dbStatus!='frozen' :
      self.pbookCore.kill(job.JobID)
    self.performRetry(job)
    job.prunAttemptCount += 1

  def getOutput(self, job) :
    '''Start process of dq2-getting output from completed job'''
    if job.isCurrentlyDownloading :
      status = job.checkDownloadProgress()
      if status == 'incomplete' :
        print "Job",job.JobID,"currently downloading, status==incomplete"
        return
      else :
        print "Job",job.JobID,"done downloading!"
        del self.currentJobs[job.JobID]
        self.nDQ2StreamsInProgress -= 1
        if job.isDownloadFailed == True :
          self.failedJobs.append(job.JobID)
        else :
          self.successfulJobs.append(job.JobID)
    elif self.nDQ2StreamsInProgress >= self.maxDQ2Streams :
      print "Holding job",job.JobID,"until dq2 stream available."
      return
    else :
      print "Starting download of job",job.JobID
      self.nDQ2StreamsInProgress += 1
      job.retrieveOutput()

  def checkCurrentStatus(self,job) :
    '''I use combinations of panda job statuses to define
    jobs as finished, failed, stuck, or running.'''

    # Panda job status options are:
    # [defined, assigned, activated, running
    #  holding, transferring, finished, failed]
    # https://www.gridpp.ac.uk/wiki/ATLAS_Monitoring_For_Sites

    jobInfo = PdbUtils.readJobDB(job.JobID)
    statusstring = jobInfo.jobStatus
    status = statusstring.split(",")
    print "Status of job", job.JobID, "is", status
    if sorted(status)!=sorted(job.status) :
      job.status = status
      job.statusSince = time.time() # seconds
    timeInThisState = time.time() - job.statusSince

    # If anything has failed, need to retry
#    if status.count('failed') :
#      return 'failed'

    # If everything finished, so is this job.
#    elif status.count('finished')==len(status) \
#           and jobInfo.dbStatus=='frozen':
#      return 'finished'

    if jobInfo.dbStatus=='frozen':
      if status.count('failed') :
        return 'failed'
      elif status.count('finished')==len(status) :
        return 'finished'
      else :
        print "Unrecognized option!"
        return 'none'

    # Don't currently count holding at the end as stuck
    elif status.count('defined') or status.count('assigned') \
           or status.count('activated') or status.count('transferring') \
           or status.count('starting') :
      if timeInThisState < self.defineJobAsStuck :
        return 'stillRunning'
      else :
        return 'stuck'
  
    # If not stuck or failed but something is still running, keep waiting
    elif status.count('running') or status.count('holding') \
           or status.count('sent') :
      return 'stillRunning'

    else :
      print "Unrecognized status!"
      return 'none'

