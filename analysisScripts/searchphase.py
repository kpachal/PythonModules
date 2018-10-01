#!/usr/bin/env python

import os
import ROOT
from art.morisot import Morisot
from array import array
import sys
import numpy as np
import math

def GetKeyNames( self, dir = "" ):
        self.cd(dir)
        return [key.GetName() for key in ROOT.gDirectory.GetListOfKeys()]
ROOT.TFile.GetKeyNames = GetKeyNames

class searchFileData :

  def __init__(self,filename,permitWindow=False) :

    self.permitWindow = permitWindow

    #print "opening",filename
    searchInputFile = ROOT.TFile.Open(filename,"READ")

    # Check if readable and return None if not.
    bit = searchInputFile.TestBit(ROOT.TFile.kRecovered)
    if bit :
      print "File",file,"not closed, skipping..."
      raise ValueError
    
    keys = searchInputFile.GetKeyNames()

    # Retrieve search phase inputs
    self.basicData = searchInputFile.Get("basicData")
    self.basicData.SetDirectory(0)
    self.basicBkgFromFit = searchInputFile.Get("basicBkgFrom4ParamFit")
    self.basicBkgFromFit.SetDirectory(0)
    self.residualHist = searchInputFile.Get("residualHist")
    self.residualHist.SetDirectory(0)
    self.relativeDiffHist = searchInputFile.Get("relativeDiffHist")
    self.relativeDiffHist.SetDirectory(0)
    self.sigOfDiffHist = searchInputFile.Get("sigOfDiffHist")
    self.sigOfDiffHist.SetDirectory(0)
    self.logLikelihoodPseudoStatHist = searchInputFile.Get("logLikelihoodStatHistNullCase")
    self.logLikelihoodPseudoStatHist.SetDirectory(0)
    self.chi2PseudoStatHist = searchInputFile.Get("chi2StatHistNullCase")
    self.chi2PseudoStatHist.SetDirectory(0)
    self.bumpHunterStatHist = searchInputFile.Get("bumpHunterStatHistNullCase")
    self.bumpHunterStatHist.SetDirectory(0)
    self.bumpHunterTomographyPlot = searchInputFile.Get('bumpHunterTomographyFromPseudoexperiments')

    fitRange = searchInputFile.Get("FitRange")
    self.fitLow = fitRange[0]
    self.fitHigh = fitRange[1]

    bumpHunterStatOfFitToData = searchInputFile.Get('bumpHunterStatOfFitToData')
    bumpHunterStatOfFitToDataInitial = searchInputFile.Get('bumpHunterStatOfFitToDataInitial')
    bumpHunterStatOfFitToDataRefined = searchInputFile.Get('bumpHunterStatOfFitToDataRefined')
    logLOfFitToDataVec = searchInputFile.Get('logLOfFitToData')
    chi2OfFitToDataVec = searchInputFile.Get('chi2OfFitToData')
    statOfFitToData = searchInputFile.Get('bumpHunterPLowHigh')
    self.logLOfFitToData = logLOfFitToDataVec[0]
    self.logLPVal = logLOfFitToDataVec[1]
    self.chi2OfFitToData = chi2OfFitToDataVec[0]
    self.chi2PVal = chi2OfFitToDataVec[1]
    self.bumpHunterStatFitToData = statOfFitToData[0]
    self.bumpHunterPVal = bumpHunterStatOfFitToData[1]
    self.bumpLowEdge = statOfFitToData[1]
    self.bumpHighEdge = statOfFitToData[2]
    self.bumpHunterStatFitToDataInitial = bumpHunterStatOfFitToDataInitial[0]
    self.bumpHunterPValInitial = bumpHunterStatOfFitToDataInitial[1]
    self.bumpHunterStatFitToDataRefined = bumpHunterStatOfFitToDataRefined[0]
    self.bumpHunterPValRefined = bumpHunterStatOfFitToDataRefined[1]

    self.NDF = searchInputFile.Get('NDF')[0]

    excludeWindowNums = searchInputFile.Get('excludeWindowNums')
    self.excludeWindow = int(excludeWindowNums[0]+0.5)
    self.bottomWindowEdge = excludeWindowNums[1]
    self.topWindowEdge = excludeWindowNums[2]

    if (self.excludeWindow and self.permitWindow) :
      statsOfRemainingSpectrum = searchInputFile.Get("BHLogLAndChi2OfRemainderAfterWindow")
      self.BHPValRemainder = statsOfRemainingSpectrum[0]
      self.LogLPValRemainder = statsOfRemainingSpectrum[1]
      self.Chi2PValRemainder = statsOfRemainingSpectrum[2]

    # Stat uncertainty on fit: only available sometimes
    if "nominalBkgFromFit_plus1Sigma" in keys :
      self.nominalPlus1Stat = searchInputFile.Get('nominalBkgFromFit_plus1Sigma')
      self.nominalPlus1Stat.SetDirectory(0)
      self.nominalMinus1Stat = searchInputFile.Get('nominalBkgFromFit_minus1Sigma')
      self.nominalMinus1Stat.SetDirectory(0)
    
      # Make ratio plots for stat uncertainty
      self.plusNomRatio = self.nominalPlus1Stat.Clone()
      self.plusNomRatio.SetName("ratioPlot_fitUncertainty_plus1Sigma")
      self.plusNomRatio.Reset(); self.plusNomRatio.SetDirectory(0)
      self.minusNomRatio = self.nominalMinus1Stat.Clone()
      self.minusNomRatio.SetName("ratioPlot_fitUncertainty_minus1Sigma")
      self.minusNomRatio.Reset(); self.minusNomRatio.SetDirectory(0)
      for bin in range(0,self.basicBkgFromFit.GetNbinsX()+1) :
        if self.basicBkgFromFit.GetBinContent(bin) == 0 :
          self.plusNomRatio.SetBinContent(bin,0)
          self.minusNomRatio.SetBinContent(bin,0)
        else :
          upVal = (self.nominalPlus1Stat.GetBinContent(bin)-self.basicBkgFromFit.GetBinContent(bin))/self.basicBkgFromFit.GetBinContent(bin)
          downVal = (self.nominalMinus1Stat.GetBinContent(bin)-self.basicBkgFromFit.GetBinContent(bin))/self.basicBkgFromFit.GetBinContent(bin)
          self.plusNomRatio.SetBinContent(bin,upVal)
          self.minusNomRatio.SetBinContent(bin,downVal)

    # Alternate function things: only available sometimes
    if "alternateFitOnRealData" in keys :
      self.alternateFit = searchInputFile.Get('alternateFitOnRealData')
      self.alternateFit.SetDirectory(0)
      self.nomFit_symmetricFitChoiceErr = searchInputFile.Get('nomOnDataWithSymmetricRMSScaleFuncChoiceErr')
      self.nomFit_symmetricFitChoiceErr.SetDirectory(0)
      if "nomOnDataWithDirectedRMSScaleFuncChoiceErr" in keys :
        self.asymmFitUncertainty = searchInputFile.Get('nomOnDataWithDirectedRMSScaleFuncChoiceErr')
      else :
        self.asymmFitUncertainty = searchInputFile.Get('nomPlusDirectedRMSScaleFuncChoiceErr')
      self.asymmFitUncertainty.SetDirectory(0)
      if "nomPlusRMSScaleFuncChoiceErr_averageDirection" in keys :
        self.asymmFitUncertainty_averageDir = searchInputFile.Get('nomPlusRMSScaleFuncChoiceErr_averageDirection')
        self.asymmFitUncertainty_averageDir.SetDirectory(0)
      else :
        self.asymmFitUncertainty_averageDir = None

      # Make ratio plot for alternate function.
      self.altFuncRatio = self.alternateFit.Clone()
      self.altFuncRatio.SetName("ratioPlot_functionChoiceUncertainty")
      self.altFuncRatio.Reset(); self.altFuncRatio.SetDirectory(0)
      for bin in range(0,self.basicBkgFromFit.GetNbinsX()+1) :
        if self.basicBkgFromFit.GetBinContent(bin) == 0 :
          self.altFuncRatio.SetBinContent(bin,0)
        else :
          val = (self.asymmFitUncertainty.GetBinContent(bin)-self.basicBkgFromFit.GetBinContent(bin))/self.basicBkgFromFit.GetBinContent(bin)
          self.altFuncRatio.SetBinContent(bin,val)

    else :
      self.alternateFit = None
      self.nomFit_symmetricFitChoiceErr = None
      self.asymmFitUncertainty = None
      self.asymmFitUncertainty_averageDir = None
    
    

    searchInputFile.Close()

  def getPValErrs(self) :

    # (DeltaX/X)^2 = (1/DeltaX)^2 = 1/X: set errors
    nRightBH = self.bumpHunterStatHist.Integral(self.bumpHunterStatHist.FindBin(self.bumpHunterStatFitToData),self.bumpHunterStatHist.GetNbinsX())
    nLeftBH = self.bumpHunterStatHist.Integral() - nRightBH
    if nRightBH > 0 and nLeftBH > 0 : deltaPvalBH = self.bumpHunterPVal * math.sqrt(1/nRightBH + 1/nLeftBH)
    else : deltaPvalBH = 0
    nRightChi2 = self.chi2PseudoStatHist.Integral(self.chi2PseudoStatHist.FindBin(self.chi2OfFitToData),self.chi2PseudoStatHist.GetNbinsX())
    nLeftChi2 = self.chi2PseudoStatHist.Integral() - nRightChi2
    if nRightChi2 > 0 and nLeftChi2 > 0 : deltaPvalChi2 = self.chi2PVal * math.sqrt(1/nRightChi2 + 1/nLeftChi2)
    else : deltaPvalChi2 = 0
    nRightLogL = self.logLikelihoodPseudoStatHist.Integral(self.logLikelihoodPseudoStatHist.FindBin(self.logLOfFitToData),self.logLikelihoodPseudoStatHist.GetNbinsX())
    nLeftLogL = self.logLikelihoodPseudoStatHist.Integral() - nRightLogL
    if nRightLogL > 0 and nLeftLogL > 0 : deltaPvalLogL = self.logLPVal * math.sqrt(1/nRightLogL + 1/nLeftLogL)
    else : deltaPvalLogL = 0

    return deltaPvalBH,deltaPvalChi2,deltaPvalLogL

  def calculateRemainingChi2(self) :

    firstBin = 0
    for bin in range(1,self.basicBkgFromFit.GetNbinsX()+2) :
      firstBin = bin
      if self.basicBkgFromFit.GetBinContent(bin) > 0 :
        break
    lastBin = 0
    for bin in range(self.basicBkgFromFit.GetNbinsX()+1,0,-1) :
      lastBin = bin
      if self.basicBkgFromFit.GetBinContent(bin) > 0 :
        break
    firstWindowBin = 0
    lastWindowBin = 0
    if self.excludeWindow :
      for bin in range(1,self.basicBkgFromFit.GetNbinsX()+2) :
        if math.fabs(self.basicBkgFromFit.GetBinLowEdge(bin) - self.bottomWindowEdge) < 0.1 :
          firstWindowBin = bin
        if math.fabs(self.basicBkgFromFit.GetBinLowEdge(bin)+self.basicBkgFromFit.GetBinWidth(bin) - self.topWindowEdge) < 0.1 :
          lastWindowBin = bin

    answer = 0
    for bin in range(firstBin,lastBin+1) :

      if self.excludeWindow and bin >= firstWindowBin and bin <= lastWindowBin : continue

      d = self.basicData.GetBinContent(bin)
      if (d==0) : continue
      b = self.basicBkgFromFit.GetBinContent(bin)
      deltaB = self.basicBkgFromFit.GetBinError(bin)

      term = (d - b) / math.sqrt(b+deltaB*deltaB)
      answer = answer + (term*term)

    nRightChi2 = self.chi2PseudoStatHist.Integral(self.chi2PseudoStatHist.FindBin(answer),self.chi2PseudoStatHist.GetNbinsX())
    nTotal = self.chi2PseudoStatHist.Integral()
    return float(nRightChi2)/float(nTotal)


  def makeSearchPhasePlots(self,myPainter,lowX,highX,luminosity,folder,ext,extraLegendLines=[],suppressText=False) :
 
    firstBin = self.basicData.FindBin(lowX)
    lastBin = self.basicData.FindBin(highX)

    myPainter.dodrawUsersText = not suppressText

    if self.excludeWindow and self.permitWindow :
      myPainter.drawDataAndFitOverSignificanceHist(self.basicData,self.basicBkgFromFit,self.residualHist,\
         'm_{jj} [GeV]','Events','Significance','{0}/figure1'.format(folder)+ext,\
         luminosity,13,lowX,highX,firstBin,lastBin,True,self.bumpLowEdge,self.bumpHighEdge,doWindowLimits=self.excludeWindow,windowLow=self.bottomWindowEdge,windowHigh=self.topWindowEdge,extraLegendLines=extraLegendLines,writeOnpval=True,pval=self.bumpHunterPVal)
    else :
      myPainter.drawDataAndFitOverSignificanceHist(self.basicData,self.basicBkgFromFit,self.residualHist,\
         'm_{jj} [GeV]','Events','Significance','{0}/figure1'.format(folder)+ext,\
         luminosity,13,lowX,highX,firstBin,lastBin,True,self.bumpLowEdge,self.bumpHighEdge,extraLegendLines=extraLegendLines,writeOnpval=True,pval=self.bumpHunterPVal)
    
    myPainter.drawPseudoExperimentsWithObservedStat(self.logLikelihoodPseudoStatHist,float(self.logLOfFitToData),self.logLPVal,0,luminosity,13,\
       'logL statistic','Pseudo-exeperiments',"{0}/logLStatPlot".format(folder)+ext)
    myPainter.drawPseudoExperimentsWithObservedStat(self.chi2PseudoStatHist,float(self.chi2OfFitToData),self.chi2PVal,0,luminosity,13,\
       "#chi^{2}",'Pseudo-exeperiments',"{0}/chi2StatPlot".format(folder)+ext)
    myPainter.drawPseudoExperimentsWithObservedStat(self.bumpHunterStatHist,float(self.bumpHunterStatFitToData),self.bumpHunterPVal,0,luminosity,13,\
       'BumpHunter','Pseudo-exeperiments',"{0}/bumpHunterStatPlot".format(folder)+ext)
    myPainter.drawBumpHunterTomographyPlot(self.bumpHunterTomographyPlot,"{0}/bumpHunterTomographyPlot".format(folder)+ext)

    # Plots with alternate function
    if self.alternateFit :
    
      # Make lines I can use for function choice uncertainty plot
      # Need this to round out a pair.
      placeHolderNom = self.basicBkgFromFit.Clone()
      placeHolderNom.SetName("placeHolderNom")
    
      # Find y range to use for residuals
      ylow = 0.0
      yhigh = 0.0
      lowPoint = 1e10
      highPoint = -1e10
      for residual in [self.altFuncRatio,self.minusNomRatio,self.plusNomRatio] :
        for bin in range(firstBin,lastBin+1) :
          val = residual.GetBinContent(bin)
          if val < lowPoint :
            lowPoint = val
          if val > highPoint :
            highPoint = val
      if lowPoint < 0 :
        ylow = lowPoint*1.2
        yhigh = highPoint*(1.2)
      else :
        ylow = lowPoint - 0.9*(highPoint - lowPoint)
        yhigh = highPoint + 0.9*(highPoint - lowPoint)
      # symmetrise
      if abs(ylow) < yhigh :
        ylow = -1 * yhigh
      else :
        yhigh = -1 * ylow
    
      myPainter.drawDataWithFitAsHistogramAndResidualPaper(self.basicData,self.basicBkgFromFit,luminosity,13,\
         "m_{jj} [GeV]","Events",["Data","Fit","Statistical uncertainty on fit","Function choice"],\
         "{0}/compareFitQualityAndFitChoice_Asymm_WithRatio".format(folder)+ext,drawError=True,\
         errors = [[self.nominalPlus1Stat,self.nominalMinus1Stat],[placeHolderNom,self.asymmFitUncertainty]],\
         residualList = [self.altFuncRatio,self.minusNomRatio,self.plusNomRatio],residYRange = [ylow,yhigh],\
         binlow = firstBin, binhigh = lastBin, doLogY = True, doLogX = True, drawAsSmoothCurve = True,\
         doRectangular = False, doLegTopRight = False)

