#!/usr/bin/env python

import os
import ROOT
from art.morisot import Morisot
from array import array
import sys
import numpy as np
import math

class searchFileData :

  def __init__(self,filename,permitWindow=False) :

    self.permitWindow = permitWindow

    searchInputFile = ROOT.TFile.Open(filename,"READ")

    # Retrieve search phase inputs
    self.basicData = searchInputFile.Get("basicData")
    self.basicData.SetDirectory(0)
    self.basicBkgFrom4ParamFit = searchInputFile.Get("basicBkgFrom4ParamFit")
    self.basicBkgFrom4ParamFit.SetDirectory(0)
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

    self.NDF = searchInputFile.Get('NDF')[0]

    excludeWindowNums = searchInputFile.Get('excludeWindowNums')
    self.excludeWindow = int(excludeWindowNums[0]+0.5)
    self.bottomWindowEdge = excludeWindowNums[1]
    self.topWindowEdge = excludeWindowNums[2]

    if (self.excludeWindow and self.permitWindow) :
      statsOfRemainingSpectrum = searchInputFile.Get("BHLogLAndChi2OfRemainderAfterWindow")
      self.BHPValRemainder = statsOfRemainingSpectrum[0]
      self.LogLPValRemainder = statsOfRemainingSpectrum[1]
      self.Chi2PValRemainder = self.calculateRemainingChi2() #statsOfRemainingSpectrum[2]

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
    for bin in range(1,self.basicBkgFrom4ParamFit.GetNbinsX()+2) :
      firstBin = bin
      if self.basicBkgFrom4ParamFit.GetBinContent(bin) > 0 :
        break
    lastBin = 0
    for bin in range(self.basicBkgFrom4ParamFit.GetNbinsX()+1,0,-1) :
      lastBin = bin
      if self.basicBkgFrom4ParamFit.GetBinContent(bin) > 0 :
        break
    firstWindowBin = 0
    lastWindowBin = 0
    if self.excludeWindow :
      for bin in range(1,self.basicBkgFrom4ParamFit.GetNbinsX()+2) :
        if math.fabs(self.basicBkgFrom4ParamFit.GetBinLowEdge(bin) - self.bottomWindowEdge) < 0.1 :
          firstWindowBin = bin
        if math.fabs(self.basicBkgFrom4ParamFit.GetBinLowEdge(bin)+self.basicBkgFrom4ParamFit.GetBinWidth(bin) - self.topWindowEdge) < 0.1 :
          lastWindowBin = bin

    answer = 0
    for bin in range(firstBin,lastBin+1) :

      if self.excludeWindow and bin >= firstWindowBin and bin <= lastWindowBin : continue

      d = self.basicData.GetBinContent(bin)
      if (d==0) : continue
      b = self.basicBkgFrom4ParamFit.GetBinContent(bin)
      deltaB = self.basicBkgFrom4ParamFit.GetBinError(bin)

      term = (d - b) / math.sqrt(b+deltaB*deltaB)
      answer = answer + (term*term)

    nRightChi2 = self.chi2PseudoStatHist.Integral(self.chi2PseudoStatHist.FindBin(answer),self.chi2PseudoStatHist.GetNbinsX())
    nTotal = self.chi2PseudoStatHist.Integral()
    return float(nRightChi2)/float(nTotal)


  def makeSearchPhasePlots(self,myPainter,lowX,highX,luminosity,folder,ext,extraLegendLines=[]) :
 
    firstBin = self.basicData.FindBin(lowX)
    lastBin = self.basicData.FindBin(highX)

    if self.excludeWindow and self.permitWindow :
      myPainter.drawDataAndFitOverSignificanceHist(self.basicData,self.basicBkgFrom4ParamFit,self.residualHist,\
         'm_{jj} [TeV]','Events','Significance','{0}/figure1'.format(folder)+ext,\
         luminosity,13,lowX,highX,firstBin,lastBin,True,self.bumpLowEdge,self.bumpHighEdge,doWindowLimits=self.excludeWindow,windowLow=self.bottomWindowEdge,windowHigh=self.topWindowEdge,extraLegendLines=extraLegendLines)
    else :
      myPainter.drawDataAndFitOverSignificanceHist(self.basicData,self.basicBkgFrom4ParamFit,self.residualHist,\
         'm_{jj} [TeV]','Events','Significance','{0}/figure1'.format(folder)+ext,\
         luminosity,13,lowX,highX,firstBin,lastBin,True,self.bumpLowEdge,self.bumpHighEdge,extraLegendLines=extraLegendLines)
    myPainter.drawPseudoExperimentsWithObservedStat(self.logLikelihoodPseudoStatHist,float(self.logLOfFitToData),self.logLPVal,0,luminosity,13,\
       'logL statistic','Pseudo-exeperiments',"{0}/logLStatPlot".format(folder)+ext)
    myPainter.drawPseudoExperimentsWithObservedStat(self.chi2PseudoStatHist,float(self.chi2OfFitToData),self.chi2PVal,0,luminosity,13,\
       "#chi^{2}",'Pseudo-exeperiments',"{0}/chi2StatPlot".format(folder)+ext)
    myPainter.drawPseudoExperimentsWithObservedStat(self.bumpHunterStatHist,float(self.bumpHunterStatFitToData),self.bumpHunterPVal,0,luminosity,13,\
       'BumpHunter','Pseudo-exeperiments',"{0}/bumpHunterStatPlot".format(folder)+ext)
    myPainter.drawBumpHunterTomographyPlot(self.bumpHunterTomographyPlot,"{0}/bumpHunterTomographyPlot".format(folder)+ext)

