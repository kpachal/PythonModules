# Kate Pachal, March 2013
#
# Pyroot plotting library for importing to individual
# analyses. Will assist in maintaining consistent style
# and comparable plots.

import sys
import ROOT
import AtlasStyle
#import AtlasUtils
import math
import time
from array import array
from colourPalette import ColourPalette
from collections import Iterable

def colorInterpolate(col1, col2,  w = 0.5):
    c1 = ROOT.gROOT.GetColor(col1);
    c2 = ROOT.gROOT.GetColor(col2);
    r = c1.GetRed()  * (1 - w) + c2.GetRed()  * w;
    g = c1.GetGreen()* (1 - w) + c2.GetGreen()* w;
    b = c1.GetBlue() * (1 - w) + c2.GetBlue() * w;
    return ROOT.TColor.GetColor(r, g, b);


class Morisot(object) :

  ## ----------------------------------------------------
  ## Initialisers

  def __init__(self) :

    # Set up style
    AtlasStyle.SetAtlasStyle()
    ROOT.gROOT.ForceStyle()
    global nsigfigs
    nsigfigs = 1 # for setting number of significant figures in lumi

    # Lydia configurable
    global epsorpdf
    epsorpdf = '.eps'
    global doLumiInPb
    doLumiInPb = False
    # Add text that the user would like to the plot
    self.dodrawUsersText = False # Turn to false for official plots
    global saveCfile # Option to save a .C version of all the plots
    saveCfile = False
    global saveRfile # Option to save a .root version of all the plots
    saveRfile = False
    global saveEfile # Option to save a .eps version of all the plots
    saveEfile = False

    self.cutstring = "|y*| < 0.6"

    self.colourpalette = ColourPalette()
#    self.colourpalette.setColourPalette("Teals")
    self.colourpalette.setColourPalette("Tropical")

    self.shortGoodColours = [1001,1002,1003]
    self.defaultGoodColours = [1001,1002,1003,1004,1000]
    self.mediumGoodColours = [ROOT.kCyan+4,ROOT.kCyan+2,ROOT.kCyan,\
                              ROOT.kBlue,ROOT.kBlue+2,\
                              ROOT.kMagenta+2,ROOT.kMagenta,\
                              ROOT.kRed,ROOT.kRed+2,ROOT.kOrange+10,\
                              ROOT.kOrange,ROOT.kYellow]
    self.longGoodColours = [ROOT.kCyan+4,ROOT.kCyan+3,ROOT.kCyan+2,ROOT.kCyan+1,ROOT.kCyan,\
                     ROOT.kBlue,ROOT.kBlue+1,ROOT.kBlue+2,ROOT.kBlue+3,ROOT.kBlue+4,\
                     ROOT.kMagenta+4,ROOT.kMagenta+3,ROOT.kMagenta+2,ROOT.kMagenta+1,ROOT.kMagenta,\
                     ROOT.kRed,ROOT.kRed+1,ROOT.kRed+2,ROOT.kOrange+9,ROOT.kOrange+10,\
                     ROOT.kOrange+7,ROOT.kOrange,ROOT.kYellow]

    self.myLatex = ROOT.TLatex()
    self.myLatex.SetTextColor(ROOT.kBlack)
    self.myLatex.SetNDC()

    self.myLatex2 = ROOT.TLatex()
    self.myLatex2.SetTextColor(ROOT.kBlack)
    self.myLatex2.SetNDC()

    self.whitebox = ROOT.TPaveText()
    self.whitebox.SetFillColor(0)
    self.whitebox.SetFillStyle(1001)
    self.whitebox.SetTextColor(ROOT.kBlack)
    self.whitebox.SetTextFont(42)
    self.whitebox.SetTextAlign(11)
    self.whitebox.SetBorderSize(0)

    self.line = ROOT.TLine()
    self.line2 = ROOT.TLine()

    self.labeltype = 2 # ATLAS internal
    
    self.labelDict = {
      0 : "", # ATLAS public
      1 : "Preliminary",
      2 : "Internal",
      3 : "Simulation Preliminary",
      4 : "Simulation Internal",
      5 : "Simulation",
      6 : "Work in Progress"
    }

#    self.set2DPalette()

    # 1 "Preliminary"
    # 2 "Internal"
    # 3 "Simulation Preliminary"
    # 4 "Simulation Internal"
    # 5 "Simulation"
    # 6 "Work in Progress"

    self.doAxisTeV = False # Use histogram's natural units

  def setColourPalette(self,palette) :
    self.colourpalette.setColourPalette(palette)

  def setLabelType(self,type) :
    self.labeltype = type

  def useDoAxisTeV(self,dotev=True) :
    self.doAxisTeV = dotev

  def setEPS(self,doEPS) :
    if doEPS :
      epsorpdf = ".eps"

  def set2DPalette(self,name="palette", ncontours=4):

    """Set a color palette from a given RGB list
    stops, red, green and blue should all be lists of the same length
    see set_decent_colors for an example"""

    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.14, 0.00]
        # elif name == "whatever":
        # (define more palettes)
    else:
        # default palette, looks cool
        #jamaica
        # stops = [0.00, 0.20, 0.61, 0.84, 1.00]
        # red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        # green = [0.00, 0.81, 1.00, 0.0, 0.00]
        # blue  = [1.00, 0.20, 0.12, 0.00, 0.00]
        stops = [0.00, 0.20, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]

    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)

    npoints = len(s)
    ROOT.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)

    # For older ROOT versions
    #gStyle.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    ROOT.gStyle.SetNumberContours(ncontours)


  ## ----------------------------------------------------
  ## User-accessible functions

  def drawBasicDataPlot(self,dataHist,luminosity,CME,xname,yname,legendlines,name,binlow=-1,binhigh=-1,doLogY=False,doLogX=False,doRectangular=False) :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Eoutputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,doRectangular)
    c.SetLogx(doLogX)
    c.SetLogy(doLogY)

    if (binlow==-1 and binhigh==-1) :
      firstBin,lastBin = self.getAxisRangeFromHist(dataHist)
    else :
      firstBin = binlow
      lastBin = binhigh
    #dataHist.GetXaxis().SetMoreLogLabels(ROOT.kTRUE)
    if dataHist.GetBinLowEdge(firstBin) > 0.001 and dataHist.GetBinLowEdge(firstBin) < 1 :
      dataHist.GetXaxis().SetNoExponent(ROOT.kTRUE)

    self.drawDataHist(dataHist,firstBin,lastBin,xname,yname,False,1,False)

    lumiInFb = round(float(luminosity)/float(1000),nsigfigs)

    legendsize = 0.04*len(legendlines)
    if legendsize > 0 :
      if doLogX :
        leftOfLegend = 0.20
        legendbottom = 0.20
        legendtop = legendbottom + legendsize
        legend = self.makeLegend(leftOfLegend,legendbottom,0.50,legendtop)
      else :
        leftOfLegend = 0.55
        rightOfLegend = 0.95
        legendtop = 0.90
        legendbottom = legendtop - legendsize
        legend = self.makeLegend(leftOfLegend,legendbottom,rightOfLegend,legendtop)
      legend.AddEntry(dataHist,legendlines[0],"PL")
      legend.Draw()

    if legendsize > 0 and doLogX :
      self.drawLumiAndCMEVert(0.2,legendtop+0.03,lumiInFb,CME)
      self.drawATLASLabels(0.42,0.88)
    else :
      self.drawATLASLabels(0.2, 0.2)
      self.drawLumiAndCMEVert(0.22,0.28,lumiInFb,CME)


    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)

  def drawBasicHistogram(self,hist,binlow,binhigh,xname,yname,name="",makeCanvas=True,doLogY=False,doLogX=False,doErrors=False,fillColour = ROOT.kRed,doRectangular = False) :

    if makeCanvas :
      canvasname = name+'_cv'
      outputname = name+epsorpdf
      if saveCfile:
        Coutputname = name+'.C'
      if saveRfile:
        Routputname = name+'.root'
      if saveEfile:
        Eoutputname = name+'.eps'
      c = self.makeCanvas(canvasname,doRectangular)
      c.SetLogx(doLogX)
      c.SetLogy(doLogY)

    if (binlow==-1 and binhigh==-1) :
      firstBin,lastBin = self.getAxisRangeFromHist(hist)
    else :
      firstBin = binlow
      lastBin = binhigh
    if hist.GetBinLowEdge(firstBin) > 0.001 and hist.GetBinLowEdge(firstBin) < 1 :
      hist.GetXaxis().SetNoExponent(ROOT.kTRUE)

    if (binlow==-1 and binhigh==-1) :
      firstBin,lastBin = self.getAxisRangeFromHist(hist)
    else :
      firstBin = binlow
      lastBin = binhigh
    if hist.GetBinLowEdge(firstBin) > 0.001 and hist.GetBinLowEdge(firstBin) < 1 :
      hist.GetXaxis().SetNoExponent(ROOT.kTRUE)

    hist.GetXaxis().SetRange(firstBin,lastBin)

    hist.SetLineColor(ROOT.kBlack)
    hist.SetFillColor(fillColour)
    hist.GetXaxis().SetTitle(xname);
    hist.GetYaxis().SetTitle(yname);

    hist.Draw("HIST")

    if makeCanvas :
      c.RedrawAxis()
      c.Update()
      c.SaveAs(outputname)
      if saveCfile:
        c.SaveSource(Coutputname)
      if saveRfile:
        c.SaveSource(Routputname)
      if saveEfile:
        c.SaveAs(Eoutputname)

  def drawBasicMatrix(self,matrix,xname,yname,name) :
    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,False)
    c.SetLogx(1)
    c.SetLogy(1)

#    matrix.SetLineColor(ROOT.kBlue+2)

    ROOT.gStyle.SetPalette(1)

    matrix.Draw("COL")
    matrix.GetXaxis().SetTitle(xname)
    matrix.GetYaxis().SetTitle(yname)

    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)


  def drawDataWithFitAsFunction(self,dataHist,function,luminosity,CME,xname,yname,legendlines,name,binlow=-1,binhigh=-1,doLogY=False,doLogX=False,doRectangular=False) :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,doRectangular)
    c.SetLogx(doLogX)
    c.SetLogy(doLogY)

    if (binlow==-1 and binhigh==-1) :
      firstBin,lastBin = self.getAxisRangeFromHist(dataHist)
    else :
      firstBin = binlow
      lastBin = binhigh
    if dataHist.GetBinLowEdge(firstBin) > 0.001 and dataHist.GetBinLowEdge(firstBin) < 1 :
      dataHist.GetXaxis().SetNoExponent(ROOT.kTRUE)

    self.drawDataHist(dataHist,firstBin,lastBin,xname,yname)
    function.SetLineColor(colorpalette.fitLineColor)
    function.Draw("SAME")

    lumiInFb = round(float(luminosity)/float(1000),nsigfigs)

    legendsize = 0.04*len(legendlines)
    if legendsize > 0 :
      if doLogX :
        leftOfLegend = 0.20
        legendbottom = 0.20
        legendtop = legendbottom + legendsize
        legend = self.makeLegend(leftOfLegend,legendbottom,0.50,legendtop)
      else :
        leftOfLegend = 0.55
        rightOfLegend = 0.95
        legendtop = 0.90
        legendbottom = legendtop - legendsize
        legend = self.makeLegend(leftOfLegend,legendbottom,rightOfLegend,legendtop)
      legend.AddEntry(dataHist,legendlines[0],"PL")
      legend.AddEntry(fitHist,legendlines[1],"PL")
      legend.Draw()

    if legendsize > 0 and doLogX :
      self.drawLumiAndCMEVert(0.2,legendtop+0.03,lumiInFb,CME)
      self.drawATLASLabels(0.42,0.88)
    else :
      self.drawATLASLabels(0.2, 0.2)
      self.drawLumiAndCMEVert(0.22,0.28,lumiInFb,CME)

    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)

  def drawDataWithFitAsHistogram(self,dataHist,fitHist,luminosity,CME,xname,yname,legendlines,name,drawError=False,errors = [],binlow=-1,binhigh=-1,doLogY=False,doLogX=False,drawAsSmoothCurve=False,doRectangular=False,doLegTopRight=True,doLabels=True,doEndLines=False) :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,doRectangular)
    c.SetLogx(doLogX)
    c.SetLogy(doLogY)

    if (binlow==-1 and binhigh==-1) :
      firstBin,lastBin = self.getAxisRangeFromHist(dataHist)
    else :
      firstBin = binlow
      lastBin = binhigh
    if dataHist.GetBinLowEdge(firstBin) > 0.001 and dataHist.GetBinLowEdge(firstBin) < 1 :
      dataHist.GetXaxis().SetNoExponent(ROOT.kTRUE)

#    self.drawDataHist(dataHist,firstBin,lastBin,xname,yname,False,1)
#    temp = self.drawFitHist(fitHist,firstBin,lastBin,"","",True,False,drawError, errors,drawAsSmoothCurve,-1,1,doEndLines)
    temp = self.drawFitHist(fitHist,firstBin,lastBin,xname,yname,False,False,drawError,errors,drawAsSmoothCurve,-1,1,doEndLines)
    self.drawDataHist(dataHist,firstBin,lastBin,"","",True,1)
    # Lydia adding analysis cuts values to plot
    if self.dodrawUsersText :
      self.drawUsersText(0.605,0.72,self.cutstring,0.05)

    lumiInFb = round(float(luminosity)/float(1000),nsigfigs)

    legendsize = 0.04*len(legendlines)
    if legendsize > 0 :
      if doLogX and (not doLegTopRight) :
        leftOfLegend = 0.20
        legendbottom = 0.20
        legendtop = legendbottom + legendsize
        legend = self.makeLegend(leftOfLegend,legendbottom,0.50,legendtop)
      else :
        leftOfLegend = 0.55
        rightOfLegend = 0.95
        legendtop = 0.90
        legendbottom = legendtop - legendsize
        legend = self.makeLegend(leftOfLegend,legendbottom,rightOfLegend,legendtop)
      legend.AddEntry(dataHist,legendlines[0],"PL")
      if drawError :
        legend.AddEntry(fitHist,legendlines[1],"L")
        if errors != [] :
          for errhist in errors :
            legend.AddEntry(errhist[1],legendlines[2+errors.index(errhist)],"L")
      else :
        legend.AddEntry(fitHist,legendlines[1],"PL")
      legend.Draw()

    if (doLabels) :
      if legendsize > 0 and doLogX and not doLegTopRight:
        self.drawCMEAndLumi(0.51,0.78,CME,lumiInFb,0.04)
        self.drawATLASLabels(0.53,0.84,True)
      else :
        self.drawATLASLabels(0.2, 0.2)
        self.drawLumiAndCMEVert(0.22,0.28,lumiInFb,CME,0.04)

    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)

  def drawDataOnPrediction(self,sigHist,predHist,xname,yname,legendlines,name,binlow=-1,binhigh=-1,ylow=-1,yhigh=-1,doLabels=False,luminosity=-1,CME=-1,doLogX=False,doLogY=False,doRectangular=False) :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,doRectangular)
    c.SetLogx(doLogX)
    c.SetLogy(doLogY)

    if (binlow==-1 and binhigh==-1) :
      firstBin,lastBin = self.getAxisRangeFromHist(predHist)
    else :
      firstBin = binlow
      lastBin = binhigh
    if predHist.GetBinLowEdge(firstBin) > 0.001 and predHist.GetBinLowEdge(firstBin) < 1 :
      predHist.GetXaxis().SetNoExponent(ROOT.kTRUE)

    if (ylow > 0 and yhigh > 0) :

      predHist.GetYaxis().SetRangeUser(ylow,yhigh)

#    self.drawSignificanceHist(predHist,firstBin,lastBin,xname,yname,False,False,False,False,self.colourpalette.statisticalTestFillColour)
    self.drawBasicHistogram(predHist,firstBin,lastBin,xname,yname,"",False,True,False,False,self.colourpalette.statisticalTestFillColour)
    self.drawDataHist(sigHist,firstBin,lastBin,"","",True,1)

    legendsize = 0.04*len(legendlines)
    if legendsize > 0 :
      if doLogX and (not doLegTopRight) :
        leftOfLegend = 0.20
        legendbottom = 0.20
        legendtop = legendbottom + legendsize
        legend = self.makeLegend(leftOfLegend,legendbottom,0.50,legendtop)
      else :
        leftOfLegend = 0.55
        rightOfLegend = 0.95
        legendtop = 0.90
        legendbottom = legendtop - legendsize
        legend = self.makeLegend(leftOfLegend,legendbottom,rightOfLegend,legendtop)
      legend.AddEntry(sigHist,legendlines[0],"PL")
      legend.AddEntry(predHist,legendlines[1],"PL")
      legend.Draw()

    if (doLabels) :
      if legendsize > 0 and doLogX and not doLegTopRight:
        self.drawCMEAndLumi(0.51,0.78,CME,lumiInFb,0.04)
        self.drawATLASLabels(0.53,0.84,True)
      else :
        self.drawATLASLabels(0.2, 0.2)
        self.drawLumiAndCMEVert(0.22,0.28,lumiInFb,CME,0.04)

    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)


  def drawDataWithFitAsHistogramAndResidual(self,dataHist,fitHist,luminosity,CME,xname,yname,legendlines,name,drawError=False,errors = [],residualList = [],binlow=-1,binhigh=-1,doLogY=False,doLogX=False,drawAsSmoothCurve=False,doRectangular=False,doLegTopRight=True,doLabels=True,doEndLines=False,writeOnpval = False, pval = -999, writeOnFit = False, FitMin =-999,FitMax =-999) :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,doRectangular)
    c.SetLogx(doLogX)
    c.SetLogy(doLogY)

    outpad,pad1,pad2 = self.setStandardTwoPads()
    pad1.SetLogx(doLogX)
    pad2.SetLogx(doLogX)

    # Draw data and fit and uncertainty histograms
    pad1.cd()
    if (binlow==-1 and binhigh==-1) :
      firstBin,lastBin = self.getAxisRangeFromHist(dataHist)
    else :
      firstBin = binlow
      lastBin = binhigh
    if dataHist.GetBinLowEdge(firstBin) > 0.001 and dataHist.GetBinLowEdge(firstBin) < 1 :
      dataHist.GetXaxis().SetNoExponent(ROOT.kTRUE)

#    self.drawDataHist(dataHist,firstBin,lastBin,xname,yname,False,1)
#    temp = self.drawFitHist(fitHist,firstBin,lastBin,"","",True,False,drawError, errors,drawAsSmoothCurve,-1,1,doEndLines)
    self.drawFitHist(fitHist,firstBin,lastBin,xname,yname,False,False,drawError,errors,drawAsSmoothCurve,-1,1,doEndLines)
    self.drawDataHist(dataHist,firstBin,lastBin,"","",True,1)
    pad1.Update()

    outpad.cd()
    # Lydia EOYE adding cuts to plots
    # Lydia adding observedStat value to plot
    print "dodrawUsersText = ",self.dodrawUsersText
    if self.dodrawUsersText:
      if writeOnFit:
        if writeOnpval:
          self.drawUsersText(0.17,0.42,"#splitline{#it{p}-value = "+str(round(pval,2))+"}{#splitline{Fit Range: "+str(round(FitMin/1000,1))+" - "+str(round(FitMax/1000,1))+" TeV}{"+"{0}}".format(self.cutstring),0.039)
        else:
          self.drawUsersText(0.17,0.42,"#splitline{Fit Range: "+str(round(FitMin/1000,1))+" - "+str(round(FitMax/1000,1))+" TeV}{"+"{0}}".format(self.cutstring),0.039)
      else:
        if writeOnpval:
          self.drawUsersText(0.17,0.42,"#it{p}-value = "+str(round(pval,2))+"",0.039)
        elif UserScaleText != "":
          self.drawUsersText(0.17,0.42,self.cutstring,0.039)
    lumiInFb = round(float(luminosity)/float(1000),nsigfigs)

    legendsize = 0.048*len(legendlines)
    if legendsize > 0 :
      if doLogX and (not doLegTopRight) :
        # Lydia EOYE
        leftOfLegend = 0.44
        legendtop = 0.9
        legendbottom = legendtop - legendsize
        legend = self.makeLegend(leftOfLegend,legendbottom,0.75,legendtop)
      else :
        leftOfLegend = 0.55
        rightOfLegend = 0.95
        legendtop = 0.90
        legendbottom = legendtop - legendsize
        legend = self.makeLegend(leftOfLegend,legendbottom,rightOfLegend,legendtop)
      legend.AddEntry(dataHist,legendlines[0],"PL")
      if drawError :
        legend.AddEntry(fitHist,legendlines[1],"L")
        if errors != [] :
          for errhist in errors :
            legend.AddEntry(errhist[1],legendlines[2+errors.index(errhist)],"L")
      else :
        legend.AddEntry(fitHist,legendlines[1],"PL")
      #legend.SetTextSize(0.045)
      legend.Draw()

    if (doLabels) :
      if legendsize > 0 and doLogX and not doLegTopRight:
        # Lydia EOYE
        self.drawCMEAndLumi(0.08,0.49,CME,lumiInFb,0.039)
        self.drawATLASLabels(0.47,0.91,False)#,True)
      else :
        self.drawATLASLabels(0.2, 0.2)
        self.drawLumiAndCMEVert(0.22,0.28,lumiInFb,CME,0.039)

    pad1.Update()

    # Draw residual histograms
    # Use bin range within which bkgPlot has entries,
    # plus one empty on either side if available
    pad2.cd()

    goodcolours = self.getGoodColours(len(residualList))

    for index in range(len(residualList)) :
      residual = residualList[index]

      if residual.GetBinLowEdge(firstBin) > 0.001 and residual.GetBinLowEdge(firstBin) < 1 :
        residual.GetXaxis().SetNoExponent(ROOT.kTRUE)
      #self.drawSignificanceHist(residual,firstBin,lastBin,"x","yresidyList[index]",True,True,True,False,goodcolours[index])
      fixYAxis = True
      inLargerPlot = True
#      doLogX = True
#      doErrors = False

      if index == 0:
        fillColour = goodcolours[1]
        lineStyle = 2
      if index == 1 or index == 2:
        fillColour = goodcolours[0]
        lineStyle = 9
      residual.SetLineColor(fillColour)
      residual.SetLineStyle(lineStyle)
      residual.SetLineWidth(2)
      residual.GetXaxis().SetRange(firstBin,lastBin)
      if index ==0:
        lowPoint = residual.GetMaximum()
        highPoint = residual.GetMinimum()
        ylow = 0.0
        yhigh = 0.0
        for bin in range(firstBin,lastBin+1) :
          val = residual.GetBinContent(bin)
          if val < lowPoint :
            lowPoint = val
          if val > highPoint :
            highPoint = val
        if highPoint == 20 :
          highPoint = 7
        if fixYAxis==False :
          if lowPoint < 0 :
            ylow = lowPoint*1.2
            yhigh = max(highPoint*(1.2),0.15)
          else :
            ylow = lowPoint - 0.9*(highPoint - lowPoint)
            yhigh = highPoint + 0.9*(highPoint - lowPoint)
        else :
          if abs(residual.GetBinContent(residual.GetMaximumBin())) < 1.5 :
            ylow = -0.5
            yhigh = 0.5
          else :
            ylow = -3.7
            yhigh = 3.7
        residual.GetYaxis().SetRangeUser(ylow,yhigh)
        if inLargerPlot :
          residual.GetYaxis().SetTickLength(0.055)
        residual.GetXaxis().SetNdivisions(805,ROOT.kTRUE)

      ###residual.GetYaxis().SetTitle(yname)
      ###residual.GetXaxis().SetTitle(xname)

      #residual.GetYaxis().SetTitleSize(0.14)
      #residual.GetYaxis().SetTitleOffset(0.45) # 1.2 = 20% larger
      #residual.GetYaxis().SetLabelSize(0.115)
      # Lydia EOYE
      residual.GetYaxis().SetTitleSize(0.12)
      residual.GetYaxis().SetTitleOffset(0.42) # 1.2 = 20% larger
      residual.GetYaxis().SetLabelSize(0.115)

      residual.GetYaxis().SetTitle("Rel. Uncert.")
      residual.GetXaxis().SetLabelSize(0.15)
      residual.GetXaxis().SetTitleSize(0.17)
      residual.GetXaxis().SetTitleOffset(1.2)
      residual.GetXaxis().SetTitle(xname)
      residual.GetYaxis().SetNdivisions(805)#5,10,0)

      if drawError :
        residual.Draw("E")
      else :
        residual.Draw("L SAME")

      self.fixTheBloodyTickMarks(ROOT.gPad, residual, residual.GetBinLowEdge(firstBin), residual.GetBinLowEdge(lastBin+1),ylow,yhigh)


    pad2.Update()
    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)


  def drawDataWithFitAsHistogramAndResidualPaper(self,dataHist,fitHist,luminosity,CME,xname,yname,legendlines,name,drawError=False,errors = [],residualList = [],residYRange = [],binlow=-1,binhigh=-1,doLogY=False,doLogX=False,drawAsSmoothCurve=False,doRectangular=False,doLegTopRight=True,doLabels=True,doEndLines=False) :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,doRectangular)
    c.SetLogx(doLogX)
    c.SetLogy(doLogY)

    outpad,pad1,pad2 = self.setStandardTwoPads()
    pad1.SetLogx(doLogX)
    pad2.SetLogx(doLogX)

    # Draw data and fit and uncertainty histograms
    pad1.cd()
    if (binlow==-1 and binhigh==-1) :
      firstBin,lastBin = self.getAxisRangeFromHist(dataHist)
    else :
      firstBin = binlow
      lastBin = binhigh
    if dataHist.GetBinLowEdge(firstBin) > 0.001 and dataHist.GetBinLowEdge(firstBin) < 1 :
      dataHist.GetXaxis().SetNoExponent(ROOT.kTRUE)
    
    self.drawDataHist(dataHist,firstBin,lastBin,xname,yname,same=False,nPads=2)
    temp = self.drawFitHist(fitHist,firstBin,lastBin,xname,yname,True,True,drawError,errors,drawAsSmoothCurve,-1,1,doEndLines)
    pad1.Update()

    # Draw residual histograms
    # Use bin range within which bkgPlot has entries,
    # plus one empty on either side if available
    pad2.cd()
    goodcolours = self.getGoodColours(len(residualList))
    for index in range(len(residualList)) :
      residual = residualList[index]

      if residual.GetBinLowEdge(firstBin) > 0.001 and residual.GetBinLowEdge(firstBin) < 1 :
        residual.GetXaxis().SetNoExponent(ROOT.kTRUE)
        
      if index == 0:
        fillColour = goodcolours[1]
      if index == 1 or index == 2:
        fillColour = goodcolours[0]
      
      # Add some extra formatting
      residual.GetYaxis().SetTitleSize(0.156)
      residual.GetXaxis().SetTitleSize(0.156)
      
      if index==0 : drawSame = False
      else : drawSame = True
      
      self.drawSignificanceHist(residual,firstBin,lastBin,xname,"Rel. Uncert.",fixYAxis = True,\
              inLargerPlot = True, doLogX = True, doErrors = False, fillColour = fillColour, drawAsLine=True, drawSame = drawSame, yRange = residYRange)

    pad2.Update()
    
    # Go to outer pad for labelling
    outpad.cd()
    
    if self.dodrawUsersText:
      self.drawUsersText(0.55,0.76,self.cutstring,0.04)

    legendsize = 0.04*len(legendlines)
    if legendsize > 0 :
      if doLogX and (not doLegTopRight) :
        leftOfLegend = 0.20
        legendbottom = 0.3
        legendtop = legendbottom + legendsize
        legend = self.makeLegend(leftOfLegend,legendbottom,0.50,legendtop)
      else :
        leftOfLegend = 0.55
        rightOfLegend = 0.95
        legendtop = 0.90
        legendbottom = legendtop - legendsize
        legend = self.makeLegend(leftOfLegend,legendbottom,rightOfLegend,legendtop)
      legend.AddEntry(dataHist,legendlines[0],"PL")
      if drawError :
        legend.AddEntry(fitHist,legendlines[1],"L")
        if errors != [] :
          for errhist in errors :
            legend.AddEntry(errhist[1],legendlines[2+errors.index(errhist)],"L")
      else :
        legend.AddEntry(fitHist,legendlines[1],"PL")
      legend.SetTextSize(0.04)
      legend.Draw()

    lumiInFb = round(float(luminosity)/float(1000),nsigfigs)
    if (doLabels) :
      if legendsize > 0 and doLogX and not doLegTopRight:
        self.drawCMEAndLumi(0.47,0.81,CME,lumiInFb,0.04)
        self.drawATLASLabels(0.48, 0.87, True)
      else :
        self.drawATLASLabels(0.2, 0.2)
        self.drawLumiAndCMEVert(0.22,0.28,lumiInFb,CME,0.04)
    
    outpad.Update()
    
    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)

  def drawSignificanceHistAlone(self,significance,xname,yname,name,doLogX=False,doErrors=False,doRectangular=False,firstBin=None,lastBin=None) :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,doRectangular)
    c.SetLogx(doLogX)
    c.SetLogy(0)

    if not firstBin and not lastBin :
      firstBin,lastBin = self.getAxisRangeFromHist(significance)
    print "Using firstBin,lastBin",firstBin,lastBin
    significance.GetXaxis().SetTitleSize(0.04)
    significance.GetXaxis().SetLabelSize(0.04)
    significance.GetYaxis().SetTitleSize(0.04)
    significance.GetYaxis().SetLabelSize(0.04)
    significance.GetXaxis().SetTitleOffset(1.4)
    significance.GetYaxis().SetTitleOffset(1.4)
    self.drawSignificanceHist(significance,firstBin,lastBin,xname,yname,False,False,doLogX,doErrors)

    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)

  def drawStackedHistograms(self,histograms,names,xname,yname,name,xmin,xmax,ymin,ymax,doRectangular=False) :
    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,doRectangular)
    c.SetLogx(0)
    c.SetLogy(1)

    legend = self.makeLegend(0.5,0.55,0.95,0.90)

    goodcolours = self.getGoodColours(len(histograms))

    stack = ROOT.THStack("stack","stacked histograms")
    for histogram in histograms :
      index = histograms.index(histogram)
      histogram.SetLineColor(goodcolours[index])
      histogram.SetLineWidth(2)
      histogram.SetFillColor(goodcolours[index])
      legend.AddEntry(histogram,names[index],"F")
      histogram.SetTitle("")
      stack.Add(histogram,"hist")

    stack.Draw()
    stack.GetXaxis().SetRangeUser(xmin,xmax)
    stack.SetMaximum(ymax)
    stack.SetMinimum(ymin)

    stack.Draw()
    legend.Draw()

    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)


  def drawManyOverlaidHistograms(self,histograms,names,xname,yname,name,xmin,xmax,ymin,ymax,extraLegendLines = [],doLogX=False,doLogY=True,doErrors=False,doRectangular=False,doLegend=True,doLegendLow=True,doLegendLocation="Left",doLegendOutsidePlot=False,doATLASLabel="Low",pairNeighbouringLines=False,dotLines = [],addHorizontalLines=[]) :
    # Left, Right, or Wide

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    if doLegendOutsidePlot :
      if len(histograms)> 12 :
        c = self.makeCanvas(canvasname,doRectangular,2.0,1.0)
      else :
        c = self.makeCanvas(canvasname,doRectangular,1.5,1.0)
    else :
      c = self.makeCanvas(canvasname,doRectangular)
    c.SetLogx(doLogX)
    c.SetLogy(doLogY)

    if doLegendOutsidePlot :
      outpad = ROOT.TPad("extpad","extpad",0,0,1,1) # For marking outermost dimensions
      if len(histograms) > 12 :
        pad1 = ROOT.TPad("pad1","pad1",0,0,0.5,1) # For main histo
        pad2 = ROOT.TPad("pad2","pad2",0.5,0,1,1) # For signal significance histo
      else :
        pad1 = ROOT.TPad("pad1","pad1",0,0,0.66,1) # For main histo
        pad2 = ROOT.TPad("pad2","pad2",0.66,0,1,1) # For signal significance histo

      # Set up to draw in right orientations
      outpad.SetFillStyle(4000) #transparent
      pad1.SetBorderMode(0)
      pad1.SetLogy(doLogY)
      pad1.SetLogx(doLogX)
      pad2.SetBorderMode(0)
      pad1.Draw()
      pad2.Draw()
      outpad.Draw()

      pad1.cd()

    lowxvals = []
    lowyvals = []
    lownonzeros = []
    highxvals = []
    highyvals = []
    flat_histograms = self.flatten_list(histograms)
    for histogram in flat_histograms : # flattens in case sub-lists
      lowx,highx = self.getAxisRangeFromHist(histogram)
      lowy,lownonzero,highy = self.getYRangeFromHist(histogram)
      lowxvals.append(lowx)
      highxvals.append(highx)
      lownonzeros.append(lownonzero)
      lowyvals.append(lowy)
      highyvals.append(highy)
    lowxvals.sort()
    lowyvals.sort()
    lownonzeros.sort()
    highxvals.sort()
    highyvals.sort()
    if xmin == 'automatic':
      minX = lowxvals[0]
    else :
      minX = xmin
    if xmax == 'automatic':
      maxX = highxvals[-1]
    else :
      maxX = xmax
    if ymin == 'automatic':
      if doLogY :
        minY = lownonzeros[0]/2.0
      else :
        minY = lowyvals[0]
    else :
      minY = ymin
    if ymax == 'automatic':
      if doLogY :
        maxY = highyvals[-1]*100
      else :
        maxY = highyvals[-1]*1.5
    else :
      maxY = ymax

    # Create legend
    legendsize = 0.04*len(histograms)
    if doLegendOutsidePlot :
      pad2.cd()
      leftOfLegend = 0
      legendbottom = 0.2
      legendtop = 0.95
      legend = self.makeLegend(leftOfLegend,legendbottom,1.0,legendtop)
      pad1.cd()
    else :
      if doLegendLow :
        legendbottom = 0.20
        legendtop = legendbottom + legendsize
      else :
        legendtop = 0.90 - (0.05)*len(extraLegendLines)
        if doATLASLabel != "Low" and doATLASLabel != "None":
          legendtop = legendtop - 0.05
        legendbottom = legendtop - legendsize
      if doLegendLocation == "Left" :
        leftOfLegend = 0.20
        rightOfLegend = 0.60
      elif doLegendLocation == "Right" :
        leftOfLegend = 0.54
        rightOfLegend = 0.95
      elif doLegendLocation == "Wide" :
        leftOfLegend = 0.20
        rightOfLegend = 0.95
      legend = self.makeLegend(leftOfLegend,legendbottom,rightOfLegend,legendtop)
      if doLegendLocation == "Wide" and len(histograms) > 5 :
        legend.SetNColumns(2)

    if pairNeighbouringLines :
      goodcolours = self.getGoodColours(len(histograms)/2+1)
    else :
      goodcolours = self.getGoodColours(len(histograms))

    for item in histograms :
      index = histograms.index(item)
      # If this is a list, I want only one legend entry
      # and for each histogram to have the same coding.
      if isinstance(item, Iterable) :
        for histogram in item :
          histogram.SetLineColor(goodcolours[index])
          histogram.SetMarkerColor(goodcolours[index])
          histogram.SetLineStyle(1)
          if dotLines != [] and dotLines[index] == True:
            histogram.SetLineStyle(index+1)
        legend.AddEntry(item[0],names[index],"PL")
      elif pairNeighbouringLines :
        print "Pairing!"
        item.SetLineColor(goodcolours[int(index/2.0)])
        item.SetMarkerColor(goodcolours[int(index/2.0)])
        legend.AddEntry(item,names[index],"PL")
      else :
        item.SetLineColor(goodcolours[index])
        item.SetMarkerColor(goodcolours[index])
        item.SetLineStyle(1)
        if dotLines != [] and dotLines[index] == True:
          item.SetLineStyle(index+1)
        legend.AddEntry(item,names[index],"PL")

    # Drawing happens for each
    for histogram in flat_histograms : # flattens in case sub-lists
      index = flat_histograms.index(histogram)
      histogram.SetLineWidth(2)
      histogram.SetFillStyle(0)
      histogram.SetTitle("")
      histogram.GetXaxis().SetRange(minX,maxX)#+5)
      histogram.GetYaxis().SetRangeUser(minY,maxY)
      histogram.GetYaxis().SetTitleSize(0.06)
      histogram.GetYaxis().SetTitleOffset(1.3) # 1.2
      histogram.GetYaxis().SetLabelSize(0.06)
      
      histogram.GetXaxis().SetTitleSize(0.06)
      histogram.GetXaxis().SetTitleOffset(1.2)
      histogram.GetXaxis().SetLabelSize(0.06)
      histogram.GetXaxis().SetNdivisions(605,ROOT.kTRUE)
      
      # Try to tidy up ugly y axes
      histogram.GetYaxis().SetNdivisions(605,ROOT.kTRUE)
      
      if (index==0) :
        histogram.GetXaxis().SetTitle(xname)
        histogram.GetYaxis().SetTitle(yname)
        if not doErrors :
          histogram.Draw("HIST")
        else :
          histogram.Draw("E")
      else :
        histogram.GetXaxis().SetTitle("")
        histogram.GetYaxis().SetTitle("")
        if not doErrors :
          histogram.Draw("HIST SAME")
        else :
          histogram.Draw("E SAME")
      if doLogX :
        if doLegendOutsidePlot :
          self.fixTheBloodyTickMarks(pad1, histogram,histogram.GetBinLowEdge(minX),histogram.GetBinLowEdge(maxX+6),minY,maxY)
        else :
          self.fixTheBloodyTickMarks(ROOT.gPad,histogram, histogram.GetBinLowEdge(minX),histogram.GetBinLowEdge(maxX+6),minY,maxY)

    if (doLegend) :

      if doLegendOutsidePlot :
        index = 0
        if len(extraLegendLines) > 0 :
          for item in range(len(extraLegendLines)) :
            toplocation = legendtop +0.02 + (0.01+0.05)*(index)
            item = self.myLatex.DrawLatex(leftOfLegend+0.02,toplocation,extraLegendLines[index])
            index = index+1
        pad2.cd()
        if len(histograms) > 15 :
          legend.SetNColumns(2)
        legend.Draw()
        pad1.cd()

      else :
        index = 0
        if len(extraLegendLines) > 0 :
          for line in extraLegendLines :
            toplocation = legendtop +0.02 + (0.01+0.04)*(index)
            if doLegendLocation == "Left" :
              self.myLatex.SetTextAlign(11)
              xLocation = leftOfLegend+0.01
            else :
              self.myLatex.SetTextAlign(31)
              xLocation = 0.90
            item = self.myLatex.DrawLatex(xLocation,toplocation,line)
            index = index+1
        legend.Draw()

    # Fix latex
    self.myLatex.SetTextAlign(11)
    if doATLASLabel == "Low" :
      if doLegendLow :
        self.drawATLASLabels(0.22,0.88,False,doRectangular)
      else :
        self.drawATLASLabels(0.2, 0.2,False,doRectangular)
    elif doATLASLabel == "None" :
      pass
    else :
      if doLegendLocation == "Left" :
        self.drawATLASLabels(0.2,0.88)
      else :
        self.drawATLASLabels(0.90,0.88,True)

    if addHorizontalLines != [] :
      for val in addHorizontalLines :
        line = ROOT.TLine(histograms[0].GetBinLowEdge(minX), val, histograms[0].GetBinLowEdge(maxX+6), val)
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(2)
        line.Draw("SAME")

    if doLegendOutsidePlot :
      pad1.RedrawAxis()
    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)

  def drawHistsAndTF1Fits(self,histograms,fits,names,xname,yname,name,xmin,xmax,ymin,ymax,doLogX=True, doRectangular=False,doLogY=False,doErrors=True) :
  
    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,doRectangular)
    c.SetLogx(doLogX)
    c.SetLogy(doLogY)

    lowxvals = []
    lowyvals = []
    lownonzeros = []
    highxvals = []
    highyvals = []
    for histogram in histograms :
      lowx,highx = self.getAxisRangeFromHist(histogram)
      lowy,lownonzero,highy = self.getYRangeFromHist(histogram)
      lowxvals.append(lowx)
      highxvals.append(highx)
      lownonzeros.append(lownonzero)
      lowyvals.append(lowy)
      highyvals.append(highy)
    lowxvals.sort()
    lowyvals.sort()
    lownonzeros.sort()
    highxvals.sort()
    highyvals.sort()
    if xmin == 'automatic':
      minX = lowxvals[0]
    else :
      minX = xmin
    if xmax == 'automatic':
      maxX = highxvals[-1]
    else :
      maxX = xmax
    if ymin == 'automatic':
      if doLogY :
        minY = lownonzeros[0]/2.0
      else :
        minY = lowyvals[0]
    else :
      minY = ymin
    if ymax == 'automatic':
      if doLogY :
        maxY = highyvals[-1]*100
      else :
        maxY = highyvals[-1]*1.5
    else :
      maxY = ymax

    # Create legend
    legendsize = 0.04*len(histograms)
    legendtop = 0.90
    legendbottom = legendtop - legendsize
    leftOfLegend = 0.20
    rightOfLegend = 0.60
    legend = self.makeLegend(leftOfLegend,legendbottom,rightOfLegend,legendtop)

    goodcolours = self.getGoodColours(len(histograms))

    # Plot histograms and fits
    for histogram,fit in zip(histograms,fits) :
      index = histograms.index(histogram)
      histogram.SetLineColor(goodcolours[index])
      histogram.SetMarkerColor(goodcolours[index])
      histogram.SetLineStyle(1)

      histogram.SetLineWidth(2)
      histogram.SetFillStyle(0)
      histogram.SetTitle("")
      histogram.GetXaxis().SetRange(minX,maxX+5)
      histogram.GetYaxis().SetRangeUser(minY,maxY)
      histogram.GetYaxis().SetTitleSize(0.06)
      histogram.GetYaxis().SetTitleOffset(1.2)
      histogram.GetYaxis().SetLabelSize(0.06)
      histogram.GetXaxis().SetTitleSize(0.06)
      histogram.GetXaxis().SetTitleOffset(1.2)
      histogram.GetXaxis().SetLabelSize(0.06)
      histogram.GetXaxis().SetNdivisions(605,ROOT.kTRUE)
      legend.AddEntry(histogram,names[index],"PL")
      if (index==0) :
        histogram.GetXaxis().SetTitle(xname)
        histogram.GetYaxis().SetTitle(yname)
        if not doErrors :
          histogram.Draw("HIST")
        else :
          histogram.Draw("E")
      else :
        histogram.GetXaxis().SetTitle("")
        histogram.GetYaxis().SetTitle("")
        if not doErrors :
          histogram.Draw("HIST SAME")
        else :
          histogram.Draw("E SAME")
          
      fit.SetLineColor(goodcolours[index])
      fit.SetLineStyle(1)
      fit.SetLineWidth(2)
      fit.SetFillStyle(0)
      fit.SetTitle("")
      fit.Draw("C SAME")
      
      if doLogX :
        self.fixTheBloodyTickMarks(ROOT.gPad,histogram, histogram.GetBinLowEdge(minX),histogram.GetBinLowEdge(maxX+6),minY,maxY)

    # Draw legend & ATLAS label
    legend.Draw()
    self.drawATLASLabels(0.2, 0.2,False,doRectangular)

    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)


  def drawTwoHistsDifferentYAxes(self,hist1,hist2,xname,yname1,yname2,name,doRectangular=False) :
    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,doRectangular)
    c.SetLogx(0)
    c.SetLogy(1)

    firstBin1,lastBin1 = self.getAxisRangeFromHist(hist1)
    firstBin2,lastBin2 = self.getAxisRangeFromHist(hist2)
    firstBin = min(firstBin1,firstBin2)
    lastBin = max(lastBin1,lastBin2)

    hists = [hist1,hist2]

    legend = self.makeLegend(0.5,0.7,0.95,0.90)
    goodcolours = self.getGoodColours(2)

    topval = 0
    for histogram in hists :
      index = hists.index(histogram)
      histogram.SetLineColor(goodcolours[index])
      histogram.SetLineWidth(2)
      histogram.SetFillStyle(0)
      histogram.SetTitle("")
      legend.AddEntry(histogram,histogram.GetName(),"PL")
      if (index==0) :
        firstbin,lastbin = self.getAxisRangeFromHist(histogram)
        firsttouse = firstbin-5
        if firsttouse<1 :
          firsttouse = 1
        lasttouse = lastbin+5
        if lasttouse > histogram.GetNbinsX() :
          lasttouse = histogram.GetNbinsX()
        histogram.GetXaxis().SetTitle(xname)
        histogram.GetXaxis().SetRange(firsttouse,lasttouse)
        histogram.GetYaxis().SetTitle(yname1)
        histogram.Draw("HIST")
        firstbin,lastbin
        topval = histogram.GetMaximum()

      else :
        scale = topval/histogram.GetMaximum()
        histogram.Scale(scale)
        histogram.Draw("HIST SAME")

      legend.Draw()

    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)


  def drawPseudoExperimentsWithObservedStat(self,pseudoStatHist,observedStat,pval,pvalerr,luminosity,CME,xname,yname,name,doRectangular=False) :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,doRectangular)
    c.SetLogx(0)
    c.SetLogy(1)

    # Draw stats from pseudoexperiments
    pseudoStatHist.SetLineColor(ROOT.kBlack)
    pseudoStatHist.SetFillColor(self.colourpalette.statisticalTestFillColour) # ROOT.kYellow
    pseudoStatHist.SetFillStyle(1001)
    pseudoStatHist.Draw("HIST")
    pseudoStatHist.GetYaxis().SetRangeUser(0.5,(pseudoStatHist.GetBinContent(pseudoStatHist.GetMaximumBin()))*50)
    pseudoStatHist.GetYaxis().SetTitleOffset(1.4)
    pseudoStatHist.GetXaxis().SetTitle(xname)
    pseudoStatHist.GetYaxis().SetTitle(yname)

    # Draw arrow to observed stat
    arrow = ROOT.TArrow()
    arrow.SetLineColor(self.colourpalette.statisticalTestArrowColour)
    arrow.SetFillColor(self.colourpalette.statisticalTestArrowColour)
    arrow.SetLineWidth(2)
    arrow.DrawArrow(observedStat,1,observedStat,0)
    #print "Drew arrow with colour",self.colourpalette.statisticalTestArrowColour

    # Create legend
    legend = self.makeLegend(0.21,0.68,0.75,0.78)
    lumInFb = round(float(luminosity)/float(1000),nsigfigs)
    legend.AddEntry(pseudoStatHist,"Pseudo-experiments","LF")
    legend.AddEntry(arrow,"Value in Data","L")
    legend.Draw()

    self.drawATLASLabels(0.21,0.88)
    self.drawCMEAndLumi(0.14,0.82,CME,lumInFb,0.04)

    # Lydia adding observedStat value to plot
    if self.dodrawUsersText:
      self.drawUsersText(0.22,0.62,"#it{p}-value = "+str(round(pval,2)),0.04)

    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)

  def drawBumpHunterTomographyPlot(self,tomographyGraph,name) :
    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,False)
    c.SetLogx(1)
    c.SetLogy(1)
    c.SetGridx(0)
    c.SetGridy(0)

    tomographyGraph.SetLineColor(self.colourpalette.tomographyGraphColour)
    tomographyGraph.SetTitle("");
    tomographyGraph.GetXaxis().SetTitle("Dijet Mass [GeV]")
    tomographyGraph.GetXaxis().SetNdivisions(805,ROOT.kTRUE)
    tomographyGraph.GetXaxis().SetMoreLogLabels(ROOT.kTRUE)
    tomographyGraph.GetYaxis().SetTitle("Poisson PVal of Interval")
    tomographyGraph.GetYaxis().SetTitleOffset(1.2);
    tomographyGraph.GetYaxis().SetMoreLogLabels(ROOT.kTRUE)

    tomographyGraph.SetMarkerColor(self.colourpalette.tomographyGraphColour)
    tomographyGraph.SetMarkerSize(0.2)

    tomographyGraph.Draw("AP");


    self.drawATLASLabels(0.55, 0.20, True)

    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)

  def drawDataAndFitOverSignificanceHist(self,dataHist,fitHist,significance,x,datay,sigy,name,luminosity,CME,FitMin,FitMax,firstBin=-1,lastBin=-1,doBumpLimits=False,bumpLow=0,bumpHigh=0,extraLegendLines=[],doLogX=True,doRectangular=False,setYRange=[],writeOnpval = False, pval = -999,doWindowLimits=False,windowLow=0,windowHigh=0) :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,False)
    c.SetLogx(1)
    c.SetLogy(doLogX)
    c.SetGridx(0)
    c.SetGridy(0)

    outpad,pad1,pad2 = self.setStandardTwoPads()
    pad1.SetLogy(1)
    pad1.SetLogx(doLogX)
    pad2.SetLogx(doLogX)

    # Draw data and fit histograms
    pad1.cd()

    # Use bin range within which bkgPlot has entries,
    # plus one empty on either side if available
    lowbin,highbin = self.getAxisRangeFromHist(dataHist)
    if (firstBin>0) :
      lowbin=firstBin
    if (lastBin>0 and lastBin>=firstBin) :
      highbin = lastBin

    self.drawDataHist(dataHist,lowbin,highbin,x,datay,same=False,nPads=2)
    self.drawFitHist(fitHist,lowbin,highbin,"","",True,True,False,[],False,ROOT.kRed,1)
    firstBinWithData,lastBinWithData = self.getAxisRangeFromHist(dataHist)

    # Draw significance histogram
    pad2.cd()
    significance.GetYaxis().SetTitleSize(0.1)
    significance.GetYaxis().SetTitleOffset(0.42) # 1.2 = 20% larger
    significance.GetYaxis().SetLabelSize(0.1)
    significance.GetXaxis().SetLabelSize(0.1)
    significance.GetXaxis().SetTitleSize(0.1)
    significance.GetXaxis().SetTitleOffset(1.2)
    self.drawSignificanceHist(significance,firstBin,lastBin,x,sigy,True)
    c.Update()

    # in place of ROOT.TLine()
    line1 = self.line.Clone("line1"); line1lims = []
    line2 = self.line.Clone("line2"); line2lims = []
    line3 = self.line.Clone("line3"); line3lims = []
    line4 = self.line.Clone("line4"); line4lims = []
    line5 = self.line.Clone("line5"); line1lims = []
    line6 = self.line.Clone("line6"); line2lims = []
    line7 = self.line.Clone("line7"); line3lims = []
    line8 = self.line.Clone("line8"); line4lims = []

    if doBumpLimits :
      heightLowEdge=0
      heightHighEdge=0
      minYvalue = dataHist.GetMinimum()
      for i in range(dataHist.GetNbinsX()) :
        locationOfTallEdge = dataHist.GetBinLowEdge(i)
        height = dataHist.GetBinContent(i)
        if locationOfTallEdge == bumpLow :
          heightLowEdge = dataHist.GetBinContent(i)
        if locationOfTallEdge == bumpHigh:
          heightHighEdge = dataHist.GetBinContent(i-1)

      lowYVal = significance.GetMinimum()#-0.2
      highYVal = significance.GetMaximum()#+0.2

      line1lims = [bumpLow,minYvalue,bumpLow,heightLowEdge]
      line2lims = [bumpHigh,minYvalue,bumpHigh,heightHighEdge]
      line3lims = [bumpLow,lowYVal,bumpLow,highYVal]
      line4lims = [bumpHigh,lowYVal,bumpHigh,highYVal]

      # Draw blue lines
      pad1.cd()
      line1.SetLineColor(ROOT.kBlue)
      line1.SetX1(line1lims[0]); line1.SetY1(line1lims[1]); line1.SetX2(line1lims[2]); line1.SetY2(line1lims[3])
      line1.Draw()
      line2.SetLineColor(ROOT.kBlue)
      line2.SetX1(line2lims[0]); line2.SetY1(line2lims[1]); line2.SetX2(line2lims[2]); line2.SetY2(line2lims[3])
      line2.Draw()
      pad2.cd()
      line3.SetLineColor(ROOT.kBlue)
      line3.SetX1(line3lims[0]); line3.SetY1(line3lims[1]); line3.SetX2(line3lims[2]); line3.SetY2(line3lims[3])
      line3.Draw()
      line4.SetLineColor(ROOT.kBlue)
      line4.SetX1(line4lims[0]); line4.SetY1(line4lims[1]); line4.SetX2(line4lims[2]); line4.SetY2(line4lims[3])
      line4.Draw()

    if doWindowLimits :
      windowHLowEdge=0
      windowHHighEdge=0
      for i in range(dataHist.GetNbinsX()) :
        locationOfTallEdge = dataHist.GetBinLowEdge(i)
        height = dataHist.GetBinContent(i)
        if locationOfTallEdge == windowLow :
          heightLowEdge = dataHist.GetBinContent(i)
        if locationOfTallEdge == windowHigh:
          heightHighEdge = dataHist.GetBinContent(i-1)

      windowlowYVal = significance.GetMinimum()#-0.2
      windowhighYVal = significance.GetMaximum()#+0.2

      line5lims = [windowLow,minYvalue,windowLow,heightLowEdge]
      line6lims = [windowHigh,minYvalue,windowHigh,heightHighEdge]
      line7lims = [windowLow,lowYVal,windowLow,highYVal]
      line8lims = [windowHigh,lowYVal,windowHigh,highYVal]

      # Draw green dashed lines
      pad1.cd()
      line5.SetLineColor(ROOT.kTeal-1)
      line5.SetLineStyle(7)
      line5.SetX1(line5lims[0]); line5.SetY1(line5lims[1]); line5.SetX2(line5lims[2]); line5.SetY2(line5lims[3])
      line5.Draw()
      line6.SetLineColor(ROOT.kTeal-1)
      line6.SetLineStyle(7)
      line6.SetX1(line6lims[0]); line6.SetY1(line6lims[1]); line6.SetX2(line6lims[2]); line6.SetY2(line6lims[3])
      line6.Draw()
      pad2.cd()
      line7.SetLineColor(ROOT.kTeal-1)
      line7.SetLineStyle(7)
      line7.SetX1(line7lims[0]); line7.SetY1(line7lims[1]); line7.SetX2(line7lims[2]); line7.SetY2(line7lims[3])
      line7.Draw()
      line8.SetLineColor(ROOT.kTeal-1)
      line8.SetLineStyle(7)
      line8.SetX1(line8lims[0]); line8.SetY1(line8lims[1]); line8.SetX2(line8lims[2]); line8.SetY2(line8lims[3])
      line8.Draw()

    c.Update()

    outpad.cd()
    leftOfLegend = 0.48
    widthOfRow = 0.04
    lumInFb = round(float(luminosity)/float(1000),nsigfigs)
    if (doLogX and not doBumpLimits) :
      self.drawATLASLabels(0.2, 0.35)
      self.drawCMEAndLumi(0.51,0.90,CME,lumInFb,0.04)
      bottomOfLegend = 0.78
      legend = self.makeLegend(leftOfLegend,bottomOfLegend,0.9,0.87)
    else :
      self.drawATLASLabels(0.45, 0.87, True)
      self.drawCMEAndLumi(0.41,0.82,CME,lumInFb,0.04)
      bottomOfLegend = 0.80 - widthOfRow*(2.0+float(doWindowLimits)+float(doBumpLimits))
      legend = self.makeLegend(leftOfLegend,bottomOfLegend,0.9,0.80,0.04)

    c.Update()

    self.myLatex.SetTextFont(42)
    self.myLatex.SetTextSize(0.04)
    index = 0
    persistent = []
    toplocation = bottomOfLegend
    if len(extraLegendLines) > 0 :
      for line in extraLegendLines :
        toplocation = bottomOfLegend - (widthOfRow)*(index+1) #topOfAll - (0.03+2*widthOfRow) - (0.01+widthOfRow)*(index)
        persistent.append(self.myLatex.DrawLatex(leftOfLegend+0.01,toplocation,line))
        index = index+1

    # Go to outer pad to fill and draw legend
    # Create legend
    outpad.cd()
    legend.AddEntry(dataHist,"Data","LFP")
    legend.AddEntry(fitHist,"Background fit","LF")
    if doBumpLimits :
      legend.AddEntry(line4,"BumpHunter interval","L")
    if doWindowLimits :
      legend.AddEntry(line8,"Excluded window","L")
    legend.Draw()
    c.Update()

    # Lydia adding observedStat value to plot
    if self.dodrawUsersText:
      if doBumpLimits:
        if not writeOnpval:
          self.drawUsersText(0.5,toplocation - 0.06 - len(extraLegendLines)*(widthOfRow+0.01),"#splitline{Fit Range: "+str(FitMin)+" - "+str(FitMax)+" GeV}{"+"{0}}}".format(self.cutstring),0.033)
        else:
          textline = "#splitline{#it{p}-value = "+str(round(pval,2))+"}{#splitline{Fit Range: "+str(FitMin)+" - "+str(FitMax)+" TeV}{"+self.cutstring+"}}"
          self.drawUsersText(0.21,0.42,textline,0.04)
      else:
        if not writeOnpval:
          self.drawUsersText(0.5,toplocation - 0.06 - len(extraLegendLines)*(widthOfRow+0.01),"#splitline{Fit Range: "+str(FitMin)+" - "+str(FitMax)+" TeV}{"+"{0}}}".format(self.cutstring),0.033)
        else:
          self.drawUsersText(0.56,0.7,"#splitline{#it{p}-value = "+str(round(pval,2))+"}{#splitline{Fit Range: "+str(FitMin)+" - "+str(FitMax)+" TeV}{"+"{0}}}".format(self.cutstring),0.04)

    # Save.
    pad1.RedrawAxis()
    pad2.RedrawAxis()
    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)

  def drawDataAndFitsOverSignificanceHists_TwoSpectra(self,dataHist,fitHist,significance,dataHist2,fitHist2,significance2,signal1,signal2,scale1,scale2,x,datay,sigy,name,lumi1,lumi2,datastring1,datastring2,CME,FitMin,FitMax,firstBin=-1,lastBin=-1,doBumpLimits=False,bumpLow=0,bumpHigh=0,bumpLow2=0,bumpHigh2=0,extraLegendLines=[],doLogX=True,doRectangular=False,setYRange=[],writeOnpval = True, pval = -999, chi2pval = -999,pval2 = -999, chi2pval2 = -999,doWindowLimits=False,windowLow=0,windowHigh=0,dataPointsOption=0,fancinessOption=0) :
  
    # Options are for alternate versions of the plot.
    # dataPointsOption:
    # 0 = solid squares for longer data spectrum, open squares for its 1 signal,
    #     solid circles for shorter spectrum, open circles for its 1 signal
    # 1 = same shapes as 0 but points are colour coded to match their spectra
    # 2 = all solid shapes for longer data spectrum, all open shapes for shorter one
    
    # fancinessOption:
    # 0 = default: data, fit, one signal per spectrum
    # 1 = least fancy: only data
    # 2 = data + fit, no signals
    # 3 = default with an extra signal per mass

    # Make canvas
    canvasname = name+'_cv'
    outputname = name+epsorpdf
    EPSoutputname = name+'.eps'
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.pdf'
      EPSoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,True) #Hanno: Was False
    c.SetLogx(1)
    c.SetLogy(doLogX)
    c.SetGridx(0)
    c.SetGridy(0)

    # Make pads
    # Dimensions: xlow, ylow, xup, yup
    outpad = ROOT.TPad("extpad","extpad",0,0,1,1) # For marking outermost dimensions
    pad1 = ROOT.TPad("pad1","pad1",0,0.4,1,1) # For main histo
    pad2 = ROOT.TPad("pad2","pad2",0,0.25,1,0.4) # For residuals histo
    pad3 = ROOT.TPad("pad2","pad2",0,0,1,0.25) # For residuals histo    

    # Set up pads to draw with right spacings
    outpad.SetFillStyle(4000) #transparent
    pad1.SetBorderMode(0)
    pad1.SetLogy(1)
    pad1.SetLogx(doLogX)
    pad1.SetLeftMargin(0.1)
    pad1.SetBottomMargin(0.00001)
    pad2.SetLeftMargin(0.1)
    pad2.SetTopMargin(0.00001)
    pad2.SetBottomMargin(0.00001)
    pad2.SetBorderMode(0)
    pad2.SetLogx(doLogX)    
    pad3.SetLeftMargin(0.1)
    pad3.SetTopMargin(0.00001)
    pad3.SetBottomMargin(0.4)
    pad3.SetBorderMode(0)
    pad3.SetLogx(doLogX)
    pad1.Draw()
    pad2.Draw()
    pad3.Draw()
    outpad.Draw()
    
    # Define a few common items
    leftOfLegend1 = 0.13
    widthOfRow = 0.04
    lumInFb1 = round(float(lumi1)/float(1000),nsigfigs)
    lumInFb2 = round(float(lumi2)/float(1000),nsigfigs)
    atlasLabelXLocation = 0.137
    atlasLabelYLocation = 0.92

    # Draw data and fit histograms
    if fancinessOption !=1 :
      pad1.cd()
    else :
      outpad.SetBorderMode(0)
      outpad.SetLogy(1)
      outpad.SetLogx(doLogX)
      outpad.SetLeftMargin(0.1)
      outpad.cd()

    # Get bin range for first spectrum
    lowbin,highbin = self.getAxisRangeFromHist(dataHist)
    if (firstBin>0) :
      lowbin=firstBin
    if (lastBin>0 and lastBin>=firstBin) :
      highbin = lastBin

    # For second spectrum just want to start drawing where we fitted it.
    lowbin2,highbin2 = self.getAxisRangeFromHist(fitHist2)
    # Gives one bin lower by default and we don't want that.
    lowbin2 = lowbin2+1
    
    # Use colours from colour scheme. Need 4 of them.
    colours = self.getGoodColours(4)

    # Draw both data hists
    self.drawDataHist(dataHist,lowbin,highbin,x,datay,True,3)
    self.drawDataHist(dataHist2,lowbin2,highbin,x,datay,True,3)
    
    # Post hoc data plot formatting.
    dataHist.SetMarkerSize(0.9)
    dataHist2.SetMarkerSize(0.9)
    if dataPointsOption == 0 :
      dataHist.SetMarkerStyle(21)
    elif dataPointsOption == 1 :
      dataHist.SetMarkerStyle(21)
      dataHist.SetMarkerColor(colours[3])
      dataHist2.SetMarkerColor(colours[0])
    elif dataPointsOption == 2 :
      dataHist.SetMarkerStyle(25)
    
    # Add signals
    if fancinessOption < 1 or fancinessOption > 2 :
      self.drawDataHist(signal1,lowbin,highbin,x,datay,True,3)
      self.drawDataHist(signal2,lowbin2,highbin,x,datay,True,3)

    # Update sig formatting
    if dataPointsOption == 2 :
      signal2.SetMarkerStyle(20)
    else :
      signal2.SetMarkerStyle(24)
    signal1.SetMarkerStyle(25)

    signal2.SetMarkerColor(colours[1])
    signal2.SetLineColor(colours[1])
    signal2.SetMarkerSize(0.9)
    signal1.SetMarkerSize(0.9)
    signal1.SetMarkerColor(colours[2])
    signal1.SetLineColor(colours[2])
    
    # Draw fits
    if fancinessOption != 1 :
      self.drawFitHist(fitHist,lowbin,highbin,"","",same=True,twoPads=True,useError=False,errors=[],drawCurve=False,lineColor=colours[3],lineStyle=1)
      self.drawFitHist(fitHist2,lowbin2,highbin,"","",same=True,twoPads=True,useError=False,errors=[],drawCurve=False,lineColor=colours[0],lineStyle=1)

      # Update label formatting
      dataHist.GetYaxis().SetTitleFont(43)
      dataHist.GetYaxis().SetTitleSize(25)
      dataHist.GetYaxis().SetTitleOffset(1)
      dataHist.GetYaxis().SetLabelFont(43)
      dataHist.GetYaxis().SetLabelSize(25)
      
    dataHist.GetYaxis().SetTitle("Events / Bin")
    dataHist.GetXaxis().SetLabelSize(0)
    
    # Update axis range
    actualMin = 1E10
    actualMax = 0
    for bin in range(firstBin,lastBin+1) :
      if dataHist.GetBinContent(bin) > actualMax : actualMax = dataHist.GetBinContent(bin)
      elif dataHist.GetBinContent(bin) < actualMin : actualMin = dataHist.GetBinContent(bin)
    minYaxis = 0.7 if actualMin < 2.5 else float(actualMin)/2.0
    dataHist.GetYaxis().SetRangeUser(minYaxis,10*actualMax)

    dataHist.GetXaxis().SetMoreLogLabels()

    # If we only want data, we're done now.
    if fancinessOption == 1 :

      # Draw x axis labels that aren't a mess
      self.fixTheBloodyLabels(outpad,significance.GetBinLowEdge(lowbin),significance.GetBinLowEdge(highbin+1),fontSize=0.05,overrideY=0.13)
      outpad.cd()
    
      # Do ATLAS and channel labels
      atlasLabelYLocation = 0.88 # Looks squashed with same dimensions as other plot
      self.drawATLASLabels(atlasLabelXLocation, atlasLabelYLocation, isRectangular = True)

      self.myLatex.SetTextFont(42)
      self.myLatex.SetTextSize(0.04)
      index = 0
      persistent = []
      if len(extraLegendLines) > 0 :
        for line in extraLegendLines :
          toplocation = atlasLabelYLocation - (widthOfRow)*(index+1)-0.02
          persistent.append(self.myLatex.DrawLatex(atlasLabelXLocation,toplocation,line))
          index = index+1
    
      # Do legend: only one for this plot
      legend = self.makeLegend(0.14,0.2,0.4,0.35,fontSize = 0.04*0.93)
      if datastring1 :
        thisline = "Data, {0} fb^{1}, {2}".format(lumInFb1,"{-1}",datastring1)
      else :
        thisline = "Data, {0} fb^{1}".format(lumInFb1,"{-1}")
      legend.AddEntry(dataHist,thisline,"LFP")

      if datastring2 :
        thisline2 = "Data, {0} fb^{1}, {2}".format(lumInFb2,"{-1}",datastring2)
      else :
        thisline2 = "Data, {0} fb^{1}".format(lumInFb2,"{-1}")
      legend.AddEntry(dataHist2,thisline2,"LFP")
      legend.Draw()
      c.Update()
    
      c.RedrawAxis()
      c.Update()
      c.SaveAs(outputname)
      if saveCfile:
        c.SaveSource(Coutputname)
      if saveRfile:
        c.SaveSource(Routputname)
      if saveEfile:
        c.SaveAs(Eoutputname)
        c.SaveAs(EPSoutputname)
      
      return
      
    # Otherwise, carry on and go to two-legend structure.

    # Draw significance histogram
    pad2.cd()
    self.drawSignificanceHist(significance2,firstBin,lastBin,"","",fixYAxis=True,inLargerPlot = False, doLogX=False, doErrors=False, fillColour = colours[0])
    
    # Update label formatting
    significance2.GetYaxis().SetTitleFont(43)
    significance2.GetYaxis().SetTitleSize(25)
    significance2.GetYaxis().SetTitleOffset(1) # 1.2 = 20% larger
    significance2.GetYaxis().SetLabelFont(43)
    significance2.GetYaxis().SetLabelSize(25)

    pad3.cd()
    self.drawSignificanceHist(significance,firstBin,lastBin,x,"",fixYAxis=True,inLargerPlot = True, doLogX=False, doErrors=False, fillColour = colours[3])
    
    # Update label formatting
    significance.GetXaxis().SetNoExponent(1)
    significance.GetYaxis().SetTitleFont(43)
    significance.GetYaxis().SetTitleSize(25)
    significance.GetYaxis().SetTitleOffset(1) # 1.2 = 20% larger
    significance.GetYaxis().SetLabelFont(43)
    significance.GetYaxis().SetLabelSize(25)
    significance.GetXaxis().SetLabelFont(43)
    significance.GetXaxis().SetLabelSize(25)
    # This allows custom label formatting.
    significance.GetXaxis().SetLabelSize(0)
    significance.GetXaxis().SetTitleFont(43)
    significance.GetXaxis().SetTitleSize(25)
    significance.GetXaxis().SetTitleOffset(4) # 1.2
    significance.GetXaxis().SetNdivisions(802,ROOT.kTRUE)
    
    # Redraw y axis label for significances
    outpad.cd()
    self.myLatex2.SetTextAngle(90)
    self.myLatex2.SetTextFont(43)
    self.myLatex2.SetTextSize(24)
    self.myLatex2.DrawLatex(0.039, 0.155, "Significance")

    c.Update()

    # Make a bunch of lines
    self.line.SetLineStyle(2)
    line1 = self.line.Clone("line1"); line1lims = []
    line2 = self.line.Clone("line2"); line2lims = []
    line3 = self.line.Clone("line3"); line3lims = []
    line4 = self.line.Clone("line4"); line4lims = []
    line5 = self.line.Clone("line5"); line1lims = []
    line6 = self.line.Clone("line6"); line2lims = []
    line7 = self.line.Clone("line7"); line3lims = []
    line8 = self.line.Clone("line8"); line4lims = []
    line3_pad2 = self.line.Clone("line3_pad2"); line3lims = []
    line4_pad2 = self.line.Clone("line4_pad2"); line4lims = []
    line1_2 = self.line.Clone("line1"); line1lims_2 = []
    line2_2 = self.line.Clone("line2"); line2lims_2 = []
    line3_2 = self.line.Clone("line3"); line3lims_2 = []
    line4_2 = self.line.Clone("line4"); line4lims_2 = []
    line5_2 = self.line.Clone("line5"); line1lims_2 = []
    line6_2 = self.line.Clone("line6"); line2lims_2 = []
    line7_2 = self.line.Clone("line7"); line3lims_2 = []
    line8_2 = self.line.Clone("line8"); line4lims_2 = []
    
    # Draw bump lines
    if doBumpLimits :
      heightLowEdge=0
      heightHighEdge=0
      minYvalue = dataHist.GetMinimum()
      for i in range(dataHist.GetNbinsX()) :
        locationOfTallEdge = dataHist.GetBinLowEdge(i)
        height = dataHist.GetBinContent(i)
        if locationOfTallEdge == bumpLow :
          heightLowEdge = dataHist.GetBinContent(i)
        if locationOfTallEdge == bumpHigh:
          heightHighEdge = dataHist.GetBinContent(i-1)

      lowYVal = significance.GetMinimum()#-0.2
      highYVal = significance.GetMaximum()#+0.2

      line1lims = [bumpLow,0.5*minYvalue,bumpLow,heightLowEdge]
      line2lims = [bumpHigh,0.5*minYvalue,bumpHigh,heightHighEdge]
      line3lims = [bumpLow,lowYVal,bumpLow,highYVal]
      line4lims = [bumpHigh,lowYVal,bumpHigh,highYVal]
      helperlinelims = [850,-4.68,943,-4.68]
      
      heightLowEdge_2=0
      heightHighEdge_2=0
      minYvalue_2 = dataHist2.GetMinimum()
      for i in range(dataHist2.GetNbinsX()) :
        locationOfTallEdge_2 = dataHist2.GetBinLowEdge(i)
        height_2 = dataHist2.GetBinContent(i)
        if locationOfTallEdge_2 == bumpLow2 :
          heightLowEdge_2 = dataHist2.GetBinContent(i)
        if locationOfTallEdge_2 == bumpHigh2:
          heightHighEdge_2 = dataHist2.GetBinContent(i-1)

      lowYVal_2 = significance2.GetMinimum()#-0.2
      highYVal_2 = significance2.GetMaximum()#+0.2
      line1lims_2 = [bumpLow2,0.5*minYvalue_2,bumpLow2,heightLowEdge_2]
      line2lims_2 = [bumpHigh2,0.5*minYvalue_2,bumpHigh2,heightHighEdge_2]
      line3lims_2 = [bumpLow2,lowYVal_2,bumpLow2,highYVal_2]
      line4lims_2 = [bumpHigh2,lowYVal_2,bumpHigh2,highYVal_2]      

      # Draw lines
      pad1.cd()
      line1.SetLineColor(colours[3])
      line1.SetX1(line1lims[0]); line1.SetY1(line1lims[1]); line1.SetX2(line1lims[2]); line1.SetY2(line1lims[3])
      line1.Draw()
      line2.SetLineColor(colours[3])
      line2.SetX1(line2lims[0]); line2.SetY1(line2lims[1]); line2.SetX2(line2lims[2]); line2.SetY2(line2lims[3])
      line2.Draw()
      pad2.cd()
      line3_pad2.SetLineColor(colours[3])
      line3_pad2.SetX1(line3lims[0]); line3_pad2.SetY1(line3lims[1]); line3_pad2.SetX2(line3lims[2]); line3_pad2.SetY2(line3lims[3])
      line3_pad2.Draw()
      line4_pad2.SetLineColor(colours[3])
      line4_pad2.SetX1(line4lims[0]); line4_pad2.SetY1(line4lims[1]); line4_pad2.SetX2(line4lims[2]); line4_pad2.SetY2(line4lims[3])
      line4_pad2.Draw()
      pad3.cd()
      line3.SetLineColor(colours[3])
      line3.SetX1(line3lims[0]); line3.SetY1(line3lims[1]); line3.SetX2(line3lims[2]); line3.SetY2(line3lims[3])
      line3.Draw()
      line4.SetLineColor(colours[3])
      line4.SetX1(line4lims[0]); line4.SetY1(line4lims[1]); line4.SetX2(line4lims[2]); line4.SetY2(line4lims[3])
      line4.Draw()
      
      pad1.cd()
      line1_2.SetLineColor(colours[0])
      line1_2.SetX1(line1lims_2[0]); line1_2.SetY1(line1lims_2[1]); line1_2.SetX2(line1lims_2[2]); line1_2.SetY2(line1lims_2[3])
      line1_2.Draw()
      line2_2.SetLineColor(colours[0])
      line2_2.SetX1(line2lims_2[0]); line2_2.SetY1(line2lims_2[1]); line2_2.SetX2(line2lims_2[2]); line2_2.SetY2(line2lims_2[3])
      line2_2.Draw()
      pad2.cd()
      line3_2.SetLineColor(colours[0])
      line3_2.SetX1(line3lims_2[0]); line3_2.SetY1(line3lims_2[1]); line3_2.SetX2(line3lims_2[2]); line3_2.SetY2(line3lims_2[3])
      line3_2.Draw()
      line4_2.SetLineColor(colours[0])
      line4_2.SetX1(line4lims_2[0]); line4_2.SetY1(line4lims_2[1]); line4_2.SetX2(line4lims_2[2]); line4_2.SetY2(line4lims_2[3])
      line4_2.Draw()

    if doWindowLimits :
      windowHLowEdge=0
      windowHHighEdge=0
      for i in range(dataHist.GetNbinsX()) :
        locationOfTallEdge = dataHist.GetBinLowEdge(i)
        height = dataHist.GetBinContent(i)
        if locationOfTallEdge == windowLow :
          heightLowEdge = dataHist.GetBinContent(i)
        if locationOfTallEdge == windowHigh:
          heightHighEdge = dataHist.GetBinContent(i-1)

      windowlowYVal = significance.GetMinimum()#-0.2
      windowhighYVal = significance.GetMaximum()#+0.2

      line5lims = [windowLow,minYvalue,windowLow,heightLowEdge]
      line6lims = [windowHigh,minYvalue,windowHigh,heightHighEdge]
      line7lims = [windowLow,lowYVal,windowLow,highYVal]
      line8lims = [windowHigh,lowYVal,windowHigh,highYVal]

      # Draw green dashed lines
      pad1.cd()
      line5.SetLineColor(ROOT.kTeal-1)
      line5.SetLineStyle(7)
      line5.SetX1(line5lims[0]); line5.SetY1(line5lims[1]); line5.SetX2(line5lims[2]); line5.SetY2(line5lims[3])
      line5.Draw()
      line6.SetLineColor(ROOT.kTeal-1)
      line6.SetLineStyle(7)
      line6.SetX1(line6lims[0]); line6.SetY1(line6lims[1]); line6.SetX2(line6lims[2]); line6.SetY2(line6lims[3])
      line6.Draw()
      pad2.cd()
      line7.SetLineColor(ROOT.kTeal-1)
      line7.SetLineStyle(7)
      line7.SetX1(line7lims[0]); line7.SetY1(line7lims[1]); line7.SetX2(line7lims[2]); line7.SetY2(line7lims[3])
      line7.Draw()
      line8.SetLineColor(ROOT.kTeal-1)
      line8.SetLineStyle(7)
      line8.SetX1(line8lims[0]); line8.SetY1(line8lims[1]); line8.SetX2(line8lims[2]); line8.SetY2(line8lims[3])
      line8.Draw()

    c.Update()
    outpad.cd()

    # ATLAS labels and lines near them
    self.drawATLASLabels(atlasLabelXLocation, atlasLabelYLocation, isRectangular = True)
    self.myLatex.SetTextFont(42)
    self.myLatex.SetTextSize(0.04)
    index = 0
    persistent = []
    if len(extraLegendLines) > 0 :
      for line in extraLegendLines :
        toplocation = atlasLabelYLocation - (widthOfRow)*(index+1)
        persistent.append(self.myLatex.DrawLatex(atlasLabelXLocation,toplocation,line))
        index = index+1

    # Legends with appropriate dimensions
    bottomOfLegend1 = 0.415
    topOfLegend2 = 0.95
    depthOfLegend = min(0.04*(2+(2 if (fancinessOption!=1 and fancinessOption!=2) else 0)),0.1) # data, fit, any signal
    depthOfLegend = depthOfLegend + 0.07*(float(doWindowLimits)+float(doBumpLimits)) # window and bump limits
    topOfLegend1 = bottomOfLegend1 + depthOfLegend
    bottomOfLegend2 = topOfLegend2 - depthOfLegend
    legend = self.makeLegend(leftOfLegend1,bottomOfLegend1,0.4,topOfLegend1,fontSize = 0.04*0.93)
    legend2 = self.makeLegend(0.56,bottomOfLegend2,0.84,topOfLegend2,fontSize = 0.04*0.93)

    c.Update()

    # Create legend
    if datastring1 :
      thisline = "Data, {0} fb^{1}, {2}".format(lumInFb1,"{-1}",datastring1)
    else :
      thisline = "Data, {0} fb^{1}".format(lumInFb1,"{-1}")
    legend.AddEntry(dataHist,thisline,"LFP")
    legend.AddEntry(fitHist,"Background fit","LF")

    if datastring2 :
      thisline2 = "Data, {0} fb^{1}, {2}".format(lumInFb2,"{-1}",datastring2)
    else :
      thisline2 = "Data, {0} fb^{1}".format(lumInFb2,"{-1}")
    legend2.AddEntry(dataHist2,thisline2,"LFP")
    legend2.AddEntry(fitHist2,"Background fit","LF")

    if doBumpLimits :
      legend.AddEntry(line4,"BumpHunter interval","L")
      legend2.AddEntry(line4_2,"BumpHunter interval","L")
    if doWindowLimits :
      legend.AddEntry(line8,"Excluded window","L")
    if fancinessOption == 0 or fancinessOption == 3 :
      legend.AddEntry(signal1,"Z', #sigma x {0}".format(scale1))
      legend2.AddEntry(signal2,"Z', #sigma x {0}".format(scale2))
      helperHist = dataHist2.Clone()
      helperHist.SetMarkerColor(0)
      helperHist.SetLineColor(0)
      legend.AddEntry(helperHist,"m_{Z'} = 250 GeV, g_{q} = 0.1","LFP")
      legend2.AddEntry(helperHist,"m_{Z'} = 550 GeV, g_{q} = 0.1","LFP")
    legend.Draw()
    legend2.Draw()
    c.Update()

    # Adding p-values to plot:
    # want them adjascent to legends
    if self.dodrawUsersText:
      if doBumpLimits:
        self.drawUsersTextLeftAligned(leftOfLegend1,topOfLegend1+0.02,["","","BH #it{p}-value = "+str(round(pval,2)), "#chi^{2} #it{p}-value = "+str(round(chi2pval,2))],0.033)
        self.drawUsersTextRightAligned(0.68,bottomOfLegend2-0.08,["","","BH #it{p}-value = "+str(round(pval2,2)), "#chi^{2} #it{p}-value = "+str(round(chi2pval2,2))],0.033)

    # Save.
    pad1.RedrawAxis()
    pad2.RedrawAxis()
    pad3.RedrawAxis()

    # Draw x axis labels that aren't a mess
    self.fixTheBloodyLabels(pad2,significance.GetBinLowEdge(lowbin),significance.GetBinLowEdge(highbin+1),fontSize=0.17)
    
    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)
      c.SaveAs(EPSoutputname)


  def compareDataToLimit(self,dataHist,fitHist,significance,observedLimit,x,datay,sigy,name,luminosity,CME,firstBin=-1,lastBin=-1,doBumpLimits=False,bumpLow=0,bumpHigh=0,extraLegendLines=[],doLogX=True,doRectangular=False,setYRange=[],writeOnpval = False, pval = -999) :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,False)
    c.SetLogx(1)
    c.SetLogy(doLogX)
    c.SetGridx(0)
    c.SetGridy(0)

    outpad,pad1,pad2 = self.setStandardTwoPads()
    pad1.SetLogx(doLogX)
    pad2.SetLogx(doLogX)

    # Draw data and fit histograms
    pad1.cd()

    # Use bin range within which bkgPlot has entries,
    # plus one empty on either side if available
    lowbin,highbin = self.getAxisRangeFromHist(dataHist)
    if (firstBin>0) :
      lowbin=firstBin
    if (lastBin>0 and lastBin>=firstBin) :
      highbin = lastBin

    fitHist.GetYaxis().SetTitleSize(0.05)
    fitHist.GetYaxis().SetTitleOffset(1.0)
    fitHist.GetYaxis().SetLabelSize(0.05)
    self.drawDataHist(dataHist,lowbin,highbin,x,datay,False,2)
    self.drawFitHist(fitHist,lowbin,highbin,"","",True,True,False,[],False,ROOT.kRed,1)
    firstBinWithData,lastBinWithData = self.getAxisRangeFromHist(dataHist)

    observedLimit.Draw("PL SAME")

    # Draw significance histogram
    pad2.cd()
    significance.GetYaxis().SetTitleSize(0.1)
    significance.GetYaxis().SetTitleOffset(0.42) # 1.2 = 20% larger
    significance.GetYaxis().SetLabelSize(0.1)
    significance.GetXaxis().SetLabelSize(0.1)
    significance.GetXaxis().SetTitleSize(0.1)
    significance.GetXaxis().SetTitleOffset(1.2)
    self.drawSignificanceHist(significance,firstBin,lastBin,x,sigy,True)
    c.Update()

    # in place of ROOT.TLine()
    line1 = self.line.Clone("line1"); line1lims = []
    line2 = self.line.Clone("line2"); line2lims = []
    line3 = self.line.Clone("line3"); line3lims = []
    line4 = self.line.Clone("line4"); line4lims = []
    if doBumpLimits :
      heightLowEdge=0
      heightHighEdge=0
      minYvalue = dataHist.GetMinimum()
      for i in range(dataHist.GetNbinsX()) :
        locationOfTallEdge = dataHist.GetBinLowEdge(i)
        height = dataHist.GetBinContent(i)
        if locationOfTallEdge == bumpLow :
          heightLowEdge = dataHist.GetBinContent(i)
        if locationOfTallEdge == bumpHigh:
          heightHighEdge = dataHist.GetBinContent(i-1)

      lowYVal = significance.GetMinimum()#-0.2
      highYVal = significance.GetMaximum()#+0.2

      line1lims = [bumpLow,minYvalue,bumpLow,heightLowEdge]
      line2lims = [bumpHigh,minYvalue,bumpHigh,heightHighEdge]
      line3lims = [bumpLow,lowYVal,bumpLow,highYVal]
      line4lims = [bumpHigh,lowYVal,bumpHigh,highYVal]

      # Draw blue lines
      pad1.cd()
      line1.SetLineColor(ROOT.kBlue)
      line1.SetX1(line1lims[0]); line1.SetY1(line1lims[1]); line1.SetX2(line1lims[2]); line1.SetY2(line1lims[3])
      line1.Draw()
      line2.SetLineColor(ROOT.kBlue)
      line2.SetX1(line2lims[0]); line2.SetY1(line2lims[1]); line2.SetX2(line2lims[2]); line2.SetY2(line2lims[3])
      line2.Draw()
      pad2.cd()
      line3.SetLineColor(ROOT.kBlue)
      line3.SetX1(line3lims[0]); line3.SetY1(line3lims[1]); line3.SetX2(line3lims[2]); line3.SetY2(line3lims[3])
      line3.Draw()
      line4.SetLineColor(ROOT.kBlue)
      line4.SetX1(line4lims[0]); line4.SetY1(line4lims[1]); line4.SetX2(line4lims[2]); line4.SetY2(line4lims[3])
      line4.Draw()

    c.Update()

    outpad.cd()
    leftOfLegend = 0.48
    widthOfRow = 0.04
    lumInFb = round(float(luminosity)/float(1000),nsigfigs)
    if (doLogX and not doBumpLimits) :
      self.drawATLASLabels(0.2, 0.35)
      self.drawCMEAndLumi(0.51,0.90,CME,lumInFb,0.04)
      bottomOfLegend = 0.78
      legend = self.makeLegend(leftOfLegend,bottomOfLegend,0.9,0.87)
    else :
      self.drawATLASLabels(0.5, 0.87, True)
      self.drawCMEAndLumi(0.41,0.82,CME,lumInFb,0.04)
      bottomOfLegend = 0.70
      legend = self.makeLegend(leftOfLegend,bottomOfLegend,0.9,0.805)

    c.Update()

    self.myLatex.SetTextFont(42)
    self.myLatex.SetTextSize(0.04)
    index = 0
    persistent = []
    if len(extraLegendLines) > 0 :
#      toplocation = bottomOfLegend - (0.01+widthOfRow)*(index) #topOfAll - (0.03+2*widthOfRow) - (0.01+widthOfRow)*(index)
      for line in extraLegendLines :
        toplocation = bottomOfLegend - (0.01+widthOfRow)*(index) #topOfAll - (0.03+2*widthOfRow) - (0.01+widthOfRow)*(index)
        persistent.append(self.myLatex.DrawLatex(leftOfLegend+0.01,toplocation,line))
        index = index+1

    # Go to outer pad to fill and draw legend
    # Create legend
    outpad.cd()
    legend.AddEntry(dataHist,"Data","LFP")
    legend.AddEntry(fitHist,"Background fit","LF")
    if doBumpLimits :
      legend.AddEntry(line4,"BumpHunter interval","L")
    legend.Draw()
    c.Update()

    # Save.
    pad1.RedrawAxis()
    pad2.RedrawAxis()
    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)


  def drawHistsOverSignificanceHists(self,histograms,names,significancehists,xname,yname,sigy,name,luminosity,CME,xmin,xmax,ymin,ymax,doLogX=True,doLogY=True,doRectangular=False,doErrMain=False,doErrSig=False,sigHistRange=[]) :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,False)
    c.SetLogx(1)
    c.SetLogy(doLogX)
    c.SetGridx(0)
    c.SetGridy(0)

    outpad,pad1,pad2 = self.setStandardTwoPads()
    pad1.SetLogx(doLogX)
    pad2.SetLogx(doLogX)

    # Draw data and fit histograms
    pad1.cd()

    # Use bin range within which are all plot entries,
    # plus one empty on either side if available
    minmaxes = []
    index=-1
    for histlist in [histograms,significancehists] :
      index = index+1
      lowxvals = []
      lowyvals = []
      lownonzeros = []
      highxvals = []
      highyvals = []
      for histogram in histlist :
        lowx,highx = self.getAxisRangeFromHist(histogram)
        lowy,lownonzero,highy = self.getYRangeFromHist(histogram)
        lowxvals.append(lowx)
        highxvals.append(highx)
        lownonzeros.append(lownonzero)
        lowyvals.append(lowy)
        highyvals.append(highy)
      lowxvals.sort()
      lowyvals.sort()
      lownonzeros.sort()
      highxvals.sort()
      highyvals.sort()
      if xmin == 'automatic':
        minX = lowxvals[0]
      else :
        minX = xmin
      if xmax == 'automatic':
        maxX = highxvals[-1]
      else :
        maxX = xmax
      if ymin == 'automatic':
        if doLogY and index==0 :
          minY = lownonzeros[0]/2.0
        else:
          minY = lowyvals[0]
      else :
        minY = ymin
      if ymax == 'automatic':
        if doLogY and index==0:
          maxY = highyvals[-1]*100
        else :
          maxY = highyvals[-1]*1.5
      else :
        maxY = ymax
      minmaxes.append([minX,maxX,minY,maxY])

    goodcolours = self.getGoodColours(len(histograms))

    range = minmaxes[0]
    minX = range[0]; maxX = range[1]; minY = range[2]; maxY = range[3]
    for histogram in histograms :
      index = histograms.index(histogram)
      histogram.SetLineColor(goodcolours[index])
      histogram.SetMarkerColor(goodcolours[index])
      histogram.SetLineStyle(1)
      histogram.SetLineWidth(2)
      histogram.SetFillStyle(0)
      histogram.SetTitle("")
      histogram.GetXaxis().SetRange(minX,maxX+5)
      histogram.GetYaxis().SetRangeUser(minY,maxY)
      histogram.GetXaxis().SetNdivisions(605,ROOT.kTRUE)
      if (index==0) :
        histogram.GetYaxis().SetTitleSize(0.05)
        histogram.GetYaxis().SetTitleOffset(1.2) #1.0
        histogram.GetYaxis().SetLabelSize(0.05)

        histogram.GetXaxis().SetTitle(xname)
        histogram.GetYaxis().SetTitle(yname)
        if doErrMain :
          histogram.Draw("E")
        else :
          histogram.Draw("HIST")
      else :
        histogram.GetXaxis().SetTitle("")
        histogram.GetYaxis().SetTitle("")
        if doErrMain :
          histogram.Draw("E SAME")
        else :
          histogram.Draw("HIST SAME")

    # Draw significance histograms
    pad2.cd()
    range = minmaxes[1]
    minY = range[2]; maxY = range[3]
    if sigHistRange != [] :
      minY = sigHistRange[0]
      maxY = sigHistRange[1]
    # find nearest 0.25 to maxY
    for histogram in significancehists :
      index = significancehists.index(histogram)
      histogram.SetLineColor(goodcolours[index])
      histogram.SetMarkerColor(goodcolours[index])
      histogram.SetLineStyle(1)
      histogram.SetLineWidth(2)
      histogram.SetFillStyle(0)
      histogram.SetTitle("")
      histogram.GetXaxis().SetRange(minX,maxX+5)
      histogram.GetYaxis().SetRangeUser(minY,maxY)
      histogram.GetXaxis().SetNdivisions(605,ROOT.kTRUE)
      if (index==0) :
        histogram.GetYaxis().SetTitleSize(0.1)
        histogram.GetYaxis().SetTitleOffset(0.6) #0.42 # 1.2 = 20% larger
        histogram.GetYaxis().SetLabelSize(0.1)
        histogram.GetXaxis().SetLabelSize(0.1)
        histogram.GetXaxis().SetTitleSize(0.1)
        histogram.GetXaxis().SetTitleOffset(1.2)
        histogram.GetXaxis().SetNdivisions(805,ROOT.kTRUE)
        histogram.GetXaxis().SetTitle(xname)
        histogram.GetYaxis().SetTitle(sigy)
        if doErrSig :
          histogram.Draw("E")
        else :
          histogram.Draw("HIST")
      else :
        histogram.GetXaxis().SetTitle("")
        histogram.GetYaxis().SetTitle("")
        if doErrSig :
          histogram.Draw("E SAME")
        else :
          histogram.Draw("HIST SAME")
      if doLogX :
        self.fixTheBloodyTickMarks(ROOT.pad2, histogram, minX, maxX,minY,maxY)

    outpad.cd()
    persistent = []
    lumInFb = round(float(luminosity)/float(1000),nsigfigs)

    lshift = 0
    maxlen = 0
    for name in names :
      if len(name) > 12 and len(name) > maxlen :
        lshift = 0.0 - 0.01*(len(name)-12)
        maxlen = len(name)

    if (doLogX) :
      legend = self.makeLegend(0.60+lshift,0.75,0.9,0.87)
    else :
      legend = self.makeLegend(0.60+lshift,0.71,0.9,0.82)
    for histogram in histograms :
      legend.AddEntry(histogram,names[histograms.index(histogram)],"PL")
    if (doLogX) :
      self.drawATLASLabels(0.2, 0.35)
      self.drawCMEAndLumi(0.51,0.90,CME,lumInFb,0.04)
    else :
      self.drawATLASLabels(0.53, 0.88, True)
      self.drawCMEAndLumi(0.51,0.83,CME,lumInFb,0.04)

    # Go to outer pad to fill and draw legend
    # Create legend
    legend.Draw()

    # Save.
    pad1.RedrawAxis()
    pad2.RedrawAxis()
    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)

  def test(self,dataHist,fitHist,signalsignificance,residual,signalsForSpec,signalsForSig,signalmasses,legendlist,x,datay,sigy,residy,name,luminosity,CME,firstBin=-1,lastBin=-1,doBumpLimits=False,bumpLow=0,bumpHigh=0,doLogX=True,doRectangular=False,rightLegend=False, UserScaleText = "",writeOnpval = False, pval = -999, writeOnFit = False, FitMin =-999,FitMax =-999,mcHist=None,mcratioHist=None,mcupratioHist=None,mcdownratioHist=None) :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    #c = self.makeCanvas(canvasname,False)
    c = self.makeCanvas(canvasname,False,1,1.15)

    c.SetLogx(1)
    c.SetGridx(0)
    c.SetGridy(0)

    outpad = ROOT.TPad("extpad","extpad",0,0,1,1) # For marking outermost dimensions
    pad1 = ROOT.TPad("pad1","pad1",0,0.36,1,1) # For main histo
    pad2 = ROOT.TPad("pad2","pad2",0,0.23,1,0.36) # For residuals histo
    pad3 = ROOT.TPad("pad3","pad3",0,0,1,0.23) # For MC comparison histo

    # Set up to draw in right orientations
    outpad.SetFillStyle(4000) #transparent
    pad1.SetBottomMargin(0.00001)
    pad1.SetBorderMode(0)
    pad1.SetLogy(1)
    pad1.SetLogx(doLogX)

    pad2.SetTopMargin(0.00001)
    pad2.SetBottomMargin(0.00001)
    pad2.SetBorderMode(0)
    pad2.SetLogx(doLogX)

    pad3.SetTopMargin(0.00001)
    pad3.SetBottomMargin(0.43)
    pad3.SetBorderMode(0)
    pad3.SetLogx(doLogX)
    pad1.Draw()
    pad2.Draw()
    pad3.Draw()

    # Publication-friendly margins
    pad1.SetLeftMargin(0.1)
    pad2.SetLeftMargin(0.1)
    pad3.SetLeftMargin(0.1)

    pad1.SetTopMargin(0.02)

    pad1.SetRightMargin(0.02)
    pad2.SetRightMargin(0.02)
    pad3.SetRightMargin(0.02)

    outpad.Draw()

    # Use bin range within which bkgPlot has entries,
    # plus one empty on either side if available
    lowbin,highbin = self.getAxisRangeFromHist(dataHist)
    if (firstBin>0) :
      lowbin=firstBin
    if (lastBin>0 and lastBin>=firstBin) :
      highbin = lastBin

    ## Add a few more bins on high end if we need that extra legend space
    highbin = highbin+17# SWITCH

    # Draw data and fit histograms (and MC if applicable)
    pad1.cd()
    fitHist.GetYaxis().SetTitleSize(0.06)
    fitHist.GetYaxis().SetTitleOffset(0.8)
    fitHist.GetYaxis().SetLabelSize(0.05)

    self.drawSignalOverlaidOnDataAndFit(dataHist,fitHist,signalsForSpec,signalmasses,[],luminosity,CME,datay,"",firstBin,lastBin,doLogX,True,False,False,3)

    # Lydia adding observedStat value to plot
    if self.dodrawUsersText:
      if writeOnFit:
        if writeOnpval:
          if doBumpLimits:
            if "it{q}" in UserScaleText and "BM" in UserScaleText:
              self.drawUsersText(0.3,0.2,"#splitline{#splitline"+UserScaleText+"}{#splitline{#it{p}-value = "+str(round(pval,2))+"}{#splitline{Fit Range: "+str(round(FitMin/1000,1))+" - "+str(round(FitMax/1000,1))+" TeV}{"+"{0}}}}".format(self.cutstring),0.045)
            else:
              self.drawUsersText(0.3,0.2,"#splitline{"+UserScaleText+"}{#splitline{#it{p}-value = "+str(round(pval,2))+"}{#splitline{Fit Range: "+str(round(FitMin/1000,1))+" - "+str(round(FitMax/1000,1))+" TeV}{"+"{0}}}}".format(self.cutstring),0.045)
          else:
            if "it{q}" in UserScaleText and "BM" in UserScaleText:
              self.drawUsersText(0.15,0.2,"#splitline{#splitline"+UserScaleText+"}{#splitline{#it{p}-value = "+str(round(pval,2))+"}{#splitline{Fit Range: "+str(round(FitMin/1000,1))+" - "+str(round(FitMax/1000,1))+" TeV}{"+"{0}}}}".format(self.cutstring),0.045)
            else:
              self.drawUsersText(0.15,0.2,"#splitline{"+UserScaleText+"}{#splitline{#it{p}-value = "+str(round(pval,2))+"}{#splitline{Fit Range: "+str(round(FitMin/1000,1))+" - "+str(round(FitMax/1000,1))+" TeV}{"+"{0}}}}".format(self.cutstring),0.045)
        else:
          if doBumpLimits:
            self.drawUsersText(0.15,0.2,"#splitline{"+UserScaleText+"}{#splitline{Fit Range: "+str(round(FitMin/1000,1))+" - "+str(round(FitMax/1000,1))+" TeV}{"+"{0}}}}".format(self.cutstring),0.045)
          else:
            self.drawUsersText(0.15,0.2,"#splitline{"+UserScaleText+"}{#splitline{Fit Range: "+str(round(FitMin/1000,1))+" - "+str(round(FitMax/1000,1))+" TeV}{"+"{0}}}}".format(self.cutstring),0.045)
      else:
        if writeOnpval:
          self.drawUsersText(0.15,0.2,"#splitline{"+UserScaleText+"}{#it{p}-value = "+str(round(pval,2))+"}",0.06)
        elif UserScaleText != "":
          #if "{" in UserScaleText:
          #  self.drawUsersText(0.23,0.81,"#splitline{#splitline"+UserScaleText+"}{|y*| < 0.6}",0.055)
            #self.drawUsersText(0.32,0.12,"#splitline{#splitline"+UserScaleText+"}{|y*| < 0.6}",0.055)
          #else:
          self.drawUsersText(0.23,0.84,"#splitline{"+UserScaleText+"}{|y*| < 0.6}",0.055)
          #self.drawUsersText(0.32,0.12,"#splitline{"+UserScaleText+"}{|y*| < 0.6}",0.055)
          # Lydia adding text to say user how much signal scaled by
#      self.drawUsersText(0.15,0.25, UserScaleText,0.04)
    pad1.Update()

    line1 = self.line.Clone("line1"); line1lims = []
    line2 = self.line.Clone("line2"); line2lims = []
    line3 = self.line.Clone("line3"); line3lims = []
    line4 = self.line.Clone("line4"); line4lims = []

    if doBumpLimits :
      heightLowEdge=0
      heightHighEdge=0
      minYvalue = dataHist.GetMinimum()
      for i in range(dataHist.GetNbinsX()) :
        locationOfTallEdge = dataHist.GetBinLowEdge(i)
        height = dataHist.GetBinContent(i)
        if locationOfTallEdge == bumpLow :
          heightLowEdge = dataHist.GetBinContent(i)
        if locationOfTallEdge == bumpHigh:
          heightHighEdge = dataHist.GetBinContent(i-1)

      lowYVal = signalsignificance.GetMinimum()#-0.2
      highYVal = signalsignificance.GetMaximum()+0.5

      line1lims = [bumpLow,minYvalue,bumpLow,heightLowEdge]
      line2lims = [bumpHigh,minYvalue,bumpHigh,heightHighEdge]
      line3lims = [bumpLow,lowYVal,bumpLow,highYVal]
      line4lims = [bumpHigh,lowYVal,bumpHigh,highYVal]

      # Draw blue lines
      pad1.cd()
      line1.SetLineColor(ROOT.kBlue)
      line1.SetLineWidth(2)
      line1.SetX1(line1lims[0]); line1.SetY1(line1lims[1]); line1.SetX2(line1lims[2]); line1.SetY2(line1lims[3])
      line1.Draw()
      line2.SetLineColor(ROOT.kBlue)
      line2.SetLineWidth(2)
      line2.SetX1(line2lims[0]); line2.SetY1(line2lims[1]); line2.SetX2(line2lims[2]); line2.SetY2(line2lims[3])
      line2.Draw()
      pad2.cd()
      line3.SetLineColor(ROOT.kBlue)
      line3.SetLineWidth(2)
      line3.SetX1(line3lims[0]); line3.SetY1(line3lims[1]); line3.SetX2(line3lims[2]); line3.SetY2(line3lims[3])
      line4.SetLineColor(ROOT.kBlue)
      line4.SetLineWidth(2)
      line4.SetX1(line4lims[0]); line4.SetY1(line4lims[1]); line4.SetX2(line4lims[2]); line4.SetY2(line4lims[3])

    c.Update()

    # Draw residual histogram
    pad2.cd()
    residual.GetYaxis().SetTitleSize(0.19)
    residual.GetYaxis().SetTitleOffset(0.2) # 1.2 = 20% larger
    residual.GetYaxis().SetLabelSize(0.16)
    #residual.GetYaxis().SetNdivisions(5,10,0)
    residual.GetYaxis().SetNdivisions(604)
    residual.GetYaxis().CenterTitle()

    residual.GetXaxis().SetLabelSize(0.13)
    residual.GetXaxis().SetTitleSize(0.17)
    residual.GetXaxis().SetTitleOffset(1.2)

    if residual.GetBinLowEdge(firstBin) > 0.001 and residual.GetBinLowEdge(firstBin) < 1 :
      residual.GetXaxis().SetNoExponent(ROOT.kTRUE)
    self.drawSignificanceHist(residual,firstBin,lastBin,x,residy,True,True)
    pad2.Update()

    pad2.cd()
    line3.Draw()
    line4.Draw()

    # Draw mc comparison residual histogram
    pad3.cd()
    mcratioHist.GetYaxis().SetTitleSize(0.1)
    mcratioHist.GetYaxis().SetTitleOffset(0.4) # 1.2 = 20% larger
    mcratioHist.GetYaxis().CenterTitle()
    mcratioHist.GetYaxis().SetLabelSize(0.09)

    # TEST
    #mcratioHist.GetYaxis().SetNdivisions(5,10,0)
    mcratioHist.GetXaxis().SetLabelSize(0.15)
    mcratioHist.GetXaxis().SetTitleSize(0.17)
    mcratioHist.GetXaxis().SetTitleOffset(1.2)

    if mcratioHist.GetBinLowEdge(firstBin) > 0.001 and mcratioHist.GetBinLowEdge(firstBin) < 1 :
      mcratioHist.GetXaxis().SetNoExponent(ROOT.kTRUE)
    mcratioHist.GetYaxis().SetRangeUser(-0.5,0.5)
    mcratioHist.GetYaxis().SetNdivisions(5,10,0)
    self.drawSignificanceHistWithJESBands(mcratioHist,mcupratioHist,mcdownratioHist,firstBin,lastBin,x,sigy,True,True)
    jeslegend = self.persistentlegend
    jeslegend.Draw()
    pad3.Update()

    # Go to outer pad to fill and draw legend
    # Create legend
    outpad.cd()

    lumInFb = round(float(luminosity)/float(1000),nsigfigs)
    widthOfRow = 0.037# SWITCH 0.036

    if (doLogX and not rightLegend) :
      topOfLegend = 0.878#0.865 SWITCH
      leftOfLegend = 0.475 # 0.47 SWITCH
      if self.labeltype == 0:
        self.drawATLASLabels(leftOfLegend+0.02, topOfLegend+0.055, True)
      else:
        self.drawATLASLabels(leftOfLegend-0.05, topOfLegend+0.055, True)
      self.drawCMEAndLumi(leftOfLegend-0.07,topOfLegend+0.015,CME,lumInFb,0.04)
      #self.drawCMEAndLumi(0.51,0.82,CME,lumInFb,0.04)

      bottomOfLegend = topOfLegend - (widthOfRow*(len(self.saveplots)+2))
      if doBumpLimits :
        bottomOfLegend = topOfLegend - (widthOfRow*(len(self.saveplots)+3))
    else :
      #self.drawATLASLabels(0.2, 0.4)
      #self.drawATLASLabels(0.5, 0.4)
      topOfLegend = 0.878 # 0.865 SWITCH
      leftOfLegend = 0.515 #0.51 SWITCH
      if self.labeltype == 0:
        self.drawATLASLabels(leftOfLegend+0.02, topOfLegend+0.055, True)
      else:
        self.drawATLASLabels(leftOfLegend-0.05, topOfLegend+0.055, True)
      self.drawCMEAndLumi(leftOfLegend-0.07,topOfLegend+0.015,CME,lumInFb,0.04)
      bottomOfLegend = topOfLegend-(widthOfRow*(len(self.saveplots)+2))
      if doBumpLimits :
        bottomOfLegend = topOfLegend-(widthOfRow*(len(self.saveplots)+3))

    rightOfLegend = leftOfLegend+0.29
    legend = self.makeLegend(leftOfLegend,bottomOfLegend,rightOfLegend,topOfLegend,0.038)#,0.0375)
    legend.SetFillStyle(0)
    legend.AddEntry(dataHist,"Data","P")
    legend.AddEntry(fitHist,"Background fit","LF")
    if doBumpLimits :
      legend.AddEntry(line4,"BumpHunter interval","L")
    for plot in self.saveplots :
      index = self.saveplots.index(plot)
      legend.AddEntry(plot,legendlist[index],"LP")
    legend.SetEntrySeparation(0.5)
    legend.Draw()

    # Save.
    c.Update()
    ROOT.gPad.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)


  def drawDataAndFitWithSignalsOverSignificances(self,dataHist,fitHist,signalsignificance,residual,signalsForSpec,signalsForSig,signalmasses,legendlist,x,datay,sigy,residy,name,luminosity,CME,firstBin=-1,lastBin=-1,doBumpLimits=False,bumpLow=0,bumpHigh=0,doLogX=True,doRectangular=False,rightLegend=False, UserScaleText = "",writeOnpval = False, pval = -999, writeOnFit = False, FitMin =-999,FitMax =-999,mcHist=None) :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,False)
    c.SetLogx(1)
    c.SetGridx(0)
    c.SetGridy(0)

    drawMC = False
    if not mcHist==None :
      drawMC=True

    # Dimensions: xlow, ylow, xup, yup
    outpad = ROOT.TPad("extpad","extpad",0,0,1,1) # For marking outermost dimensions
    pad1 = ROOT.TPad("pad1","pad1",0,0.3,1,1) # For main histo
    #pad2 = ROOT.TPad("pad2","pad2",0,0.23,1,0.36) # For signal significance histo
    #pad3 = ROOT.TPad("pad3","pad3",0,0,1,0.23) # For residuals histo
    pad3 = ROOT.TPad("pad3","pad3",0,0,1,0.30) # For residuals histo

    # Set up to draw in right orientations
    outpad.SetFillStyle(4000) #transparent
    pad1.SetBottomMargin(0.00001)
    pad1.SetBorderMode(0)
    pad1.SetLogy(1)
    pad1.SetLogx(doLogX)
    #pad2.SetTopMargin(0.00001)
    #pad2.SetBottomMargin(0.00001)
    #pad2.SetBorderMode(0)
    #pad2.SetLogx(doLogX)
    pad3.SetTopMargin(0.00001)
    pad3.SetBottomMargin(0.43)
    pad3.SetBorderMode(0)
    pad3.SetLogx(doLogX)
    pad1.Draw()
    #pad2.Draw()
    pad3.Draw()

    # Publication-friendly margins
    pad1.SetLeftMargin(0.1)
#    pad2.SetLeftMargin(0.2)
    pad3.SetLeftMargin(0.1)
    pad1.SetTopMargin(0.02)
    pad1.SetRightMargin(0.02)
    pad3.SetRightMargin(0.02)
    outpad.Draw()

    # Use bin range within which bkgPlot has entries,
    # plus one empty on either side if available
    lowbin,highbin = self.getAxisRangeFromHist(dataHist)
    if (firstBin>0) :
      lowbin=firstBin
    if (lastBin>0 and lastBin>=firstBin) :
      highbin = lastBin

    ## Add a few more bins on high end if we need that extra legend line
    if drawMC :
      highbin = highbin+3

    # Draw data and fit histograms (and MC if applicable)
    pad1.cd()
    fitHist.GetYaxis().SetTitleSize(0.06)
    fitHist.GetYaxis().SetTitleOffset(0.8)
    fitHist.GetYaxis().SetLabelSize(0.05)

    if drawMC :
      self.drawSignalOverlaidOnDataAndFit(dataHist,fitHist,signalsForSpec,signalmasses,[],luminosity,CME,datay,"",firstBin,lastBin,doLogX,True,False,False,3,mcHist)
    else :
      self.drawSignalOverlaidOnDataAndFit(dataHist,fitHist,signalsForSpec,signalmasses,[],luminosity,CME,datay,"",firstBin,lastBin,doLogX,True,False,False,3)

    # Lydia adding observedStat value to plot
    if self.dodrawUsersText:
      if writeOnFit:
        if writeOnpval:
          if doBumpLimits:
            if "it{q}" in UserScaleText and "BM" in UserScaleText:
              self.drawUsersText(0.3,0.2,"#splitline{#splitline"+UserScaleText+"}{#splitline{#it{p}-value = "+str(round(pval,2))+"}{#splitline{Fit Range: "+str(round(FitMin/1000,1))+" - "+str(round(FitMax/1000,1))+" TeV}{"+"{0}}}}".format(self.cutstring),0.045)
            else:
              self.drawUsersText(0.3,0.2,"#splitline{"+UserScaleText+"}{#splitline{#it{p}-value = "+str(round(pval,2))+"}{#splitline{Fit Range: "+str(round(FitMin/1000,1))+" - "+str(round(FitMax/1000,1))+" TeV}{"+"{0}}}}".format(self.cutstring),0.045)
          else:
            if "it{q}" in UserScaleText and "BM" in UserScaleText:
              self.drawUsersText(0.15,0.2,"#splitline{#splitline"+UserScaleText+"}{#splitline{#it{p}-value = "+str(round(pval,2))+"}{#splitline{Fit Range: "+str(round(FitMin/1000,1))+" - "+str(round(FitMax/1000,1))+" TeV}{"+"{0}}}}".format(self.cutstring),0.045)
            else:
              self.drawUsersText(0.15,0.2,"#splitline{"+UserScaleText+"}{#splitline{#it{p}-value = "+str(round(pval,2))+"}{#splitline{Fit Range: "+str(round(FitMin/1000,1))+" - "+str(round(FitMax/1000,1))+" TeV}{"+"{0}}}}".format(self.cutstring),0.045)
        else:
          if doBumpLimits:
            self.drawUsersText(0.15,0.2,"#splitline{"+UserScaleText+"}{#splitline{Fit Range: "+str(round(FitMin/1000,1))+" - "+str(round(FitMax/1000,1))+" TeV}{"+"{0}}}}".format(self.cutstring),0.045)
          else:
            self.drawUsersText(0.15,0.2,"#splitline{"+UserScaleText+"}{#splitline{Fit Range: "+str(round(FitMin/1000,1))+" - "+str(round(FitMax/1000,1))+" TeV}{"+"{0}}}}".format(self.cutstring),0.045)
      else:
        if writeOnpval:
          self.drawUsersText(0.15,0.2,"#splitline{"+UserScaleText+"}{#it{p}-value = "+str(round(pval,2))+"}",0.06)
        elif UserScaleText != "":
          #if "{" in UserScaleText:
          #  self.drawUsersText(0.23,0.81,"#splitline{#splitline"+UserScaleText+"}{|y*| < 0.6}",0.055)
            #self.drawUsersText(0.32,0.12,"#splitline{#splitline"+UserScaleText+"}{|y*| < 0.6}",0.055)
          #else:
          self.drawUsersText(0.23,0.84,"#splitline{"+UserScaleText+"}{"+"{0}}".format(self.cutstring),0.055)
          #self.drawUsersText(0.32,0.12,"#splitline{"+UserScaleText+"}{|y*| < 0.6}",0.055)
          # Lydia adding text to say user how much signal scaled by
#      self.drawUsersText(0.15,0.25, UserScaleText,0.04)
    pad1.Update()

    # Draw signals significance histogram
    #pad2.cd()
    #signalsignificance.GetYaxis().SetTitleSize(0.21)
    #signalsignificance.GetYaxis().SetTitleOffset(0.18) # 1.2 = 20% larger
    #signalsignificance.GetYaxis().SetLabelSize(0.2)
    #self.drawSignalOverlaidOnBkgPlot(signalsignificance,signalsForSig,signalmasses,[],luminosity,CME,sigy,"",firstBin,lastBin,False,True,False)

    #pad2.Update()

    line1 = self.line.Clone("line1"); line1lims = []
    line2 = self.line.Clone("line2"); line2lims = []
    line3 = self.line.Clone("line3"); line3lims = []
    line4 = self.line.Clone("line4"); line4lims = []

    if doBumpLimits :
      heightLowEdge=0
      heightHighEdge=0
      minYvalue = dataHist.GetMinimum()
      for i in range(dataHist.GetNbinsX()) :
        locationOfTallEdge = dataHist.GetBinLowEdge(i)
        height = dataHist.GetBinContent(i)
        if locationOfTallEdge == bumpLow :
          heightLowEdge = dataHist.GetBinContent(i)
        if locationOfTallEdge == bumpHigh:
          heightHighEdge = dataHist.GetBinContent(i-1)

      lowYVal = signalsignificance.GetMinimum()#-0.2
      highYVal = signalsignificance.GetMaximum()+0.5

      line1lims = [bumpLow,minYvalue,bumpLow,heightLowEdge]
      line2lims = [bumpHigh,minYvalue,bumpHigh,heightHighEdge]
      line3lims = [bumpLow,lowYVal,bumpLow,highYVal]
      line4lims = [bumpHigh,lowYVal,bumpHigh,highYVal]

      # Draw blue lines
      pad1.cd()
      line1.SetLineColor(ROOT.kBlue)
      line1.SetLineWidth(2)
      line1.SetX1(line1lims[0]); line1.SetY1(line1lims[1]); line1.SetX2(line1lims[2]); line1.SetY2(line1lims[3])
      line1.Draw()
      line2.SetLineColor(ROOT.kBlue)
      line2.SetLineWidth(2)
      line2.SetX1(line2lims[0]); line2.SetY1(line2lims[1]); line2.SetX2(line2lims[2]); line2.SetY2(line2lims[3])
      line2.Draw()
      pad3.cd()
      line3.SetLineColor(ROOT.kBlue)
      line3.SetLineWidth(2)
      line3.SetX1(line3lims[0]); line3.SetY1(line3lims[1]); line3.SetX2(line3lims[2]); line3.SetY2(line3lims[3])
      line4.SetLineColor(ROOT.kBlue)
      line4.SetLineWidth(2)
      line4.SetX1(line4lims[0]); line4.SetY1(line4lims[1]); line4.SetX2(line4lims[2]); line4.SetY2(line4lims[3])

    c.Update()

    # Draw residual histogram
    pad3.cd()
    residual.GetYaxis().SetTitleSize(0.12)
    residual.GetYaxis().SetTitleOffset(0.32) # 1.2 = 20% larger
    residual.GetYaxis().SetLabelSize(0.115)

    # TEST
    #residual.GetYaxis().SetNdivisions(5,10,0)
    residual.GetYaxis().SetNdivisions(604)

    residual.GetXaxis().SetLabelSize(0.15)
    residual.GetXaxis().SetTitleSize(0.17)
    residual.GetXaxis().SetTitleOffset(1.2)
    if residual.GetBinLowEdge(firstBin) > 0.001 and residual.GetBinLowEdge(firstBin) < 1 :
      residual.GetXaxis().SetNoExponent(ROOT.kTRUE)
    self.drawSignificanceHist(residual,firstBin,lastBin,x,residy,True,True)
    pad3.Update()

    pad3.cd()
    line3.Draw()
    line4.Draw()
    # Go to outer pad to fill and draw legend
    # Create legend
    outpad.cd()

    lumInFb = round(float(luminosity)/float(1000),nsigfigs)
    widthOfRow = 0.0415

    if (doLogX and not rightLegend) :
      topOfLegend = 0.85
      leftOfLegend = 0.46
      if len(legendlist) != 0:
        if "(QBH)" in legendlist[0]:
          leftOfLegend = 0.44
      if self.labeltype == 0:
        self.drawATLASLabels(leftOfLegend+0.02, topOfLegend+0.065, True)
      else:
        self.drawATLASLabels(leftOfLegend-0.05, topOfLegend+0.065, True)
      self.drawCMEAndLumi(leftOfLegend-0.07,topOfLegend+0.015,CME,lumInFb,0.04)
      #self.drawCMEAndLumi(0.51,0.82,CME,lumInFb,0.04)

      bottomOfLegend = topOfLegend - (widthOfRow*(len(self.saveplots)+2))
      if doBumpLimits :
        bottomOfLegend = topOfLegend - (widthOfRow*(len(self.saveplots)+3))
    else :
      #self.drawATLASLabels(0.2, 0.4)
      #self.drawATLASLabels(0.5, 0.4)
      topOfLegend = 0.85
      leftOfLegend = 0.53
      if "Z'" in legendlist[0]:
        leftOfLegend = 0.51
      if self.labeltype == 0:
        self.drawATLASLabels(leftOfLegend+0.02, topOfLegend+0.065, True)
      else:
        self.drawATLASLabels(leftOfLegend-0.05, topOfLegend+0.065, True)
      self.drawCMEAndLumi(leftOfLegend-0.07,topOfLegend+0.015,CME,lumInFb,0.04)
      bottomOfLegend = topOfLegend-(widthOfRow*(len(self.saveplots)+2))
      if doBumpLimits :
        bottomOfLegend = topOfLegend-(widthOfRow*(len(self.saveplots)+3))

    if drawMC :
      bottomOfLegend = bottomOfLegend - widthOfRow

    rightOfLegend = leftOfLegend+0.36
    legend = self.makeLegend(leftOfLegend,bottomOfLegend,rightOfLegend,topOfLegend)
    legend.SetFillStyle(0)
    legend.AddEntry(dataHist,"Data","P")
    legend.AddEntry(fitHist,"Background fit","LF")
    if drawMC :
      legend.AddEntry(mcHist,"SM MC","LF")
    if doBumpLimits :
      legend.AddEntry(line4,"BumpHunter interval","L")
    for plot in self.saveplots :
      index = self.saveplots.index(plot)
      legend.AddEntry(plot,legendlist[index],"LP")
    legend.SetEntrySeparation(0.5)
    legend.Draw()

    # Save.
    c.Update()
    ROOT.gPad.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)

  def drawMultipleFitsAndResiduals(self,dataHist,fitHistList,residualList,legendlist,x,datay,residyList,name,luminosity,CME,firstBin=-1,lastBin=-1,doBumpLimits=False,bumpLow=0,bumpHigh=0,doLogX=True,doRectangular=False,notLogY=False,lowY=11,highY=-1) :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,False)
    c.SetLogx(1)

    if notLogY :
      c.SetLogy(0)
    c.SetGridx(0)
    c.SetGridy(0)

    # Dimensions: xlow, ylow, xup, yup
    outpad = ROOT.TPad("extpad","extpad",0,0,1,1) # For marking outermost dimensions
    pads = []
    if len(residualList) == 1 :
      padsize = 0.2
    elif len(residualList) == 2 :
      padsize = 0.13
    elif len(residualList) == 3 :
      padsize = 0.10
    else :
      padsize = 0.4/len(residualList)
    topOfSubplots = 0.1 + padsize * len(residualList)
    for ipad in range(len(residualList)+1) :
      padname = "pad_{0}".format(ipad)
      if ipad == 0 :
        pad = ROOT.TPad(padname,padname,0,topOfSubplots,1,1) # for main histo
      elif ipad!= len(residualList) :
        pad = ROOT.TPad(padname,padname,0,topOfSubplots - ipad*padsize, 1, topOfSubplots - (ipad-1)*padsize)
      else :
        pad = ROOT.TPad(padname,padname,0, 0, 1, topOfSubplots - (ipad-1)*padsize)
      pads.append(pad)

    # Set up to draw in right orientations
    outpad.SetFillStyle(4000) #transparent
    for pad in pads :
      pad.SetBorderMode(0)
      pad.SetLogx(doLogX)
      if pads.index(pad)==0 :
        pad.SetBottomMargin(0.00001)
        pad.SetLogy(1)
        if (notLogY) :
          pad.SetLogy(0)
      elif pads.index(pad)==len(pads)-1 :
        pad.SetTopMargin(0.00001)
        pad.SetBottomMargin(0.1/(0.1+padsize))
      else :
        pad.SetTopMargin(0.00001)
        pad.SetBottomMargin(0.00001)
      pad.Draw()
    outpad.Draw()

    # Use bin range within which bkgPlot has entries,
    # plus one empty on either side if available
    lowbin,highbin = self.getAxisRangeFromHist(dataHist)
    if (firstBin>0) :
      lowbin=firstBin
    if (lastBin>0 and lastBin>=firstBin) :
      highbin = lastBin

    goodcolours = self.getGoodColours(len(fitHistList))

    # Draw data and fit histograms
    pads[0].cd()

    self.drawDataHist(dataHist,lowbin,highbin,"",datay,False,3) # was True
    fitHistList[0].GetYaxis().SetTitleSize(0.06)
    fitHistList[0].GetYaxis().SetTitleOffset(0.8)
    fitHistList[0].GetYaxis().SetLabelSize(0.05)
    for hist in fitHistList :
      self.drawFitHist(hist, lowbin, highbin, "", datay, True, True, False,[],True,\
              goodcolours[fitHistList.index(hist)], 1)
    pads[0].Update()

    # Draw residual histograms
    for index in range(len(residualList)) :
      pad = pads[index+1]
      pad.cd()
      residual = residualList[index]
      if index != len(residualList)-1 :
        residual.GetYaxis().SetTitleSize(0.21)
        residual.GetYaxis().SetTitleOffset(0.18) # 1.2 = 20% larger
        residual.GetYaxis().SetLabelSize(0.2)
      else :
        if (len(residualList) > 1) :
          residual.GetYaxis().SetTitleSize(0.12)
          residual.GetYaxis().SetTitleOffset(0.32) # 1.2 = 20% larger
          residual.GetYaxis().SetLabelSize(0.115)
          residual.GetXaxis().SetLabelSize(0.15)
          residual.GetXaxis().SetTitleSize(0.17)
          residual.GetXaxis().SetTitleOffset(1.2)
        else :
          residual.GetYaxis().SetTitleSize(0.1)
          residual.GetYaxis().SetTitleOffset(0.42) # 1.2 = 20% larger
          residual.GetYaxis().SetLabelSize(0.1)
          residual.GetXaxis().SetLabelSize(0.1)
          residual.GetXaxis().SetTitleSize(0.1)
          residual.GetXaxis().SetTitleOffset(1.2)

      if residual.GetBinLowEdge(firstBin) > 0.001 and residual.GetBinLowEdge(firstBin) < 1 :
        residual.GetXaxis().SetNoExponent(ROOT.kTRUE)
        #significance,firstBin,lastBin,xname,yname,fixYAxis=False,inLargerPlot=False,doLogX=False,doErrors=False,fillColour = ROOT.kRed
      self.drawSignificanceHist(residual,firstBin,lastBin,x,residyList[index],True,True,True,False,goodcolours[index])
      pad.Update()

    # Go to outer pad to fill and draw legend
    # Create legend
    outpad.cd()

    lumInFb = round(float(luminosity)/float(1000),nsigfigs)
    widthOfRow = 0.05
    if (doLogX) :
      self.drawATLASLabels(0.53, 0.88, True)
      bottomOfLegend = topOfSubplots + 0.02
      leftOfLegend = 0.2
      self.drawCMEAndLumi(0.51,0.82,CME,lumInFb,0.04)
      topOfLegend = bottomOfLegend + (widthOfRow*(len(fitHistList)+1))
    else :
      self.drawATLASLabels(0.2, 0.4)
      topOfLegend = 0.87
      leftOfLegend = 0.5
      self.drawCMEAndLumi(leftOfLegend,topOfLegend+0.02,CME,lumInFb,0.04)
      bottomOfLegend = topOfLegend-(widthOfRow*(len(fitHistList)+1))
    rightOfLegend = leftOfLegend+0.4
    legend = self.makeLegend(leftOfLegend,bottomOfLegend,rightOfLegend,topOfLegend)

    legend.AddEntry(dataHist,"Data","LFP")
    for fit in fitHistList :
      legend.AddEntry(fit,legendlist[fitHistList.index(fit)],"LF")
    legend.Draw()

    # Save.
    c.Update()
    for pad in pads :
      pad.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)

  # Emma's version ;)
  def drawLimitSettingPlotObservedExpected(self,observed,expected, expected1sigma,expected2sigma,signals,signalslegend,name,nameX,nameY,luminosity,CME,xmin,xmax,ymin,ymax,doRectangular=False,drawExistingLimit = False, ExistingLimit = 0, ExistingLimitName = "",doCanvas = True,lineWidth = 3) :

    if type(signals) is not list :
      signals = [signals]
      signalslegend = [signalslegend]

    if doCanvas :
      canvasname = name+'_cv'
      outputname = name+epsorpdf
      if saveCfile: Coutputname = name+'.C'
      if saveRfile: Routputname = name+'.root'
      if saveEfile: Eoutputname = name+'.eps'
      c = self.makeCanvas(canvasname,doRectangular)
      (c.SetLogx(0),c.SetLogy(1),c.SetGridx(0),c.SetGridy(0))

    # Set automatic axis range from graphs.
    # X axis range will be exactly ends of graphs
    xVals = [observed.GetX()[i] for i in range(observed.GetN())]
    xVals.sort()

    minX = xVals[0] if xmin == 'automatic' else xmin
    maxX = xVals[-1]  if xmax == 'automatic' else xmax

    # Y axis range will be 3 orders of magnitude
    # above highest point of observed, because want space for legend
    # Lowest point should be 2 orders of magnitude below lowest point in observed

    minY = observed.GetMinimum()/100  if ymin == 'automatic' else ymin
    maxY = observed.GetMaximum()*1000 if ymax == 'automatic' else ymax

    print "minX, maxX are:",minX,maxX

    # Set axis names.
    # So far, should always be thus so don't pass as parameters.

    # Create legend At Bottom Left
    if doCanvas :
      leftOfLegend = 0.18
      bottomOfLegend = 0.18#20
      topOfLegend = bottomOfLegend+ 0.055*(3+len(signals)) # 0.64
      self.persistentlegend = self.makeLegend(leftOfLegend,bottomOfLegend,0.5,topOfLegend) # 0.88
    else :     # Create legend At Top Right
      leftOfLegend = 0.8
      rightOfLegend = 1
      topOfLegend = 0.85

      if len(signalslegend) ==1:
        lensigleg = len(signalslegend[0]) # Count letters in legend, used below to calculate legend position
        if "it" in signalslegend[0]:
          lensigleg = lensigleg-5
        leftOfLegend = leftOfLegend - 0.03*lensigleg
      else:
        leftOfLegend = 0.57
        rightOfLegend = 0.9
        topOfLegend = 0.9
      bottomOfLegend = topOfLegend-0.1*(len(signals)-1) # 0.64
      self.persistentlegend = self.makeLegend(leftOfLegend,bottomOfLegend,rightOfLegend,topOfLegend) # 0.88
      self.persistentlegend.SetTextSize(0.075)

    # Set up display for expectations
    for graph,colour in [[expected2sigma,self.colourpalette.twoSigmaBandColour],[expected1sigma,self.colourpalette.oneSigmaBandColour],[expected,self.colourpalette.oneSigmaBandColour]] :
      graph.SetMarkerColor(1)
      graph.SetMarkerSize(1)
      graph.SetMarkerStyle(20)
      graph.SetLineColor(1)
      graph.SetLineWidth(lineWidth)
      graph.SetLineStyle(3)
      graph.SetFillColor(colour)
      graph.GetXaxis().SetTitle(nameX)
      graph.GetYaxis().SetTitle(nameY)
      graph.GetXaxis().SetLimits(minX,maxX)
      graph.GetXaxis().SetNdivisions(705,ROOT.kTRUE)
      graph.GetYaxis().SetRangeUser(minY,maxY)
      graph.GetXaxis().SetTitleOffset(1.3)
      if minX > 0.001 and minX < 1 :
        graph.GetXaxis().SetNoExponent(ROOT.kTRUE)

    # Set up display for signal
    for signal in signals :
      thiscolour = self.colourpalette.signalLineColours[signals.index(signal)]
      thiserrorcolour = self.colourpalette.signalErrorColours[signals.index(signal)]
      signal.SetMarkerColor(4)
      signal.SetMarkerSize(1)
      signal.SetMarkerStyle(24)
      signal.SetLineColor(thiscolour)
      signal.SetLineWidth(lineWidth)
      signal.GetXaxis().SetTitleOffset(1.3)
      if lineWidth > 2 :
        signal.SetLineStyle(9 - signals.index(signal))
      else :
        signal.SetLineStyle(7 - signals.index(signal))
      signal.SetFillColor(0) #thiserrorcolour)
      signal.GetXaxis().SetTitle(nameX)
      signal.GetXaxis().SetNdivisions(705,ROOT.kTRUE)
      signal.GetYaxis().SetTitle(nameY)
      signal.GetYaxis().SetRangeUser(minY,maxY)
      signal.GetYaxis().SetLimits(minY,maxY)
      if minX > 0.001 and minX < 1 :
        signal.GetXaxis().SetNoExponent(ROOT.kTRUE)

    # Set up display for observations
    observed.SetMarkerColor(1)
    observed.SetMarkerSize(1)
    observed.SetMarkerStyle(20)
    observed.SetLineColor(1)
    observed.SetLineWidth(lineWidth)
    observed.SetLineStyle(1)
    observed.SetFillColor(0)
    observed.GetXaxis().SetTitle(nameX)
    observed.GetYaxis().SetTitle(nameY)
    observed.GetYaxis().SetRangeUser(minY,maxY)
    observed.GetXaxis().SetLimits(minX,maxX)
    observed.GetXaxis().SetTitleOffset(1.3)
    if minX > 0.001 and minX < 1 :
      observed.GetXaxis().SetNoExponent(ROOT.kTRUE)

    # First one has to include axes or everything comes out blank
    # Rest have to NOT include axes or each successive one overwrites
    # previous. "SAME option does not exist for TGraph classes.
    #observed.Draw("APL") # Data points of measurement
    expected2sigma.Draw("AF") # 2-sigma expectation error bands
    expected1sigma.Draw("F") # 1-sigma expectation error bands
    expected.Draw("LX") # Center of expectation
    for signal in signals :
      signal.Draw("03") #L03
    for signal in signals :
      signal.Draw("LX") # was CX
    observed.Draw("PL") # Data points of measurement

    # Draw arrow to existing limit
    if drawExistingLimit:
      arrow = ROOT.TArrow()
      arrow.SetLineColor(ROOT.kRed)
      arrow.SetFillColor(ROOT.kRed)
      arrow.SetLineWidth(2)
      #arrow.DrawArrow(ExistingLimit,0.0032,ExistingLimit,0)

      arrow.DrawArrow(ExistingLimit,0.002,ExistingLimit,0) # Matches to 1E-3
      #arrow.DrawArrow(ExistingLimit,0.00025,ExistingLimit,0) # Matches to 1E-4

      self.persistentlegend.AddEntry(arrow,ExistingLimitName,"L")

    # Fill and draw legend
    for signal in signals :
      index = signals.index(signal)
      if signalslegend != [] and signalslegend !=[[]] :
        self.persistentlegend.AddEntry(signal,signalslegend[index],"LF")#"L")
    ## When not doing canvas
    if not doCanvas :
      self.persistentlegend.Draw()
    ## When doing canvas
    else :
      self.persistentlegend.AddEntry(observed,"Observed 95% CL upper limit","PL")
      self.persistentlegend.AddEntry(expected1sigma, "Expected 95% CL upper limit","L")
      self.persistentlegend.AddEntry( "NULL" , "68% and 95% bands","")

      self.drawATLASLabels(0.58,0.88)
      self.persistentlegend.Draw()

      shadeBox = ROOT.TBox()

      # Legend in bottom left-hand corner
      boxX1 = minX + (maxX - minX)*0.045#25#3#2 #21
      boxX2 = boxX1 + (maxX - minX)*0.0735 # 135

      if doCanvas:
        if c.GetLogy() :
          boxY1 = math.exp(math.log(minY) + (math.log(maxY) - math.log(minY))*(0.035))#5))#35-0.06*(len(signals)-1))) #0.66
          boxY2 = math.exp(math.log(boxY1) + (math.log(maxY) - math.log(minY))*0.0155)
          boxY3 = math.exp(math.log(boxY2) + (math.log(maxY) - math.log(minY))*0.025)
          boxY4 = math.exp(math.log(boxY3) + (math.log(maxY) - math.log(minY))*0.0155)
      else :
        boxY1 = minY + (maxY - minY)*(0.205-0.06*(len(signals)-1)) #0.66
        boxY2 = boxY1 + (maxY - minY)*0.0155
        boxY3 = boxY2 + (maxY - minY)*0.025
        boxY4 = boxY3 + (maxY - minY)*0.0155

      shadeBox.SetFillColor(self.colourpalette.twoSigmaBandColour)
      shadeBox.DrawBox(boxX1,boxY1,boxX2,boxY4)
      shadeBox.SetFillColor(self.colourpalette.oneSigmaBandColour)
      shadeBox.DrawBox(boxX1,boxY2,boxX2,boxY3)

      lumInFb = round(float(luminosity)/float(1000),nsigfigs)
      self.drawCMEAndLumi(0.5,0.825,CME,lumInFb,0.04)

    # Lydia adding analysis cuts values to plot
    if self.dodrawUsersText and doCanvas:
      self.drawUsersText(0.585,0.775,self.cutstring,0.04)

    if doCanvas :
      c.RedrawAxis()
      c.Update()
      c.SaveAs(outputname)
      c.SaveAs(name+".png")
      if saveCfile:
        c.SaveSource(Coutputname)
      if saveRfile:
        c.SaveSource(Routputname)
      if saveEfile:
        c.SaveAs(Eoutputname)

    if len(signals)==0:
      return [None, None]
    elif len(signals)==1:
      signal = signals[0]
      obsLimits = self.calculateIntersectionOfGraphs(signal,observed,True,True)
      expLimits = self.calculateIntersectionOfGraphs(signal,expected1sigma,True,True)
      return [obsLimits,expLimits]
    else :
      output = []
      for signal in signals :
        obsLimits = self.calculateIntersectionOfGraphs(signal,observed,True,True)
        expLimits = self.calculateIntersectionOfGraphs(signal,expected1sigma,True,True)
        output.append([obsLimits,expLimits])
      return output

  def drawLimitSettingPlot2Sigma(self,observed,expected1sigma,expected2sigma,signals,signalslegend,name,nameX,nameY,luminosity,CME,xmin,xmax,ymin,ymax,doRectangular=False,drawExistingLimit = False, ExistingLimit = 0, ExistingLimitName = "",doCanvas = True,lineWidth = 3) :

    if type(signals) is not list :
      signals = [signals]
      signalslegend = [signalslegend]

    if doCanvas :
      canvasname = name+'_cv'
      outputname = name+epsorpdf
      if saveCfile:
        Coutputname = name+'.C'
      if saveRfile:
        Routputname = name+'.root'
      if saveEfile:
        Eoutputname = name+'.eps'
      c = self.makeCanvas(canvasname,doRectangular)
      c.SetLogx(0)
      c.SetLogy(1)
      c.SetGridx(0)
      c.SetGridy(0)

    # Set automatic axis range from graphs.
    # X axis range will be exactly ends of graphs
    xVals = []
    for i in range(observed.GetN()) :
      xVals.append(observed.GetX()[i])
    xVals.sort()
    if xmin == 'automatic':
      minX = xVals[-1]
    else :
      print "setting minx = ",xmin
      minX = xmin
    if xmax == 'automatic':
      maxX = xVals[0]
    else :
      maxX = xmax
    # Y axis range will be 3 orders of magnitude
    # above highest point of observed, because want space for legend
    # Lowest point should be 2 orders of magnitude below lowest point in observed
    if ymin == 'automatic' :
      minY = observed.GetMinimum()/100
    else :
      minY = ymin
    if ymax == 'automatic' :
      maxY = observed.GetMaximum()*1000
    else :
      maxY = ymax

    print "minX, maxX are:",minX,maxX

    # Set axis names.
    # So far, should always be thus so don't pass as parameters.

    # Create legend At Bottom Left
    if doCanvas :
      leftOfLegend = 0.18
      bottomOfLegend = 0.18#20
      topOfLegend = bottomOfLegend+ 0.055*(3+len(signals)) # 0.64
      self.persistentlegend = self.makeLegend(leftOfLegend,bottomOfLegend,0.5,topOfLegend) # 0.88
    else :     # Create legend At Top Right
      leftOfLegend = 0.8
      rightOfLegend = 1
      topOfLegend = 0.85

      if len(signalslegend) ==1:
        lensigleg = len(signalslegend[0]) # Count letters in legend, used below to calculate legend position
        if "it" in signalslegend[0]:
          lensigleg = lensigleg-5
        leftOfLegend = leftOfLegend - 0.03*lensigleg
      else:
        leftOfLegend = 0.57
        rightOfLegend = 0.9
        topOfLegend = 0.9
      bottomOfLegend = topOfLegend-0.1*(len(signals)-1) # 0.64
      self.persistentlegend = self.makeLegend(leftOfLegend,bottomOfLegend,rightOfLegend,topOfLegend) # 0.88
      self.persistentlegend.SetTextSize(0.075)
#      leftOfLegend = 0.65
#      topOfLegend = 0.90
#      bottomOfLegend = topOfLegend - 0.055*(len(signals)) # 0.64
#      self.persistentlegend = self.makeLegend(leftOfLegend,bottomOfLegend,1.0,topOfLegend) # 0.88
      #leftOfLegend = 0.06 + ROOT.gPad.GetLeftMargin()
      #bottomOfLegend = 0.08 + ROOT.gPad.GetBottomMargin()
      #topOfLegend = bottomOfLegend+ 0.055*(len(signals)) # 0.64
      #self.persistentlegend = self.makeLegend(leftOfLegend,bottomOfLegend,0.0,topOfLegend,0.05) # 0.88

    # Set up display for expectations
    for graph,colour in [[expected2sigma,self.colourpalette.twoSigmaBandColour],[expected1sigma,self.colourpalette.oneSigmaBandColour]] :
      graph.SetMarkerColor(1)
      graph.SetMarkerSize(1)
      graph.SetMarkerStyle(20)
      graph.SetLineColor(1)
      graph.SetLineWidth(lineWidth)
      graph.SetLineStyle(3)
      graph.SetFillColor(colour)
      graph.GetXaxis().SetTitle(nameX)
      graph.GetYaxis().SetTitle(nameY)
      graph.GetXaxis().SetLimits(minX,maxX)
      graph.GetXaxis().SetNdivisions(705,ROOT.kTRUE)
      graph.GetYaxis().SetRangeUser(minY,maxY)
      graph.GetXaxis().SetTitleOffset(1.3)
      if minX > 0.001 and minX < 1 :
        graph.GetXaxis().SetNoExponent(ROOT.kTRUE)

    # Set up display for signal
    for signal in signals :
      thiscolour = self.colourpalette.signalLineColours[signals.index(signal)]
      thiserrorcolour = self.colourpalette.signalErrorColours[signals.index(signal)]
      signal.SetMarkerColor(4)
      signal.SetMarkerSize(1)
      signal.SetMarkerStyle(24)
      signal.SetLineColor(thiscolour)
      signal.SetLineWidth(lineWidth)
      signal.GetXaxis().SetTitleOffset(1.3)
      if lineWidth > 2 :
        signal.SetLineStyle(9 - signals.index(signal))
      else :
        signal.SetLineStyle(7 - signals.index(signal))
      signal.SetFillColor(0) #thiserrorcolour)
      signal.GetXaxis().SetTitle(nameX)
      signal.GetXaxis().SetNdivisions(705,ROOT.kTRUE)
      signal.GetYaxis().SetTitle(nameY)
      signal.GetYaxis().SetRangeUser(minY,maxY)
      signal.GetYaxis().SetLimits(minY,maxY)
      if minX > 0.001 and minX < 1 :
        signal.GetXaxis().SetNoExponent(ROOT.kTRUE)

    # Set up display for observations
    observed.SetMarkerColor(1)
    observed.SetMarkerSize(1)
    observed.SetMarkerStyle(20)
    observed.SetLineColor(1)
    observed.SetLineWidth(lineWidth)
    observed.SetLineStyle(1)
    observed.SetFillColor(0)
    observed.GetXaxis().SetTitle(nameX)
    observed.GetYaxis().SetTitle(nameY)
    observed.GetYaxis().SetRangeUser(minY,maxY)
    observed.GetXaxis().SetLimits(minX,maxX)
    observed.GetXaxis().SetTitleOffset(1.3)
    if minX > 0.001 and minX < 1 :
      observed.GetXaxis().SetNoExponent(ROOT.kTRUE)

    # First one has to include axes or everything comes out blank
    # Rest have to NOT include axes or each successive one overwrites
    # previous. "SAME option does not exist for TGraph classes.
    expected2sigma.Draw("A3") # 2-sigma expectation error bands
    expected1sigma.Draw("L3") # 1-sigma expectation error bands
    expected1sigma.Draw("LX") # Center of expectation
    for signal in signals :
      signal.Draw("03") #L03
    for signal in signals :
      signal.Draw("LX") # was CX
    observed.Draw("PL") # Data points of measurement

    # Draw arrow to existing limit
    if drawExistingLimit:
      arrow = ROOT.TArrow()
      arrow.SetLineColor(ROOT.kRed)
      arrow.SetFillColor(ROOT.kRed)
      arrow.SetLineWidth(2)
      #arrow.DrawArrow(ExistingLimit,0.0032,ExistingLimit,0)

      arrow.DrawArrow(ExistingLimit,0.002,ExistingLimit,0) # Matches to 1E-3
      #arrow.DrawArrow(ExistingLimit,0.00025,ExistingLimit,0) # Matches to 1E-4

      self.persistentlegend.AddEntry(arrow,ExistingLimitName,"L")

    # Fill and draw legend
    for signal in signals :
      index = signals.index(signal)
      if signalslegend != [] and signalslegend !=[[]] :
        self.persistentlegend.AddEntry(signal,signalslegend[index],"LF")#"L")
    ## When not doing canvas
    if not doCanvas :
      self.persistentlegend.Draw()
    ## When doing canvas
    else :
      self.persistentlegend.AddEntry(observed,"Observed 95% CL upper limit","PL")
      self.persistentlegend.AddEntry(expected1sigma, "Expected 95% CL upper limit","L")
      self.persistentlegend.AddEntry( "NULL" , "68% and 95% bands","")

      self.drawATLASLabels(0.58,0.88)
      self.persistentlegend.Draw()

      shadeBox = ROOT.TBox()

      # Legend in bottom left-hand corner
      boxX1 = minX + (maxX - minX)*0.045#25#3#2 #21
      boxX2 = boxX1 + (maxX - minX)*0.0735 # 135

      if doCanvas:
        if c.GetLogy() :
          boxY1 = math.exp(math.log(minY) + (math.log(maxY) - math.log(minY))*(0.035))#5))#35-0.06*(len(signals)-1))) #0.66
          boxY2 = math.exp(math.log(boxY1) + (math.log(maxY) - math.log(minY))*0.0155)
          boxY3 = math.exp(math.log(boxY2) + (math.log(maxY) - math.log(minY))*0.025)
          boxY4 = math.exp(math.log(boxY3) + (math.log(maxY) - math.log(minY))*0.0155)
      else :
        boxY1 = minY + (maxY - minY)*(0.205-0.06*(len(signals)-1)) #0.66
        boxY2 = boxY1 + (maxY - minY)*0.0155
        boxY3 = boxY2 + (maxY - minY)*0.025
        boxY4 = boxY3 + (maxY - minY)*0.0155


    # Legend in top right-hand corner
    #boxX1 = minX + (maxX - minX)*0.19 #21
    #boxX2 = boxX1 + (maxX - minX)*0.125 # 135

#    if c.GetLogy() :
#      boxY1 = math.exp(math.log(minY) + (math.log(maxY) - math.log(minY))*(0.665-0.06*(len(signals)-1))) #0.66
#      boxY2 = math.exp(math.log(boxY1) + (math.log(maxY) - math.log(minY))*0.0155)
#      boxY3 = math.exp(math.log(boxY2) + (math.log(maxY) - math.log(minY))*0.025)
#      boxY4 = math.exp(math.log(boxY3) + (math.log(maxY) - math.log(minY))*0.0155)
#    else :
#      boxY1 = minY + (maxY - minY)*(0.665-0.06*(len(signals)-1)) #0.66
#      boxY2 = boxY1 + (maxY - minY)*0.0155
#      boxY3 = boxY2 + (maxY - minY)*0.025
#      boxY4 = boxY3 + (maxY - minY)*0.0155

      shadeBox.SetFillColor(self.colourpalette.twoSigmaBandColour)
      shadeBox.DrawBox(boxX1,boxY1,boxX2,boxY4)
      shadeBox.SetFillColor(self.colourpalette.oneSigmaBandColour)
      shadeBox.DrawBox(boxX1,boxY2,boxX2,boxY3)

      lumInFb = round(float(luminosity)/float(1000),nsigfigs)
      self.drawCMEAndLumi(0.5,0.825,CME,lumInFb,0.04)

    # Lydia adding analysis cuts values to plot
    if self.dodrawUsersText and doCanvas:
      self.drawUsersText(0.585,0.775,self.cutstring,0.04)

    if doCanvas :
      c.RedrawAxis()
      c.Update()
      c.SaveAs(outputname)
      if saveCfile:
        c.SaveSource(Coutputname)
      if saveRfile:
        c.SaveSource(Routputname)
      if saveEfile:
        c.SaveAs(Eoutputname)

    if len(signals)==0:
      return
    elif len(signals)==1:
      signal = signals[0]
      obsLimits = self.calculateIntersectionOfGraphs(signal,observed,True,True)
      expLimits = self.calculateIntersectionOfGraphs(signal,expected1sigma,True,True)
      return [obsLimits,expLimits]
    else :
      output = []
      for signal in signals :
        obsLimits = self.calculateIntersectionOfGraphs(signal,observed,True,True)
        expLimits = self.calculateIntersectionOfGraphs(signal,expected1sigma,True,True)
        output.append([obsLimits,expLimits])
      return output


  def drawFourLimitPlots_Grid(self,plot1Materials,plot2Materials,plot3Materials,plot4Materials,name,nameX,nameY,luminosity,CME,xmin1,xmax1,xmin2,xmax2,ymin1,ymax1,ymin2,ymax2) :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    #c = self.makeCanvas(canvasname,False,2,2.5)
    c = self.makeCanvas(canvasname,False,2,2)
    c.SetLogx(0)
    c.SetLogy(1)
    c.SetGridx(0)
    c.SetGridy(0)

    # Dimensions: xlow, ylow, xup, yup
    outpad = ROOT.TPad("extpad","extpad",0,0,1,1) # For marking outermost dimensions
    outpad = ROOT.TPad("extpad","extpad",0,0,1,1) # For marking outermost dimensions
    pad1 = ROOT.TPad("pad1","pad1",0,0.52,0.525,0.97) # For first histo
    pad2 = ROOT.TPad("pad2","pad2",0.525,0.52,1,0.97) # For second histo
    pad3 = ROOT.TPad("pad3","pad3",0,0,0.525,0.52) # For third histo
    pad4 = ROOT.TPad("pad4","pad4",0.525,0,1,0.52) # For fourth histo
    #pad1 = ROOT.TPad("pad1","pad1",0,0.57,0.525,0.97) # For first histo
    #pad2 = ROOT.TPad("pad2","pad2",0.525,0.57,1,0.97) # For second histo
    #pad3 = ROOT.TPad("pad3","pad3",0,0.15,0.525,0.57) # For third histo
    #pad4 = ROOT.TPad("pad4","pad4",0.525,0.15,1,0.57) # For fourth histo

    # Set up to draw in right orientations
    outpad.SetFillStyle(4000) #transparent
    pad1.SetBottomMargin(0.00001)
    pad1.SetRightMargin(0.00001)
    pad1.SetBorderMode(0)
    pad1.SetLogy(1)
    pad1.SetLogx(0)
    pad2.SetLeftMargin(0.00001)
    pad2.SetBottomMargin(0.00001)
    pad2.SetBorderMode(0)
    pad2.SetLogy(1)
    pad2.SetLogx(0)
    pad3.SetTopMargin(0.00001)
    pad3.SetRightMargin(0.00001)
    pad3.SetBorderMode(0)
    pad3.SetLogy(1)
    pad3.SetLogx(0)
    pad4.SetTopMargin(0.00001)
    pad4.SetLeftMargin(0.00001)
    pad4.SetBorderMode(0)
    pad4.SetLogy(1)
    pad4.SetLogx(0)

    pad1.Draw()
    pad2.Draw()
    pad3.Draw()
    pad4.Draw()
    outpad.Draw()

    pad1.cd()
    self.drawLimitSettingPlot2Sigma(plot1Materials["observed"],plot1Materials["expected_1sigma"],\
          plot1Materials["expected_2sigma"],[plot1Materials["signal"],plot1Materials["extrasignal"],plot1Materials["extraextrasignal"]],[plot1Materials["signalLabel"],plot1Materials["extrasignalLabel"],plot1Materials["extraextrasignalLabel"]],"test",nameX,nameY,luminosity,CME,xmin1,xmax1,ymin1,ymax1,False,False,0,"",False,2)
    savelegend1 = self.persistentlegend
    savelegend1.Draw()
    c.Update()

    pad2.cd()
    # Removed plot2Materials["signalLabel"]
    self.drawLimitSettingPlot2Sigma(plot2Materials["observed"],plot2Materials["expected_1sigma"],\
          plot2Materials["expected_2sigma"],plot2Materials["signal"],plot2Materials["signalLabel"],"",nameX,nameY,luminosity,CME,xmin2,xmax2,ymin1,ymax1,False,False,0,"",False,2)
    savelegend2 = self.persistentlegend
    savelegend2.Draw()
    c.Update()

    pad3.cd()
    self.drawLimitSettingPlot2Sigma(plot3Materials["observed"],plot3Materials["expected_1sigma"],\
          plot3Materials["expected_2sigma"],plot3Materials["signal"],plot3Materials["signalLabel"],"",nameX,nameY,luminosity,CME,xmin1,xmax1,ymin2,ymax2,False,False,0,"",False,2)
    savelegend3 = self.persistentlegend
    savelegend3.Draw()
    c.Update()

    pad4.cd()
    self.drawLimitSettingPlot2Sigma(plot4Materials["observed"],plot4Materials["expected_1sigma"],\
          plot4Materials["expected_2sigma"],plot4Materials["signal"],plot4Materials["signalLabel"],"",nameX,nameY,luminosity,CME,xmin2,xmax2,ymin2,ymax2,False,False,0,"",False,2)
    savelegend4 = self.persistentlegend
    savelegend4.Draw()
    c.Update()

    # Make new legend for stripes etc
    #outpad.cd()
    pad2.cd()
    leftOfLegend = 0.05#05
    topOfLegend = 0.25#90
    bottomOfLegend = topOfLegend - 0.2 # 0.64
    savelegend5 = self.makeLegend(leftOfLegend,bottomOfLegend,0.92,topOfLegend,0.065)#04) # 0.88
#    savelegend5.AddEntry(plot2Materials["signal"],plot2Materials["signalLabel"],"LF")#"L")
    savelegend5.AddEntry(plot2Materials["observed"],"Observed 95% CL upper limit","PL")
    savelegend5.AddEntry(plot2Materials["expected_1sigma"], "Expected 95% CL upper limit","L")
    savelegend5.AddEntry( "NULL" , "68% and 95% bands","")
    savelegend5.SetMargin(0.15)
    savelegend5.Draw()


    # Legend in centre bottom of pad2 set coordinates using outpad though!
    shadeBox = ROOT.TBox()
    outpad.cd()

    boxX1 = 0.56#leftOfLegend#0.207 #21
    boxX2 = boxX1 + 0.042 # 135

    boxY1 = 0.54#bottomOfLegend #0.66
    boxY2 = boxY1 + 0.009
    boxY3 = boxY2 + 0.0139
    boxY4 = boxY3 + 0.009

    shadeBox.SetFillColor(self.colourpalette.twoSigmaBandColour)
    shadeBox.DrawBox(boxX1,boxY1,boxX2,boxY4)
    shadeBox.SetFillColor(self.colourpalette.oneSigmaBandColour)
    shadeBox.DrawBox(boxX1,boxY2,boxX2,boxY3)

    c.Update()

    # Add ATLAS labels
    # Add lumi etc
    outpad.cd()
    lumInFb = round(float(luminosity)/float(1000),nsigfigs)
    #self.drawATLASLabels(0.1,0.96,False,False,0.033)
    #self.drawCMEAndLumi(0.65,0.96,CME,lumInFb,0.028)

    # ATLAS labels and CME and lumi outside plot on top!
    #self.drawATLASLabels(0.08,0.96,False,False,0.036)
    #self.drawCMEAndLumi(0.65,0.96,CME,lumInFb,0.032)
    # ATLAS labels and CME and lumi right hand plot!
    #self.drawATLASLabels(0.72,0.85,False,False,0.036)
    #self.drawLumiAndCMEVert(0.72,0.73,lumInFb,CME,0.032)
    # ATLAS labels and CME and lumi left hand plot!
    self.drawATLASLabels(0.105,0.9,False,False,0.036)
    #self.drawLumiAndCMEVert(0.11,0.78,lumInFb,CME,0.032)
    p1 = self.drawCME(0.10,0.86,CME,0.032)
    p2 = self.drawLumi(0.105,0.82,lumInFb,0.032)
    if self.dodrawUsersText :
      self.drawUsersText(0.1,0.78,self.cutstring,0.032)
    c.Update()

    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)


  def drawThreePlusOneLimitPlots_Grid(self,plot1Materials,plot2Materials,plot3Materials,TwoDPlotMaterials,name,nameX,nameY,luminosity,CME,xmin1,xmax1,xmin2,xmax2,ymin1,ymax1,ymin2,ymax2) :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    # For right aligned
    c = self.makeCanvas(canvasname,False,2,2)
    # For not right aligned
#    c = self.makeCanvas(canvasname,False,2.2,2)
    c.SetLogx(0)
    c.SetLogy(1)
    c.SetGridx(0)
    c.SetGridy(0)

    # Dimensions: xlow, ylow, xup, yup for right-hand stuff aligned
    outpad = ROOT.TPad("extpad","extpad",0,0,1,1) # For marking outermost dimensions
    pad1 = ROOT.TPad("pad1","pad1",0,0.5,0.525,1) # For first histo
    pad2 = ROOT.TPad("pad2","pad2",0.525,0.445,1,1) # For second histo
    pad3 = ROOT.TPad("pad3","pad3",0,0,0.525,0.5) # For third histo
    pad4 = ROOT.TPad("pad4","pad4",0.525,0,1,0.445) # For fourth histo
     # For not right aligned
#    outpad = ROOT.TPad("extpad","extpad",0,0,1,1) # For marking outermost dimensions
#    pad1 = ROOT.TPad("pad1","pad1",0,0.5,0.475,1) # For first histo
#    pad2 = ROOT.TPad("pad2","pad2",0.475,0.5,1,1) # For second histo
#    pad3 = ROOT.TPad("pad3","pad3",0,0,0.475,0.5) # For third histo
#    pad4 = ROOT.TPad("pad4","pad4",0.475,0,1,0.5) # For fourth histo

    # Set up to draw in right orientations
    outpad.SetFillStyle(4000) #transparent
    pad1.SetBottomMargin(0.00001)
    pad1.SetRightMargin(0.00001)
    pad1.SetBorderMode(0)
    pad1.SetLogy(1)
    pad1.SetLogx(0)
    pad2.SetLeftMargin(0.00001)

    # Below: for when we want right hand stuff not aligned
#    pad2.SetRightMargin(0.2)
    # Below: for when we want right hand stuff aligned
    pad2.SetBottomMargin(0.1)
    pad2.SetTopMargin(0.045)

    pad2.SetBorderMode(0)
    pad2.SetLogy(1)
    pad2.SetLogx(0)
    pad3.SetTopMargin(0.00001)
    pad3.SetRightMargin(0.00001)
    pad3.SetBorderMode(0)
    pad3.SetLogy(1)
    pad3.SetLogx(0)

    # Below: for when we want right hand stuff aligned
    pad4.SetRightMargin(0.23)
    pad4.SetBottomMargin(0.19)

    # Below: for when we want right hand stuff not aligned
#    pad4.SetRightMargin(0.2)
    pad4.SetBorderMode(0)
    pad4.SetLogy(0)
    pad4.SetLogx(0)

    pad1.Draw()
    pad2.Draw()
    pad3.Draw()
    pad4.Draw()
    outpad.Draw()

    # Tiny square version
    plot2Materials["expected_2sigma"].GetXaxis().SetTitleOffset(0.92)
    plot3Materials["expected_2sigma"].GetXaxis().SetTitleOffset(0.92)
    TwoDPlotMaterials["hist"].GetXaxis().SetTitleOffset(0.92)

    pad1.cd()
    self.drawLimitSettingPlot2Sigma(plot1Materials["observed"],plot1Materials["expected_1sigma"],\
          plot1Materials["expected_2sigma"],plot1Materials["signal"],plot1Materials["signalLabel"],"test",nameX,nameY,luminosity,CME,xmin1,xmax1,ymin1,ymax1,False,False,0,"",False,2)
    savelegend1 = self.persistentlegend
    savelegend1.Draw()
    c.Update()

    pad2.cd()
    # Removed plot2Materials["signalLabel"]
    self.drawLimitSettingPlot2Sigma(plot2Materials["observed"],plot2Materials["expected_1sigma"],\
          plot2Materials["expected_2sigma"],plot2Materials["signal"],plot2Materials["signalLabel"],"",nameX,nameY,luminosity,CME,xmin2,xmax2,ymin1,ymax1,False,False,0,"",False,2)
    savelegend2 = self.persistentlegend
    savelegend2.Draw()
    c.Update()

    pad3.cd()
    self.drawLimitSettingPlot2Sigma(plot3Materials["observed"],plot3Materials["expected_1sigma"],\
          plot3Materials["expected_2sigma"],plot3Materials["signal"],plot3Materials["signalLabel"],"",nameX,nameY,luminosity,CME,xmin1,xmax1,ymin2,ymax2,False,False,0,"",False,2)
    savelegend3 = self.persistentlegend
    savelegend3.Draw()
    c.Update()


    pad4.cd()
    TwoDPlotMaterials["hist"].Draw("colz")
    TwoDPlotMaterials["hist"].GetZaxis().SetRangeUser(0,4)
    TwoDPlotMaterials["hist"].GetZaxis().SetTitle(TwoDPlotMaterials["zAxisName"])
    TwoDPlotMaterials["hist"].GetZaxis().SetTitleOffset(1.40)
    TwoDPlotMaterials["hist"].GetXaxis().SetRangeUser(1.25,3.75)
    TwoDPlotMaterials["hist"].GetXaxis().SetTitleOffset(1.40)
    TwoDPlotMaterials["hist"].GetYaxis().SetRangeUser(0,0.5)
    TwoDPlotMaterials["hist"].GetYaxis().SetTitle(TwoDPlotMaterials["yAxisName"])
    TwoDPlotMaterials["hist"].GetXaxis().SetTitle(TwoDPlotMaterials["xAxisName"])
    TwoDPlotMaterials["hist"].GetYaxis().SetTitleOffset(1.40)
    TwoDPlotMaterials["hist"].GetYaxis().SetLabelSize(0.075)
    TwoDPlotMaterials["hist"].Draw("textsame")
    c.Update()

    # Make new legend for stripes etc
    pad2.cd()
    leftOfLegend = 0.1#0.05
    topOfLegend = 0.90
    bottomOfLegend = topOfLegend - 0.15 # 0.64
    savelegend5 = self.makeLegend(leftOfLegend,bottomOfLegend,1.0,topOfLegend,0.05) # 0.88
    savelegend5.AddEntry(plot2Materials["observed"],"Observed 95% CL upper limit","PL")
    savelegend5.AddEntry(plot2Materials["expected_1sigma"], "Expected 95% CL upper limit","L")
    savelegend5.AddEntry( "NULL" , "68% and 95% bands","")
    savelegend5.SetMargin(0.15)
    savelegend5.Draw()

    shadeBox = ROOT.TBox()

    # Legend in top right-hand corner
    boxX1 = xmin2 + (xmax2 - xmin2)*leftOfLegend + 0.15 #0.35
    boxX2 = boxX1 + 0.5 # 0.60
    if c.GetLogy() :
      boxY1 = math.exp(math.log(ymin2) + (math.log(ymax2) - math.log(ymin2))*(topOfLegend-0.15)) #0.66
      boxY2 = math.exp(math.log(boxY1) + (math.log(ymax2) - math.log(ymin2))*0.0155)
      boxY3 = math.exp(math.log(boxY2) + (math.log(ymax2) - math.log(ymin2))*0.025)
      boxY4 = math.exp(math.log(boxY3) + (math.log(ymax2) - math.log(ymin2))*0.0155)
    else :
      boxY1 = ymin2 + (ymax2 - ymin2)*(topOfLegend-0.06) #0.66
      boxY2 = boxY1 + (ymax2 - ymin2)*0.0155
      boxY3 = boxY2 + (ymax2 - ymin2)*0.025
      boxY4 = boxY3 + (ymax2 - ymin2)*0.0155

    shadeBox.SetFillColor(self.colourpalette.twoSigmaBandColour)
    shadeBox.DrawBox(boxX1,boxY1,boxX2,boxY4)
    shadeBox.SetFillColor(self.colourpalette.oneSigmaBandColour)
    shadeBox.DrawBox(boxX1,boxY2,boxX2,boxY3)

    c.Update()

    pad1.cd()
    lumInFb = round(float(luminosity)/float(1000),nsigfigs)
    self.drawATLASLabels(0.55,0.85,False,False,0.06)
    self.drawCMEAndLumi(0.47,0.80,CME,lumInFb,0.05)
    c.Update()

    pad2.RedrawAxis()
    pad2.RedrawAxis()
    pad3.RedrawAxis()
    c.Update()

    # Add ATLAS labels
    # Add lumi etc
    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)


  def drawFourLimitPlots_Line(self,plot1Materials,plot2Materials,plot3Materials,plot4Materials,name,nameY,luminosity,CME,xmin1,xmax1,xmin2,xmax2,xmin3,xmax3,xmin4,xmax4,ymin,ymax) :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,False,4,1.2)
    c.SetLogx(0)
    c.SetLogy(1)
    c.SetGridx(0)
    c.SetGridy(0)

    # Dimensions: xlow, ylow, xup, yup
    outpad = ROOT.TPad("extpad","extpad",0,0,1,1) # For marking outermost dimensions
    pad1 = ROOT.TPad("pad1","pad1",0,0,0.27,1) # For first histo
    pad2 = ROOT.TPad("pad2","pad2",0.27,0,0.27+0.73/3.0,1) # For second histo
    pad3 = ROOT.TPad("pad3","pad3",0.27+0.73/3.0,0,0.27+(2.0*0.73)/3.0,1) # For third histo
    pad4 = ROOT.TPad("pad4","pad4",0.27+(2.0*0.73)/3.0,0,1,1) # For fourth histo

    # Set up to draw in right orientations
    outpad.SetFillStyle(4000) #transparent
    pad1.SetRightMargin(0.00001)
    pad1.SetBorderMode(0)
    pad1.SetLogy(1)
    pad1.SetLogx(0)
    pad2.SetLeftMargin(0.00001)
    pad2.SetRightMargin(0.00001)
    pad2.SetBorderMode(0)
    pad2.SetLogy(1)
    pad2.SetLogx(0)
    pad3.SetLeftMargin(0.00001)
    pad3.SetRightMargin(0.00001)
    pad3.SetBorderMode(0)
    pad3.SetLogy(1)
    pad3.SetLogx(0)
    pad4.SetLeftMargin(0.00001)
    pad4.SetBorderMode(0)
    pad4.SetLogy(1)
    pad4.SetLogx(0)

    pad1.Draw()
    pad2.Draw()
    pad3.Draw()
    pad4.Draw()
    outpad.Draw()

    pad1.cd()
    self.drawLimitSettingPlot2Sigma(plot1Materials["observed"],plot1Materials["expected_1sigma"],\
          plot1Materials["expected_2sigma"],plot1Materials["signal"],plot1Materials["signalLabel"],"test",plot1Materials["xAxisName"],nameY,luminosity,CME,xmin1,xmax1,ymin,ymax,False,False,0,"",False,1)
    savelegend1 = self.persistentlegend
    savelegend1.Draw()
    c.Update()

    pad2.cd()
    self.drawLimitSettingPlot2Sigma(plot2Materials["observed"],plot2Materials["expected_1sigma"],\
          plot2Materials["expected_2sigma"],plot2Materials["signal"],plot2Materials["signalLabel"],"",plot2Materials["xAxisName"],nameY,luminosity,CME,xmin2,xmax2,ymin,ymax,False,False,0,"",False,1)
    savelegend2 = self.persistentlegend
    savelegend2.Draw()
    c.Update()

    pad3.cd()
    self.drawLimitSettingPlot2Sigma(plot3Materials["observed"],plot3Materials["expected_1sigma"],\
          plot3Materials["expected_2sigma"],plot3Materials["signal"],plot3Materials["signalLabel"],"",plot3Materials["xAxisName"],nameY,luminosity,CME,xmin3,xmax3,ymin,ymax,False,False,0,"",False,1)
    savelegend3 = self.persistentlegend
    savelegend3.Draw()
    c.Update()

    pad4.cd()
    self.drawLimitSettingPlot2Sigma(plot4Materials["observed"],plot4Materials["expected_1sigma"],plot4Materials["expected_2sigma"],plot4Materials["signal"],plot4Materials["signalLabel"],"",plot4Materials["xAxisName"],nameY,luminosity,CME,xmin4,xmax4,ymin,ymax,False,False,0,"",False,1)
    savelegend4 = self.persistentlegend
    savelegend4.Draw()
    c.Update()

    # Make new legend for stripes etc
    # Make new legend for stripes etc
    pad4.cd()
    leftOfLegend = 0.05 # 0.25
    topOfLegend = 0.90
    bottomOfLegend = topOfLegend - 0.15 # 0.64
    savelegend5 = self.makeLegend(leftOfLegend,bottomOfLegend,1.0,topOfLegend,0.05) # 0.88
#    savelegend5.AddEntry(plot4Materials["signal"],plot4Materials["signalLabel"],"LF")#"L")
    savelegend5.AddEntry(plot4Materials["observed"],"Observed 95% CL upper limit","PL")
    savelegend5.AddEntry(plot4Materials["expected_1sigma"], "Expected 95% CL upper limit","L")
    savelegend5.AddEntry( "NULL" , "68% and 95% bands","")
    savelegend5.Draw()

    shadeBox = ROOT.TBox()

    # Legend in top right-hand corner
    boxX1 = xmin4 + (xmax4 - xmin4)*leftOfLegend + 0.18 #0.18
    boxX2 = boxX1 + 0.45 # 0.45
    if c.GetLogy() :
      boxY1 = math.exp(math.log(ymin) + (math.log(ymax) - math.log(ymin))*(topOfLegend-0.15)) #0.66
      boxY2 = math.exp(math.log(boxY1) + (math.log(ymax) - math.log(ymin))*0.0155)
      boxY3 = math.exp(math.log(boxY2) + (math.log(ymax) - math.log(ymin))*0.025)
      boxY4 = math.exp(math.log(boxY3) + (math.log(ymax) - math.log(ymin))*0.0155)
    else :
      boxY1 = ymin + (ymax - ymin)*(topOfLegend-0.06) #0.66
      boxY2 = boxY1 + (ymax - ymin)*0.0155
      boxY3 = boxY2 + (ymax - ymin)*0.025
      boxY4 = boxY3 + (ymax - ymin)*0.0155

    shadeBox.SetFillColor(self.colourpalette.twoSigmaBandColour)
    shadeBox.DrawBox(boxX1,boxY1,boxX2,boxY4)
    shadeBox.SetFillColor(self.colourpalette.oneSigmaBandColour)
    shadeBox.DrawBox(boxX1,boxY2,boxX2,boxY3)

    c.Update()

    pad1.cd()
    lumInFb = round(float(luminosity)/float(1000),nsigfigs)

    self.drawATLASLabels(0.50,0.88,False,False,0.06)
    self.drawCMEAndLumi(0.42,0.83,CME,lumInFb,0.05)
#    self.drawATLASLabels(0.18,0.20)
#    self.drawCMEAndLumi(0.1,0.26,CME,lumInFb,0.04)
    c.Update()

    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)

  def draw2DLimit(self,hist,name,xAxisName,xlow,xhigh,yAxisName,ylow,yhigh,zAxisName,luminosity=-1,CME=-1,doRectangular=False) :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,doRectangular,1.2)
    c.SetLogx(0)
    c.SetLogy(0)
    c.SetGridx(0)
    c.SetGridy(0)

    c.SetRightMargin(0.2)

    hist.Draw("colz")
    hist.GetZaxis().SetRangeUser(0,4)
    hist.GetZaxis().SetTitle(zAxisName)
    hist.GetZaxis().SetTitleOffset(1.40)
    hist.GetXaxis().SetRangeUser(xlow,xhigh)
    hist.GetXaxis().SetTitleOffset(1.40)
    hist.GetYaxis().SetRangeUser(ylow,yhigh)
    hist.GetYaxis().SetTitle(yAxisName)
    hist.GetXaxis().SetTitle(xAxisName)
    hist.GetYaxis().SetTitleOffset(1.40)
    hist.GetYaxis().SetLabelSize(0.075)
    hist.Draw("textsame")

    lumInFb = round(float(luminosity)/float(1000),nsigfigs)
    self.drawATLASLabels(0.17,0.88,False,True,0.05)
    #self.drawCMEAndLumi(0.08,0.82,CME,lumInFb,0.04)

    p1 = self.drawCME(0.165,0.81,CME,0.05)
    p2 = self.drawLumi(0.17,0.74,lumInFb,0.05)
    if self.dodrawUsersText :
      self.drawUsersText(0.165,0.695,self.cutstring,0.039)

    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)

  def drawOverlaid2DPlots(self,histBase,histsTop,name,xAxisName,xlow,xhigh,yAxisName,ylow,yhigh,zAxisName,luminosity=-1,CME=-1,doRectangular=False) :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,doRectangular,1.2)
    c.SetLogx(0)
    c.SetLogy(0)
    c.SetGridx(0)
    c.SetGridy(0)

    c.SetRightMargin(0.1)

    histBase.GetZaxis().SetRangeUser(0.0,2.0)
    histBase.GetZaxis().SetTitle(zAxisName)
    histBase.GetZaxis().SetTitleOffset(1.40)
    histBase.GetXaxis().SetRangeUser(xlow,xhigh)
    histBase.GetXaxis().SetTitleOffset(1.40)
    histBase.GetYaxis().SetRangeUser(ylow,yhigh)
    histBase.GetYaxis().SetTitle(yAxisName)
    histBase.GetXaxis().SetTitle(xAxisName)
    histBase.GetYaxis().SetTitleOffset(1.80)
    histBase.GetYaxis().SetLabelSize(0.075)
    histBase.Draw("colz")

#    hist

    for histContour in histsTop :
      histContour.SetLineColor(ROOT.kBlack)
      histContour.Draw("CONT1 SAME")

    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)

  def drawSeveralObservedAndExpected(self,observeds,expecteds1sigma,expecteds2sigma,legendnames,name,nameX,nameY,luminosity,CME,xmin,xmax,ymin,ymax,doRectangular=False) :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,doRectangular)
    c.SetLogx(0)
    c.SetLogy(1)
    c.SetGridx(0)
    c.SetGridy(0)

    # Set automatic axis range from graphs.
    # X axis range will be exactly ends of graphs
    xVals = []
    for i in range(observeds[0].GetN()) :
      xVals.append(observeds[0].GetX()[i])
    xVals.sort()
    if xmin == 'automatic':
      minX = xVals[-1]
    else :
      minX = xmin
    if xmax == 'automatic':
      maxX = xVals[0]
    else :
      maxX = xmax

    # Y axis range will be 3 orders of magnitude
    # above highest point of observed, because want space for legend
    # Lowest point should be 2 orders of magnitude below lowest point in observed
    if ymin == 'automatic' :
      minY = observeds[0].GetMinimum()/100
    else :
      minY = ymin
    if ymax == 'automatic' :
      maxY = observeds[0].GetMaximum()*1000
    else :
      maxY = ymax

    # Set axis names.
    # So far, should always be thus so don't pass as parameters.

    # Create legend
    leftOfLegend = 0.28
    bottomOfLegend = 0.68-0.05*(len(observeds)-1) # 0.64
    legend = self.makeLegend(leftOfLegend,bottomOfLegend,0.95,0.92) # 0.88

    # Set up display for expectations
    allerrorgraphs = expecteds1sigma+expecteds2sigma
    for graph in allerrorgraphs :
      graph.SetMarkerColor(1)
      graph.SetMarkerSize(1)
      graph.SetMarkerStyle(20)
      graph.SetLineColor(1)
      graph.SetLineWidth(3)
      graph.SetLineStyle(3)
      graph.GetXaxis().SetTitle(nameX)
      graph.GetYaxis().SetTitle(nameY)
      graph.GetXaxis().SetRangeUser(minX,maxX)
      graph.GetXaxis().SetNdivisions(705,ROOT.kTRUE)
      graph.GetYaxis().SetRangeUser(minY,maxY)
      if minX > 0.001 and minX < 1 :
        graph.GetXaxis().SetNoExponent(ROOT.kTRUE)

    for graph in expecteds1sigma :
      graph.SetFillColor(self.colourpalette.oneSigmaBandColour)
    for graph in expecteds2sigma :
      graph.SetFillColor(self.colourpalette.twoSigmaBandColour)

    # Set up display for observations
    for observed in observeds :
      index = observeds.index(observed)
      observed.SetMarkerColor(1)
      observed.SetMarkerSize(1)
#      observed.SetMarkerStyle(20)
      observed.SetMarkerStyle(24+index)
      observed.SetLineColor(1)
#      observed.SetLineWidth(3)
      observed.SetLineWidth(2)
      observed.SetLineStyle(1)
      observed.SetFillColor(0)
      observed.GetXaxis().SetTitle(nameX)
      observed.GetYaxis().SetTitle(nameY)
      observed.GetXaxis().SetRangeUser(minX,maxX)
      observed.GetYaxis().SetRangeUser(minY,maxY)
      if minX > 0.001 and minX < 1 :
        observed.GetXaxis().SetNoExponent(ROOT.kTRUE)

    # First one has to include axes or everything comes out blank
    # Rest have to NOT include axes or each successive one overwrites
    # previous. "SAME option does not exist for TGraph classes.
    expecteds1sigma[-1].Draw("A3")
    for graph in expecteds2sigma :
      graph.Draw("3") # no axis
    for graph in expecteds1sigma :
      graph.Draw("L3") # 1-sigma expectation error bands
#    for graph in expecteds1sigma :
#      graph.Draw("LX") # Center of expectation
    for observed in observeds :
      observed.Draw("PL") # Data points of measurement

    # Fill and draw legend
    for observed in observeds :
      index = observeds.index(observed)
      legend.AddEntry(observed,legendnames[index],"PL")#"L")
    legend.AddEntry(expecteds1sigma[0], "Expected 95% CL upper limits","L")
    legend.AddEntry( "NULL" , "68% and 95% bands","")

    self.drawATLASLabels(0.20,0.20)
    legend.Draw()

    shadeBox = ROOT.TBox()
    boxX1 = minX + (maxX - minX)*0.19 #21
    boxX2 = boxX1 + (maxX - minX)*0.125 # 135
    boxY1 = math.exp(math.log(minY) + (math.log(maxY) - math.log(minY))*(0.665-0.06*(len(observeds)-1))) #0.66
    boxY2 = math.exp(math.log(boxY1) + (math.log(maxY) - math.log(minY))*0.0155)
    boxY3 = math.exp(math.log(boxY2) + (math.log(maxY) - math.log(minY))*0.025)
    boxY4 = math.exp(math.log(boxY3) + (math.log(maxY) - math.log(minY))*0.0155)
    shadeBox.SetFillColor(self.colourpalette.twoSigmaBandColour)
    shadeBox.DrawBox(boxX1,boxY1,boxX2,boxY4)
    shadeBox.SetFillColor(self.colourpalette.oneSigmaBandColour)
    shadeBox.DrawBox(boxX1,boxY2,boxX2,boxY3)

    lumInFb = round(float(luminosity)/float(1000),nsigfigs)
    lumiloc = min(0.45,bottomOfLegend-0.12)
    self.drawLumiAndCMEVert(0.65,lumiloc,lumInFb,CME,0.04)

    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)

  def drawSignalOverlaidOnBkgPlot(self,bkgPlot,signalPlots,signalMasses,legendlist,lumi,CME,yname,name,firstBin=-1,lastBin=-1,doLogY=False,FixY=False,printCanvas=True,doRectangular=False) :

    if (printCanvas) :
      canvasname = name+'_cv'
      outputname = name+epsorpdf
      if saveCfile:
        Coutputname = name+'.C'
      if saveRfile:
        Routputname = name+'.root'
      if saveEfile:
        Eoutputname = name+'.eps'
      c = self.makeCanvas(canvasname,doRectangular)
      c.SetLogx(1)
      c.SetLogy(doLogY)
      c.SetGridx(0)
      c.SetGridy(0)

    # Use bin range within which bkgPlot has entries,
    # plus one empty on either side if available
    lowbin,highbin = self.getAxisRangeFromHist(bkgPlot)
    if (firstBin>0) :
      lowbin=firstBin
    if (lastBin>0 and lastBin>=firstBin) :
      highbin = lastBin

    # Create legend
    if (printCanvas) :
      lumInFb = round(float(lumi)/float(1000),nsigfigs)
      topOfLegend = 0.81
      widthOfRow = 0.05
      bottomOfLegend = topOfLegend-(widthOfRow * len(signalPlots))
      legend = self.makeLegend(0.2,bottomOfLegend,0.9,topOfLegend)

    # Calculate y axis range
    maxval = -10
    minval = 10
    for bin in range(lowbin,highbin) :
      if bkgPlot.GetBinContent(bin) + bkgPlot.GetBinError(bin) > maxval :
        maxval = bkgPlot.GetBinContent(bin) + bkgPlot.GetBinError(bin)
        if len(signalPlots) > 0 :
          maxval = maxval + signalPlots[0].GetBinContent(bin)
      if bkgPlot.GetBinContent(bin) - bkgPlot.GetBinError(bin) < minval :
        minval = bkgPlot.GetBinContent(bin) - bkgPlot.GetBinError(bin)
    locationzero = (0.0 - minval)/(maxval - minval)
    cutoff = 0.50
    if (printCanvas) :
      cutoff = 0.35
    if locationzero > cutoff :
      maxval = (0.0 - minval)/cutoff + minval

    # If part of bigger plot, want axis labels to be near to some value such that labels easy to see
    print 'FixY is',FixY
    if printCanvas :
      bkgPlot.GetYaxis().SetRangeUser(minval - 0.1,maxval+0.1)
    elif FixY :
      print "val is", bkgPlot.GetBinContent(bkgPlot.GetMinimumBin())
      if bkgPlot.GetBinContent(bkgPlot.GetMinimumBin()) >= -1.0 :
        bkgPlot.GetYaxis().SetRangeUser(-1.3,3.7)
      else :
        bkgPlot.GetYaxis().SetRangeUser(-3.7,3.7)
      print "setting narrow axis range"
    else :
      lowerint = math.floor(minval)
      if minval - lowerint > 0.7 :
        minval = lowerint+0.5
      else :
        minval = lowerint - 0.3
      upperint = math.ceil(maxval)
      if upperint - maxval > 0.7 :
        maxval = upperint - 0.5
      else :
        maxval = upperint + 0.3
      # Lydia example: To modify middle panel of Fancy Figure can comment out line below and replace with bkgPlot.GetYaxis().SetRangeUser(-4,4)
      bkgPlot.GetYaxis().SetRangeUser(minval,maxval)

    xname = "Reconstructed m_{jj} [TeV]"
    if (printCanvas) :
      if FixY :
        self.drawDataHist(bkgPlot,lowbin,highbin,xname,yname,False,1,False,True)
      else : self.drawDataHist(bkgPlot,lowbin,highbin,xname,yname,False,1,False)
    else :
      if FixY :
        self.drawDataHist(bkgPlot,lowbin,highbin,xname,yname,False,3,False,True)
      else :
        self.drawDataHist(bkgPlot,lowbin,highbin,xname,yname,False,3,False)
    if (printCanvas) :
      legend.AddEntry(bkgPlot,yname,"LP");

    goodcolours = self.getGoodColours(len(signalPlots))

    self.savegraphs = []
    for observed in signalPlots :
      index = signalPlots.index(observed)
      mass = signalMasses[index]
      newname = observed.GetName()+"_graph"
      plot = ROOT.TGraph()
      plot.SetName(newname)
      pointn = 0

      # Trim to only desired bins
      for bin in range(observed.GetNbinsX()+2) :
        if ((observed.GetBinLowEdge(bin)+observed.GetBinWidth(bin) < 0.68*mass) or\
           (observed.GetBinLowEdge(bin) > 1.85*mass)) :
          continue
        plot.SetPoint(pointn,observed.GetBinCenter(bin),observed.GetBinContent(bin))
        pointn = pointn+1

      plot.SetMarkerSize(1)
      plot.SetMarkerColor(goodcolours[index])
      plot.SetMarkerStyle(24+index) # was 20 when things were nice
      #plot.SetMarkerStyle(ROOT.kOpenCircle)

      plot.SetLineColor(goodcolours[index])
      plot.SetLineStyle(2);
      plot.GetXaxis().SetTitle("")
      plot.GetYaxis().SetTitle("")
      plot.SetTitle("")

      self.savegraphs.append(plot)

      if (printCanvas) :
        legend.AddEntry(self.savegraphs[index],legendlist[index],"LP");

      self.savegraphs[index].Draw("SAME CP")

    # Draw ratio plot once more to make sure it's on top
    bkgPlot.Draw("E SAME")

    if (printCanvas) :
      self.drawATLASLabels(0.22,0.88)
      self.drawCMEAndLumi(0.22,0.83,CME,lumInFb,0.04)
      legend.Draw()

      c.RedrawAxis()
      c.Update()
      c.SaveAs(outputname)
      if saveCfile:
        c.SaveSource(Coutputname)
      if saveRfile:
        c.SaveSource(Routputname)
      if saveEfile:
        c.SaveAs(Eoutputname)

  def drawSignalOverlaidOnDataAndFit(self,dataHist,fitHist,signalPlots,signalMasses,legendlist,lumi,CME,yname,name,firstBin=-1,lastBin=-1,doLogX=False,doLogY=True,printCanvas=True,doRectangular=False,nPads = 1,mcHist=None) :

    if (printCanvas) :
      canvasname = name+'_cv'
      outputname = name+epsorpdf
      if saveCfile:
        Coutputname = name+'.C'
      if saveRfile:
        Routputname = name+'.root'
      if saveEfile:
        Eoutputname = name+'.eps'
      c = self.makeCanvas(canvasname,doRectangular)
      c.SetLogx(doLogX)
      c.SetLogy(doLogY)

    drawMC = False
    if not mcHist==None :
      drawMC=True

    # Use bin range within which bkgPlot has entries,
    # plus one empty on either side if available
    lowbin,highbin = self.getAxisRangeFromHist(dataHist)
    if (firstBin>0) :
      lowbin=firstBin
    if (lastBin>0 and lastBin>=firstBin) :
      highbin = lastBin

    # Create legend
    if (printCanvas) :
      lumInFb = round(float(lumi)/float(1000),nsigfigs)
      if (doLogX) :
        topOfLegend = 0.43
        leftOfLegend = 0.18
      else :
        topOfLegend = 0.69
        leftOfLegend = 0.45
      widthOfRow = 0.05
      bottomOfLegend = topOfLegend-(widthOfRow*(len(signalPlots)+2))
      if drawMC :
        bottomOfLegend -= widthOfRow
      rightOfLegend = leftOfLegend+0.4
      legend = self.makeLegend(leftOfLegend,bottomOfLegend,rightOfLegend,topOfLegend)

    xname = "Reconstructed m_{jj} [TeV]"

    if drawMC :
      self.drawFitHist(mcHist,lowbin,highbin,xname,yname,False,True,False,[],False,ROOT.kGreen+2,1,False,drawMC)
      self.drawFitHist(fitHist,lowbin,highbin,"","",True,True,False,[],False,ROOT.kRed,1,False,drawMC)

    else :
      self.drawFitHist(fitHist,lowbin,highbin,xname,yname,False,True,False,[],False,ROOT.kRed,1,False,drawMC)

    goodcolours = self.getGoodColours(len(signalPlots))

    self.saveplots = []
    for observed in signalPlots :

      index = signalPlots.index(observed)

      mass = signalMasses[index]
      newname = observed.GetName()+"_graph_2"
      plot = ROOT.TGraph()
      plot.SetName(newname)
      pointn = 0

      # Trim to only desired bins
      for bin in range(observed.GetNbinsX()+2) :
        if ((observed.GetBinLowEdge(bin)+observed.GetBinWidth(bin) < 0.68*mass) or\
           (observed.GetBinLowEdge(bin) > 1.85*mass)) :
          continue
        plot.SetPoint(pointn,observed.GetBinCenter(bin),observed.GetBinContent(bin)+fitHist.GetBinContent(bin))
        pointn = pointn+1

      #plot.SetMarkerStyle(ROOT.kOpenCircle)
      plot.SetMarkerStyle(24+index) # was 20 when things were nice
      plot.SetMarkerSize(1)
      plot.SetMarkerColor(goodcolours[index])

      plot.SetMarkerColor(goodcolours[index])
      plot.SetLineColor(goodcolours[index])
      plot.SetLineStyle(2);
      plot.GetXaxis().SetTitle("")
      plot.GetYaxis().SetTitle("")
      plot.SetTitle("")

      self.saveplots.append(plot)
      if (printCanvas) :
        legend.AddEntry(self.saveplots[index],legendlist[index],"LP")

      self.saveplots[index].Draw("SAME CP")

    self.drawDataHist(dataHist,lowbin,highbin,"","",True,nPads,False,False,drawMC)
    if drawMC :
      self.drawFitHist(mcHist,lowbin,highbin,xname,yname,True,True,False,[],False,ROOT.kGreen+2,1,False,drawMC)
    self.drawFitHist(fitHist,lowbin,highbin,xname,yname,True,True,False,[],False,ROOT.kRed,1,False,drawMC)

    if (printCanvas) :
      legend.AddEntry(dataHist,"Data","LP")
      legend.AddEntry(fitHist,"Fit","LF")
      if drawMC :
        legend.AddEntry(mcHist,"QCD MC","LF")
      if (doLogX) :
        self.drawATLASLabels(0.57, 0.85)
      else :
        self.drawATLASLabels(0.52, 0.88)
      # Lydia Moved labels self.drawCMEAndLumi(leftOfLegend+0.02,topOfLegend+0.02,CME,lumInFb,0.04)
      self.drawLumiAndCMEVert(0.57,0.73,lumInFb,CME,0.04)
      legend.SetFillStyle(0)
      legend.Draw()

      c.RedrawAxis()
      c.Update()
      c.SaveAs(outputname)
      if saveCfile:
        c.SaveSource(Coutputname)
      if saveRfile:
        c.SaveSource(Routputname)
      if saveEfile:
        c.SaveAs(Eoutputname)

  def drawSeveralObservedLimits(self,observedlist,signallegendlist,name,nameX,nameY,luminosity,CME,xmin,xmax,ymin,ymax,extraLegendLines = [], doLogY=True,doLogX=False,doRectangular=False,doLegendLocation="Right",ATLASLabelLocation="BottomL",isTomBeingDumb=False,addHorizontalLines=[],pairNeighbouringLines=False,cutLocation="Right") :

    # LegendLocation should be "Right","Left", or "Wide"
    # ATLASLabelLocation should be "BottomL", "BottomR", "byLegend", or "None"
    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,doRectangular)
    c.SetLogx(doLogX)
    c.SetLogy(doLogY)
    c.SetGridx(0)
    c.SetGridy(0)

    # Set automatic axis range from graphs.
    # X axis range will be exactly ends of graphs
    xVals = []
    for list in observedlist :
      for i in range(list.GetN()) :
        xVals.append(list.GetX()[i])
    xVals.sort()
    if xmin == 'automatic':
      minX = xVals[0]
    else :
      minX = xmin
    if xmax == 'automatic':
      maxX = xVals[-1]
    else :
      maxX = xmax
    # Y axis range will be 2 orders of magnitude
    # above highest point of signal, because want space for legend
    # Lowest point should be 2 orders of magnitude below lowest point in observed
    minval = 1E10
    maxval = -1E10
    for graph in observedlist :
      testMax = ROOT.TMath.MaxElement(graph.GetN(), graph.GetY());
      testMin = ROOT.TMath.MinElement(graph.GetN(), graph.GetY());
      if testMin < minval :
        minval = testMin
      if testMax > maxval :
        maxval = testMax
    if ymin == 'automatic' :
      if doLogY :
        minY = minval/100
      else :
        if minval < 0 :
          minY = 1.2*minval
        else :
          minY = 0.8*minval
    else :
      minY = ymin
    if ymax == 'automatic' :
      if doLogY :
        maxY = maxval*100
      else :
        if maxval < 0 :
          maxY = 0.8*maxval
        else :
          maxY = 1.2*maxval
    else :
      maxY = ymax

    # Create legend
    if doLegendLocation=="Right" :
#      leftOfAll = 0.60
      leftOfAll = 0.55
    elif doLegendLocation=="Center" :
      leftOfAll = 0.35
    else :
       leftOfAll = 0.20
    topOfAll = 0.88
    widthOfRow = 0.05

    leftOfLegend = leftOfAll
    if doLegendLocation=="Left" :
      rightOfLegend = 0.50
    else :
      rightOfLegend = 0.90
    topOfLegend = topOfAll - (widthOfRow+0.05)*len(extraLegendLines)#-0.03# 0.75
    if cutLocation != "Left" :
      topOfLegend = topOfLegend - 0.04 # was 0.09
    else :
      topOfLegend = topOfLegend - 0.01
    bottomOfLegend = topOfLegend-(widthOfRow * len(observedlist)) 

    legend = self.makeLegend(leftOfLegend,bottomOfLegend,rightOfLegend,topOfLegend)

    # Set up display for expectations
    if pairNeighbouringLines :
      goodcolours = self.getGoodColours(len(observedlist)/2+1)
    else :
      goodcolours = self.getGoodColours(len(observedlist))

    for observed in observedlist :

      index = observedlist.index(observed)
      if pairNeighbouringLines :
        if len(observedlist) < 3 :
          colour = goodcolours[int(index/2.0)+1]
        else :
          colour = goodcolours[int(index/2.0)]
      else :
        if len(observedlist) < 3 :
          colour = goodcolours[index+1]
        else :
          colour = goodcolours[index]

      # Set up display for observations
      observed.SetMarkerColor(colour)
      if not isTomBeingDumb :
        observed.SetMarkerSize(0.7)  # was 0.5 back when things were nice
        observed.SetMarkerStyle(24+index) # was 20 when things were nice
      else :
        observed.SetMarkerSize(1.0)  # was 0.5 back when things were nice
        observed.SetMarkerStyle(20+index) # was 20 when things were nice
      observed.SetLineColor(colour)
      observed.SetLineWidth(2)
      observed.SetLineStyle(1)
      if pairNeighbouringLines :
        print "Pairing!"
        if index % 2 == 1 :
          observed.SetLineStyle(2)
      observed.SetFillColor(0)
      observed.GetXaxis().SetTitle(nameX)
      observed.GetYaxis().SetTitle(nameY)
      observed.GetXaxis().SetLimits(minX,maxX)
      observed.GetYaxis().SetRangeUser(minY,maxY)
      observed.GetXaxis().SetNdivisions(605,ROOT.kTRUE)

      # First one has to include axes or everything comes out blank
      # Rest have to NOT include axes or each successive one overwrites
      # previous. "SAME option does not exist for TGraph classes.
      if index==0 :
        observed.Draw("APL") # Data points of measurement
      else :
        observed.Draw("PL")

      # Fill and draw legend
      legend.AddEntry(observed,signallegendlist[index],"PL")

    if addHorizontalLines != [] :
      for val in addHorizontalLines :
        line = ROOT.TLine(minX, val, maxX, val)
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(2)
        line.Draw("SAME")

    if (ATLASLabelLocation=="BottomL") :
      self.drawATLASLabels(0.20,0.20)
    elif (ATLASLabelLocation=="byLegend") :
      self.drawATLASLabels(leftOfLegend,bottomOfLegend-0.08,True)
    elif (ATLASLabelLocation=="BottomR") :
      self.drawATLASLabels(0.53,0.20,True)

    persistent = []

    if luminosity != [] and CME != [] :
      lumInFb = round(float(luminosity)/float(1000),nsigfigs)
      # From when we had luminosity sign and L
#      persistent.append(self.drawCME(leftOfAll,topOfAll,CME,0.04))
#      c.Update()
#      persistent.append(self.drawLumi(leftOfAll,topOfAll-0.02-widthOfRow,lumInFb,0.04))
      persistent.append(self.drawCMEAndLumi(leftOfAll-0.08,topOfAll,CME,lumInFb,0.04))

    if self.dodrawUsersText :
      if not cutLocation == "Left" :
        self.drawUsersText(leftOfAll+0.01,topOfAll-0.05,self.cutstring,0.04)
      else :
        self.drawUsersText(0.20,topOfAll,self.cutstring,0.04)

    c.Update()

    self.myLatex.SetTextFont(42)
    self.myLatex.SetTextSize(0.04)
    index = 0
    if len(extraLegendLines) > 0 :
      toplocation = topOfAll - (0.03+widthOfRow) - (0.01+widthOfRow)*(index) # first one was 2*widthOfRow when had lumi and cme separately
      for line in extraLegendLines :
        toplocation = topOfAll - (0.03+widthOfRow) - (0.01+widthOfRow)*(index) # first one was 2*widthOfRow when had lumi and cme separately
        persistent.append(self.myLatex.DrawLatex(leftOfLegend+0.01,toplocation,line))
        index = index+1

    if doLegendLocation=="Wide" :
      if len(observedlist) > 5 :
        legend.SetNColumns(2)
    legend.Draw()

    # Should have draw-box for the bands
    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)

  # Emma's version ;)
  def drawSeveralObservedExpectedLimits(self,observedlist,expectedlist,expected1list,expected2list,signallegendlist,name,nameX,nameY,luminosity,CME,xmin,xmax,ymin,ymax,extraLegendLines = [], doLogY=True,doLogX=False,doRectangular=False,doLegendLocation="Left",ATLASLabelLocation="BottomL",addHorizontalLines=[],cutLocation="Right",labels=[]) :

    # LegendLocation should be "Right","Left", or "Wide"
    # ATLASLabelLocation should be "BottomL", "BottomR", "byLegend", or "None"
    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,doRectangular)
    c.SetLogx(doLogX)
    c.SetLogy(doLogY)
    c.SetGridx(0)
    c.SetGridy(0)

    doSeparatedLegend = True

    # Set automatic axis range from graphs.
    # X axis range will be exactly ends of graphs
    xVals = []
    for obs in observedlist :
      if isinstance(obs, (list, tuple)): 
        for g in obs:
          if g == None: continue
          for i in range(g.GetN()) :  xVals.append(g.GetX()[i])
      else:
        for i in range(obs.GetN()) : xVals.append(obs.GetX()[i])
    xVals.sort()
    minX = xVals[0]  if xmin == 'automatic' else xmin
    maxX = xVals[-1] if xmax == 'automatic' else xmax

    # Y axis range will be 2 orders of magnitude
    # above highest point of signal, because want space for legend
    # Lowest point should be 2 orders of magnitude below lowest point in observed

    (minval,maxval) = (1E10, -1E10)
    
    for e in observedlist :
      if isinstance(e, (list, tuple)):
        for g in e:
          if g == None: continue
          minval = min(minval, ROOT.TMath.MinElement(g.GetN(), g.GetY()))
          maxval = max(maxval, ROOT.TMath.MaxElement(g.GetN(), g.GetY()))
      else:
        minval = min(minval, ROOT.TMath.MinElement(e.GetN(), e.GetY()))
        maxval = max(maxval, ROOT.TMath.MaxElement(e.GetN(), e.GetY()))
      
    if ymin == 'automatic' :
      if doLogY: minY = minval/100
      else:      minY = 1.2*minval  if minval < 0 else 0.8*minval
    else :
      minY = ymin
    if ymax == 'automatic' :
      if doLogY: maxY = maxval*100
      else:      maxY = 0.8*maxval if maxval < 0 else 1.2*maxval
    else :
      maxY = ymax

    # Create legend
    if doLegendLocation=="Right" :
#      leftOfAll = 0.60
      leftOfAll = 0.55
    elif doLegendLocation=="Center" :
      leftOfAll = 0.55
    else :
       leftOfAll = 0.20 
    topOfAll = 0.88
    widthOfRow = 0.05

    leftOfLegend = leftOfAll #- (0.1 if len(observedlist) == 1 else 0)
    if doLegendLocation=="Left" :
      rightOfLegend = 0.50
    else :
      rightOfLegend = 0.90
    topOfLegend = topOfAll - (widthOfRow+0.05)*len(extraLegendLines)#-0.03# 0.75
    if cutLocation != "Left" :
      topOfLegend = topOfLegend - 0.04 # was 0.09 
    else :
      topOfLegend = topOfLegend - 0.01
    legendLength = self.cutstring.count(';') + (3 if len(observedlist) == 1 else len(observedlist)) 
    bottomOfLegend = topOfLegend-(widthOfRow * (legendLength)) + (-0.25 if len(observedlist) > 1 else 0)

    legend = self.makeLegend(leftOfLegend,bottomOfLegend,rightOfLegend,topOfLegend,fontSize = 0.035)
    if ';' in self.cutstring: [legend.AddEntry( "NULL" , "#kern[-0.8]{%s}"%substr, "") for substr in self.cutstring.split(';')]

    else: legend.AddEntry( "NULL" , self.cutstring, "")
    #legend->SetNColumns(len(observed[0]))

    # Set up display for expectations

    def GetContourColor(index,offset = 0):
      if index == 0: return ROOT.TColor.GetColor("#003EFF")
      if index == 1: return ROOT.kViolet+1
      if index == 2: return ROOT.kGreen+2
      goodcolours = self.getGoodColours(len(observedlist))
      return goodcolours[index]
      
    persistent = []
      
    def drawHelper(g,opt):
      if g == None: return
      if isinstance(g, (list, tuple)):
        [drawHelper(i,opt) for i in g if g != None]
      else:
        persistent.append(g)
        g.Draw(opt)

    for index, observed in enumerate(observedlist) :
      expected =  expectedlist[index]
      expected1 = expected1list[index]
      expected2 = expected2list[index]


      # Set up display for observations
      
      def formatObservedGraph(obs,i=0):
        if isinstance(obs, (list, tuple)): [formatObservedGraph(g,j) for j,g in enumerate(obs) if g != None]
        else:
          colour = GetContourColor(index,i)
          if len(observedlist) == 1: colour = ROOT.kBlack
          obs.SetMarkerColor(colour)
          obs.SetMarkerSize(0.7)  
          obs.SetMarkerStyle(24+index) 
          if len(observedlist) == 1:
            obs.SetMarkerSize(0.7)
            obs.SetMarkerStyle(20)

          obs.SetLineColor(colour)
          obs.SetLineWidth(2)
          obs.SetLineStyle(1)
          obs.SetFillColor(0)
          
          obs.GetXaxis().SetTitle(nameX)
          obs.GetYaxis().SetTitle(nameY)
          obs.GetXaxis().SetLimits(minX,maxX)
          obs.GetYaxis().SetRangeUser(minY,maxY)
          obs.GetXaxis().SetNdivisions(605,ROOT.kTRUE)
          
      def formatExpectedGraph(exp,i=0):
        if isinstance(exp, (list, tuple)): [formatExpectedGraph(g,j) for j,g in enumerate(exp) if g != None ]
        else:
          colour = GetContourColor(index,i)
          if len(observedlist) == 1: colour = ROOT.kBlack
          exp.SetMarkerColor(0)
          exp.SetLineColor(colour)
          exp.SetFillColor(0)
          exp.SetLineWidth(2)
          exp.SetLineStyle(ROOT.kDashed)
        
      def formatBands(band, bandcolorfactor, i=0):
        if isinstance(band, (list, tuple)): [formatBands(g, bandcolorfactor,j) for j,g in enumerate(band) if g != None]
        else:
          colour = GetContourColor(index,i)
          bandcolor = colorInterpolate(colour,ROOT.kWhite,bandcolorfactor)
          if len(observedlist) == 1:
            bandcolor = self.colourpalette.oneSigmaBandColour
            if bandcolorfactor > 0.5: bandcolor = self.colourpalette.twoSigmaBandColour
          #else: bandcolor = 0
          band.SetFillColor(bandcolor)
        
      formatObservedGraph(observed)
      formatExpectedGraph(expected)
      formatBands(expected1, 0.5)
      formatBands(expected2, 0.8)

      # First one has to include axes or everything comes out blank
      # Rest have to NOT include axes or each successive one overwrites
      # previous. "SAME option does not exist for TGraph classes.
      
      if index == 0:
        if isinstance(observed, (list, tuple)):
          for g in observed:
            if g == None: continue
            g.Draw( "APL" )
            break
        else:
          observed.Draw( "APL" )
      #expected.Draw("L")

      if isinstance(observed, (list, tuple)):
        for i,g in enumerate(observed):
          if g == None: continue
          if len(observedlist) == 1:
            legend.AddEntry(g,"Observed 95% CL Upper Limit","PL")
          break

      else:
        legend_str = "Observed 95% CL Upper Limit" if len(observedlist) == 1 else signallegendlist[index]
        legend.AddEntry(observed,legend_str,"PL")
        #legend.AddEntry(observed,signallegendlist[index],"PL")

    for i,g in enumerate(expectedlist):
      if g == None: continue
      if i == 0: legend.AddEntry("NULL", "#kern[-0.5]{Expected 95% CL Upper Limits for: }","")
      if isinstance(g, (list, tuple)):
        if g[0] == None: g = g[1]
        else: g = g[0]
      legend.AddEntry(g, signallegendlist[i]  + (" (#pm 1-2#sigma)" if i ==0 else ''),"L") 
    for i,g in enumerate(observedlist):
      if g == None: continue
      if i == 0: legend.AddEntry("NULL", "#kern[-0.5]{Observed 95% CL Upper Limits for: }","")
      if isinstance(g, (list, tuple)):
        if g[0] == None: g = g[1]
        else: g = g[0]
      legend.AddEntry(g, signallegendlist[i],"PL") 
      #legend.AddEntry(expected[i], signallegendlist[i] + " Expected 95% CL Upper Limit " + ("(#pm 1-2#sigma)" if i ==0 else ''),"L")

    if len(observedlist) == 1:
      if isinstance(expected, (list, tuple)):
        legend.AddEntry(expected[0], "Expected 95% CL Upper Limit (#pm 1-2#sigma)","L")
      else:
        legend.AddEntry(expected, "Expected 95% CL Upper Limit (#pm 1-2#sigma)","L")
      #legend.AddEntry( "NULL" , "68% and 95% bands","")

      fill_line = ROOT.TLine()
      fill_line.SetNDC()
      fill_line.SetLineColor(self.colourpalette.twoSigmaBandColour)
      fill_line.SetLineWidth(20)
      offset = (3-legendLength)*0.05
      fill_line.DrawLineNDC(leftOfAll+0.01, 0.745 + offset, leftOfAll + 0.065, 0.745+offset)
      fill_line.SetLineColor(self.colourpalette.oneSigmaBandColour)
      fill_line.SetLineWidth(10)
      fill_line.DrawLineNDC(leftOfAll+0.01, 0.745 + offset, leftOfAll + 0.065, 0.745+offset)
    else:
      fill_line = ROOT.TLine()
      fill_line.SetNDC()
      fill_line.SetLineColor(ROOT.TColor.GetColor("#D1DCFF"))
      fill_line.SetLineWidth(20)
      offset = -0.02
      fill_line.DrawLineNDC(leftOfAll+0.01, 0.745 + offset, leftOfAll + 0.075, 0.745+offset)
      fill_line.SetLineColor(ROOT.TColor.GetColor("#A2B9FF"))
      fill_line.SetLineWidth(10)
      fill_line.DrawLineNDC(leftOfAll+0.01, 0.745 + offset, leftOfAll + 0.075, 0.745+offset)



    for index, observed in enumerate(observedlist) :
      drawHelper(expected2list[index],"F")
      drawHelper(expected1list[index],"F")

    for index, observed in enumerate(observedlist) :
      expected = expectedlist[index]
      drawHelper(expected,"L")
      drawHelper(observed,"PL")

    if addHorizontalLines != [] :
      for val in addHorizontalLines :
        line = ROOT.TLine(minX, val, maxX, val)
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(2)
        line.Draw("SAME")

    if (ATLASLabelLocation=="BottomL") :
      self.drawATLASLabels(0.20,0.20)
    elif (ATLASLabelLocation=="byLegend") :
      self.drawATLASLabels(leftOfLegend-0.05,bottomOfLegend-0.08,True)
    elif (ATLASLabelLocation=="BottomR") :
      self.drawATLASLabels(0.53,0.20,True)


    if luminosity != [] and CME != []:
      if isinstance(luminosity, (list,tuple)):
        lumInFb = [round(float(l)/float(1000),nsigfigs) for l in luminosity]
      else: lumInFb = round(float(luminosity)/float(1000),nsigfigs)
      if luminosity > 0 or len(luminosity) > 1:
        persistent.append(self.drawCMEAndLumi(leftOfAll+(-0.25 if doLegendLocation=="Center" else 0),topOfAll,CME,lumInFb,0.04))
      else: persistent.append(self.drawCME(leftOfAll+(-0.25 if doLegendLocation=="Center" else 0),topOfAll,CME,0.04))

    c.Update()

    self.myLatex.SetTextFont(42)
    self.myLatex.SetTextSize(0.04)
    index = 0
    if len(extraLegendLines) > 0 :
      toplocation = topOfAll - (0.03+widthOfRow) - (0.01+widthOfRow)*(index) # first one was 2*widthOfRow when had lumi and cme separately
      for line in extraLegendLines :
        toplocation = topOfAll - (0.03+widthOfRow) - (0.01+widthOfRow)*(index) # first one was 2*widthOfRow when had lumi and cme separately
        persistent.append(self.myLatex.DrawLatex(leftOfLegend+0.01,toplocation,line))
        index = index+1

    if doLegendLocation=="Wide" :
      if len(observedlist) > 5 :
        legend.SetNColumns(2)
    legend.Draw()

    # Should have draw-box for the bands
    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    c.SaveAs(name + '.png')
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)

  def drawPosteriorsWithCLs(self,posteriorsandclslist,legendlist,luminosity,CME,name,align=2,central=True,drawAsHist=False,addlinestolegend=False,doLogY=False,drawPriors=[],doRectangular=False,newxname="",newyname="") :

    canvasname = name+'_cv'
    outputname = name+epsorpdf
    if saveCfile:
      Coutputname = name+'.C'
    if saveRfile:
      Routputname = name+'.root'
    if saveEfile:
      Eoutputname = name+'.eps'
    c = self.makeCanvas(canvasname,doRectangular)
    c.SetLogx(0)
    c.SetLogy(doLogY)
    c.SetGridx(0)
    c.SetGridy(0)

    # Set axis names.
    # So far, should always be thus so don't pass as parameters.
    if newxname == "" :
      nameX = "Number of signal events"
    else :
      nameX = newxname
    if newyname == "" :
      nameY = "p(signal | data)"
    else :
      nameY = newyname

    # Create legend
    ncomp = len(legendlist)
    y1forleg = 0.65 - 0.1*(ncomp-6)
    if y1forleg < 0 :
      y1forleg = 0

    legendsize = 0.04*(len(legendlist)+len(drawPriors))
    if addlinestolegend:
      legendsize = 2*legendsize
    if align == 0:
      leftOfLegend = 0.21
      rightOfLegend = 0.65
    else :
      leftOfLegend = 0.51
      rightOfLegend = 0.95
    if align == 1 :
      legendtop = 0.83
    else :
      legendtop = 0.76
    legendbottom = legendtop - legendsize
    legend = self.makeLegend(leftOfLegend,legendbottom,rightOfLegend,legendtop)

    # Set up display for expectations
    goodcolours = self.getGoodColours(len(posteriorsandclslist)+len(drawPriors))

    for observed in posteriorsandclslist :

      index = posteriorsandclslist.index(observed)
      colour = goodcolours[index]

      posterior = observed[0]
      cl = observed[1]

      # Set up display for observations
      posterior.SetMarkerColor(colour)
      posterior.SetMarkerSize(0.5)
      posterior.SetMarkerStyle(20)
      posterior.SetLineColor(colour)
      posterior.SetLineWidth(3)
      posterior.SetLineStyle(1)
      posterior.SetFillColor(0)
      posterior.GetXaxis().SetTitle(nameX)
      posterior.GetYaxis().SetTitle(nameY)

      if len(posteriorsandclslist) + len(drawPriors) < 3 :
        posterior.SetMarkerColor(1000)
        posterior.SetLineColor(1000)

      # Standard draw option:
      if index==0 :
        if posterior.GetMaximum() < 0.5 or posterior.GetMaximum > 1E4 :
          c.SetLeftMargin(0.175)
          posterior.GetYaxis().SetTitleOffset(1.8)
        if central == True:
          posterior.GetYaxis().SetRangeUser(0,1.8*posterior.GetMaximum())
        posterior.Draw("C")
      else :
        posterior.Draw("C SAME")

      # Draw line for CL in matching colour :
      clLine = ROOT.TLine()
      clLine.SetLineColor(colour)
      if len(posteriorsandclslist) + len(drawPriors) < 3 :
       clLine.SetLineColor(1000)
      lineMinY = 0
      lineMaxY = posterior.GetMaximum()/2
      if type(cl) is float :
        clLine.DrawLine(cl,lineMinY,cl,lineMaxY)
      else :
        clLine.DrawLine(cl[0],lineMinY,cl[0],lineMaxY)

      # Fill and draw legend
      legend.AddEntry(posterior,legendlist[index],"P")
      if addlinestolegend :
        legend.AddEntry(clLine,"95% quantile","L")

    colourindex = len(posteriorsandclslist)
    for prior in drawPriors :
      prior.SetLineWidth(3)
      prior.SetLineStyle(1)
      colour = goodcolours[colourindex]
      if len(posteriorsandclslist) + len(drawPriors) == 2 :
        colour = 1002
      prior.SetLineColor(colour)
      prior.SetFillColor(0)
      colourindex = colourindex+1
      legend.AddEntry(prior,"Prior","L")
      prior.Draw("C SAME")

    legend.Draw()
    lumInFb = round(float(luminosity)/float(1000),nsigfigs)
    if align == 0 :
      self.drawATLASLabels(0.2, 0.87)
      self.drawCMEAndLumi(leftOfLegend+0.01,0.805,CME,lumInFb,0.04)
    elif align == 2 :
      self.drawATLASLabels(0.53, 0.87, True)
      self.drawCMEAndLumi(leftOfLegend+0.01,0.805,CME,lumInFb,0.04)
    else :
      self.drawATLASLabels(0.2,0.2)
      self.drawCMEAndLumi(leftOfLegend+0.01,0.87,CME,lumInFb,0.04)

    # Should have draw-box for the bands
    c.RedrawAxis()
    c.Update()
    c.SaveAs(outputname)
    if saveCfile:
      c.SaveSource(Coutputname)
    if saveRfile:
      c.SaveSource(Routputname)
    if saveEfile:
      c.SaveAs(Eoutputname)


  ## ----------------------------------------------------
  ## Internal functions

  def makeCanvas(self,canvasname,doRectangular,scaleX=1.0,scaleY=1.0) :

    if doRectangular :
      dim = int(800*scaleX),int(600*scaleY)
    else :
      dim = int(600*scaleX),int(600*scaleY)
    return ROOT.TCanvas(canvasname,'',0,0,dim[0],dim[1])

  def makeLegend(self,legX1,legY1,legX2,legY2,fontSize = 0.04) :

    legend = ROOT.TLegend(legX1,legY1,legX2,legY2)
    legend.SetTextFont(42)
    legend.SetTextSize(fontSize)
    legend.SetBorderSize(0)
    legend.SetLineColor(0)
    legend.SetLineStyle(1)
    legend.SetLineWidth(1)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)#1001)
    return legend

  def getAxisRangeFromHist(self,hist) :
    # Axis range should be decided by data hist
    firstBin =0
    while (hist.GetBinContent(firstBin+1)==0 and firstBin < hist.GetNbinsX()) :
      firstBin+=1
    lastBin = hist.GetNbinsX()+1
    while (hist.GetBinContent(lastBin-1)==0 and lastBin > 0) :
      lastBin-=1
    if (firstBin > lastBin) :
      firstBin=1
      lastBin = hist.GetNbinsX()
    return firstBin,lastBin

  def getYRangeFromHist(self,hist) :
    lowyval = 1E10
    lownonzero = 1E10
    highyval = -1E10
    lowbin,highbin = self.getAxisRangeFromHist(hist)
    for bin in range(lowbin,highbin) :
      if hist.GetBinContent(bin) < lowyval :
        lowyval = hist.GetBinContent(bin)
      if hist.GetBinContent(bin) < lownonzero and hist.GetBinContent(bin) > 0 :
        lownonzero = hist.GetBinContent(bin)
      if hist.GetBinContent(bin) > highyval :
        highyval = hist.GetBinContent(bin)
    return lowyval,lownonzero,highyval

  def getGoodColours(self, ncolours) :
    if ncolours < 4 :
      return self.colourpalette.shortGoodColours
    elif ncolours < 6 :
      return self.colourpalette.defaultGoodColours
    elif ncolours < 13 :
      return self.colourpalette.mediumGoodColours
    else :
      return self.colourpalette.longGoodColours


  def drawDataHist(self, dataHist,firstBin,lastBin,xname,yname,same=False,nPads=1,FixYAxis=False,LeaveAxisAlone=False,extraRoom=False) :

    # Data hist must be in data points with weighted error bars
    dataHist.SetMarkerColor(ROOT.kBlack)
    dataHist.SetLineColor(ROOT.kBlack)
    dataHist.GetXaxis().SetTitle(xname)
    dataHist.GetYaxis().SetTitle(yname)
    dataHist.GetXaxis().SetRange(firstBin,lastBin)
    dataHist.SetMarkerSize(0.75) # Lydia EOYE, was 1.2 before!
    
    actualMin = 1E10
    actualMax = 0
    for bin in range(firstBin,lastBin+1) :
      if dataHist.GetBinContent(bin) > actualMax : actualMax = dataHist.GetBinContent(bin)
      elif dataHist.GetBinContent(bin) < actualMin : actualMin = dataHist.GetBinContent(bin)
    
    if not LeaveAxisAlone :
      if not ROOT.gPad.GetLogy() :
        y1  = 0
        y2 = actualMax*1.2
      else :
        y1 = max(0.3, actualMin/5.0)
        y2 = actualMax*5.0
      if FixYAxis :
        y1 = max(0.5, actualMin/5.0)
        y2 = dataHist.GetBinContent(firstBin)*5
      # moving
      if extraRoom :
        y2 = y2 * 3
      dataHist.GetYaxis().SetRangeUser(y1,y2)
    else :
      # Get max and min within the available bins.
      y1 = dataHist.GetMinimum()
      y2 = dataHist.GetMaximum()

    if nPads==2 :
      dataHist.GetYaxis().SetNdivisions(805,ROOT.kTRUE)
    elif nPads==1 :
      dataHist.GetYaxis().SetNdivisions(605,ROOT.kTRUE)
    elif nPads==3 :
      dataHist.GetYaxis().SetNdivisions(605,ROOT.kTRUE)
    dataHist.GetXaxis().SetMoreLogLabels(ROOT.kTRUE)
    if same :
      dataHist.Draw("E SAME")
    else :
      dataHist.Draw("E")
      self.fixTheBloodyTickMarks(ROOT.gPad, dataHist, dataHist.GetBinLowEdge(firstBin), dataHist.GetBinLowEdge(lastBin+1),y1,y2)

  def drawFitHist(self,fitHist,firstBin,lastBin,xname,yname,same=False,twoPads=False,useError=False,errors = [],\
        drawCurve=False, lineColor = -1, lineStyle = 1, doEndLines = False,extraRoom=False) :

    if lineColor == -1 :
      lineColor = self.colourpalette.signalLineColours[0]
    # Fit hist must be expressed as a histo in a red line
    if (useError) :
      if errors == [] :
        fitHist.SetLineColor(self.colourpalette.signalLineColours[0])
        fitHist.SetFillColor(self.colourpalette.signalErrorColours[0])
        fitHist.SetFillStyle(1001)
        fitHist.SetMarkerStyle(20)
        fitHist.SetMarkerSize(0.0)
      else :
        fitHist.SetLineColor(self.colourpalette.signalLineColours[0])
        fitHist.SetLineStyle(lineStyle)
        fitHist.SetFillStyle(0)
        goodcolours = self.getGoodColours(len(errors))
        for errhistpair in errors :
          for thishist in errhistpair :
            thishist.SetLineColor(goodcolours[errors.index(errhistpair)])
            thishist.SetFillStyle(0)
            if errors.index(errhistpair) == 0:
              thishist.SetLineStyle(9)
            elif errors.index(errhistpair) == 1:
              thishist.SetLineStyle(2)
            else:
              thishist.SetLineStyle(lineStyle)
    else :
      fitHist.SetLineColor(lineColor)
      fitHist.SetLineStyle(lineStyle)
      fitHist.SetFillStyle(0)
    fitHist.SetLineWidth(2)
    fitHist.SetTitle("")
    fitHist.GetXaxis().SetTitle(xname)
    fitHist.GetYaxis().SetTitle(yname)
    fitHist.GetXaxis().SetRange(firstBin,lastBin)
    if not ROOT.gPad.GetLogy() :
      y1  = 0
      y2 = fitHist.GetBinContent(fitHist.GetMaximumBin())*1.2
      if extraRoom :
        y2 = y2*1.2
    else :
      y1 = 0.3
      y2 = fitHist.GetBinContent(fitHist.GetMaximumBin())*5
      if extraRoom :
        y2 = y2 * 3
    fitHist.GetYaxis().SetRangeUser(y1,y2)

    for errhistpair in errors:
      for thishist in errhistpair :
        thishist.GetXaxis().SetRange(firstBin,lastBin)
        thishist.GetYaxis().SetRangeUser(y1,y2)
        thishist.GetXaxis().SetTitle(xname)
        thishist.GetYaxis().SetTitleSize(0.06)
        thishist.GetYaxis().SetTitleOffset(1.1) # 1.2 = 20% larger
        thishist.GetYaxis().SetTitle(yname)
        thishist.GetXaxis().SetMoreLogLabels()
    if twoPads :
      fitHist.GetXaxis().SetNdivisions(805,ROOT.kTRUE)
    else :
      fitHist.GetXaxis().SetNdivisions(605,ROOT.kTRUE)

    drawOption = "HIST ]["
    if useError and errors == [] :
      drawOption = "][ E3"
    if same :
      drawOption = drawOption + " SAME"

    if errors != [] :
      count = 0
      for histpair in errors :
        for hist in histpair :
          # Lydia EOYE
          hist.GetYaxis().SetTitleSize(0.06)
          hist.GetYaxis().SetTitleOffset(0.8)
          hist.GetYaxis().SetLabelSize(0.05)
          drawOption = "HIST ][ SAME" # was L
          if drawCurve :
            drawOption = drawOption+" L"
          hist.Draw(drawOption)
          count = count+1

    if (doEndLines) :
      drawOption = drawOption.replace("]["," ")

    if drawCurve :
      drawOption = drawOption+" L"
    fitHist.Draw(drawOption)

    if not same :
      self.fixTheBloodyTickMarks(ROOT.gPad, fitHist, fitHist.GetBinLowEdge(firstBin), fitHist.GetBinLowEdge(lastBin+1),y1,y2)


    # To get line on top need to redraw
    persistent = fitHist.Clone("newone")
    persistent.SetFillStyle(0)
    if useError :
      persistent.SetLineColor(lineColor)
      persistent.Draw("HIST L SAME") # was L

    return persistent

  def setStandardTwoPads(self) :
  
    # Dimensions: xlow, ylow, xup, yup
    outpad = ROOT.TPad("extpad","extpad",0,0,1,1) # For marking outermost dimensions
    pad1 = ROOT.TPad("pad1","pad1",0,0.27,1,1) # For main histo
    pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.27) # For residuals histo

    # Set up to draw in right orientations
    outpad.SetFillStyle(4000) #transparent
    pad1.SetBottomMargin(0.00001)
    pad1.SetBorderMode(0)
    pad1.SetLogy(1)
    pad2.SetTopMargin(0.00001)
    pad2.SetBottomMargin(0.35) # 0.3
    pad2.SetBorderMode(0)
    pad1.Draw()
    pad2.Draw()
    outpad.Draw()
    return outpad,pad1,pad2

  def drawSignificanceHist(self,significance,firstBin,lastBin,xname,yname,fixYAxis=False,\
        inLargerPlot=False,doLogX=False,doErrors=False,fillColour = ROOT.kRed, drawAsLine=False, drawSame=False, yRange=[]) :

    if drawAsLine :
      significance.SetLineColor(fillColour)
      significance.SetFillStyle(0)
      # For "draw as line" we also want a simple dotted line at 0
      # which should be at the back.
      if not drawSame :
        line = ROOT.TLine(significance.GetBinLowEdge(firstBin), 0.0, significance.GetBinLowEdge(lastBin), 0.0)
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(2)
        line.Draw("SAME")
    else :
      significance.SetLineColor(ROOT.kBlack)
      significance.SetLineWidth(2)
      significance.SetFillColor(fillColour)
      significance.SetFillStyle(1001)

    significance.GetXaxis().SetTitle(xname)
    
    significance.GetYaxis().SetTitleSize(0.12) # 0.1
    significance.GetYaxis().SetTitleOffset(0.45) # 0.42
    significance.GetYaxis().SetLabelSize(0.1) # 0.1
    significance.GetXaxis().SetLabelSize(0.1) # 0.1
    significance.GetXaxis().SetTitleSize(0.12) # 0.1
    significance.GetXaxis().SetTitleOffset(1.3) # 1.2

    significance.GetXaxis().SetRange(firstBin,lastBin)
    lowPoint = significance.GetMaximum(significance.GetMaximumBin())
    highPoint = significance.GetMinimum()
    ylow = 0.0
    yhigh = 0.0
    for bin in range(firstBin,lastBin+1) :
      val = significance.GetBinContent(bin)
      if val < lowPoint :
        lowPoint = val
      if val > highPoint :
        highPoint = val
    if highPoint == 20 :
      highPoint = 7
    if fixYAxis==False :
      if lowPoint < 0 :
        ylow = lowPoint*1.2
        yhigh = highPoint*(1.2)
      else :
        ylow = lowPoint - 0.9*(highPoint - lowPoint)
        yhigh = highPoint + 0.9*(highPoint - lowPoint)
    else :
      if yRange :
        ylow = yRange[0]
        yhigh = yRange[1]
      else :
        if abs(significance.GetMaximum()) < 1.5 :
          ylow = -1.7
          yhigh = 1.7
        else :
          ylow = -3.7
          yhigh = 3.7
    significance.GetYaxis().SetRangeUser(ylow,yhigh)
    if inLargerPlot :
      significance.GetYaxis().SetTickLength(0.055)
    else :
      significance.GetYaxis().SetTickLength(0.035)
    significance.GetXaxis().SetTickLength(0.08)
    significance.GetXaxis().SetNdivisions(805,ROOT.kTRUE)
    significance.GetYaxis().SetNdivisions(805,ROOT.kTRUE)

    significance.GetYaxis().SetTitle(yname)
    significance.GetXaxis().SetTitle(xname)
    drawOption = ""
    if doErrors :
      drawOption = "E"
    else :
      if drawAsLine :
        drawOption = "L"
      else :
        drawOption = "HIST"
    if drawSame :
      drawOption = drawOption + " SAME"
    significance.Draw(drawOption)

    self.fixTheBloodyTickMarks(ROOT.gPad, significance, significance.GetBinLowEdge(firstBin), significance.GetBinLowEdge(lastBin+1),ylow,yhigh)
    
  def drawSignificanceHistWithJESBands(self,significance,upsignificance,downsignificance,firstBin,lastBin,xname,yname,fixYAxis=False,\
        inLargerPlot=False,doLogX=False,doErrors=False) :
    #significance.SetMarkerColor(ROOT.kBlack)
    significance.SetTitle("")

    significance.GetXaxis().SetRange(firstBin,lastBin)
    upsignificance.GetXaxis().SetRange(firstBin,lastBin)
    downsignificance.GetXaxis().SetRange(firstBin,lastBin)
    ylow = -1.2
    yhigh = 1.2
    significance.GetYaxis().SetRangeUser(ylow,yhigh)
    upsignificance.GetYaxis().SetRangeUser(ylow,yhigh)
    downsignificance.GetYaxis().SetRangeUser(ylow,yhigh)
    #significance.GetYaxis().SetNdivisions(5,10,0)
    #upsignificance.GetYaxis().SetNdivisions(5,10,0)
    #downsignificance.GetYaxis().SetNdivisions(5,10,0)
    significance.GetYaxis().SetNdivisions(604)
    upsignificance.GetYaxis().SetNdivisions(604)
    downsignificance.GetYaxis().SetNdivisions(604)

    if inLargerPlot :
      significance.GetYaxis().SetTickLength(0.055)
    significance.GetXaxis().SetNdivisions(805,ROOT.kTRUE)

    significance.GetYaxis().SetTitle(yname)
    significance.GetXaxis().SetTitle(xname)
    upsignificance.GetYaxis().SetTitle(yname)
    upsignificance.GetXaxis().SetTitle(xname)
    downsignificance.GetYaxis().SetTitle(yname)
    downsignificance.GetXaxis().SetTitle(xname)
    self.persistentlegend = self.makeLegend(0.135,0.78,0.32,0.9,0.12) # 0.88
    self.persistentlegend.AddEntry(upsignificance,"JES Uncertainty","F")

    if doErrors : # FIXME ok?
      upsignificance.Draw("EFHISTSAME")
      significance.Draw("ESAME")
      downsignificance.Draw("EFHISTSAME")
    else :
      #ratioHists[i].SetStats(0)
      #if i == 0: ratioHists[i].DrawCopy("p") #e0
      #if i == 0: ratioHists[i].SetMarkerSize(0) #e0
      #if i == 0: ratioHists[i].DrawCopy("same e0") #e0
      #else: ratioHists[i].Draw( "fhistsame" )
      significance.SetStats(0)
      significance.DrawCopy("p") #e0
      significance.SetMarkerSize(0) #e0
      significance.DrawCopy("same e0") #e0
      #significance.Draw("SAME")
      upsignificance.Draw("FHISTSAME")
      downsignificance.Draw("FHISTSAME")

    significance.GetXaxis().SetMoreLogLabels(ROOT.kTRUE)
    upsignificance.GetXaxis().SetMoreLogLabels(ROOT.kTRUE)
    downsignificance.GetXaxis().SetMoreLogLabels(ROOT.kTRUE)

    self.fixTheBloodyTickMarks(ROOT.gPad, significance, significance.GetBinLowEdge(firstBin), significance.GetBinLowEdge(lastBin+1),ylow,yhigh)
    self.fixTheBloodyTickMarks(ROOT.gPad, upsignificance, significance.GetBinLowEdge(firstBin), significance.GetBinLowEdge(lastBin+1),ylow,yhigh)
    self.fixTheBloodyTickMarks(ROOT.gPad, downsignificance, significance.GetBinLowEdge(firstBin), significance.GetBinLowEdge(lastBin+1),ylow,yhigh)

  def drawATLASLabels(self,xstart,ystart,rightalign=False,isRectangular=False,fontSize=0.05) :
    if self.labeltype < 0 :
      return
    # If we have set "rightalign" = true, we will take the x location as the right side of the label.
    self.myLatex.SetTextSize(fontSize)
    string = "#font[72]{0} #font[42]{1}".format("{ATLAS}","{"+self.labelDict[self.labeltype]+"}")
    if rightalign :
      self.myLatex.SetTextAlign(31)
    self.myLatex.DrawLatex(xstart, ystart, string)
    self.myLatex.SetTextAlign(11)

  # Text height required is around 0.04 for size 0.04, 0.05 for size 0.05, etc

  def drawLumi(self,xstart,ystart,lumiInFb,fontsize=0.05) :
    self.whitebox.Clear()
    self.whitebox.SetTextSize(fontsize)
    self.whitebox.SetX1NDC(xstart-0.01)
    self.whitebox.SetY1NDC(ystart-0.01)
    self.whitebox.SetX2NDC(xstart+0.25)
    self.whitebox.SetY2NDC(ystart+0.06)
    self.whitebox.SetTextFont(42)
    mystring = "#scale[0.7]{#int}L dt"
    myfb = "fb^{-1}"
    mypb = "pb^{-1}"
    #self.whitebox.AddText(0.04,1.0/8.0,"{0} = {1} {2}".format(mystring,lumiInFb,myfb))
    if doLumiInPb:
      self.whitebox.AddText(0.04,1.0/8.0,"{0} {1}".format(int(lumiInFb*1000),mypb))
    else:
      self.whitebox.AddText(0.04,1.0/8.0,"{0} {1}".format(lumiInFb,myfb))
    persistent = self.whitebox.DrawClone()
    return persistent

  def drawCME(self,xstart,ystart,CME,fontsize=0.05) :
    self.whitebox.Clear()
    self.whitebox.SetTextSize(fontsize)
    self.whitebox.SetX1NDC(xstart-0.01)
    self.whitebox.SetY1NDC(ystart-0.01)
    self.whitebox.SetX2NDC(xstart+0.25)
    self.whitebox.SetY2NDC(ystart+0.05)
    mysqrt = "#sqrt{s}"
    self.whitebox.AddText(0.04,1.0/6.0,"{0}={1} TeV".format(mysqrt,CME))
    persistent = self.whitebox.DrawClone()
    return persistent

  def drawLumiAndCMEVert(self,xstart,ystart,lumiInFb,CME,fontsize=0.05) :
    if lumiInFb < 0 and CME < 0 :
      return
    self.whitebox.Clear()
    self.whitebox.SetFillColor(0)
    #self.whitebox.SetFillStyle(1001)
    self.whitebox.SetFillStyle(3000) # make box transparent
    self.whitebox.SetTextColor(ROOT.kBlack)
    self.whitebox.SetTextFont(42)
    self.whitebox.SetTextAlign(11)
    self.whitebox.SetBorderSize(0)
    self.whitebox.SetTextSize(fontsize)
    self.whitebox.SetX1NDC(xstart-0.02)
    self.whitebox.SetY1NDC(ystart-0.01)
    self.whitebox.SetX2NDC(xstart+0.25)
    self.whitebox.SetY2NDC(ystart+0.11)
    mystring = "#scale[0.7]{#int}L dt"
    myfb = "fb^{-1}"
    mypb = "pb^{-1}"
#    inputstring1 = "{0} = {1} {2}".format(mystring,lumiInFb,myfb)
    if doLumiInPb:
      inputstring1 = "{0} {1}".format(int(lumiInFb*1000),mypb)
    else:
      inputstring1 = "{0} {1}".format(lumiInFb,myfb)
    #self.whitebox.AddText(0.0,0.55,inputstring1)
    self.whitebox.AddText(0.0,0.4,inputstring1)
    mysqrt = "#sqrt{s}"
    inputstring2 = "{0}={1} TeV".format(mysqrt,CME)
    self.whitebox.AddText(0.0,0.7,inputstring2)
    self.whitebox.Draw()
    return

  def drawCMEAndLumi(self,xstart,ystart,CME,lumiInFb,fontsize=0.05) :
    if lumiInFb < 0 and CME < 0 :
      return
    # Readjust x location after removing integral and lumi sign
    xstart = xstart + 0.08
    '''self.whitebox.Clear()
    self.whitebox.SetFillStyle(3000) # make box transparent
    self.whitebox.SetTextSize(fontsize)
    self.whitebox.SetX1NDC(xstart-0.01)
    self.whitebox.SetY1NDC(ystart-0.01)
    self.whitebox.SetX2NDC(xstart+0.32) # was 0.4 before
    self.whitebox.SetY2NDC(ystart+0.05)'''
    mysqrt = "#sqrt{s}"
    mystring = "#scale[0.7]{#int}L dt"
    myfb = "fb^{-1}"
    mypb = "pb^{-1}"

    if doLumiInPb:
      mytext = "{0}={1} TeV, {2} {3}".format(mysqrt,CME,int(lumiInFb*1000),mypb)
      #self.whitebox.AddText(0.04,1.0/8.0,"{0}={1} TeV, {2} {3}".format(mysqrt,CME,int(lumiInFb*1000),mypb))
    else:
      #lumitext = ''
      if isinstance(lumiInFb, (list,tuple)):
        lumitext = '-'.join(["{0}".format(l) for l in lumiInFb])
      else: lumitext = lumiInFb
      mytext = "{0}={1} TeV, {2} {3}".format(mysqrt,CME,lumitext,myfb)
      #self.whitebox.AddText(0.04,1.0/8.0,"{0}={1} TeV, {2} {3}".format(mysqrt,CME,lumiInFb,myfb))
    #self.whitebox.Draw()
    newtext = ROOT.TLatex()
    newtext.SetNDC()
    newtext.SetTextSize(fontsize)
    newtext.SetTextFont(42)
    newtext.SetTextAlign(11)
    newtext.DrawLatex(xstart,ystart,"{0}".format(mytext)) 
    return

  def drawXAxisInTeV(self,xmin,xmax,ymin,ymax,ndiv=510) :
    axisfunc = ROOT.TF1("axisfunc","x",float(xmin)/1000.0,float(xmax)/1000.0)
    newaxis = ROOT.TGaxis(xmin,ymin,xmax,ymax,"axisfunc",ndiv)
    return newaxis,axisfunc

  # Generic function to add text to plot e.g. to write a value on it
  def drawUsersText(self,xstart,ystart,text,fontsize=0.06) :

    newtext = ROOT.TLatex()
    newtext.SetNDC()
    newtext.SetTextSize(fontsize)
    newtext.SetTextFont(42)
    newtext.SetTextAlign(11)
    newtext.DrawLatex(xstart,ystart,"{0}".format(text)) 
    return
  
  def drawUsersTextLeftAligned(self,xstart,ystart,text,fontsize=0.06) :
    self.newwhitebox = ROOT.TPaveText()
    self.newwhitebox.SetFillColor(0)
    #self.newwhitebox.SetFillStyle(1001)
    self.newwhitebox.SetFillStyle(3000) # make box transparent
    self.newwhitebox.SetTextColor(ROOT.kBlack)
    self.newwhitebox.SetTextFont(42)
    self.newwhitebox.SetTextAlign(12)
    self.newwhitebox.SetBorderSize(0)

    self.newwhitebox.SetTextSize(fontsize)
    self.newwhitebox.SetX1NDC(xstart-0.01)
    self.newwhitebox.SetY1NDC(ystart-0.01)
    self.newwhitebox.SetX2NDC(xstart+0.25)
    self.newwhitebox.SetY2NDC(ystart+0.15)
    for textitem in text:
      self.newwhitebox.AddText("{0}".format(textitem))
    self.newwhitebox.Draw()
    return

  def drawUsersTextRightAligned(self,xstart,ystart,text,fontsize=0.06) :
    self.newwhitebox2 = ROOT.TPaveText()
    self.newwhitebox2.SetFillColor(0)
    #self.newwhitebox.SetFillStyle(1001)
    self.newwhitebox2.SetFillStyle(3000) # make box transparent
    self.newwhitebox2.SetTextColor(ROOT.kBlack)
    self.newwhitebox2.SetTextFont(42)
    self.newwhitebox2.SetTextAlign(32)
    self.newwhitebox2.SetBorderSize(0)

    self.newwhitebox2.SetTextSize(fontsize)
    self.newwhitebox2.SetX1NDC(xstart-0.01)
    self.newwhitebox2.SetY1NDC(ystart-0.01)
    self.newwhitebox2.SetX2NDC(xstart+0.25)
    self.newwhitebox2.SetY2NDC(ystart+0.15)
    for textitem in text:
      self.newwhitebox2.AddText("{0}".format(textitem))
    self.newwhitebox2.Draw()
    return

  def flatten_list(self,list) :
    flat_list = []
    for item in list:
      if isinstance(item,Iterable) :
        for subitem in item:
          flat_list.append(subitem)
      else :
        flat_list.append(item)
    return flat_list

  def calculateIntersectionOfGraphs(self, graph1, graph2, doLogGraph1=False, doLogGraph2=False) :

    crossings = []

    for point in range(graph1.GetN()-1) :
      graph1HigherAtLeft=False
      graph1HigherAtRight=False
      thisX1 = graph1.GetX()[point]
      thisX2 = graph1.GetX()[point+1]

      if graph1.Eval(thisX1) <= 0 or graph2.Eval(thisX1) <= 0 \
           or graph2.GetX()[0] > thisX1 or graph2.GetX()[graph2.GetN()-1] < thisX2 :
        continue

      graph1y1 = graph1.GetY()[point]
      graph1y2 = graph1.GetY()[point+1]
      if doLogGraph2 :
        graph2y1 = self.getGraphAtXWithLog(graph2,thisX1)
        graph2y2 = self.getGraphAtXWithLog(graph2,thisX2)
      else :
        graph2y1 = graph2.Eval(thisX1)
        graph2y2 = graph2.Eval(thisX2)

      if (graph1y1 > graph2y1) :
        graph1HigherAtLeft = True
      if (graph1y2 > graph2y2) :
        graph1HigherAtRight = True

      # If no crossing in this interval carry on
      if (graph1HigherAtLeft == graph1HigherAtRight) :
        continue

      # Otherwise, figure out where and keep it
      xtest = 0.5*(thisX1+thisX2)
      while(abs(thisX1-thisX2) > 0.001) :

        if doLogGraph1 :
          graph1y1 = self.getGraphAtXWithLog(graph1,thisX1)
          graph1y2 = self.getGraphAtXWithLog(graph1,thisX2)
          graph1ytest = self.getGraphAtXWithLog(graph1,xtest)
        else :
          graph1y1 = graph1.Eval(thisX1)
          graph1y2 = graph1.Eval(thisX2)
          graph1ytest = graph1.Eval(xtest)
        if doLogGraph2 :
          graph2y1 = self.getGraphAtXWithLog(graph2,thisX1)
          graph2y2 = self.getGraphAtXWithLog(graph2,thisX2)
          graph2ytest = self.getGraphAtXWithLog(graph2,xtest)
        else :
          graph2y1 = graph2.Eval(thisX1)
          graph2y2 = graph2.Eval(thisX2)
          graph2ytest = graph2.Eval(xtest)

        if ((graph1y1 >= graph2y1) and (graph1y2 <= graph2y2)) :
          if graph1ytest > graph2ytest :
            thisX1 = xtest
          else :
            thisX2 = xtest

        elif ((graph1y1 <= graph2y1) and (graph1y2 >= graph2y2)) :
          if graph1ytest > graph2ytest :
            thisX2 = xtest
          else :
            thisX1 = xtest
        xtest = 0.5*(thisX1+thisX2)

      crossings.append(xtest)

    return crossings

  def getGraphAtXWithLog(self, graph, x) :

    for point in range(graph.GetN()-1) :
      thisX1 = graph.GetX()[point]
      thisX2 = graph.GetX()[point+1]
      if thisX1 > x or thisX2 < x :
        continue
      thisY1 = graph.GetY()[point]
      thisY2 = graph.GetY()[point+1]
      m = (math.log(thisY2) - math.log(thisY1))/(thisX2 - thisX1)
      lny = math.log(thisY1) + m*(x-thisX1)
    return math.exp(lny)

  def fixTheBloodyTickMarks(self, pad, hist, x1, x2, y1, y2, override = False) :

    if not pad.GetLogx() :
      return

    hist.GetXaxis().SetMoreLogLabels()

    if (x2/x1 > 100 or x2/x1 < 10) and not override :
      return

    tick = ROOT.TLine()

    tick.SetLineWidth(1)
    tick.SetLineColor(1)

    pad.Update()
    tickminy = y1
    tickmaxy = y2

    if pad.GetLogy() :
      length = 0.5 - tickminy
      shortlength = 0.5*length
      if tickminy > 0.5 :
        length = 1.5
        tickminy = 1.0
    else :
      if y2 == -y1 and abs(y2) < 5 :
        length = (0.2/1.7)*y2
      else :
        length = 0.2

    TeVScale = True
    if x2 > 1000 :
      TeVScale = False

    if TeVScale :
      xlow = int(math.floor(x1/1000)*1000)
      xup = int(math.ceil(x2/1000)*1000)
    else :
      xlow = int(math.floor(x1))
      xup = int(math.ceil(x2))

    ylatex =  tickminy - 1.3*length if pad.GetLogy() else tickminy - 2*length +10

    if TeVScale :
      startRange = int(xlow*1000)
      stopRange = int(xup*1000)
    else :
      startRange = xlow
      stopRange = xup

    print startRange,stopRange

    for i in xrange(startRange, stopRange, 100) :
      if TeVScale :
        xx = float(i)/1000
      else :
        xx = float(i)
      if xx < x1 or xx > x2 :
        continue
      if pad.GetLogy() :
          order = 0
          result = 1
          while result != 0 :
            order = order+1
            result = int(hist.GetBinContent(hist.GetMaximumBin())/pow(10,order))
          tick.DrawLine(xx, tickminy, xx, tickminy + (length if i%1000 == 0 else shortlength))
          factor = math.ceil(tickmaxy/pow(10,order))+1
          tick.DrawLine(xx, tickmaxy - factor*pow(10,order)*(length if i%1000 == 0 else shortlength), xx, tickmaxy)

      else :
          tick.DrawLine(xx, tickminy, xx, tickminy + (2*length if i%1000 == 0 else length))
          tick.DrawLine(xx, tickmaxy - (2*length if i%1000 == 0 else length), xx, tickmaxy)

    ROOT.gPad.Update()

  def fixTheBloodyLabels(self,pad,x1,x2,fontSize=0.04,nLabels=7,overrideY=0) :

    if not pad.GetLogx() :
      return

    # First point to mark will be first above x1
    firstTick = self.roundUpOrderOfMagnitude(x1)

    # Last will be last below x2
    lastTick = self.roundDownOrderOfMagnitude(x2)

    # Set of possible labels are spaced by order of magnitude.
    possibleLabels = []
    val = firstTick
    for order in range(self.magnitude(firstTick),self.magnitude(lastTick)+1) :
      order_list = []
      for integer in range(1,10) :
        possibleVal = integer * pow(10,order)
        if possibleVal >= firstTick and possibleVal <= lastTick :
          order_list.append(possibleVal)
      possibleLabels.append(order_list)

    # Optimal number of labels is about 7, unless otherwise specified.
    # Take them off the top end of the lowest order of magnitude
    # until there are more left at the next order up, then go there.
    while self.countNestedList(possibleLabels) > nLabels :
      for order in range(len(possibleLabels)-1) :
        if len(possibleLabels[order]) > len(possibleLabels[order+1]) :
          possibleLabels[order].pop()
    finalLabels = [item for sublist in possibleLabels for item in sublist]

    print "Chose",nLabels,"optimal labels:",finalLabels
    
    # Now decide where to put them and draw them on.
    labelmaker = ROOT.TLatex()
    labelmaker.SetTextColor(ROOT.kBlack)
    labelmaker.SetTextFont(42)
    labelmaker.SetTextSize(fontSize)
    labelmaker.SetTextAlign(22)
    
    ymin = pad.GetUymin()
    ymax = pad.GetUymax()
    
    # If it's a log y axis, we need to go to a place
    # which probably doesn't exist in user coordinates for this.
    # We have to switch to pad coordinates, but also keep the
    # x locations.
    if pad.GetLogy() :
      
      ylocation = 0.7*pad.GetBottomMargin()
      if overrideY :
        ylocation = overrideY

      # The following then sets the label text to turn it into [0,1] coordinates.
      labelmaker.SetNDC()
      # And overwrite the x locations.
      for label in finalLabels :
        xndc = self.convertXUserToNDC(pad,label,x1,x2)
        labeltext = "{0}".format(label)
        # Center horizontally and vertically around specified location
        labelmaker.DrawLatex(xndc,ylocation,labeltext)

    else :

      ylocation = ymin - abs(float(ymax) - float(ymin))/6.0
      if overrideY :
        ylocation = overrideY
    
      for label in finalLabels :
        labeltext = "{0}".format(label)
        # Center horizontally and vertically around specified location
        labelmaker.DrawLatex(label,ylocation,labeltext)

  def magnitude (self,value):
    value = float(value)
    if (value == 0): return 0
    return int(math.floor(math.log10(abs(value))))

  def roundUpOrderOfMagnitude(self, value) :
    value = float(value)
    magnitude = self.magnitude(value)
    xhigh = int(math.ceil(value/pow(10,magnitude)))*pow(10,magnitude)
    return xhigh

  def roundDownOrderOfMagnitude(self,value) :
    value = float(value)
    magnitude = self.magnitude(value)
    xlow = int(math.floor(value/pow(10,magnitude)))*pow(10,magnitude)
    return xlow

  def countNestedList(self,testlist) :
    length = 0
    for item in testlist :
      if isinstance(item,list)  :
        length = length+self.countNestedList(item)
      else :
        length = length+1
    return length

  def convertXUserToNDC(self,pad,val,minAxis,maxAxis) :

    pad.Update() #this is necessary!

    xpLeft = pad.GetLeftMargin()
    xpRight = pad.GetRightMargin()

    if pad.GetLogx() :
      xPercentageAxis = (math.log(float(val))-math.log(minAxis))/(math.log(maxAxis)-math.log(minAxis))
    else :
      xPercentageAxis = (float(val)-minAxis)/(maxAxis-minAxis)

    # Left margin + right margin + axis = 1.0
    axisWidth = 1.0-xpLeft-xpRight
    xndc = xpLeft + xPercentageAxis*axisWidth
    return xndc

