#!/usr/bin/env python3
import ROOT

class ColourPalette(object) :

  ## ----------------------------------------------------
  ## Initialisers

  def __init__(self) :

    self.dataLineColor = ROOT.kBlack
    self.dataPointColor = ROOT.kBlack
    self.dataFillColor = ROOT.kWhite
    self.palette = None

  ## ----------------------------------------------------
  ## User-accessible functions

  def setColourPalette(self, colourPalette) :
    self.palette = colourPalette
    if self.palette == "ATLAS" :
      self.setATLASColours()
    elif self.palette == "Oxford" :
      self.setOxfordColours()
    elif self.palette == "Teals" :
      self.setTealColours()
    elif self.palette == "Tropical" :
      self.setTropicalColours()
    elif self.palette == "notSynthwave" : 
      self.setnotSynthwaveColours()

  def getColourPalette(self) :
    return self.palette

  ## ----------------------------------------------------
  ## Individual palette definitions

  def setnotSynthwaveColours(self) : 

    # Plot fill colours         
    self.oneFillColour = ROOT.TColor(500,72,145,220)
    self.oneFitColour = ROOT.TColor(501,190,15,52)
    self.statisticaTestFillColour = ROOT.TColor(502,158,206,235)
    self.originalStatMarkerColour = ROOT.TColor(503,0,33,71)

    self.fitLineColor = ROOT.kRed

    # Limit-setting plots                   
    self.oneSigmaBandColour = ROOT.kAzure+5
    self.twoSigmaBandColour = ROOT.kAzure+5
    self.signalLineColours = [ROOT.kAzure+4,ROOT.kAzure+3]
    self.signalErrorColours = [2013,2013]


    # Stat plots                                
    self.statisticalTestFillColour = 2000
    self.statisticalTestArrowColour = 2003
    self.tomographyGraphColour = 2003

    # Colours for lists
    self.shortLineColours = [ROOT.kAzure+5, ROOT.kBlue+2,ROOT.kRed-3,ROOT.kOrange+1]
    self.defaultLineColours = [ROOT.kAzure+5, ROOT.kBlue+2,ROOT.kRed-3,ROOT.kOrange+1,ROOT.kViolet+2]
    self.mediumLineColours = [ROOT.kAzure+5, ROOT.kBlue+2,ROOT.kRed-3,ROOT.kOrange+1,ROOT.kViolet+2]#+...
    self.longLineColours = self.mediumLineColours
    
    # Paired colours for matching things
    self.pairedComplementaryColours = [[ROOT.kAzure+5, ROOT.kBlue+2],[ROOT.kRed-3,ROOT.kOrange+1]]

    # For now, use the same colours for fills, but code separately
    # to leave option for making lighter ones in future
    self.shortFillColours = [ROOT.kAzure+5, ROOT.kBlue+2,ROOT.kRed-3,ROOT.kOrange+1]
    self.defaultFillColours = [ROOT.kAzure+5, ROOT.kBlue+2,ROOT.kRed-3,ROOT.kOrange+1,ROOT.kViolet+2]
    self.mediumFillColours = [ROOT.kAzure+5, ROOT.kBlue+2,ROOT.kRed-3,ROOT.kOrange+1,ROOT.kViolet+2]#+...
    self.longFillColours = self.mediumFillColours
    
  def setATLASColours(self) :

    # Plot fill colours
    self.oneFillColour = ROOT.kYellow
    self.oneFitColour = ROOT.kRed
    self.statisticaTestFillColour = ROOT.kYellow
    self.originalStatMarkerColour = ROOT.kBlue

    self.fitLineColor = ROOT.kRed

    # Stat plots
    self.statisticaTestFillColour = ROOT.kYellow
    self.statisticalTestArrowColour = ROOT.kBlue
    self.tomographyGraphColour = ROOT.kRed

    # Limit-setting plots
    self.oneSigmaBandColour = ROOT.kGreen#ROOT.kSpring
    self.twoSigmaBandColour = ROOT.kYellow
    self.signalLineColours = [ROOT.kBlue,ROOT.kAzure+7,ROOT.kCyan+3]#ROOT.kTeal+9]
    self.signalErrorColours = [ROOT.kBlue-9,ROOT.kCyan-9,ROOT.kTeal+8]#ROOT.kAzure+6,ROOT.kTeal+8]

    self.shortLineColours = self.mediumLineColours
    self.defaultLineColours = self.mediumlineColours

  def setOxfordColours(self) :

    # Plot fill colours
    self.oneFillColour = ROOT.TColor(500,72,145,220)
    self.oneFitColour = ROOT.TColor(501,190,15,52)
    self.statisticaTestFillColour = ROOT.TColor(502,158,206,235)
    self.originalStatMarkerColour = ROOT.TColor(503,0,33,71)

    self.fitLineColor = ROOT.kRed

    # Limit-setting plots
    self.oneSigmaBandColour = ROOT.TColor(500)
    self.twoSigmaBandColour = ROOT.TColor(502)
    self.signalLineColour = ROOT.TColor(501)

  def setTealColours(self) :

    self.blue1 = ROOT.TColor(2010,0.833,1.00,0.967)
    self.blue2 = ROOT.TColor(2011,0.392,0.825,0.850)
    self.bluegreen = ROOT.TColor(2012,0.200,0.470,0.521)
    self.darkred = ROOT.TColor(2013,0.808,0.063,0.175)
    self.mauve = ROOT.TColor(2014,0.817,0.792,0.783)

    # Limit-setting plots
    self.oneSigmaBandColour = 2011
    self.twoSigmaBandColour = 2012
    self.signalLineColours = [2013,2013]
#    self.transparentColor1 = ROOT.gROOT.GetColor(1013)#GetColorTransparent(1013,0.5)
#    self.transparentColor1.SetAlpha(0.5)
#    self.signalErrorColours = [self.transparentColor1,self.transparentColor1]
    self.signalErrorColours = [2013,2013]

    # Stat plots
    self.statisticalTestFillColour = 2011
    self.statisticalTestArrowColour = 2013
    self.tomographyGraphColour = 2013

    self.defaultLineColours = [2010,2011,2012,2013,2014]
    self.shortLineColours = [2012,2011,2013]

  def setTropicalColours(self) :

    self.blue1 = ROOT.TColor(2000,13.0/255.0,159.0/255.0,216.0/255.0)
    self.mauve = ROOT.TColor(2001,133.0/255.0,105.0/255,207.0/255.0)
    self.green = ROOT.TColor(2002,124.0/255.0,194.0/255.0,66.0/255.0)
    self.darkred = ROOT.TColor(2003,248.0/255.0,14.0/255.0,39.0/255.0)
    self.mango = ROOT.TColor(2004,238.0/255.0, 206.0/255.0, 0.0/255.0)

    # Plot fill colours
    self.oneFillColour = ROOT.TColor(500,72,145,220)
    self.oneFitColour = ROOT.TColor(501,190,15,52)
    self.statisticaTestFillColour = ROOT.TColor(502,158,206,235)
    self.originalStatMarkerColour = ROOT.TColor(503,0,33,71)

    self.fitLineColor = ROOT.kRed

    # Limit-setting plots
    self.oneSigmaBandColour = 2002
    self.twoSigmaBandColour = 2004
    self.signalLineColours = [2003,2003]

#    self.transparentColor1 = ROOT.gROOT.GetColor(1003)
#    self.transparentColor1.SetAlpha(0.5)
#    self.signalErrorColours = [self.transparentColor1,self.transparentColor1]
    self.signalErrorColours = [2003,2003]

    # Stat plots
    self.statisticalTestFillColour = 2000
    self.statisticalTestArrowColour = 2003
    self.tomographyGraphColour = 2003

    # Colours for lists
    self.shortLineColours = [ROOT.kViolet+6,ROOT.kSpring-5,ROOT.kPink-1]
    self.defaultLineColours = [ROOT.kViolet+6,ROOT.kAzure+1,ROOT.kSpring-5,ROOT.kPink-1,ROOT.kRed]
    self.mediumLineColours = [ROOT.kCyan+4,ROOT.kCyan+2,ROOT.kCyan,\
                              ROOT.kBlue,ROOT.kBlue+2,\
                              ROOT.kMagenta+2,ROOT.kMagenta,\
                              ROOT.kRed,ROOT.kRed+2,ROOT.kOrange+10,\
                              ROOT.kOrange,ROOT.kYellow]
    self.longLineColours = [ROOT.kCyan+4,ROOT.kCyan+3,ROOT.kCyan+2,ROOT.kCyan+1,ROOT.kCyan,\
                     ROOT.kBlue,ROOT.kBlue+1,ROOT.kBlue+2,ROOT.kBlue+3,ROOT.kBlue+4,\
                     ROOT.kMagenta+4,ROOT.kMagenta+3,ROOT.kMagenta+2,ROOT.kMagenta+1,ROOT.kMagenta,\
                     ROOT.kRed,ROOT.kRed+1,ROOT.kRed+2,ROOT.kOrange+9,ROOT.kOrange+10,\
                     ROOT.kOrange+7,ROOT.kOrange,ROOT.kYellow,ROOT.kGreen-4,ROOT.kGreen-1,ROOT.kGreen+2,ROOT.kTeal]

