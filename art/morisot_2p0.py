#import sys
import ROOT
import AtlasStyle
import math
from array import array
from colourPalette import ColourPalette

class Morisot_2p0(object) :

  ## ----------------------------------------------------
  ## Initialisers

  def __init__(self) :

    # Set up style
    AtlasStyle.SetAtlasStyle()
    ROOT.gROOT.ForceStyle()

    # Alternate output formats if requested
    self.saveCFile=False
    self.saveRootFile=False
    self.savePDF=False

    # Colours
    self.colourpalette = ColourPalette()
    self.colourpalette.setColourPalette("Tropical")

    # Set one of these
    self.labeltype = 2 # ATLAS Internal

    # 1 "Preliminary"
    # 2 "Internal"
    # 3 "Simulation Preliminary"
    # 4 "Simulation Internal"
    # 5 "Simulation"
    # 6 "Work in Progress"

    # For internal use
    self.myLatex = ROOT.TLatex()
    self.myLatex.SetTextColor(ROOT.kBlack)
    self.myLatex.SetNDC()

  ###------------------------------------------------###
  ### Helpful utilities

  def getFileKeys(self,open_rootfile,dir="") :
    open_rootfile.cd(dir)
    return sorted([key.GetName() for key in ROOT.gDirectory.GetListOfKeys()])

  ###------------------------------------------------###
  ### Support functions

  # Make a canvas
  def makeCanvas(self,name,logx=False,logy=True,doRectangular=False,scaleX=1.0,scaleY=1.0) :
    canvasname = name+'_cv'
    if doRectangular :
      dim = int(800*scaleX),int(600*scaleY)
    else :
      dim = int(600*scaleX),int(600*scaleY)
    canvas = ROOT.TCanvas(canvasname,'',0,0,dim[0],dim[1])
    canvas.SetLogx(logx)
    canvas.SetLogy(logy)
    return canvas

  def saveCanvas(self,canvas,outputname) :
    canvas.RedrawAxis()
    canvas.Update()
    canvas.SaveAs(outputname+".eps")
    simpleName = outputname.split(".")[0]
    if self.saveCFile:
      canvas.SaveSource(simpleName+".C")
    if self.saveRootFile:
      canvas.SaveSource(simpleName+".root")
    if self.savePDF:
      canvas.SaveAs(simpleName+".pdf")

  def setStandardTwoPads(self,logx=False,logy=False) :
      
    # Dimensions: xlow, ylow, xup, yup
    outpad = ROOT.TPad("extpad","extpad",0,0,1,1) # For marking outermost dimensions
    pad1 = ROOT.TPad("pad1","pad1",0,0.27,1,1) # For main histo
    pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.27) # For residuals histo

    # Set up to draw in right orientations
    outpad.SetFillStyle(4000) #transparent
    pad1.SetBottomMargin(0.01)
    pad1.SetBorderMode(0)
    pad1.SetLogy(1)
    pad2.SetTopMargin(0.04)
    pad2.SetBottomMargin(0.4) # 0.3
    pad2.SetBorderMode(0)
    pad1.Draw()
    pad2.Draw()
    outpad.Draw()
    pad1.SetLogx(logx)
    pad2.SetLogx(logx)
    pad1.SetLogy(logy)
    pad2.SetLogy(0)
    return outpad,pad1,pad2

  def setColourPalette(self,palette) :
    self.colourpalette.setColourPalette(palette)

  def makeLegend(self,legX1,legY1,legX2,legY2,fontSize = 0.04,nColumns=1) :

    legend = ROOT.TLegend(legX1,legY1,legX2,legY2)
    legend.SetTextFont(42)
    legend.SetTextSize(fontSize)
    legend.SetBorderSize(0)
    legend.SetLineColor(0)
    legend.SetLineStyle(1)
    legend.SetLineWidth(1)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)#1001)
    if nColumns > 1 :
      legend.SetNColumns(nColumns)
    return legend

  def drawHorizontalLine(self,x1,x2,y,solid=False) :
    line = ROOT.TLine(x1,y,x2,y)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(2)
    if solid :
      line.SetLineStyle(1)
    else :
      line.SetLineStyle(2)
    line.Draw("SAME")
    return line

  def getLineColours(self, ncolours) :
    if ncolours == 0 :
      return []
    if ncolours < 5 :
      return self.colourpalette.shortLineColours
    elif ncolours < 7 :
      return self.colourpalette.defaultLineColours
    elif ncolours < 13 :
      return self.colourpalette.mediumLineColours
    else :
      return self.colourpalette.longLineColours

  def getFillColours(self,ncolours) :
    if ncolours == 0 :
      return []    
    if ncolours < 5 :
      return self.colourpalette.shortFillColours
    elif ncolours < 7 :
      return self.colourpalette.defaultFillColours
    elif ncolours < 13 :
      return self.colourpalette.mediumFillColours
    else :
      return self.colourpalette.longFillColours

  def getSignalColours(self,ncolours) :
    return self.colourpalette.signalLineColours

  def getAxisRangesFromHist(self,hist) :

    # Get x axis limits
    firstBin =0
    while (hist.GetBinContent(firstBin+1)==0 and firstBin < hist.GetNbinsX()+1) :
      firstBin+=1
    lastBin = hist.GetNbinsX()+1
    while (hist.GetBinContent(lastBin-1)==0 and lastBin > 0) :
      lastBin-=1
    if (firstBin > lastBin) :
      firstBin=1
      lastBin = hist.GetNbinsX()

    # Get y axis limits
    actualMin = 1E10
    actualMax = 0
    for bin in range(firstBin,lastBin) :
      if hist.GetBinContent(bin) > actualMax : actualMax = hist.GetBinContent(bin)
      elif hist.GetBinContent(bin) < actualMin and hist.GetBinContent(bin) != 0 : 
        actualMin = hist.GetBinContent(bin)  

    return hist.GetBinLowEdge(firstBin), hist.GetBinLowEdge(lastBin)+hist.GetBinWidth(lastBin), actualMin, actualMax

  def getAxisRangesFromGraph(self,graph) :

    xVals = []
    yVals = []
    for i in range(graph.GetN()) :
      xVals.append(graph.GetX()[i])
      yVals.append(graph.GetY()[i])
    xVals.sort()
    yVals.sort()
    return xVals[0],xVals[-1],yVals[0],yVals[-1]

  def pickNiceYLimits(self, ylow, yhigh, extraRoom ) :

    if not ROOT.gPad.GetLogy() :
      y1  = 0.01
      y2 = max(y1,ylow*1.5)
    else :
      y1 = max(0.3, ylow/5.0)
      # Amount varies depending on how big y2 is.
      # Want it to come up about 3/4 of the way.
      if yhigh > 0 :
        y2 = yhigh*5.0*math.log(yhigh,10)
      else :
        y2 = yhigh*5.0
    if extraRoom :
      y2 = y2 * 3

    return y1, y2

  def createRatio(self,hist_num,histlist_denom) :
    # Error is just bin error of hist 1 divided by hist 2
    # unless otherwise specified.
    ratio = hist_num.Clone()
    ratio.SetName(hist_num.GetName()+"_ratioplot")
    ratio.SetDirectory(0)
    ratio.Reset()
    total_hist = hist_num.Clone()
    total_hist.SetName(hist_num.GetName()+"_stackForDenominator")
    total_hist.SetDirectory(0)
    total_hist.Reset()
    for hist_denom in histlist_denom :
      total_hist.Add(hist_denom)
    for bin in range (0,ratio.GetNbinsX()+2) : 
      if total_hist.GetBinContent(bin) != 0 :
        ratio.SetBinContent(bin,hist_num.GetBinContent(bin)/total_hist.GetBinContent(bin))
        ratio.SetBinError(bin,hist_num.GetBinError(bin)/total_hist.GetBinContent(bin))
      else :
        ratio.SetBinContent(bin,0)
        ratio.SetBinError(bin,0)
    return ratio

  def formatForPad2(self,hist,low=None,high=None) :
    hist.GetYaxis().SetTitleSize(0.15)
    hist.GetYaxis().SetTitleOffset(0.4) # 1.2 = 20% larger
    hist.GetYaxis().SetLabelSize(0.14)
    hist.GetXaxis().SetLabelSize(0.14)
    hist.GetXaxis().SetTitleSize(0.15)
    hist.GetXaxis().SetTitleOffset(1.2)
    if not low and not high :
      hist.GetYaxis().SetRangeUser(0.45,1.55)
    else :
       hist.GetYaxis().SetRangeUser(low,high)

    # Nice tick marks
    hist.GetYaxis().SetTickLength(0.055)
    hist.GetXaxis().SetTickLength(0.055)

    return hist

  def drawATLASLabel(self,xval,yval,rightalign=False,isRectangular=False,fontSize=0.05) :
    if self.labeltype < 0 :
      return
    self.myLatex.SetTextSize(fontSize)
    self.myLatex.SetTextFont(72)
    if rightalign :
      self.myLatex.SetTextAlign(31)
    else :
      self.myLatex.SetTextAlign(11)
    self.myLatex.DrawLatex(xval, yval, "ATLAS")

    if self.labeltype==0 :
      return
    spacing = 0.17*(fontSize/0.05)
    if (isRectangular) :
      spacing = 0.14*(fontSize/0.05)

    self.myLatex.SetTextFont(42)
    if self.labeltype==1 :
      self.myLatex.DrawLatex(xval + spacing, yval, "Preliminary")
    elif self.labeltype==2 :
      self.myLatex.DrawLatex(xval + spacing, yval, "Internal")
    elif self.labeltype==3 :
      self.myLatex.DrawLatex(xval + spacing, yval, "Simulation Preliminary")
    elif self.labeltype==4 :
      self.myLatex.DrawLatex(xval + spacing, yval, "Simulation Internal")
    elif self.labeltype==5 :
      self.myLatex.DrawLatex(xval + spacing, yval, "Simulation")
    elif self.labeltype==6 :
      self.myLatex.DrawLatex(xval + spacing, yval, "Work in Progress")
    return    

  def drawCMEAndLumi(self,xval,yval,CME,lumiInFb,fontSize=0.05) :
    if lumiInFb < 0 and CME < 0 :
      return

    mysqrt = "#sqrt{s}"
    myfb = "fb^{-1}"

    if isinstance(lumiInFb, (list,tuple)):
      lumitext = '-'.join(["{0}".format(l) for l in lumiInFb])
    else: 
      lumitext = lumiInFb
      mytext = "{0}={1} TeV, {2} {3}".format(mysqrt,CME,lumitext,myfb)

    newtext = ROOT.TLatex()
    newtext.SetNDC()
    newtext.SetTextSize(fontSize)
    newtext.SetTextFont(42)
    newtext.SetTextAlign(11)
    newtext.DrawLatex(xval,yval,"{0}".format(mytext)) 
    return

  def drawText(self,xval,yval,text,rightalign=False,fontSize=0.04) :

    usertext = ROOT.TLatex()
    usertext.SetNDC()
    usertext.SetTextFont(42)
    usertext.SetTextSize(fontSize)
    if not rightalign :
      usertext.SetTextAlign(11)
    else :
      usertext.SetTextAlign(31)
    usertext.DrawLatex(xval,yval,text) 
    return usertext

  ###------------------------------------------------###
  ### Not-quite-standalone constituent plots
  def drawDataHist(self, dataHist,xlow,xhigh,xname,yname,same=False,nPads=1,FixYAxis=False,LeaveAxisAlone=False,extraRoom=False) :

    # Data hist must be in data points with weighted error bars
    dataHist.SetMarkerColor(ROOT.kBlack)
    dataHist.SetLineColor(ROOT.kBlack)
    dataHist.GetXaxis().SetTitle(xname)
    dataHist.GetYaxis().SetTitle(yname)
    dataHist.GetYaxis().SetTitleSize(0.06)
    dataHist.GetYaxis().SetTitleOffset(1.0)
    dataHist.GetXaxis().SetRangeUser(xlow,xhigh)
    dataHist.SetMarkerSize(1.0) # was 0.75 in dijet ISR

    xmin_auto, xmax_auto, actualMin, actualMax = self.getAxisRangesFromHist(dataHist)
    
    y1, y2 = self.pickNiceYLimits(actualMin, actualMax, extraRoom)
    if not LeaveAxisAlone :
      dataHist.GetYaxis().SetRangeUser(y1,y2)

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

  def drawBackgroundHist(self, histo, colour, xlow,xhigh,xname,yname,same=False,doErrors=False,nPads=1,LeaveAxisAlone=False,extraRoom=False) :

    # Format with intention of drawing a single line with no fill.
    histo.SetLineColor(colour)
    histo.SetMarkerColor(colour)
    histo.SetLineStyle(1)
    histo.SetLineWidth(2)
    histo.SetFillStyle(0)

    # Set up labels and titles
    histo.SetTitle("")
    histo.GetYaxis().SetRangeUser()
    histo.GetYaxis().SetTitleSize(0.06)
    histo.GetYaxis().SetTitleOffset(1.3) # 1.2
    histo.GetYaxis().SetLabelSize(0.06)
      
    histogram.GetXaxis().SetTitleSize(0.06)
    histogram.GetXaxis().SetTitleOffset(1.2)
    if (doLogX) :
      histogram.GetXaxis().SetLabelSize(0)
    else :
      histogram.GetXaxis().SetLabelSize(0.06)

    # Axis ranges and divisions
    histo.GetXaxis().SetRangeUser(xlow,xhigh)   
    histo.GetXaxis().SetNdivisions(605,ROOT.kTRUE)
    histo.GetYaxis().SetNdivisions(605,ROOT.kTRUE)
      
    if not same :
      histo.GetXaxis().SetTitle(xname)
      histo.GetYaxis().SetTitle(yname)
      if not doErrors :
        histo.Draw("HIST")
      else :
        histo.Draw("E")
    else :
      histo.GetXaxis().SetTitle("")
      histo.GetYaxis().SetTitle("")
      if not doErrors :
        histo.Draw("HIST SAME")
      else :
        histo.Draw("E SAME")

    # if doLogX :
    #   self.fixTheBloodyTickMarks(ROOT.gPad,histogram, minX,maxX,minY,maxY)
    #   self.fixTheBloodyLabels(ROOT.gPad,minX,maxX,fontSize=0.06,nLabels=7,overrideY=0.15,suppressFirstOrder=True)    

  ###------------------------------------------------###
  ### Primary functions


  def drawOverlaidHistos(self,histos_list,data=None,signal_list=None,histos_labels=[],data_label="",xlabel="",ylabel="",plotname="",doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,luminosity=None,CME=13,nLegendColumns=1,extraLines=[],xLow=None,xHigh=None,yLow=None,yHigh=None) :

    # Make a canvas
    outpad,pad1,pad2 = None,None,None
    c = self.makeCanvas(plotname,logx,logy)

    # Make a legend
    legendDepth = 0.04*len(histos_list)/float(nLegendColumns) + 0.04*(histos_list is not None)
    legendTop = 0.90 if not doRatio else 0.95
    legend = self.makeLegend(0.5,legendTop-legendDepth,0.92,legendTop,nColumns=nLegendColumns,fontSize = 0.04)

    if doRatio :
      outpad,pad1,pad2 = self.setStandardTwoPads(logx,logy)
      pad1.cd()

    if data :
      xlow, xhigh, actualMin, actualMax = self.getAxisRangesFromHist(data)
      legend.AddEntry(data,data_label,"LP")

    # Get colours
    goodcolours = self.getLineColours(len(histos_list))

    # Draw histograms
    for histogram in histos_list :
      index = histos_list.index(histogram)
      if index == 0 :
        self.drawBackgroundHist()
      else :
        self.drawBackgroundHist()

    # Draw data, if there is data
    if data : self.drawDataHist(data,xlow,xhigh,xlabel,ylabel,same=True,nPads=1,extraRoom=False)

    return


  def drawStackedHistos(self,stack_list,data=None,signal_list=None,stack_labels=[],data_label="",xlabel="",ylabel="",plotname="",doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,luminosity=None,CME=13,nLegendColumns=2,extraLines=[],xLow=None,xHigh=None,yLow=None,yHigh=None) :

    # Make a canvas
    outpad,pad1,pad2 = None,None,None
    c = self.makeCanvas(plotname,logx,logy)

    # Make a legend
    legendDepth = 0.04*len(stack_labels)/float(nLegendColumns) + 0.04*(signal_list is not None)
    legendTop = 0.90 if not doRatio else 0.95
    legend = self.makeLegend(0.5,legendTop-legendDepth,0.92,legendTop,nColumns=nLegendColumns,fontSize = 0.03)

    if doRatio :
      outpad,pad1,pad2 = self.setStandardTwoPads(logx,logy)
      pad1.cd()

    if data :
      xlow, xhigh, actualMin, actualMax = self.getAxisRangesFromHist(data)
      self.drawDataHist(data,xlow,xhigh,xlabel,ylabel,same=False,nPads=1,extraRoom=False)
      legend.AddEntry(data,data_label,"LP")

    # Make a stack for items in the stack_list
    goodcolours = self.getFillColours(len(stack_list))

    stack = ROOT.THStack("stack","stacked histograms")
    for histogram in stack_list :
      index = stack_list.index(histogram)
      histogram.SetLineColor(goodcolours[index])
      histogram.SetLineWidth(2)
      histogram.SetFillColor(goodcolours[index])
      if nLegendColumns != 2 :
        legend.AddEntry(histogram,stack_labels[index],"F")
      histogram.SetTitle("")
      stack.Add(histogram,"hist")

    # Want legend entries to read naturally which means filling it
    # in a weird order.
    # Also, this goes under "Data" so need to start with second half.
    if nLegendColumns == 2 :
      split_1 = stack_list[0:int((len(stack_list)+1)/2.0)]
      split_2 = stack_list[int((len(stack_list)+1)/2.0):]
      names_1 = stack_labels[0:int((len(stack_list)+1)/2.0)]
      names_2 = stack_labels[int((len(stack_list)+1)/2.0):]
      reordered_hists = []
      reordered_names = []
      for i in reversed(range(max(len(split_1),len(split_2)))) :
        if len(split_2) > i :
          reordered_hists.append(split_2[i])
          reordered_names.append(names_2[i])
        if len(split_1) > i :
          reordered_hists.append(split_1[i])
          reordered_names.append(names_1[i])
      for hist, name in zip(reordered_hists,reordered_names) :
        legend.AddEntry(hist,name,"F")

    # Without data, stack sets plot ranges and formats
    if not data :
      stack.Draw()
      stack.SetMaximum(stack.GetMaximum()*5.0 if logy else stack.GetMaximum()*1.4)
      stack.SetMinimum(0.5 if logy else 0)
      # A sensible stack should have the biggest contribution (i.e. longest tails) on top
      mc_xlow, mc_xhigh, false_ymin, false_ymax = self.getAxisRangesFromHist(stack_list[-1])
      stack.GetXaxis().SetRangeUser(mc_xlow,mc_xhigh)      
      stack.GetXaxis().SetTitle(xlabel)
      stack.GetYaxis().SetTitle(ylabel)
      c.Update()

    # With data, put it on top
    else :
      stack.Draw("SAME")
      # Put the data hist on top
      self.drawDataHist(data,xlow,xhigh,xlabel,ylabel,same=True,nPads=1,extraRoom=False)

    # Add signals if requested
    if signal_list :
      sigcolours = self.getSignalColours(len(signal_list))
      for signal in signal_list :
        signal.SetLineColor(sigcolours[signal_list.index(signal)])
        signal.SetLineWidth(3)
        signal.SetLineStyle(2)
        signal.Draw("HIST SAME")
        legend.AddEntry(signal,"Signal, #sigma x 100","L")

    # Format and add ratio plot
    if doRatio :
      pad1.RedrawAxis()
      ratioHist = self.createRatio(data,stack_list)
      pad2.cd()
      updatedHist = self.formatForPad2(ratioHist)
      updatedHist.GetYaxis().SetTitle(ratioName)
      updatedHist.Draw("E")
      dotted_line = self.drawHorizontalLine(xlow,xhigh,1)
      outpad.cd()

    c.Update()

    legend.Draw()

    # Finally, add some labels
    self.drawATLASLabel(0.2,0.9,fontSize=0.04)
    if luminosity :
      self.drawCMEAndLumi(0.2,0.85,CME,luminosity,fontSize=0.035)
    if extraLines :
      for line in extraLines :
        text = self.drawText(0.2,0.80-0.04*extraLines.index(line),line,fontSize=0.035)      

    self.saveCanvas(c,plotname)

  def drawOverlaidTGraphs(self,graphs_list,names_list,xlabel="",ylabel="",plotname="",logx=False,logy=True,luminosity=None,CME=13,xmin=None,xmax=None,ymin=None,ymax=None,addHorizontalLines=[],extraLines=[]) :

    c = self.makeCanvas(plotname,logx,logy)    

    # Set automatic axis range from graphs.
    xVals = []
    yVals = []
    if not graphs_list :
      return

    for thisgraph in graphs_list :
      thisxmin,thisxmax,thisymin,thisymax = self.getAxisRangesFromGraph(thisgraph)
      for i in range(thisgraph.GetN()) :
        xVals.append(thisxmin)
        xVals.append(thisxmax)
        yVals.append(thisymin)
        yVals.append(thisymax)
    # Don't plot if all graphs empty
    if thisymin == 0 and thisymax == 0 : return
    xVals.sort()
    if not xmin :
      xmin = 0.9 * xVals[0]
    if not xmax :
      xmax = 1.1 * xVals[-1]
    yVals.sort()
    if not ymin :
      if logy :
        ymin = yVals[0]/10.0
      else :
        ymin = yVals[0]*0.8
    if not ymax :
      if logy :
        # More space above than below, for legends
        ymax = yVals[-1]*50.0
      else :
        ymax = yVals[-1]*1.2

    # Create legend.
    # Currently putting ATLAS label + cut info at top left
    # and legend at right.
    legendDepth = 0.04*len(names_list)
    legendTop = 0.90
    legend = self.makeLegend(0.5,legendTop-legendDepth,0.92,legendTop,nColumns=1,fontSize = 0.03)    

    goodcolours = self.getLineColours(len(graphs_list))

    for graph in graphs_list :

      index = graphs_list.index(graph)
      if len(graphs_list) < 4 :
        colour = goodcolours[index+1]
      else :
        colour = goodcolours[index]

      # Set up display for observations
      graph.SetMarkerColor(colour)
      graph.SetMarkerSize(1.0)  # was 0.5 back when things were nice
      graph.SetMarkerStyle(20+index) # was 20 when things were nice
      graph.SetLineColor(colour)
      graph.SetLineWidth(2)
      graph.SetLineStyle(1)
      graph.SetFillColor(0)
      graph.GetXaxis().SetTitle(xlabel)
      graph.GetYaxis().SetTitle(ylabel)
      graph.GetXaxis().SetLimits(xmin,xmax)
      graph.GetYaxis().SetRangeUser(ymin,ymax)
      graph.GetXaxis().SetNdivisions(605,ROOT.kTRUE)

      # First one has to include axes or everything comes out blank
      # Rest have to NOT include axes or each successive one overwrites
      # previous. "SAME option does not exist for TGraph classes.
      if index==0 :
        graph.Draw("APL") # Data points of measurement
      else :
        graph.Draw("PL")

      # Fill and draw legend
      legend.AddEntry(graph,names_list[index],"PL")

    if addHorizontalLines != [] :
      for val in addHorizontalLines :
        dotted_line = self.drawHorizontalLine(xmin,xmax,val,solid=False)

    c.Update()

    legend.Draw()

    # Finally, add some labels
    self.drawATLASLabel(0.2,0.85,fontSize=0.04)
    if luminosity :
      self.drawCMEAndLumi(0.2,0.80,CME,luminosity,fontSize=0.035)
    if extraLines :
      for line in extraLines :
        self.drawText(0.2,0.76-0.04*extraLines.index(line),line,fontSize=0.035)

    self.saveCanvas(c,plotname) 
