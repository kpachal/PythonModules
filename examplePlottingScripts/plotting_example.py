# An example python script for making a tidy plot with several lines on it
# I'm using TGraphs because that's what I have in a file, but you can
# do the same thing with histograms if you want to
# The only difference is for histograms you should always do
# hist.SetDirectory(0)
# after you get it out of the root file,
# and you have to do slightly different things to plot
# several graphs on the same canvas (see below)

import ROOT
import math

# Imports ATLAS style for plotting
import AtlasStyle
AtlasStyle.SetAtlasStyle()

# Load some graphs from the example file
# (Gaussian limits from 2016 TLA conf)
infile = ROOT.TFile.Open("example_graphs.root","READ")
graphs = []
legendLines = []

graph1 = infile.Get("gaussians_resolutionwidth_scaled")
graphs.append(graph1)
legendLines.append("#sigma_{G}/m_{G} = Res.")

graph2 = infile.Get("gaussians_70_scaled")
graphs.append(graph2)
legendLines.append("#sigma_{G}/m_{G} = 0.07")

graph3 = infile.Get("gaussians_100_scaled")
graphs.append(graph3)
legendLines.append("#sigma_{G}/m_{G} = 0.10")

# Make a canvas to put the plot on.
# Format it to have two log axes.
c = ROOT.TCanvas("canvas",'',0,0,600,600)
c.SetLogx(True)
c.SetLogy(False)
c.SetGridx(0)
c.SetGridy(0)

# Decide what x and y range to use in the display.
xRange = [0.4,1.3]
yRange = [0.1,10]

# Decide what colours to use.
# These ones look decent, but obviously use
# whatever you like best.
goodColours = [ROOT.kCyan+2,ROOT.kBlue+1,ROOT.kMagenta+1]

# If you want a complicated but fun example of how you can define your own
# colours based on RGB codes, email me!

# Make a legend.
# These are the locations of the left side, bottom side, right
# side, and top, as fractions of the canvas.
legend = ROOT.TLegend(0.6,0.72,0.92,0.92)
# Make the text a nice fond, and big enough
legend.SetTextFont(42)
legend.SetTextSize(0.04)
# A few more formatting things .....
legend.SetBorderSize(0)
legend.SetLineColor(0)
legend.SetLineStyle(1)
legend.SetLineWidth(1)
legend.SetFillColor(0)
legend.SetFillStyle(0)

# Draw each graph
for graph, line in zip(graphs,legendLines) :

  index = graphs.index(graph)
  colour = goodColours[index]

  # Set up marker to look nice
  graph.SetMarkerColor(colour)
  graph.SetMarkerSize(0.7)  # was 0.5 back when things were nice
  graph.SetMarkerStyle(20+index) # was 20 when things were nice

  # Set up line to look nice
  graph.SetLineColor(colour)
  graph.SetLineWidth(2)
  graph.SetLineStyle(1)

  # Make sure we don't get a fill
  graph.SetFillColor(0)

  # Label my axes!!
  graph.GetXaxis().SetTitle("m_{G} [TeV]")
  graph.GetYaxis().SetTitle("#sigma x A x BR [pb]")

  # Set the limit for the axes
  graph.GetXaxis().SetLimits(xRange[0],xRange[1])
  graph.GetYaxis().SetRangeUser(yRange[0],yRange[1])

  # This makes you get nice tick marks on the mass plot
  # It doesn't help much for this one because we have
  # a really inconvenient mass range, but should look
  # pretty good on other plots.
  graph.GetXaxis().SetNdivisions(805,ROOT.kTRUE)

  # This stuff is less annoying for histograms. For graphs,
  # the first one has to include axes or everything comes out blank
  # Rest have to NOT include axes or each successive one overwrites
  # previous. "SAME option does not exist for TGraph classes.
  if index==0 :
    graph.Draw("APL") # Data points of measurement
  else :
    graph.Draw("PL")

  # Fill entry into legend
  # "PL" means both the line and point style
  # will show up in the legend.
  legend.AddEntry(graph,line,"PL")

# Actually draw the legend
legend.Draw()

# This is one way to draw text on the plot
myLatex = ROOT.TLatex()
myLatex.SetTextColor(ROOT.kBlack)
myLatex.SetNDC()

# Put an ATLAS Internal label
# I think it has to be Helvetica
myLatex.SetTextSize(0.05)
myLatex.SetTextFont(72)
# These are the x and y coordinates of the bottom left corner of the text
# as fractions of the canvas
myLatex.DrawLatex(0.18,0.18,"ATLAS")
# Now we switch back to normal font for the "Internal"
myLatex.SetTextFont(42)
myLatex.DrawLatex(0.35,0.18,"Internal")

# Update the canvas
c.Update()

# Save the output as a .eps, a .C, and a .root
c.SaveAs("plot.eps")
c.SaveSource("plot.C")
c.SaveSource("plot.root")
