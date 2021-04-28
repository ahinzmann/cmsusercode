import os, sys
from ROOT import * 
from DataFormats.FWLite import Events,Handle
import array
from math import *
from scipy import stats

#gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
#gROOT.Reset()
gROOT.SetBatch(True)
gROOT.SetStyle("Plain")
gStyle.SetOptStat(0)
gStyle.SetOptFit(0)
gStyle.SetTitleOffset(1.2,"Y")
gStyle.SetPadLeftMargin(0.18)
gStyle.SetPadBottomMargin(0.15)
gStyle.SetPadTopMargin(0.03)
gStyle.SetPadRightMargin(0.05)
gStyle.SetMarkerSize(1.5)
gStyle.SetHistLineWidth(1)
gStyle.SetStatFontSize(0.020)
gStyle.SetTitleSize(0.06, "XYZ")
gStyle.SetLabelSize(0.05, "XYZ")
gStyle.SetNdivisions(510, "XYZ")
gStyle.SetLegendBorderSize(0)

def rebin(h1,nbins,binning):
    for b in range(h1.GetXaxis().GetNbins()):
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)*h1.GetBinWidth(b+1))
        h1.SetBinError(b+1,h1.GetBinError(b+1)*h1.GetBinWidth(b+1))
    h1=h1.Rebin(nbins,h1.GetName()+"_rebin",binning)
    for b in range(h1.GetXaxis().GetNbins()):
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)/h1.GetBinWidth(b+1))
        h1.SetBinError(b+1,h1.GetBinError(b+1)/h1.GetBinWidth(b+1))
    return h1

def cloneNormalize(h1):
    h1=h1.Clone(h1.GetName()+"clone")
    for b in range(h1.GetXaxis().GetNbins()):
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)/h1.GetBinWidth(b+1))
        h1.SetBinError(b+1,h1.GetBinError(b+1)/h1.GetBinWidth(b+1))
    return h1

def smooth(h1,f):
    fit=TF1(h1.GetName()+"smooth",f,1,16)
    h1.Fit(fit,"NQ")
    for chi_bin in range(h1.GetXaxis().GetNbins()):
      h1.SetBinContent(chi_bin+1,fit.Eval(h1.GetBinCenter(chi_bin+1)))
    return h1

def smoothChi(h1):
    for b in range(h1.GetXaxis().GetNbins()):
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)/h1.GetBinWidth(b+1))
        h1.SetBinError(b+1,h1.GetBinError(b+1)/h1.GetBinWidth(b+1))
    fit=TF1(h1.GetName()+"smooth","pol3",3,16)
    h1.Fit(fit,"RNQ")
    for chi_bin in range(2,h1.GetXaxis().GetNbins()):
      h1.SetBinContent(chi_bin+1,fit.Eval(h1.GetBinCenter(chi_bin+1)))
    for b in range(h1.GetXaxis().GetNbins()):
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)*h1.GetBinWidth(b+1))
        h1.SetBinError(b+1,h1.GetBinError(b+1)*h1.GetBinWidth(b+1))
    return h1

if __name__ == '__main__':

    useLensData=False
    useUnfoldedData=False
    injectSignal=False
    only6000=False # mass binning
    useNNLO=True # choice for QCD
    useM2=True # choice of mu-scale for QCD
    
    if useNNLO:
      pdfset="ct14nnlo"
    else:
      pdfset="ct14nlo"
    if useM2:
      muScale="m2"
      muAltScale="pt12"
    else:
      muScale="pt12"
      muAltScale="m2"

    prefixs=["versions/run2NNLOMar25/datacard_shapelimit13TeV"]
 
    # negligible jes sources removed
    jessources=["AbsoluteScale",
                "AbsoluteMPFBias",
		"SinglePionHCAL",
		"FlavorQCD",
		"TimePtEta2016","TimePtEta2017","TimePtEta2018",
		"RelativePtBB",
		"RelativeBal",
		"RelativeFSR",
		"RelativeSample2016","RelativeSample2017","RelativeSample2018",
		"RelativeStatEC2016","RelativeStatEC2017","RelativeStatEC2018",
		"RelativeJEREC12016","RelativeJEREC12017","RelativeJEREC12018",
		"PileUpDataMC",
		"PileUpPtRef",
		"PileUpPtBB",
		"PileUpPtEC1",
		"SumInQuadrature",
		]
    print len(jessources)-1,"jes sources"
    jersources=["JER12016","JER12017","JER12018",
                "JER22016","JER22017","JER22018",
		"SumInQuadrature",
		]

 
    chi_bins=[(1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,3,6,9,12,16),
              ]
    massbins=[(1200,1500),
              (1500,1900),
              (1900,2400),
              (2400,3000),
              (3000,3600),
              (3600,4200),
              (4200,4800),
              (4800,5400),
              (5400,6000),
              (6000,7000),
              (7000,13000)]
    mass_bins_nlo3={}
    mass_bins_nlo3[0]=1200
    mass_bins_nlo3[1]=1500
    mass_bins_nlo3[2]=1900
    mass_bins_nlo3[3]=2400
    mass_bins_nlo3[4]=3000
    mass_bins_nlo3[5]=3600
    mass_bins_nlo3[6]=4200
    mass_bins_nlo3[7]=4800
    mass_bins_nlo3[8]=5400
    mass_bins_nlo3[9]=6000
    mass_bins_nlo3[10]=7000
    mass_bins_nlo3[11]=13000
    mass_bins_nlo_list=[(0,),
                  (1,),
                  (2,),
                  (3,),
                  (4,),
                  (5,),
                  (6,),
              (7,),
              (8,),
              (9,),
	      (10,)
                 ]
    if only6000:
      chi_bins=[(1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,2,3,4,5,6,7,8,9,10,12,14,16),
               (1,3,6,9,12,16),
              ]
      massbins=[(1200,1500),
              (1500,1900),
              (1900,2400),
              (2400,3000),
              (3000,3600),
              (3600,4200),
              (4200,4800),
              (4800,5400),
              (5400,6000),
              (6000,13000)]
      mass_bins_nlo_list=[(0,),
                  (1,),
                  (2,),
                  (3,),
                  (4,),
                  (5,),
                  (6,),
              (7,),
              (8,),
              (9,10,)
                 ]
		 
    if useUnfoldedData:
       massbins=massbins[1:]
       mass_bins_nlo_list=mass_bins_nlo_list[1:]
       chi_bins=chi_bins[1:]
    chi_binnings=[]
    for mass_bin in chi_bins:
        chi_binnings+=[array.array('d')]
        for chi_bin in mass_bin:
            chi_binnings[-1].append(chi_bin)

    samples=[]
    samples1=[]
    samples2=[]
    samples3=[]
    samples4=[]
    samples5=[]
    samples6=[]
    samplesSkip=[]

    samples+=[("QCD",[("pythia8_ci_m1000_1500_50000_1_0_0_13TeV_Nov14",3.769e-05),
                       ("pythia8_ci_m1500_1900_50000_1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_50000_1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_50000_1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_50000_1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_50000_1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_50000_1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_50000_1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
            ]
    samples+=[("QCDAntiCIplusLL12000",[("pythia8_ci_m1500_1900_12000_1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_12000_1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_12000_1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_12000_1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_12000_1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_12000_1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_12000_1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
             ]
    samples+=[("QCDCIplusLL12000",[("pythia8_ci_m1500_1900_12000_1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_12000_1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_12000_1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_12000_1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_12000_1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_12000_1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_12000_1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
             ]
    samplesSkip=[("QCDCIplusLL8000",[("pythia8_ci_m1500_1900_8000_1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_8000_1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_8000_1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_8000_1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_8000_1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_8000_1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_8000_1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDCIplusLL9000",[("pythia8_ci_m1500_1900_9000_1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_9000_1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_9000_1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_9000_1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_9000_1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_9000_1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_9000_1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDCIplusLL10000",[("pythia8_ci_m1500_1900_10000_1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_10000_1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_10000_1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_10000_1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_10000_1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_10000_1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_10000_1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDCIplusLL11000",[("pythia8_ci_m1500_1900_11000_1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_11000_1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_11000_1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_11000_1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_11000_1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_11000_1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_11000_1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDCIplusLL12000",[("pythia8_ci_m1500_1900_12000_1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_12000_1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_12000_1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_12000_1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_12000_1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_12000_1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_12000_1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDCIplusLL13000",[("pythia8_ci_m1500_1900_13000_1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_13000_1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_13000_1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_13000_1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_13000_1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_13000_1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_13000_1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDCIplusLL14000",[("pythia8_ci_m1500_1900_14000_1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_14000_1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_14000_1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_14000_1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_14000_1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_14000_1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_14000_1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDCIplusLL16000",[("pythia8_ci_m1500_1900_16000_1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_16000_1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_16000_1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_16000_1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_16000_1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_16000_1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_16000_1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDCIplusLL18000",[("pythia8_ci_m1500_1900_18000_1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_18000_1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_18000_1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_18000_1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_18000_1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_18000_1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_18000_1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
             ]
    samplesSkip+=[("QCDCIminusLL8000",[("pythia8_ci_m1500_1900_8000_-1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_8000_-1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_8000_-1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_8000_-1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_8000_-1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_8000_-1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_8000_-1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDCIminusLL9000",[("pythia8_ci_m1500_1900_9000_-1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_9000_-1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_9000_-1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_9000_-1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_9000_-1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_9000_-1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_9000_-1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDCIminusLL10000",[("pythia8_ci_m1500_1900_10000_-1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_10000_-1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_10000_-1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_10000_-1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_10000_-1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_10000_-1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_10000_-1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDCIminusLL11000",[("pythia8_ci_m1500_1900_11000_-1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_11000_-1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_11000_-1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_11000_-1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_11000_-1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_11000_-1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_11000_-1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDCIminusLL12000",[("pythia8_ci_m1500_1900_12000_-1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_12000_-1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_12000_-1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_12000_-1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_12000_-1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_12000_-1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_12000_-1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDCIminusLL13000",[("pythia8_ci_m1500_1900_13000_-1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_13000_-1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_13000_-1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_13000_-1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_13000_-1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_13000_-1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_13000_-1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDCIminusLL14000",[("pythia8_ci_m1500_1900_14000_-1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_14000_-1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_14000_-1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_14000_-1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_14000_-1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_14000_-1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_14000_-1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDCIminusLL16000",[("pythia8_ci_m1500_1900_16000_-1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_16000_-1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_16000_-1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_16000_-1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_16000_-1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_16000_-1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_16000_-1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDCIminusLL18000",[("pythia8_ci_m1500_1900_18000_-1_0_0_13TeV_Nov14",3.307e-06),
                       ("pythia8_ci_m1900_2400_18000_-1_0_0_13TeV_Nov14",8.836e-07),
                       ("pythia8_ci_m2400_2800_18000_-1_0_0_13TeV_Nov14",1.649e-07),
                       ("pythia8_ci_m2800_3300_18000_-1_0_0_13TeV_Nov14",6.446e-08),
                       ("pythia8_ci_m3300_3800_18000_-1_0_0_13TeV_Nov14",1.863e-08),
                       ("pythia8_ci_m3800_4300_18000_-1_0_0_13TeV_Nov14",5.867e-09),
                       ("pythia8_ci_m4300_13000_18000_-1_0_0_13TeV_Nov14",3.507e-09),
                       ]),
             ]
    samples1+=[#("QCDADD6000",[("pythia8_add_m1500_1900_6000_0_0_0_1_13TeV_Nov14",3.307e-06),
              #         ("pythia8_add_m1900_2400_6000_0_0_0_1_13TeV_Nov14",8.836e-07),
              #         ("pythia8_add_m2400_2800_6000_0_0_0_1_13TeV_Nov14",1.649e-07),
              #         ("pythia8_add_m2800_3300_6000_0_0_0_1_13TeV_Nov14",6.446e-08),
              #         ("pythia8_add_m3300_3800_6000_0_0_0_1_13TeV_Nov14",1.863e-08),
              #         ("pythia8_add_m3800_4300_6000_0_0_0_1_13TeV_Nov14",5.867e-09),
              #         ("pythia8_add_m4300_13000_6000_0_0_0_1_13TeV_Nov14",3.507e-09),
              #         ]),
             #("QCDADD7000",[("pythia8_add_m1500_1900_7000_0_0_0_1_13TeV_Nov14",3.307e-06),
             #          ("pythia8_add_m1900_2400_7000_0_0_0_1_13TeV_Nov14",8.836e-07),
             #          ("pythia8_add_m2400_2800_7000_0_0_0_1_13TeV_Nov14",1.649e-07),
             #          ("pythia8_add_m2800_3300_7000_0_0_0_1_13TeV_Nov14",6.446e-08),
             #          ("pythia8_add_m3300_3800_7000_0_0_0_1_13TeV_Nov14",1.863e-08),
             #          ("pythia8_add_m3800_4300_7000_0_0_0_1_13TeV_Nov14",5.867e-09),
             #          ("pythia8_add_m4300_13000_7000_0_0_0_1_13TeV_Nov14",3.507e-09),
             #          ]),
             #("QCDADD8000",[("pythia8_add_m1500_1900_8000_0_0_0_1_13TeV_Nov14",3.307e-06),
             #          ("pythia8_add_m1900_2400_8000_0_0_0_1_13TeV_Nov14",8.836e-07),
             #          ("pythia8_add_m2400_2800_8000_0_0_0_1_13TeV_Nov14",1.649e-07),
             #          ("pythia8_add_m2800_3300_8000_0_0_0_1_13TeV_Nov14",6.446e-08),
             #          ("pythia8_add_m3300_3800_8000_0_0_0_1_13TeV_Nov14",1.863e-08),
             #          ("pythia8_add_m3800_4300_8000_0_0_0_1_13TeV_Nov14",5.867e-09),
             #          ("pythia8_add_m4300_13000_8000_0_0_0_1_13TeV_Nov14",3.507e-09),
             #          ]),
             ("QCDADD9000",[("pythia8_add_m1500_1900_9000_0_0_0_1_13TeV_Nov14",3.307e-06),
                       ("pythia8_add_m1900_2400_9000_0_0_0_1_13TeV_Nov14",8.836e-07),
                       ("pythia8_add_m2400_2800_9000_0_0_0_1_13TeV_Nov14",1.649e-07),
                       ("pythia8_add_m2800_3300_9000_0_0_0_1_13TeV_Nov14",6.446e-08),
                       ("pythia8_add_m3300_3800_9000_0_0_0_1_13TeV_Nov14",1.863e-08),
                       ("pythia8_add_m3800_4300_9000_0_0_0_1_13TeV_Nov14",5.867e-09),
                       ("pythia8_add_m4300_13000_9000_0_0_0_1_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDADD10000",[("pythia8_add_m1500_1900_10000_0_0_0_1_13TeV_Nov14",3.307e-06),
                       ("pythia8_add_m1900_2400_10000_0_0_0_1_13TeV_Nov14",8.836e-07),
                       ("pythia8_add_m2400_2800_10000_0_0_0_1_13TeV_Nov14",1.649e-07),
                       ("pythia8_add_m2800_3300_10000_0_0_0_1_13TeV_Nov14",6.446e-08),
                       ("pythia8_add_m3300_3800_10000_0_0_0_1_13TeV_Nov14",1.863e-08),
                       ("pythia8_add_m3800_4300_10000_0_0_0_1_13TeV_Nov14",5.867e-09),
                       ("pythia8_add_m4300_13000_10000_0_0_0_1_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDADD11000",[("pythia8_add_m1500_1900_11000_0_0_0_1_13TeV_Nov14",3.307e-06),
                       ("pythia8_add_m1900_2400_11000_0_0_0_1_13TeV_Nov14",8.836e-07),
                       ("pythia8_add_m2400_2800_11000_0_0_0_1_13TeV_Nov14",1.649e-07),
                       ("pythia8_add_m2800_3300_11000_0_0_0_1_13TeV_Nov14",6.446e-08),
                       ("pythia8_add_m3300_3800_11000_0_0_0_1_13TeV_Nov14",1.863e-08),
                       ("pythia8_add_m3800_4300_11000_0_0_0_1_13TeV_Nov14",5.867e-09),
                       ("pythia8_add_m4300_13000_11000_0_0_0_1_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDADD12000",[("pythia8_add_m1500_1900_12000_0_0_0_1_13TeV_Nov14",3.307e-06),
                       ("pythia8_add_m1900_2400_12000_0_0_0_1_13TeV_Nov14",8.836e-07),
                       ("pythia8_add_m2400_2800_12000_0_0_0_1_13TeV_Nov14",1.649e-07),
                       ("pythia8_add_m2800_3300_12000_0_0_0_1_13TeV_Nov14",6.446e-08),
                       ("pythia8_add_m3300_3800_12000_0_0_0_1_13TeV_Nov14",1.863e-08),
                       ("pythia8_add_m3800_4300_12000_0_0_0_1_13TeV_Nov14",5.867e-09),
                       ("pythia8_add_m4300_13000_12000_0_0_0_1_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDADD13000",[("pythia8_add_m1500_1900_13000_0_0_0_1_13TeV_Nov14",3.307e-06),
                       ("pythia8_add_m1900_2400_13000_0_0_0_1_13TeV_Nov14",8.836e-07),
                       ("pythia8_add_m2400_2800_13000_0_0_0_1_13TeV_Nov14",1.649e-07),
                       ("pythia8_add_m2800_3300_13000_0_0_0_1_13TeV_Nov14",6.446e-08),
                       ("pythia8_add_m3300_3800_13000_0_0_0_1_13TeV_Nov14",1.863e-08),
                       ("pythia8_add_m3800_4300_13000_0_0_0_1_13TeV_Nov14",5.867e-09),
                       ("pythia8_add_m4300_13000_13000_0_0_0_1_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDADD14000",[("pythia8_add_m1500_1900_14000_0_0_0_1_13TeV_Nov14",3.307e-06),
                       ("pythia8_add_m1900_2400_14000_0_0_0_1_13TeV_Nov14",8.836e-07),
                       ("pythia8_add_m2400_2800_14000_0_0_0_1_13TeV_Nov14",1.649e-07),
                       ("pythia8_add_m2800_3300_14000_0_0_0_1_13TeV_Nov14",6.446e-08),
                       ("pythia8_add_m3300_3800_14000_0_0_0_1_13TeV_Nov14",1.863e-08),
                       ("pythia8_add_m3800_4300_14000_0_0_0_1_13TeV_Nov14",5.867e-09),
                       ("pythia8_add_m4300_13000_14000_0_0_0_1_13TeV_Nov14",3.507e-09),
                       ]),
             ("QCDADD15000",[("pythia8_add_m1500_1900_15000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m1900_2400_15000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m2400_2800_15000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m2800_3300_15000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m3300_3800_15000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m3800_4300_15000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m4300_5200_15000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m5200_13000_15000_0_0_0_1_13TeV_Nov14",1),
                       ]),
             ("QCDADD16000",[("pythia8_add_m1500_1900_16000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m1900_2400_16000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m2400_2800_16000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m2800_3300_16000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m3300_3800_16000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m3800_4300_16000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m4300_5200_16000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m5200_13000_16000_0_0_0_1_13TeV_Nov14",1),
                       ]),
             ("QCDADD17000",[("pythia8_add_m1500_1900_17000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m1900_2400_17000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m2400_2800_17000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m2800_3300_17000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m3300_3800_17000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m3800_4300_17000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m4300_5200_17000_0_0_0_1_13TeV_Nov14",1),
                       ("pythia8_add_m5200_13000_17000_0_0_0_1_13TeV_Nov14",1),
                       ]),
             #("QCDADD18000",[("pythia8_add_m1500_1900_18000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m1900_2400_18000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m2400_2800_18000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m2800_3300_18000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m3300_3800_18000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m3800_4300_18000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m4300_5200_18000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m5200_13000_18000_0_0_0_1_13TeV_Nov14",1),
             #          ]),
             #("QCDADD19000",[("pythia8_add_m1500_1900_19000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m1900_2400_19000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m2400_2800_19000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m2800_3300_19000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m3300_3800_19000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m3800_4300_19000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m4300_5200_19000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m5200_13000_19000_0_0_0_1_13TeV_Nov14",1),
             #          ]),
             #("QCDADD20000",[("pythia8_add_m1500_1900_20000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m1900_2400_20000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m2400_2800_20000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m2800_3300_20000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m3300_3800_20000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m3800_4300_20000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m4300_5200_20000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m5200_13000_20000_0_0_0_1_13TeV_Nov14",1),
             #          ]),
             #("QCDADD21000",[("pythia8_add_m1500_1900_21000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m1900_2400_21000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m2400_2800_21000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m2800_3300_21000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m3300_3800_21000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m3800_4300_21000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m4300_5200_21000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m5200_13000_21000_0_0_0_1_13TeV_Nov14",1),
             #          ]),
             #("QCDADD22000",[("pythia8_add_m1500_1900_22000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m1900_2400_22000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m2400_2800_22000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m2800_3300_22000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m3300_3800_22000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m3800_4300_22000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m4300_5200_22000_0_0_0_1_13TeV_Nov14",1),
             #          ("pythia8_add_m5200_13000_22000_0_0_0_1_13TeV_Nov14",1),
             #          ]),
             ]

    for m in range(5,31):
    #for m in [13]:
       samples2+=[("cs_"+pdfset+"_"+str(m*1000)+"_LL+",[]),
               ("cs_"+pdfset+"_"+str(m*1000)+"_LL-",[]),
               ("cs_"+pdfset+"_"+str(m*1000)+"_RR+",[]),
               ("cs_"+pdfset+"_"+str(m*1000)+"_RR-",[]),
               ("cs_"+pdfset+"_"+str(m*1000)+"_VV+",[]),
               ("cs_"+pdfset+"_"+str(m*1000)+"_VV-",[]),
               ("cs_"+pdfset+"_"+str(m*1000)+"_AA+",[]),
               ("cs_"+pdfset+"_"+str(m*1000)+"_AA-",[]),
               ("cs_"+pdfset+"_"+str(m*1000)+"_V-A+",[]),
               ("cs_"+pdfset+"_"+str(m*1000)+"_V-A-",[]),
               ]

    for mass in [1700,2000,2300,2600,2900,3200,3500,3800,4100,4400,4700,5000,5300,5600,5900,6200,6500,6800,7100]:
     for xsec in [5e-5,1e-4,2e-4,5e-4,1e-3,2e-3,5e-3,1e-2,2e-2,5e-2,1e-2]:
      for width in ["kMG1439","kMG2035","kMG2493","kMG3218"]:
        for decay in ["GluonGluon","QuarkQuark"]:
          samplesSkip+=[#("wide"+str(mass)+"_"+str(xsec)+"_"+decay+"_"+width,["wide"+str(mass)+"_"+str(xsec)+"_"+decay+"_"+width,1])
	            ]

    #xsecs={}
    #for l in open("xsecs_13TeV_dm.txt").readlines():
    #  xsecs[l.split("     ")[0]]=eval(l.split("     ")[1])
    for mass in [1000,1500,1750,2000,2250,2500,3000,3500,4000,4500,5000,6000,7000]:
    #for mass in [7000]:
    #for mass in [2000,2250,2500]:
    #for mass in [3000,3500,4000]:
    #for mass in [4500,5000,6000,7000]:
    #for mass in [2000,2250,2500,3000,3500,4000,4500,5000,6000,7000]:
    #for mass in [2250]:
     #if mass==6000:
     #  mDMs=[1,2990]
     #elif mass==7000:
     #  mDMs=[1,4000]
     #elif mass==8000:
     #  mDMs=[1,3990]
     #else:
     #  mDMs=[1,3000]
     if mass==7000:
       mDMs=[4000] #1 files are broken
     else:
       mDMs=[1]
     for mDM in mDMs:
      for weight in ['gdmv_1p0_gdma_0_gv_0p01_ga_0', 'gdmv_1p0_gdma_0_gv_0p05_ga_0', 'gdmv_1p0_gdma_0_gv_0p1_ga_0', 'gdmv_1p0_gdma_0_gv_0p2_ga_0', 'gdmv_1p0_gdma_0_gv_0p25_ga_0', 'gdmv_1p0_gdma_0_gv_0p3_ga_0', 'gdmv_1p0_gdma_0_gv_0p5_ga_0', 'gdmv_1p0_gdma_0_gv_0p75_ga_0', 'gdmv_1p0_gdma_0_gv_1_ga_0', 'gdmv_1p0_gdma_0_gv_1p5_ga_0', 'gdmv_1p0_gdma_0_gv_2p0_ga_0', 'gdmv_1p0_gdma_0_gv_2p5_ga_0', 'gdmv_1p0_gdma_0_gv_3p0_ga_0']:
         samplesSkip+=[("DMVector_Dijet_LO_Mphi_"+str(mass)+"_"+str(mDM)+"_1p0_1p0_Mar5_"+weight,[("DMVector_Dijet_LO_Mphi_"+str(mass)+"_"+str(mDM)+"_1p0_1p0_Mar5_"+weight,0)]),
             ]
      #for weight in ['gdmv_0_gdma_1p0_gv_0_ga_0p01', 'gdmv_0_gdma_1p0_gv_0_ga_0p05', 'gdmv_0_gdma_1p0_gv_0_ga_0p1', 'gdmv_0_gdma_1p0_gv_0_ga_0p2', 'gdmv_0_gdma_1p0_gv_0_ga_0p25', 'gdmv_0_gdma_1p0_gv_0_ga_0p3', 'gdmv_0_gdma_1p0_gv_0_ga_0p5', 'gdmv_0_gdma_1p0_gv_0_ga_0p75', 'gdmv_0_gdma_1p0_gv_0_ga_1', 'gdmv_0_gdma_1p0_gv_0_ga_1p5', 'gdmv_0_gdma_1p0_gv_0_ga_2p0', 'gdmv_0_gdma_1p0_gv_0_ga_2p5', 'gdmv_0_gdma_1p0_gv_0_ga_3p0']:
      for weight in ['gdmv_0_gdma_1p0_gv_0_ga_0p1', 'gdmv_0_gdma_1p0_gv_0_ga_0p2', 'gdmv_0_gdma_1p0_gv_0_ga_0p25', 'gdmv_0_gdma_1p0_gv_0_ga_0p3', 'gdmv_0_gdma_1p0_gv_0_ga_0p5', 'gdmv_0_gdma_1p0_gv_0_ga_0p75', 'gdmv_0_gdma_1p0_gv_0_ga_1', 'gdmv_0_gdma_1p0_gv_0_ga_1p5']:
      #for weight in ['gdmv_0_gdma_1p0_gv_0_ga_1p5']:
         samples4+=[("DMAxial_Dijet_LO_Mphi_"+str(mass)+"_"+str(mDM)+"_1p0_1p0_Mar5_"+weight,[("DMAxial_Dijet_LO_Mphi_"+str(mass)+"_"+str(mDM)+"_1p0_1p0_Mar5_"+weight,0)]),
             ]

    for m in [[7500,0.01678352],[8000,0.004871688],[8500,0.001292072],[9000,0.0003054339],[9500,0.00006221544],[10000,0.00001040396],[10500,0.000001327101],[11000,0.0000001145733]]:
        samples3+=[("QBH_"+str(m[0])+"_6",[("QBH_"+str(m[0])+"_6",m[1])]),]

    for m in [[4500,0.05148],[5000,0.01829],[5500,0.006472],[6000,0.002250],[6500,0.0007599],[7000,0.0002461]]:
        samples3+=[("QBH_"+str(m[0])+"_RS1",[("QBH_"+str(m[0])+"_RS1",m[1])]),]

    for weight in ['fa1000','fa1500','fa2000','fa2500','fa3000','fa3500','fa4000','fa4500','fa5000','fa50000']:
         samples5+=[("alp_QCD_"+weight,[("alp_QCD_"+weight,0)]),]

    for weight in ["CG0p1","CG0p05","CG0p04","CG0p03","CG0p025","CG0p02","CG0p015","CG0p01","CG0p0075","CG0p005","CG0p0025","CG0p0"]:
         samples6+=[("tripleG_QCD_"+weight,[("tripleG_QCD_"+weight,0)]),]

    # all samples
    samples=samples+samples1+samples2+samples3+samples4+samples5+samples6
    # for alp+tripleG
    samples=samples5
    #samples=samples6
    # for postfit plots
    #samples=[("DMAxial_Dijet_LO_Mphi_7000_4000_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_1",[("DMAxial_Dijet_LO_Mphi_7000_4000_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_1",0)]), ]
    #samples=[("cs_ct14nnlo_14000_V-A-",[])]
    
    if len(sys.argv)>1:
      if int(sys.argv[1])>=len(samples):
        samples=[]
      else:
        samples=[samples[int(sys.argv[1])]]

    print samples

    dataevents={}
    datahist={}
    data=None
    dataplot={}
    qcdnorm={}
    cinorm={}
    for prefix in prefixs: 
     # signal cards
     for i in range(len(samples)):
      if samples[i][0]=="QCDCIplusLL8000":
        sample=prefix + '_GENnp-0-run2_chi.root'
      elif samples[i][0]=="QCDCIplusLL9000":
        sample=prefix + '_GENnp-1-run2_chi.root'
      elif samples[i][0]=="QCDCIplusLL10000":
        sample=prefix + '_GENnp-2-run2_chi.root'
      elif samples[i][0]=="QCDCIplusLL11000":
        sample=prefix + '_GENnp-3-run2_chi.root'
      elif samples[i][0]=="QCDCIplusLL12000":
        sample=prefix + '_GENnp-4-run2_chi.root'
      elif samples[i][0]=="QCDCIplusLL13000":
        sample=prefix + '_GENnp-5-run2_chi.root'
      elif samples[i][0]=="QCDCIplusLL14000":
        sample=prefix + '_GENnp-6-run2_chi.root'
      elif samples[i][0]=="QCDCIplusLL16000":
        sample=prefix + '_GENnp-7-run2_chi.root'
      elif samples[i][0]=="QCDCIplusLL18000":
        sample=prefix + '_GENnp-8-run2_chi.root'
      elif samples[i][0]=="QCDCIminusLL8000":
        sample=prefix + '_GENnp-9-run2_chi.root'
      elif samples[i][0]=="QCDCIminusLL9000":
        sample=prefix + '_GENnp-10-run2_chi.root'
      elif samples[i][0]=="QCDCIminusLL10000":
        sample=prefix + '_GENnp-11-run2_chi.root'
      elif samples[i][0]=="QCDCIminusLL11000":
        sample=prefix + '_GENnp-12-run2_chi.root'
      elif samples[i][0]=="QCDCIminusLL12000":
        sample=prefix + '_GENnp-13-run2_chi.root'
      elif samples[i][0]=="QCDCIminusLL13000":
        sample=prefix + '_GENnp-14-run2_chi.root'
      elif samples[i][0]=="QCDCIminusLL14000":
        sample=prefix + '_GENnp-15-run2_chi.root'
      elif samples[i][0]=="QCDCIminusLL16000":
        sample=prefix + '_GENnp-16-run2_chi.root'
      elif samples[i][0]=="QCDCIminusLL18000":
        sample=prefix + '_GENnp-17-run2_chi.root'
      elif samples[i][0]=="QCDADD6000":
        sample=prefix + '_GENnp-18-run2_chi.root'
      elif samples[i][0]=="QCDADD7000":
        sample=prefix + '_GENnp-19-run2_chi.root'
      elif samples[i][0]=="QCDADD8000":
        sample=prefix + '_GENnp-20-run2_chi.root'
      elif samples[i][0]=="QCDADD9000":
        sample=prefix + '_GENnp-21-run2_chi.root'
      elif samples[i][0]=="QCDADD10000":
        sample=prefix + '_GENnp-22-run2_chi.root'
      elif samples[i][0]=="QCDADD11000":
        sample=prefix + '_GENnp-23-run2_chi.root'
      elif samples[i][0]=="QCDADD12000":
        sample=prefix + '_GENnp-24-run2_chi.root'
      elif samples[i][0]=="QCDADD13000":
        sample=prefix + '_GENnp-25-run2_chi.root'
      elif samples[i][0]=="QCDADD14000":
        sample=prefix + '_GENnp-26-run2_chi.root'
      elif samples[i][0]=="QCDADD15000":
        sample=prefix + '_GENnp-27-run2_chi.root'
      elif samples[i][0]=="QCDADD16000":
        sample=prefix + '_GENnp-28-run2_chi.root'
      elif samples[i][0]=="QCDADD17000":
        sample=prefix + '_GENnp-29-run2_chi.root'
      elif samples[i][0]=="QCDADD18000":
        sample=prefix + '_GENnp-30-run2_chi.root'
      elif samples[i][0]=="QCDADD19000":
        sample=prefix + '_GENnp-31-run2_chi.root'
      elif samples[i][0]=="QCDADD20000":
        sample=prefix + '_GENnp-32-run2_chi.root'
      elif samples[i][0]=="QCDADD21000":
        sample=prefix + '_GENnp-33-run2_chi.root'
      elif samples[i][0]=="QCDADD22000":
        sample=prefix + '_GENnp-34-run2_chi.root'
      elif samples[i][0]=="QCDAntiCIplusLL12000":
        sample=prefix + '_GENnp-antici-run2_chi.root'
      elif samples[i][0]=="QCD":
        sample=prefix + '_GEN-QCD-run2_chi.root'
      elif "alp" in samples[i][0] or "tripleG" in samples[i][0] or "DM" in samples[i][0] or "ll" in samples[i][0] or "cs" in samples[i][0] or "wide" in samples[i][0] or "QBH" in samples[i][0]:
        sample=prefix + "_" + samples[i][0] + '-run2_chi.root'
      #if "ADD" in samples[i][0]:
      #  sample=prefix + '_GENaddv3_chi2016.root'
      #elif "CIplus" in samples[i][0]:
      #  sample=prefix + '_GENciv3_chi2016.root'
      #elif "CIminus" in samples[i][0]:
      #  sample=prefix + '_GENciminusv3_chi2016.root'
      else:
        sample=prefix + '_GEN-QCD-run2_chi.root'
      print sample
  
      out=TFile(sample,'UPDATE')
      closefiles=[out]
      
      # LO QCD file
      sample2=prefix + '_GEN-QCD-run2_chi.root'
      print sample2
      in2=TFile(sample2,'READ')

      # Madgraph QCD file
      sampleMadgraph=prefix + '_alp_QCD_fa50000-run2_chi.root'
      print sampleMadgraph
      inMadgraph=TFile(sampleMadgraph,'READ')

      # data file
      #insample='datacards/chiHist_dataReReco_v3_PFHT900.root' #2016
      #insample='datacards/datacard_shapelimit13TeV_25nsData13combi_chi.root' # buggy data
      insample="datacard_shapelimit13TeV_run2_2016_chi.root" # Aug rereco
      print insample
      infile=TFile(insample,'READ')
      #insample17='datacards/uhh2.AnalysisModuleRunner.DATA.Run2017_RunBCDEF_17Nov2017-v1.root' #2017
      insample17="datacard_shapelimit13TeV_run2_2017_chi.root"
      print insample17
      infile17=TFile(insample17,'READ')
      #insample18='datacards/uhh2.AnalysisModuleRunner.DATA.Run2018_RunABCD_RunII_102X_v1.root' #2018
      insample18="datacard_shapelimit13TeV_run2_2018_chi.root"
      print insample18
      infile18=TFile(insample18,'READ')

      # unfolded data file
      #unfoldsample='datacards/Unfolded_chiNtuple_dataReReco_v3_Coarse_PFHT900_fromCB_AK4SF_pythia8_Pt_170toInf.root'
      unfoldsample='datacards/Unfolded_chiNtuple_dataReReco_v3_Coarse_PFHT900_fromCB_AK4SF_pythia8_Pt_170toInf_MatrixInvert.root'
      print unfoldsample
      unfoldfile=TFile(unfoldsample,'READ')

      # (N)NLO correction
      #filename1nu2="fastnlo/RunII/fnl5662j_v23_fix_CT14nlo_allmu_ak4.root"
      if useNNLO:
        filename1nu2="fastnlo/NNLO/2jet.NNLO.fnl5662j_mjj_chi_ct14nnlo_cppread_mu_"+muScale+".root"
      else:
        filename1nu2="fastnlo/NNLO/2jet.NNLO.fnl5662j_mjj_chi_ct14nlo_cppread_mu_pt12.root"
      print filename1nu2
      nlofile2 = TFile.Open(filename1nu2)
      closefiles+=[nlofile2]

      # (N)NLO uncertainties
      filename1nusys="fastnlo/NNLO/fnl5662j_cs_"+pdfset+"_"+muScale+"_30000_LL+.root"
      print filename1nusys
      nlofilesys = TFile.Open(filename1nusys)
      closefiles+=[nlofilesys]
      filename1altscale="fastnlo/NNLO/fnl5662j_cs_"+pdfset+"_"+muAltScale+"_30000_LL+.root"
      print filename1altscale
      nlofilealtscale = TFile.Open(filename1altscale)
      closefiles+=[nlofilealtscale]
      
       # DM uncertainties
      filename1dmpdf="datacards/chi_dm_pdf_plots6000_13TeV_2016.root" # FIX recompute DM PDF uncertainties
      print filename1dmpdf
      dmpdffile = TFile.Open(filename1dmpdf)
      closefiles+=[dmpdffile]

      # EWK correction
      filename1ewk="fastnlo/RunII/DijetAngularCMS13_ewk.root"
      print filename1ewk
      ewkfile = TFile.Open(filename1ewk)
      closefiles+=[ewkfile]

      # JES uncertainty QCD
      jesfiles=[]
      for n in range(5):
        filename1jes="chi_systematic_plotschi_QCDJES_"+str(n)+"_13TeV_run2.root"
        print filename1jes
        jesfiles += [TFile.Open(filename1jes)]
        closefiles+=[jesfiles[-1]]
      jeshists={}
      for j in range(len(massbins)):
        for source in jessources:
          jeshists[str(j)+source]=[]
	  for jesfile in jesfiles:
	    jespad=jesfile.Get("jes")
            jes=jespad.GetListOfPrimitives()[j+useUnfoldedData*1]
            jeshists[str(j)+source]+=[a for a in jes.GetListOfPrimitives() if source in str(a)]
          if len(jeshists[str(j)+source])!=2:
	    print "JES source not found", source
	    error

      # JES uncertainty CI
      jescifiles=[]
      for n in range(5):
        filename1jesci="chi_systematic_plotschi_QCDJES_"+str(n)+"_13TeV_run2.root"
        print filename1jesci
        jescifiles += [TFile.Open(filename1jesci)]
        closefiles+=[jescifiles[-1]]
      jescihists={}
      for j in range(len(massbins)):
        for source in jessources:
          jescihists[str(j)+source]=[]
	  for jescifile in jescifiles:
	    jespad=jescifile.Get("jes")
            jes=jespad.GetListOfPrimitives()[j+useUnfoldedData*1]
            jescihists[str(j)+source]+=[a for a in jes.GetListOfPrimitives() if source in str(a)]
          if len(jescihists[str(j)+source])!=2:
	    print "JES source not found", source
	    error

      # JER uncertainty QCD
      jerfiles=[]
      for n in range(2):
        filename1jer="chi_systematic_plotschi_QCDJER_"+str(n)+"_13TeV_run2.root"
        print filename1jer
        jerfiles += [TFile.Open(filename1jer)]
        closefiles+=[jerfiles[-1]]
      jerhists={}
      for j in range(len(massbins)):
        for source in jersources:
          jerhists[str(j)+source]=[]
	  for jerfile in jerfiles:
	    jerpad=jerfile.Get("jer")
            jer=jerpad.GetListOfPrimitives()[j+useUnfoldedData*1]
            jerhists[str(j)+source]+=[a for a in jer.GetListOfPrimitives() if source in str(a)]
          if len(jerhists[str(j)+source])!=2:
	    print "JER source not found", source
	    error

      # JER uncertainty CI
      jercifiles=[]
      for n in range(2):
        filename1jerci="chi_systematic_plotschi_QCDJER_"+str(n)+"_13TeV_run2.root"
        print filename1jerci
        jercifiles += [TFile.Open(filename1jerci)]
        closefiles+=[jercifiles[-1]]
      jercihists={}
      for j in range(len(massbins)):
        for source in jersources:
          jercihists[str(j)+source]=[]
	  for jercifile in jercifiles:
	    jerpad=jercifile.Get("jer")
            jer=jerpad.GetListOfPrimitives()[j+useUnfoldedData*1]
            jercihists[str(j)+source]+=[a for a in jer.GetListOfPrimitives() if source in str(a)]
          if len(jercihists[str(j)+source])!=2:
	    print "JER source not found", source
	    error

      canvas = TCanvas("","",0,0,800,600)
      canvas.Divide(4,3)
      plots=[]
      legends=[]

      for j in range(len(massbins)):
        histname="dijet_"+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_chi"
        print histname
        if useLensData:
          if "13000" in str(massbins[j]):
            histname2="dijet_m_chi_2__projY_"+str(massbins[j]).strip("()").replace(',',"-").replace(' ',"")
          else:
            histname2="dijet_m_chi_2__projY_"+str(massbins[j]).strip("()").replace(',',"-").replace(' ',"")
          print histname2
          data = TH1D(unfoldfile.Get(histname2))
          data.SetName(histname)
          data=data.Rebin(len(chi_binnings[j])-1,data.GetName()+"_rebin1",chi_binnings[j])
        elif useUnfoldedData:
          histname2="dijet_mass_"+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_chi_unfolded"
          print histname2
          data = TH1F(unfoldfile.Get(histname2))
          data.SetName(histname)
          data=data.Rebin(len(chi_binnings[j])-1,data.GetName()+"_rebin1",chi_binnings[j])
        else:
          #histname2="dijet_"+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"").replace("1200_1500","1900_2400").replace("1500_1900","1900_2400").replace("6000_7000","6000_13000").replace("7000_13000","6000_13000")+"_chi"
          histname2="datacard_shapelimit13TeV_run2_2016#chi"+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print histname2
          data = TH1F(infile.Get(histname2))
          data.SetName(histname)
          data=data.Rebin(len(chi_binnings[j])-1,data.GetName()+"_rebin1",chi_binnings[j])
	  #if j<2: data.Scale(0) # FIX Recompute lowest mass bins for 2016
          #if j>=10: data.Scale(0) # FIX Recompute highest mass bins for 2016
          #histname17="Dijet/chi_"+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"").replace("1200_1500","1900_2400").replace("1500_1900","1900_2400").replace("6000_7000","6000_6600").replace("7000_13000","6600_13000")
          histname17="datacard_shapelimit13TeV_run2_2017#chi"+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print histname17
          data17 = TH1F(infile17.Get(histname17))
          data17=data17.Rebin(len(chi_binnings[j])-1,data.GetName()+"_rebin1",chi_binnings[j])
          data.Add(data17)
          #histname18="Dijet/chi_"+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")
          histname18="datacard_shapelimit13TeV_run2_2018#chi"+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print histname18
          data18 = TH1F(infile18.Get(histname18))
          data18=data18.Rebin(len(chi_binnings[j])-1,data.GetName()+"_rebin1",chi_binnings[j])
	  data.Add(data18)
        dataevents[j]=data.Integral()
        print dataevents[j]
        out.cd()
        if not injectSignal:
         histname='data_obs#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
         for k in range(0,200):
             out.Delete(histname+";"+str(k))
         data.Write(histname)

        # NLO
        nloqcd=None
        for k in mass_bins_nlo_list[j]:
         #histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])
         histname='qcd_chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-1.0"
         print histname
         hnlo = TH1F(nlofile2.Get(histname))
         #hnlo.Scale(float(mass_bins_nlo3[k+1]-mass_bins_nlo3[k]))
         hnlo=hnlo.Rebin(len(chi_binnings[j])-1,hnlo.GetName()+"_rebin1",chi_binnings[j])
	 #hnlo=rebin(hnlo,len(chi_binnings[j])-1,chi_binnings[j])
         if nloqcd:
            nloqcd.Add(hnlo)
         else:
            nloqcd=hnlo
        #for b in range(nloqcd.GetXaxis().GetNbins()):
        #  nloqcd.SetBinContent(b+1,nloqcd.GetBinContent(b+1)/nloqcd.GetBinWidth(b+1))
        nloqcd=smoothChi(nloqcd) # SMOOTH NNLO PREDICTION (FIX ME)
        nloqcdbackup=nloqcd.Clone(nloqcd.GetName()+"_backup")
	print "NLO integral (pb):", nloqcdbackup.Integral()

        # NLO normalized
        nloqcdnorm=None
        for k in mass_bins_nlo_list[j]:
         histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])
         print histname
         hnlo = TH1F(nlofilesys.Get(histname))
         #hnlo.Scale(float(mass_bins_nlo3[k+1]-mass_bins_nlo3[k]))
         hnlo=rebin(hnlo,len(chi_binnings[j])-1,chi_binnings[j])
         if nloqcdnorm:
            nloqcdnorm.Add(hnlo)
         else:
            nloqcdnorm=hnlo

        # EWK corrections
        histname='chi-'+str(massbins[j]).strip("()").replace(',',"-").replace(' ',"").replace("1200-1500","1900-2400").replace("1500-1900","1900-2400").replace("6000-7000","6000-6600").replace("6000-13000","6000-6600").replace("7000-13000","6600-13000") # FIX compute lower and higher mass bins
        print histname
        ewk=ewkfile.Get(histname)
        for b in range(nloqcd.GetXaxis().GetNbins()):
            low_bin=ewk.FindBin(nloqcd.GetXaxis().GetBinLowEdge(b+1))
            up_bin=ewk.FindBin(nloqcd.GetXaxis().GetBinUpEdge(b+1))
            correction=ewk.Integral(low_bin,up_bin-1)/(up_bin-low_bin)
            if not "EWK" in samples[i][0]:
               nloqcd.SetBinContent(b+1,nloqcd.GetBinContent(b+1)*correction)
        nloqcd.Scale(1./nloqcd.Integral())
        ewk.SetName("ewk-"+histname)

        # QCD (empty background, not used in limit)
        histname='QCD#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_backup"
        print histname
        if only6000:
          qcd=in2.Get(histname.replace("6000_13000","6000_7000"))
        else:
	  qcd=in2.Get(histname)
        out.cd()
        for k in range(0,200):
            out.Delete(histname.replace("_backup","")+";"+str(k))
        qcd.Write(histname.replace("_backup",""))
        qcd=qcd.Rebin(len(chi_binnings[j])-1,histname,chi_binnings[j])
        qcd.Scale(1e10*1e9) #mb -> pb and factor 1e-10 for backup
        if j in qcdnorm.keys():
           qcd.Scale(qcdnorm[j]/qcd.Integral())
        else:
           qcdnorm[j]=qcd.Integral()
        print "k-factor", nloqcdbackup.Integral()/qcd.Integral()

        # Madgraph reference sample
        histname='alp_QCD_fa50000#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_backup"
        print histname
        if only6000:
          qcdMadgraph=inMadgraph.Get(histname.replace("6000_13000","6000_7000"))
        else:
	  qcdMadgraph=inMadgraph.Get(histname)

        # CI (=LO CI+NLO QCD)
        histname=samples[i][0].replace("Anti","")+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_backup"
        print histname
        if "EWK" in samples[i][0]:
          histname=histname.replace("_backup","")
          ci=nloqcdbackup.Clone(histname)
          for b in range(ci.GetXaxis().GetNbins()):
            low_bin=ewk.FindBin(ci.GetXaxis().GetBinLowEdge(b+1))
            up_bin=ewk.FindBin(ci.GetXaxis().GetBinUpEdge(b+1))
            correction=ewk.Integral(low_bin,up_bin-1)/(up_bin-low_bin)
            ci.SetBinContent(b+1,ci.GetBinContent(b+1)*correction)
          ci.Scale(1./ci.Integral())
        elif "lo" in samples[i][0] or "cteq66" in samples[i][0] or "cteq6ll" in samples[i][0]:
          filenamecinlo="fastnlo/NNLO/fnl5662j_"+samples[i][0].replace("QCD","").replace("nnlo","nnlo_"+muScale)+".root" # calcTheoryUncert.py from Jingyu already gives normalized signals
          print filenamecinlo
          cinlofile = TFile.Open(filenamecinlo)
          closefiles+=[cinlofile]
          filenamecinloaltscale="fastnlo/NNLO/fnl5662j_"+samples[i][0].replace("QCD","").replace("nnlo","nnlo_"+muAltScale)+".root" # calcTheoryUncert.py from Jingyu already gives normalized signals
          print filenamecinloaltscale
          cinlofilealtscale = TFile.Open(filenamecinloaltscale)
          closefiles+=[cinlofilealtscale]
          ci=None
          for k in mass_bins_nlo_list[j]:
            histname2='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])
            print histname2
            hci = TH1F(cinlofile.Get(histname2))
            #hci=rebin(hci,len(chi_binnings[j])-1,chi_binnings[j])
            hci=hci.Rebin(len(chi_binnings[j])-1,histname.replace("_backup",""),chi_binnings[j])
            if ci:
               ci.Add(hci)
            else:
               ci=hci
          cibackup=ci.Clone(ci.GetName()+"_backup")
          # apply EWK corrections
          ci.Scale(nloqcdbackup.Integral()/nloqcd.Integral())
          ci.Add(nloqcdbackup,-1)
          ci.Scale(nloqcd.Integral()/nloqcdbackup.Integral())
          ci.Add(nloqcd,+1)
          histname=histname.replace("_backup","")
          #   ci.Scale(nloqcd.Integral()/ci.Integral()/5.) # fake signal size for lower mass bins
        elif "LOCI" in samples[i][0]:
          lambdamass=samples[i][0].split("I")[-1]
          if "QCDDNLO" in samples[i][0]:
            filenamecinlo="fastnlo/cidijet_DijetChi_DILHC_2012_Lambda-"+lambdamass+"_Order-1_xmu-1.root"
          elif "QCDNLO" in samples[i][0]:
            filenamecinlo="fastnlo/cidijet_DijetChi_CILHC_2012_Lambda-"+lambdamass+"_Order-1_xmu-1.root"
          elif "QCDADLO" in samples[i][0]:
            filenamecinlo="fastnlo/cidijet_DijetChi_DILHC_2012_Lambda-"+lambdamass+"_Order-0_xmu-1.root"
          elif "QCDDLO" in samples[i][0]:
            filenamecinlo="fastnlo/cidijet_DijetChi_DILHC_2012_Lambda-"+lambdamass+"_Order-0_xmu-1.root"
          else:
            filenamecinlo="fastnlo/cidijet_DijetChi_CILHC_2012_Lambda-"+lambdamass+"_Order-0_xmu-1.root"
          print filenamecinlo
          cinlofile = TFile.Open(filenamecinlo)
          closefiles+=[cinlofile]
          histname2="chi-"+str(massbinsci[0])+"-"+str(massbinsci[1])
          print histname2
          histname=histname.replace("_backup","")
          ci = TH1F(cinlofile.Get(histname2))
          ci=ci.Rebin(len(chi_binnings[j])-1,histname.replace("_backup",""),chi_binnings[j])
          if "QCDADLO" in samples[i][0]:
            ci.Scale(-1)
          if j>=2:
             ci.Scale(1./nloqcdbackup.Integral())
          else:
             ci.Scale(nloqcd.Integral()/ci.Integral()/5.) # fake signal size for lower mass bins
          ci.Add(nloqcd)
        elif "DM" in samples[i][0] or "tripleG" in samples[i][0]:
	  if only6000:
            cibackup=out.Get(histname.replace("6000_13000","6000_7000"))
          else:
	    cibackup=out.Get(histname)
          histname=histname.replace("_backup","")
          ci=cibackup.Clone(histname)
          ci=ci.Rebin(len(chi_binnings[j])-1,ci.GetName(),chi_binnings[j])
          ci.Scale(1./nloqcdbackup.Integral())
          #if not "zprime" in samples[i][0]:
          #  ci.Scale(5./4.) #to bug fix xsec from Phil
	  print histname,"signal fraction in first bin", ci.GetBinContent(1)/nloqcd.GetBinContent(1)
          ci.Add(nloqcd)
        elif "alp" in samples[i][0]:
	  if only6000:
            cibackup=out.Get(histname.replace("6000_13000","6000_7000"))
          else:
	    cibackup=out.Get(histname)
          histname=histname.replace("_backup","")
          ci=cibackup.Clone(histname)
	  print "ALP:",ci.Integral(), "QCDMG:",qcdMadgraph.Integral(), "QCDPY:",qcd.Integral(), "QCDNLO:",nloqcdbackup.Integral()
	  ci.Add(qcdMadgraph,-1)
	  ci.Scale(nloqcdbackup.Integral()/qcdMadgraph.Integral()) # CORRECT FOR MISSING FACTOR LO->NNLO k-factor of ~10 in Madgraph QCD cross sections
          ci=ci.Rebin(len(chi_binnings[j])-1,ci.GetName(),chi_binnings[j])
          ci.Scale(1./nloqcdbackup.Integral())
	  print histname,"signal fraction in first bin", ci.GetBinContent(1)/nloqcd.GetBinContent(1)
          ci.Add(nloqcd)
        elif "QBH" in samples[i][0]:
            histname=samples[i][0]+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
            if j<3:
                ci=nloqcd.Clone(histname)
            else:
                if samples[i][0].split("_")[2]=='6':
                    histnamein='QCDADD'+samples[i][0].split("_")[2]+samples[i][0].split("_")[0]+samples[i][0].split("_")[1]+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")
                else:
                    histnamein='QCD'+samples[i][0].split("_")[2]+samples[i][0].split("_")[0]+samples[i][0].split("_")[1]+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")
                cibackup=out.Get(histnamein)
                ci=cibackup.Clone(histname)
                ci=ci.Rebin(len(chi_binnings[j])-1,ci.GetName(),chi_binnings[j])
                ci.Scale(samples[i][1][0][1]/1000000)
                ci.Scale(1./nloqcdbackup.Integral())
                ci.Add(nloqcd)
        elif "wide" in samples[i][0]:
          cibackup=out.Get(histname)
          try:
              histname=cibackup.GetName().replace("_backup","")
          except:
            print "problem reading", histname
            break
          ci=cibackup.Clone(histname)
          ci=ci.Rebin(len(chi_binnings[j])-1,ci.GetName(),chi_binnings[j])
          ci.Scale(10./nloqcdbackup.Integral()) # make in units if 10pb
          ci.Add(nloqcd)
        elif "CIplus" in samples[i][0]:
          print "CREATE FAKE SIGNAL"
          if True:#j<=3: # fake signal for signficances
            histnamealt=samples[i][0].replace("Anti","")+'#chi'+str((4800,5400)).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_backup"
          else:
            histnamealt=histname
          cibackup=out.Get(histnamealt).Clone(histname)
          histname=cibackup.GetName().replace("_backup","")
          ci=cibackup.Clone(histname)
          ci=ci.Rebin(len(chi_binnings[j])-1,ci.GetName(),chi_binnings[j])
          ci.Scale(1e9) #mb -> pb
          cinorm[j]=ci.Integral()
          #print "BBB", qcdnorm[0],cinorm[0],qcdnorm[j]/dataevents[j]*26000.,cinorm[j]/dataevents[j]*26000.,nloqcdbackup.Integral()/dataevents[j]*26000.
          # CORRECT FORMULAT
          #ci.Scale(qcdnorm[0]/cinorm[0])
          # APPROXIMATE FORMULAT
          ci.Scale(qcdnorm[j]/cinorm[j])
          ci.Add(qcd,-1.)
          #ci.Scale(cinorm[0]/qcdnorm[0]) # trusting the QCD+CI LO prediction in mb (10^9) for the LO cross section
          #print "AAA",histname,ci.Integral()/qcdnorm[j], ci.Integral()/qcdnorm[0]*cinorm[0]/nloqcdbackup.Integral()*1e6, qcdnorm[j], nloqcdbackup.Integral()/1e6
          if "Anti" in samples[i][0]:
            ci.Scale(-1.)
            histname=histname.replace("CI","AntiCI")
            ci.SetName(histname)
          #if j>5:
          #  # APPROXIMATE FORMULA
          #  #ci.Scale(1./qcdnorm[j])
          #  # CORRECT FORMULA
          #   ci.Scale(1./nloqcdbackup.Integral())
          #else:
          #  # APPROXIMATE FORMULA
          #  #ci.Scale(nloqcd.Integral()/qcdnorm[j]/5.)
          #  # CORRECT FORMULA
          ci.Scale(nloqcd.Integral()/nloqcdbackup.Integral()) # fake signal for signficances
	  if data.Integral()>0 and ci.GetBinContent(1)>0:
            scalesignal=abs(3*data.GetBinError(1)*nloqcd.Integral()/data.Integral()/ci.GetBinContent(1)) # find 2 sigma deviation
            ci.Scale(scalesignal)
          ci.Add(nloqcd)
        else:
	  if only6000:
            cibackup=out.Get(histname.replace("6000_13000","6000_7000"))
          else:  
	    cibackup=out.Get(histname)
          histname=histname.replace("_backup","")
          ci=cibackup.Clone(histname)
          ci=ci.Rebin(len(chi_binnings[j])-1,ci.GetName(),chi_binnings[j])
          ci.Scale(1e9) #mb -> pb
          #cinorm[j]=ci.Integral()
          #print "BBB", qcdnorm[0],cinorm[0],qcdnorm[j]/dataevents[j]*26000.,cinorm[j]/dataevents[j]*26000.,nloqcdbackup.Integral()/dataevents[j]*26000.
          # CORRECT FORMULAT
          #ci.Scale(qcdnorm[0]/cinorm[0])
          # APPROXIMATE FORMULAT
          #ci.Scale(qcdnorm[j]/cinorm[j])
          ci.Add(qcd,-1.)
          #ci.Scale(cinorm[0]/qcdnorm[0]) # trusting the QCD+CI LO prediction in mb (10^9) for the LO cross section
          #print "AAA",histname,ci.Integral()/qcdnorm[j], ci.Integral()/qcdnorm[0]*cinorm[0]/nloqcdbackup.Integral()*1e6, qcdnorm[j], nloqcdbackup.Integral()/1e6
          #if "Anti" in samples[i][0]:
          #  ci.Scale(-1.)
          #  histname=histname.replace("CI","AntiCI")
          #if j>=2:
          #  # APPROXIMATE FORMULA
          #  #ci.Scale(1./qcdnorm[j])
          #  # CORRECT FORMULA
          ci.Scale(1./nloqcdbackup.Integral())
          #else:
          #  # APPROXIMATE FORMULA
          #  #ci.Scale(nloqcd.Integral()/qcdnorm[j]/5.)
          #  # CORRECT FORMULA
          #ci.Scale(nloqcd.Integral()/nloqcdbackup.Integral()/5.)
          ci.Add(nloqcd)
        if ci.Integral()!=0:
          ci.Scale(dataevents[j]/ci.Integral())
        for b in range(ci.GetXaxis().GetNbins()):
            ci.SetBinError(b+1,0)
        out.cd()
        for k in range(0,200):
            out.Delete(histname+";"+str(k))
        ci.Write(histname)

        # ALT (=NLO QCD)
        histname=samples[i][0].replace("Anti","")+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
        print histname
        #if "LOCI" in samples[i][0] or "CT10" in samples[i][0] or "cteq" in samples[i][0] or "EWK" in samples[i][0]:
        alt=nloqcd.Clone(histname)
        #else:
        #    alt=out.Get(histname)
        #alt=alt.Rebin(len(chi_binnings[j])-1,alt.GetName(),chi_binnings[j])
        #alt.Add(alt,-1)
        #alt.Add(nloqcd)
        alt.Scale(dataevents[j]/alt.Integral())
        for b in range(alt.GetXaxis().GetNbins()):
            alt.SetBinError(b+1,0)
        out.cd()
        histname=samples[i][0]+'_ALT#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
        alt.SetName(histname)
        alt.SetLineColor(2)
        alt.SetTitle("")
        alt.GetYaxis().SetTitle("dN/d#chi")
        for k in range(0,200):
            out.Delete(histname+";"+str(k))
        alt.Write(histname)
        
        if injectSignal:
         data.Add(ci,1)
         data.Add(alt,-1)
         histname='data_obs#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
         for k in range(0,200):
             out.Delete(histname+";"+str(k))
         data.Write(histname)

        # Model uncertainty
	for modeln in [str(mn[0]) for mn in massbins]+[""]:
     	 histname=samples[i][0]+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
     	 clone=ci.Clone(histname)
     	 clone=clone.Rebin(len(chi_binnings[j])-1,clone.GetName(),chi_binnings[j])
     	 cimodelup=clone.Clone(histname+"_model"+modeln+"Up")
     	 cimodeldown=clone.Clone(histname+"_model"+modeln+"Down")
     	 slopes={}
     	 slopes[1200]=0.03 # FIX compute
     	 slopes[1500]=0.03 # FIX compute
     	 slopes[1900]=0.018
     	 slopes[2400]=0.018 
     	 slopes[3000]=0.020 
     	 slopes[3600]=0.020
     	 slopes[4200]=0.034
     	 slopes[4800]=0.034
     	 slopes[5400]=0.026
     	 slopes[6000]=0.026
     	 slopes[7000]=0.03 # FIX compute
     	 if modeln=="" or str(massbins[j][0])==modeln:
	   for b in range(clone.GetNbinsX()):
     	     cimodelup.SetBinContent(b+1,clone.GetBinContent(b+1)*(1.+(clone.GetBinCenter(b+1)-8.5)/7.5*slopes[massbins[j][0]]))
     	     cimodeldown.SetBinContent(b+1,clone.GetBinContent(b+1)*(1.-(clone.GetBinCenter(b+1)-8.5)/7.5*slopes[massbins[j][0]]))
     	 #cimodelup.Scale(dataevents[j]/cimodelup.Integral())
     	 #cimodeldown.Scale(dataevents[j]/cimodeldown.Integral())
     	 cimodelup.SetLineColor(3)
     	 cimodelup.SetLineStyle(2)
     	 cimodeldown.SetLineColor(3)
     	 cimodeldown.SetLineStyle(3)
     	 out.cd()
     	 for k in range(0,200):
     	     out.Delete(histname+"_model"+modeln+"Up"+";"+str(k))
     	     out.Delete(histname+"_model"+modeln+"Down"+";"+str(k))
     	 cimodelup.Write()
     	 cimodeldown.Write()

     	 histname=samples[i][0]+'_ALT#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
     	 clone=alt.Clone(histname)
     	 clone=clone.Rebin(len(chi_binnings[j])-1,clone.GetName(),chi_binnings[j])
     	 modelup=clone.Clone(histname+"_model"+modeln+"Up")
     	 modeldown=clone.Clone(histname+"_model"+modeln+"Down")
     	 if modeln=="" or str(massbins[j][0])==modeln:
     	   for b in range(clone.GetNbinsX()):
     	     modelup.SetBinContent(b+1,clone.GetBinContent(b+1)*(1.+(clone.GetBinCenter(b+1)-8.5)/7.5*slopes[massbins[j][0]]))
     	     modeldown.SetBinContent(b+1,clone.GetBinContent(b+1)*(1.-(clone.GetBinCenter(b+1)-8.5)/7.5*slopes[massbins[j][0]]))
     	 #modelup.Scale(dataevents[j]/modelup.Integral())
     	 #modeldown.Scale(dataevents[j]/modeldown.Integral())
     	 modelup.SetLineColor(3)
     	 modelup.SetLineStyle(2)
     	 modeldown.SetLineColor(3)
     	 modeldown.SetLineStyle(3)
     	 out.cd()
     	 for k in range(0,200):
     	     out.Delete(histname+"_model"+modeln+"Up"+";"+str(k))
     	     out.Delete(histname+"_model"+modeln+"Down"+";"+str(k))
     	 modelup.Write()
     	 modeldown.Write()

        # jes uncertainty
        for source in jessources:
	  if source=="SumInQuadrature":
	    jesname=""
	  else:
	    jesname=str(jessources.index(source)+1)
          histname=samples[i][0]+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          clone=ci.Clone(histname)
          clone=clone.Rebin(len(chi_binnings[j])-1,clone.GetName(),chi_binnings[j])
          cijesup=clone.Clone(histname+"_jes"+jesname+"Up")
          cijesdown=clone.Clone(histname+"_jes"+jesname+"Down")
          for b in range(clone.GetNbinsX()):
              cijesup.SetBinContent(b+1,clone.GetBinContent(b+1)*jescihists[str(j)+source][0].GetBinContent(b+1))
              cijesdown.SetBinContent(b+1,clone.GetBinContent(b+1)*jescihists[str(j)+source][1].GetBinContent(b+1))
          #cijesup.Scale(dataevents[j]/cijesup.Integral())
          #cijesdown.Scale(dataevents[j]/cijesdown.Integral())
          #print histname, "jes down", cijesdown.GetBinContent(1)/ci.GetBinContent(1)
          #print histname, "jes up", cijesup.GetBinContent(1)/ci.GetBinContent(1)
          cijesup.SetLineColor(4)
          cijesup.SetLineStyle(2)
          cijesdown.SetLineColor(4)
          cijesdown.SetLineStyle(3)
          out.cd()
          for k in range(0,200):
              out.Delete(histname+"_jes"+jesname+"Up"+";"+str(k))
              out.Delete(histname+"_jes"+jesname+"Down"+";"+str(k))
          cijesup.Write()
          cijesdown.Write()

          histname=samples[i][0]+'_ALT#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          clone=alt.Clone(histname)
          clone=clone.Rebin(len(chi_binnings[j])-1,clone.GetName(),chi_binnings[j])
          jesup=clone.Clone(histname+"_jes"+jesname+"Up")
          jesdown=clone.Clone(histname+"_jes"+jesname+"Down")
	  for b in range(clone.GetNbinsX()):
              jesup.SetBinContent(b+1,clone.GetBinContent(b+1)*jeshists[str(j)+source][0].GetBinContent(b+1))
              jesdown.SetBinContent(b+1,clone.GetBinContent(b+1)*jeshists[str(j)+source][1].GetBinContent(b+1))
          #jesup.Scale(dataevents[j]/jesup.Integral())
          #jesdown.Scale(dataevents[j]/jesdown.Integral())
          jesup.SetLineColor(4)
          jesup.SetLineStyle(2)
          jesdown.SetLineColor(4)
          jesdown.SetLineStyle(3)
          out.cd()
          for k in range(0,200):
              out.Delete(histname+"_jes"+jesname+"Up"+";"+str(k))
              out.Delete(histname+"_jes"+jesname+"Down"+";"+str(k))
          jesup.Write()
          jesdown.Write()

        # jer uncertainty
        for source in jersources:
	  if source=="SumInQuadrature":
	    jername=""
	  else:
	    jername=str(jersources.index(source)+1)
          histname=samples[i][0]+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          clone=ci.Clone(histname)
          clone=clone.Rebin(len(chi_binnings[j])-1,clone.GetName(),chi_binnings[j])
          cijerup=clone.Clone(histname+"_jer"+jername+"Up")
          cijerdown=clone.Clone(histname+"_jer"+jername+"Down")
          for b in range(clone.GetNbinsX()):
              cijerup.SetBinContent(b+1,clone.GetBinContent(b+1)*jercihists[str(j)+source][0].GetBinContent(b+1))
              cijerdown.SetBinContent(b+1,clone.GetBinContent(b+1)*jercihists[str(j)+source][1].GetBinContent(b+1))
          #cijerup.Scale(dataevents[j]/cijerup.Integral())
          #cijerdown.Scale(dataevents[j]/cijerdown.Integral())
          #print histname, "jer down", cijerdown.GetBinContent(1)/ci.GetBinContent(1)
          #print histname, "jer up", cijerup.GetBinContent(1)/ci.GetBinContent(1)
          cijerup.SetLineColor(28)
          cijerup.SetLineStyle(2)
          cijerdown.SetLineColor(28)
          cijerdown.SetLineStyle(3)
          out.cd()
          for k in range(0,200):
              out.Delete(histname+"_jer"+jername+"Up"+";"+str(k))
              out.Delete(histname+"_jer"+jername+"Down"+";"+str(k))
          cijerup.Write()
          cijerdown.Write()

          histname=samples[i][0]+'_ALT#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          clone=alt.Clone(histname)
          clone=clone.Rebin(len(chi_binnings[j])-1,clone.GetName(),chi_binnings[j])
          jerup=clone.Clone(histname+"_jer"+jername+"Up")
          jerdown=clone.Clone(histname+"_jer"+jername+"Down")
	  for b in range(clone.GetNbinsX()):
              jerup.SetBinContent(b+1,clone.GetBinContent(b+1)*jerhists[str(j)+source][0].GetBinContent(b+1))
              jerdown.SetBinContent(b+1,clone.GetBinContent(b+1)*jerhists[str(j)+source][1].GetBinContent(b+1))
          #jerup.Scale(dataevents[j]/jerup.Integral())
          #jerdown.Scale(dataevents[j]/jerdown.Integral())
          jerup.SetLineColor(28)
          jerup.SetLineStyle(2)
          jerdown.SetLineColor(28)
          jerdown.SetLineStyle(3)
          out.cd()
          for k in range(0,200):
              out.Delete(histname+"_jer"+jername+"Up"+";"+str(k))
              out.Delete(histname+"_jer"+jername+"Down"+";"+str(k))
          jerup.Write()
          jerdown.Write()

        # NLO PDFup/down
        nloPDFupqcd=None
        for k in mass_bins_nlo_list[j]:
         histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"PDFUp"
         print histname
         hnloPDFup = TH1F(nlofilesys.Get(histname))
         #hnloPDFup.Scale(float(mass_bins_nlo3[k+1]-mass_bins_nlo3[k]))
         hnloPDFup=rebin(hnloPDFup,len(chi_binnings[j])-1,chi_binnings[j])
         if nloPDFupqcd:
            nloPDFupqcd.Add(hnloPDFup)
         else:
            nloPDFupqcd=hnloPDFup
        nloPDFupqcd.Add(nloqcdnorm,-1)
        nloPDFupqcd.Scale(1./nloqcdnorm.Integral())
        nloPDFupqcd=smooth(nloPDFupqcd,"pol3") # SMOOTH NNLO PREDICTION (FIX ME)

        nloPDFdownqcd=None
        for k in mass_bins_nlo_list[j]:
         histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"PDFDown"
         print histname
         hnloPDFdown = TH1F(nlofilesys.Get(histname))
         #hnloPDFdown.Scale(float(mass_bins_nlo3[k+1]-mass_bins_nlo3[k]))
         hnloPDFdown=rebin(hnloPDFdown,len(chi_binnings[j])-1,chi_binnings[j])
         if nloPDFdownqcd:
            nloPDFdownqcd.Add(hnloPDFdown)
         else:
            nloPDFdownqcd=hnloPDFdown
        nloPDFdownqcd.Add(nloqcdnorm,-1)
        nloPDFdownqcd.Scale(1./nloqcdnorm.Integral())
        nloPDFdownqcd=smooth(nloPDFdownqcd,"pol3") # SMOOTH NNLO PREDICTION (FIX ME)

        pdfup=alt.Clone(alt.GetName()+"_pdfUp")
        pdfdown=alt.Clone(alt.GetName()+"_pdfDown")
        pdfup.Add(nloPDFupqcd,dataevents[j])
        pdfdown.Add(nloPDFdownqcd,dataevents[j])
        for b in range(pdfup.GetXaxis().GetNbins()):
            pdfup.SetBinError(b+1,0)
            pdfdown.SetBinError(b+1,0)
            if pdfup.GetBinCenter(b+1)-8.5>0:
               tmp=pdfup.GetBinContent(b+1)
               pdfup.SetBinContent(b+1,pdfdown.GetBinContent(b+1))
               pdfdown.SetBinContent(b+1,tmp)
        #pdfup.Scale(dataevents[j]/pdfup.Integral())
        #pdfdown.Scale(dataevents[j]/pdfdown.Integral())
        pdfup.SetLineColor(6)
        pdfup.SetLineStyle(2)
        pdfdown.SetLineColor(6)
        pdfdown.SetLineStyle(3)
        out.cd()
        for k in range(0,200):
            out.Delete(alt.GetName()+"_pdfUp"+";"+str(k))
            out.Delete(alt.GetName()+"_pdfDown"+";"+str(k))
        pdfup.Write()
        pdfdown.Write()

        if "lo" in samples[i][0] or "cteq66" in samples[i][0] or "cteq6ll" in samples[i][0]:
           nloPDFupci=None
           for k in mass_bins_nlo_list[j]:
            histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"PDFUp"
            print histname
            hnloPDFup = TH1F(cinlofile.Get(histname))
            hnloPDFup=hnloPDFup.Rebin(len(chi_binnings[j])-1,histname.replace("_backup",""),chi_binnings[j])
            if nloPDFupci:
               nloPDFupci.Add(hnloPDFup)
            else:
               nloPDFupci=hnloPDFup
           nloPDFupci.Add(cibackup,-1)
           nloPDFupci.Scale(1./cibackup.Integral())

           nloPDFdownci=None
           for k in mass_bins_nlo_list[j]:
            histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"PDFDown"
            print histname
            hnloPDFdown = TH1F(cinlofile.Get(histname))
            hnloPDFdown=hnloPDFdown.Rebin(len(chi_binnings[j])-1,histname.replace("_backup",""),chi_binnings[j])
            if nloPDFdownci:
              nloPDFdownci.Add(hnloPDFdown)
            else:
              nloPDFdownci=hnloPDFdown
           nloPDFdownci.Add(cibackup,-1)
           nloPDFdownci.Scale(1./cibackup.Integral())
        elif "DM" in samples[i][0] and ("Mphi_6000" in samples[i][0] or "Mphi_7000" in samples[i][0]):
           nloPDFdownci=nloPDFdownqcd.Clone("DM_pdf_down")
           nloPDFupci=nloPDFupqcd.Clone("DM_pdf_up")
	   dmpdfcanvas=dmpdffile.Get("pdf")
	   dmpdfplot=dmpdfcanvas.GetListOfPrimitives()[7] # FIX compute for all massbins [j+useUnfoldedData*1]
	   dmpdfhist=[a for a in dmpdfplot.GetListOfPrimitives() if "mean" in str(a)][0]
	   for b in range(nloPDFupci.GetNbinsX()):
	     slope=(ci.GetBinContent(b+1)/dataevents[j]-nloqcd.GetBinContent(b+1))*dmpdfhist.GetBinError(1)/dmpdfhist.GetBinContent(1)
             nloPDFupci.SetBinContent(b+1,sqrt(pow((nloPDFupci.GetBinCenter(b+1)-8.5)/7.5*slope,2)+pow(nloPDFupqcd.GetBinContent(b+1),2)))
             nloPDFdownci.SetBinContent(b+1,-sqrt(pow((nloPDFdownci.GetBinCenter(b+1)-8.5)/7.5*slope,2)+pow(nloPDFdownqcd.GetBinContent(b+1),2)))
	else:
           nloPDFdownci=nloPDFdownqcd
           nloPDFupci=nloPDFupqcd

        cipdfup=ci.Clone(ci.GetName()+"_pdfUp")
        cipdfdown=ci.Clone(ci.GetName()+"_pdfDown")
        cipdfup.Add(nloPDFupci,dataevents[j])
        cipdfdown.Add(nloPDFdownci,dataevents[j])
        for b in range(cipdfup.GetXaxis().GetNbins()):
            cipdfup.SetBinError(b+1,0)
            cipdfdown.SetBinError(b+1,0)
            if cipdfup.GetBinCenter(b+1)-8.5>0:
               tmp=cipdfup.GetBinContent(b+1)
               cipdfup.SetBinContent(b+1,cipdfdown.GetBinContent(b+1))
               cipdfdown.SetBinContent(b+1,tmp)
        #cipdfup.Scale(dataevents[j]/cipdfup.Integral())
        #cipdfdown.Scale(dataevents[j]/cipdfdown.Integral())
        cipdfup.SetLineColor(6)
        cipdfup.SetLineStyle(2)
        cipdfdown.SetLineColor(6)
        cipdfdown.SetLineStyle(3)
        out.cd()
        for k in range(0,200):
          out.Delete(ci.GetName()+"_pdfUp"+";"+str(k))
          out.Delete(ci.GetName()+"_pdfDown"+";"+str(k))
        cipdfup.Write()
        cipdfdown.Write()

        # NLO Scaleup/down
	for scaleVariation in ["MuR","MuF","","Alt"]:
          nloScaleupqcd=None
          for k in mass_bins_nlo_list[j]:
           histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale"+scaleVariation+"Up"
           print histname
	   if scaleVariation=="":
             hnloScaleup = TH1F(nlofilesys.Get(histname))
	   elif scaleVariation=="Alt":
             hnloScaleup = TH1F(nlofilealtscale.Get("CIJET_fnl5662j_cs_001_ct14nlo_0_56_30000_LL+_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-1.0_addmu")).Clone(histname)
           elif scaleVariation=="MuR":
             hnloScaleup = TH1F(nlofilesys.Get("CIJET_fnl5662j_cs_001_ct14nlo_0_56_30000_LL+_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-2.0-1.0_addmu")).Clone(histname)
	   elif scaleVariation=="MuF":
             hnloScaleup = TH1F(nlofilesys.Get("CIJET_fnl5662j_cs_001_ct14nlo_0_56_30000_LL+_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-2.0_addmu")).Clone(histname)
	   #hnloScaleup.Scale(float(mass_bins_nlo3[k+1]-mass_bins_nlo3[k]))
           hnloScaleup=rebin(hnloScaleup,len(chi_binnings[j])-1,chi_binnings[j])
           if nloScaleupqcd:
              nloScaleupqcd.Add(hnloScaleup)
           else:
              nloScaleupqcd=hnloScaleup
          nloScaleupqcd.Add(nloqcdnorm,-1)
          nloScaleupqcd.Scale(1./nloqcdnorm.Integral())
          nloScaleupqcd=smooth(nloScaleupqcd,"pol3") # SMOOTH NNLO PREDICTION (FIX ME)

          nloScaledownqcd=None
          for k in mass_bins_nlo_list[j]:
           histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale"+scaleVariation+"Down"
           print histname
           if scaleVariation=="":
             hnloScaledown = TH1F(nlofilesys.Get(histname))
	   elif scaleVariation=="Alt":
             hnloScaleup = TH1F(nlofilesys.Get("CIJET_fnl5662j_cs_001_ct14nlo_0_56_30000_LL+_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-1.0_addmu")).Clone(histname)
           elif scaleVariation=="MuR":
             hnloScaledown = TH1F(nlofilesys.Get("CIJET_fnl5662j_cs_001_ct14nlo_0_56_30000_LL+_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-0.5-1.0_addmu")).Clone(histname)
	   elif scaleVariation=="MuF":
             hnloScaledown = TH1F(nlofilesys.Get("CIJET_fnl5662j_cs_001_ct14nlo_0_56_30000_LL+_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-0.5_addmu")).Clone(histname)
	   #hnloScaledown.Scale(float(mass_bins_nlo3[k+1]-mass_bins_nlo3[k]))
           hnloScaledown=rebin(hnloScaledown,len(chi_binnings[j])-1,chi_binnings[j])
           if nloScaledownqcd:
              nloScaledownqcd.Add(hnloScaledown)
           else:
              nloScaledownqcd=hnloScaledown
          nloScaledownqcd.Add(nloqcdnorm,-1)
          nloScaledownqcd.Scale(1./nloqcdnorm.Integral())
          nloScaledownqcd=smooth(nloScaledownqcd,"pol3") # SMOOTH NNLO PREDICTION (FIX ME)

          scaleup=alt.Clone(alt.GetName()+"_scale"+scaleVariation+"Up")
          scaledown=alt.Clone(alt.GetName()+"_scale"+scaleVariation+"Down")
          scaleup.Add(nloScaleupqcd,dataevents[j])
          scaledown.Add(nloScaledownqcd,dataevents[j])
          for b in range(scaleup.GetXaxis().GetNbins()):
              scaleup.SetBinError(b+1,0)
              scaledown.SetBinError(b+1,0)
              if scaleVariation=="" and scaleup.GetBinCenter(b+1)-8.5>0:
         	 tmp=scaleup.GetBinContent(b+1)
         	 scaleup.SetBinContent(b+1,scaledown.GetBinContent(b+1))
         	 scaledown.SetBinContent(b+1,tmp)
          #scaleup.Scale(dataevents[j]/scaleup.Integral())
          #scaledown.Scale(dataevents[j]/scaledown.Integral())
          scaleup.SetLineColor(7)
          scaleup.SetLineStyle(2)
          scaledown.SetLineColor(7)
          scaledown.SetLineStyle(3)
          out.cd()
          for k in range(0,200):
              out.Delete(alt.GetName()+"_scale"+scaleVariation+"Up"+";"+str(k))
              out.Delete(alt.GetName()+"_scale"+scaleVariation+"Down"+";"+str(k))
          scaleup.Write()
          scaledown.Write()
          
          if "lo" in samples[i][0] or "cteq66" in samples[i][0] or "cteq6ll" in samples[i][0]:
             signalmassname="_".join(samples[i][0].split("_")[2:4])
	     nloScaleupci=None
             for k in mass_bins_nlo_list[j]:
              histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale"+scaleVariation+"Up"
              print histname
              if scaleVariation=="":
                hnloScaleup = TH1F(cinlofile.Get(histname))
              elif scaleVariation=="Alt":
                hnloScaleup = TH1F(cinlofilealtscale.Get("CIJET_fnl5662j_cs_001_ct14nlo_0_56_"+signalmassname+"_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-1.0_addmu")).Clone(histname)
              elif scaleVariation=="MuR":
                hnloScaleup = TH1F(cinlofile.Get("CIJET_fnl5662j_cs_001_ct14nlo_0_56_"+signalmassname+"_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-2.0-1.0_addmu")).Clone(histname)
	      elif scaleVariation=="MuF":
                hnloScaleup = TH1F(cinlofile.Get("CIJET_fnl5662j_cs_001_ct14nlo_0_56_"+signalmassname+"_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-2.0_addmu")).Clone(histname)
              hnloScaleup=hnloScaleup.Rebin(len(chi_binnings[j])-1,histname.replace("_backup",""),chi_binnings[j])
              if nloScaleupci:
         	 nloScaleupci.Add(hnloScaleup)
              else:
         	 nloScaleupci=hnloScaleup
             nloScaleupci.Add(cibackup,-1)
             nloScaleupci.Scale(1./cibackup.Integral())

             nloScaledownci=None
             for k in mass_bins_nlo_list[j]:
              histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale"+scaleVariation+"Down"
              print histname
              if scaleVariation=="":
                hnloScaledown = TH1F(cinlofile.Get(histname))
              elif scaleVariation=="Alt":
                hnloScaledown = TH1F(cinlofile.Get("CIJET_fnl5662j_cs_001_ct14nlo_0_56_"+signalmassname+"_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-1.0_addmu")).Clone(histname)
              elif scaleVariation=="MuR":
                hnloScaledown = TH1F(cinlofile.Get("CIJET_fnl5662j_cs_001_ct14nlo_0_56_"+signalmassname+"_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-0.5-1.0_addmu")).Clone(histname)
	      elif scaleVariation=="MuF":
                hnloScaledown = TH1F(cinlofile.Get("CIJET_fnl5662j_cs_001_ct14nlo_0_56_"+signalmassname+"_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-0.5_addmu")).Clone(histname)
	      hnloScaledown=hnloScaledown.Rebin(len(chi_binnings[j])-1,histname.replace("_backup",""),chi_binnings[j])
              if nloScaledownci:
         	nloScaledownci.Add(hnloScaledown)
              else:
         	nloScaledownci=hnloScaledown
             nloScaledownci.Add(cibackup,-1)
             nloScaledownci.Scale(1./cibackup.Integral())
          else:
             nloScaledownci=nloScaledownqcd
             nloScaleupci=nloScaleupqcd

          ciscaleup=ci.Clone(ci.GetName()+"_scale"+scaleVariation+"Up")
          ciscaledown=ci.Clone(ci.GetName()+"_scale"+scaleVariation+"Down")
          ciscaleup.Add(nloScaleupci,dataevents[j])
          ciscaledown.Add(nloScaledownci,dataevents[j])
          for b in range(ciscaleup.GetXaxis().GetNbins()):
              ciscaleup.SetBinError(b+1,0)
              ciscaledown.SetBinError(b+1,0)
              if scaleVariation=="" and ciscaleup.GetBinCenter(b+1)-8.5>0:
                 tmp=ciscaleup.GetBinContent(b+1)
                 ciscaleup.SetBinContent(b+1,ciscaledown.GetBinContent(b+1))
                 ciscaledown.SetBinContent(b+1,tmp)
          #ciscaleup.Scale(dataevents[j]/ciscaleup.Integral())
          #ciscaledown.Scale(dataevents[j]/ciscaledown.Integral())
          ciscaleup.SetLineColor(7)
          ciscaleup.SetLineStyle(2)
          ciscaledown.SetLineColor(7)
          ciscaledown.SetLineStyle(3)
          out.cd()
          for k in range(0,200):
            out.Delete(ci.GetName()+"_scale"+scaleVariation+"Up"+";"+str(k))
            out.Delete(ci.GetName()+"_scale"+scaleVariation+"Down"+";"+str(k))
          ciscaleup.Write()
          ciscaledown.Write()
        
        # DATA BLINDED
        #data=alt.Clone("data_blinded")
        #for b in range(data.GetXaxis().GetNbins()):
        #    data.SetBinError(b+1,sqrt(data.GetBinContent(b+1)))
        #out.cd()
        #for k in range(0,200):
        #    out.Delete('data_obs#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"+";"+str(k))
        #data.Write('data_obs#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1")
      
        # FAKE SIGNAL
        #ci=alt.Clone("fake_signal")
        #out.cd()
        #for k in range(0,200):
        #    out.Delete(samples[i][0]+'chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"+";"+str(k))
        #ci.Write(samples[i][0]+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1")
      
        # PLOTS
        #if j<3:
        #   continue
        canvas.cd(j+1)#j-2
        legend1=TLegend(0.2,0.6,0.9,0.95,(str(massbins[j][0])+"<m_{jj}<"+str(massbins[j][1])+" GeV").replace("7000<m_{jj}<13000","m_{jj}>7000"))
        legends+=[legend1]
        legend1.AddEntry(data,"data","lpe")
        alt=cloneNormalize(alt)
        plots+=[alt]
        alt.Draw("he")
        legend1.AddEntry(alt,"QCD","l")
        
        # background uncertainties
        modelup=cloneNormalize(modelup)
        plots+=[modelup]
        modelup.Draw("hesame")
        legend1.AddEntry(modelup,"model","l")
        modeldown=cloneNormalize(modeldown)
        plots+=[modeldown]
        modeldown.Draw("hesame")
        jesup=cloneNormalize(jesup)
        plots+=[jesup]
        jesup.Draw("hesame")
        legend1.AddEntry(jesup,"JES","l")
        jesdown=cloneNormalize(jesdown)
        plots+=[jesdown]
        jesdown.Draw("hesame")
        jerup=cloneNormalize(jerup)
        plots+=[jerup]
        jerup.Draw("hesame")
        legend1.AddEntry(jerup,"JER","l")
        jerdown=cloneNormalize(jerdown)
        plots+=[jerdown]
        jerdown.Draw("hesame")
        pdfup=cloneNormalize(pdfup)
        plots+=[pdfup]
        pdfup.Draw("hesame")
        legend1.AddEntry(pdfup,"PDF","l")
        pdfdown=cloneNormalize(pdfdown)
        plots+=[pdfdown]
        pdfdown.Draw("hesame")
        scaleup=cloneNormalize(scaleup)
        plots+=[scaleup]
        scaleup.Draw("hesame")
        legend1.AddEntry(scaleup,"scale","l")
        scaledown=cloneNormalize(scaledown)
        plots+=[scaledown]
        scaledown.Draw("hesame")
        
	if samples[i][0]!="QCD":
         ci=cloneNormalize(ci)
         plots+=[ci]
         ci.Draw("hesame")
         legend1.AddEntry(ci,"New Physics","l")
        
         # signal uncertainties
         cimodelup=cloneNormalize(cimodelup)
         plots+=[cimodelup]
         cimodelup.Draw("hesame")
         cimodeldown=cloneNormalize(cimodeldown)
         plots+=[cimodeldown]
         cimodeldown.Draw("hesame")
         cijesup=cloneNormalize(cijesup)
         plots+=[cijesup]
         cijesup.Draw("hesame")
         cijesdown=cloneNormalize(cijesdown)
         plots+=[cijesdown]
         cijesdown.Draw("hesame")
         cijerup=cloneNormalize(cijerup)
         plots+=[cijerup]
         cijerup.Draw("hesame")
         cijerdown=cloneNormalize(cijerdown)
         plots+=[cijerdown]
         cijerdown.Draw("hesame")
         cipdfup=cloneNormalize(cipdfup)
         plots+=[cipdfup]
         cipdfup.Draw("hesame")
         cipdfdown=cloneNormalize(cipdfdown)
         plots+=[cipdfdown]
         cipdfdown.Draw("hesame")
         ciscaleup=cloneNormalize(ciscaleup)
         plots+=[ciscaleup]
         ciscaleup.Draw("hesame")
         ciscaledown=cloneNormalize(ciscaledown)
         plots+=[ciscaledown]
         ciscaledown.Draw("hesame")
 
        origdata=data
	datahist[j]=data
        dataplot[j]=TGraphAsymmErrors(cloneNormalize(data))
        plots+=[dataplot[j]]
        alpha=1.-0.6827
        for b in range(dataplot[j].GetN()):
            if useUnfoldedData:
              N=1./pow(origdata.GetBinError(b+1)/origdata.GetBinContent(b+1),2)
            else:
              N=origdata.GetBinContent(b+1)
            L=0
            if N>0:
              L=ROOT.Math.gamma_quantile(alpha/2.,N,1.)
            U=ROOT.Math.gamma_quantile_c(alpha/2.,N+1,1.)
            dataplot[j].SetPointEYlow(b,(N-L)/origdata.GetBinWidth(b+1))
            dataplot[j].SetPointEYhigh(b,(U-N)/origdata.GetBinWidth(b+1))
        dataplot[j].SetLineColor(1)
        dataplot[j].SetMarkerStyle(24)
        dataplot[j].SetMarkerSize(0.5)
        dataplot[j].Draw("pe0zsame")

	miny=min(cloneNormalize(datahist[j]).GetMinimum(),alt.GetMinimum())
	maxy=max(cloneNormalize(datahist[j]).GetMaximum(),alt.GetMaximum())
        alt.GetYaxis().SetRangeUser(miny*0.8,(maxy-miny)*2.+miny)
        
        legend1.SetTextSize(0.04)
        legend1.SetFillStyle(0)
        legend1.Draw("same")

      canvas.SaveAs(prefix + "_"+samples[i][0].replace("QCD","") + '_sys_run2.pdf')
      #canvas.SaveAs(prefix + "_"+samples[i][0].replace("QCD","") + '_sys_run2.eps')
      
      if not useUnfoldedData:
        canvas = TCanvas("","",0,0,800,600)
        canvas.Divide(4,3)
	print "Applying response matrix"
        matrix=TFile("datacards/responseMatrices_20190930.root",'READ')
        matrix1=matrix.Get("TMatrixT<double>;2") #low mass chi binning
        matrix2=matrix.Get("TMatrixT<double>;3") #highest mass chi binning
        matrixMassBins=[1000,1200,1500,1900,2400,3000,3600,4200,4800,5400,6000,7000,13000]
        #print (len(matrixMassBins)-1)*(len(chi_bins[0])-1)
        #matrix1.Print()
        #print (len(matrixMassBins)-1)*(len(chi_bins[-1])-1)
        #matrix2.Print()

        # modify histograms
        histogram_names=[
          (samples[i][0]+'_ALT#chi',"_rebin1"),
          (samples[i][0]+'#chi',"_rebin1"),
        ]
        for pre,post in histogram_names:
          althists=[]
          althistsclones=[]
          for j in range(len(massbins)):
            histname=pre+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+post
	    print histname
            althists+=[out.Get(histname)]
            althistsclones+=[althists[-1].Clone(althists[-1].GetName()+"original")]
            althists[-1].Scale(0)
          countgen=0
          for j1 in range(len(massbins)):
            for b1 in range(althists[j1].GetNbinsX()):
              for j2 in range(len(massbins)):
                for b2 in range(althists[j2].GetNbinsX()):
		  #if abs(j1-j2)>1 or j2<j1: continue #don't take elements too far from the diagonal to avoid statistical fluctuations
		  # map massbins on matrixMassBins
		  j1m=j1+1
		  j2m=j2+1
                  response=0
                  if j1==(len(massbins)-1) and j2==(len(massbins)-1):
                    # both in highest mass bin
                    response=matrix2[j2m+b2*(len(matrixMassBins)-1)][j1m+b1*(len(matrixMassBins)-1)]
                  elif j1!=(len(massbins)-1) and j2!=(len(massbins)-1):
                    # both in lower mass bins
                    response=matrix1[j2m+b2*(len(matrixMassBins)-1)][j1m+b1*(len(matrixMassBins)-1)]
                  elif j1==(len(massbins)-1):
                    # one in highest, one in lower mass bin
		    for bin1 in range(althists[0].GetNbinsX()):
                      if althists[0].GetXaxis().GetBinUpEdge(bin1+1)>althists[j1].GetXaxis().GetBinLowEdge(b1+1) and\
                         althists[0].GetXaxis().GetBinLowEdge(bin1+1)<althists[j1].GetXaxis().GetBinUpEdge(b1+1):
                        response+=matrix1[j2m+b2*(len(matrixMassBins)-1)][j1m+bin1*(len(matrixMassBins)-1)]/\
			    (althists[j1].GetXaxis().GetBinUpEdge(b1+1)-althists[j1].GetXaxis().GetBinLowEdge(b1+1))*\
			    (althists[0].GetXaxis().GetBinUpEdge(bin1+1)-althists[0].GetXaxis().GetBinLowEdge(bin1+1))
                  elif j2==(len(massbins)-1):
                    # one in highest, one in lower mass bin
		    for bin2 in range(althists[0].GetNbinsX()):
                      if althists[0].GetXaxis().GetBinUpEdge(bin2+1)>althists[j2].GetXaxis().GetBinLowEdge(b2+1) and\
                         althists[0].GetXaxis().GetBinLowEdge(bin2+1)<althists[j2].GetXaxis().GetBinUpEdge(b2+1):
                        response+=matrix1[j2m+bin2*(len(matrixMassBins)-1)][j1m+b1*(len(matrixMassBins)-1)]
		  #response=((b1==b2) and (j1==j2))
		  althists[j2].Fill(althists[j2].GetBinCenter(b2+1),althistsclones[j1].GetBinContent(b1+1)*response)
          for j in range(len(massbins)):
            print dataevents[j],althistsclones[j].Integral(),althists[j].Integral()
            althists[j].Scale(dataevents[j]/althists[j].Integral())
	    for b in range(althists[j].GetXaxis().GetNbins()):
	      althists[j].SetBinError(b+1,0)
            #print althists[j].GetBinContent(1),althistsclones[j].GetBinContent(1)
          out.cd()
          for hist in althists:
            for k in range(0,200):
              out.Delete(hist.GetName()+";"+str(k))
            hist.Write()
 	  for j in range(len(massbins)):
            canvas.cd(j+1)#j-2
            alt=cloneNormalize(althists[j])
	    if "ALT" in pre:
	      alt.Draw("he")
              miny=min(cloneNormalize(datahist[j]).GetMinimum(),alt.GetMinimum())
              maxy=max(cloneNormalize(datahist[j]).GetMaximum(),alt.GetMaximum())
              alt.GetYaxis().SetRangeUser(miny*0.8,(maxy-miny)*2.+miny)
            else:
              if samples[i][0]!="QCD":
	       alt.Draw("hesame")
            plots+=[alt]
	  syshists={}
	  for j in range(len(massbins)):
	    syss=["model","pdf"]
	    syss+=["model"+str(mn[0]) for mn in massbins]
	    skipInSum=["model"+str(mn[0]) for mn in massbins]
	    syss+=["jes"]
	    for n in range(len(jessources)-1):
	      syss+=["jes"+str(n+1)]
	      skipInSum+=["jes"+str(n+1)]
	    syss+=["jer"]
	    for n in range(len(jersources)-1):
	      syss+=["jer"+str(n+1)]
	      skipInSum+=["jer"+str(n+1)]
	    syss+=["scaleAlt"]
            syss+=["scale","scaleMuR","scaleMuF"]
	    skipInSum+=["scale","scaleMuR","scaleMuF"]
	    for sys in syss:
	      for shift in ["Up","Down"]:
	        histname=althists[j].GetName()+"_"+sys+shift
		#print histname
 	        sysHist=out.Get(histname)
	        sysNorm=sysHist.Integral()
	        for b in range(sysHist.GetNbinsX()):
		  if althistsclones[j].GetBinContent(b+1)*althists[j].GetBinContent(b+1)>0:
	            sysHist.SetBinContent(b+1,sysHist.GetBinContent(b+1)/althistsclones[j].GetBinContent(b+1)*althists[j].GetBinContent(b+1))
                if sysHist.Integral()>0:
		  sysHist.Scale(sysNorm/althistsclones[j].Integral()*althists[j].Integral()/sysHist.Integral())
	        #sysHist.Scale(dataevents[j]/sysHist.Integral())
                for k in range(0,200):
                  out.Delete(sysHist.GetName()+";"+str(k))
	        sysHist.Write()
		syshists[sys+shift]=sysHist
                canvas.cd(j+1)#j-2
                alt=cloneNormalize(sysHist)
	        if (samples[i][0]!="QCD"or "ALT" in pre) and not sys in skipInSum:
                 alt.Draw("hesame")
                plots+=[alt]
            dataplot[j].Draw("pe0zsame")
	    legend=legend1.Clone()
	    legend.SetHeader((str(massbins[j][0])+"<m_{jj}<"+str(massbins[j][1])+" GeV").replace("7000<m_{jj}<13000","m_{jj}>7000"))
	    legend.Draw("same")
	    plots+=[legend]
	    chi2=0
	    for b in range(datahist[j].GetNbinsX()):
	      if datahist[j].GetBinContent(b+1)==0: continue
	      unc2=max(pow(dataplot[j].GetErrorYlow(b),2),pow(dataplot[j].GetErrorYhigh(b),2))
	      print "stat",sqrt(unc2)/datahist[j].GetBinContent(b+1)
	      for sys in syss:
	       if not sys in skipInSum:
	        unc2+=max(pow(syshists[sys+"Up"].GetBinContent(b+1)-althists[j].GetBinContent(b+1),2),pow(syshists[sys+"Down"].GetBinContent(b+1)-althists[j].GetBinContent(b+1),2))
	      print "stat+sys",sqrt(unc2)/datahist[j].GetBinContent(b+1)
	      chi2+=pow(datahist[j].GetBinContent(b+1)-althists[j].GetBinContent(b+1),2)/unc2
	    pvalue=stats.distributions.chi2.sf(chi2, datahist[j].GetNbinsX())
            sign=-stats.norm.ppf(pvalue)
	    print "sign",sign,"chi2/ndof",chi2/datahist[j].GetNbinsX()
        canvas.SaveAs(prefix + "_"+samples[i][0].replace("QCD","") + '_sys_run2_smear.pdf')
        #canvas.SaveAs(prefix + "_"+samples[i][0].replace("QCD","") + '_sys_run2_smear.eps')

      for closefile in closefiles:
          closefile.Close()
