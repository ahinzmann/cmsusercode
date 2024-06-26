import os, sys
from ROOT import * 
import array
from math import *

#gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
gROOT.Reset()
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
gStyle.SetNdivisions(506, "XYZ")
gStyle.SetLegendBorderSize(0)
gROOT.SetBatch(True)

if __name__ == '__main__':

  files=[]
  for prefix in  ["jet_mass_uzh2_76W_","jet_mass_uzh2_76Z_","jet_mass_uzh2_76H_","jet_mass_uzh2_W_","jet_mass_uzh2_Z_","jet_mass_uzh2_H_","jet_mass_uzh2_WZH_"]:
  #for prefix in  ["jet_mass_uzh2_WZH_"]:
    if "76W" in prefix:    
     samples=[("W",[
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/ZprimeToWW_narrow_M-1800_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_ZprimeToWW_narrow_M-1800_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/ZprimeToWW_narrow_M-2000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_ZprimeToWW_narrow_M-2000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/ZprimeToWW_narrow_M-2500_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_ZprimeToWW_narrow_M-2500_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/ZprimeToWW_narrow_M-3000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_ZprimeToWW_narrow_M-3000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/ZprimeToWW_narrow_M-3500_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_ZprimeToWW_narrow_M-3500_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/ZprimeToWW_narrow_M-4000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_ZprimeToWW_narrow_M-4000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/ZprimeToWW_narrow_M-4500_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_ZprimeToWW_narrow_M-4500_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/ZprimeToWW_narrow_M-600_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_ZprimeToWW_narrow_M-600_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/ZprimeToWW_narrow_M-800_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_ZprimeToWW_narrow_M-800_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/BulkGravToWW_narrow_M-1200_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_BulkGravToWW_narrow_M-1200_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/BulkGravToWW_narrow_M-1400_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_BulkGravToWW_narrow_M-1400_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/BulkGravToWW_narrow_M-1600_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_BulkGravToWW_narrow_M-1600_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/BulkGravToWW_narrow_M-1800_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_BulkGravToWW_narrow_M-1800_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
                    ]),
            ]
     ptbins=[(250,350),
            (350,450),
            (550,650),
	    #(600,800),
            (800,1200),
	    (1400,1600),
            (1600,2400),
            #(2500,3500),
            #(3500,4500),
           ]
    elif "76Z" in prefix:    
     samples=[("Z",[
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RadionToZZ_narrow_M-1000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RadionToZZ_narrow_M-1000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RadionToZZ_narrow_M-1200_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RadionToZZ_narrow_M-1200_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RadionToZZ_narrow_M-1400_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RadionToZZ_narrow_M-1400_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RadionToZZ_narrow_M-1600_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RadionToZZ_narrow_M-1600_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RadionToZZ_narrow_M-1800_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RadionToZZ_narrow_M-1800_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RadionToZZ_narrow_M-2000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RadionToZZ_narrow_M-2000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RadionToZZ_narrow_M-2500_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RadionToZZ_narrow_M-2500_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RadionToZZ_narrow_M-3000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RadionToZZ_narrow_M-3000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RadionToZZ_narrow_M-3500_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RadionToZZ_narrow_M-3500_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RadionToZZ_narrow_M-4000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RadionToZZ_narrow_M-4000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RadionToZZ_narrow_M-4500_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RadionToZZ_narrow_M-4500_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RadionToZZ_narrow_M-600_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RadionToZZ_narrow_M-600_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RadionToZZ_narrow_M-800_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RadionToZZ_narrow_M-800_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
                    ]),
            ]
     ptbins=[(250,350),
            (350,450),
            (550,650),
	    #(600,800),
            (800,1200),
	    (1400,1600),
            (1600,2400),
            #(2500,3500),
            #(3500,4500),
           ]
    elif "76H" in prefix:    
     samples=[("H",[
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RSGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RSGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RSGravTohhTohbbhbb_narrow_M-1200_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RSGravTohhTohbbhbb_narrow_M-1200_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RSGravTohhTohbbhbb_narrow_M-1600_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RSGravTohhTohbbhbb_narrow_M-1600_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RSGravTohhTohbbhbb_narrow_M-1800_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RSGravTohhTohbbhbb_narrow_M-1800_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RSGravTohhTohbbhbb_narrow_M-2000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RSGravTohhTohbbhbb_narrow_M-2000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RSGravTohhTohbbhbb_narrow_M-2500_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RSGravTohhTohbbhbb_narrow_M-2500_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RSGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RSGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RSGravTohhTohbbhbb_narrow_M-3500_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RSGravTohhTohbbhbb_narrow_M-3500_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RSGravTohhTohbbhbb_narrow_M-4000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RSGravTohhTohbbhbb_narrow_M-4000_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall15/RSGravTohhTohbbhbb_narrow_M-4500_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/EXOVVTree_RSGravTohhTohbbhbb_narrow_M-4500_13TeV-madgraph_RunIIFall15MiniAODv2_PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1_1.root",
                    ]),
            ]
     ptbins=[#(250,350),
            #(350,450),
            (550,650),
	    #(600,800),
            (800,1200),
	    (1400,1600),
            (1600,2400),
            #(2500,3500),
            #(3500,4500),
           ]

    elif "WZH" in prefix:    
     samples=[("W",["dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1000_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1200_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1400_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1600_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1800_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-2000_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-2500_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-3000_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-3500_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-4000_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-4500_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-600_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-800_13TeV-madgraph_1.root",
                    ]),
            ("Z",["dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-1000_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-1200_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-1400_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-1600_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-1800_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-2000_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-2500_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-3000_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-3500_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-4000_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-4500_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-600_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-800_13TeV-madgraph_1.root",
                    ]),
            ("H",["dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-1200_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-1400_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-1600_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-1800_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-2000_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-2500_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-3500_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-4000_13TeV-madgraph_1.root",
                    #"dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-4500_13TeV-madgraph_1.root",
                    ]),
            ("q/g",["dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_2.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_3.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_4.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_5.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_6.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_7.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_8.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_9.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_10.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_11.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_12.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_13.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_14.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_15.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_16.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_17.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_18.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_19.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_20.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_21.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_22.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_23.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_24.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_25.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_26.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_27.root",
                    ]),
            ]
     ptbins=[#(250,350),
            #(350,450),
            (450,550),
	    #(600,800),
	    #(700,1100),
            (800,1200),
	    #(1400,1600),
            #(1600,2400),
            #(2500,3500),
            #(3500,4500),
           ]

    elif "W" in prefix:    
     samples=[("W",["dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1000_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1200_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1400_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1600_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1800_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-2000_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-2500_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-3000_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-3500_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-4000_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-4500_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-600_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-800_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1000_23TeV-madgraph_2.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1200_23TeV-madgraph_2.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1400_23TeV-madgraph_2.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1600_23TeV-madgraph_2.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1800_23TeV-madgraph_2.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-2000_23TeV-madgraph_2.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-2500_23TeV-madgraph_2.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-3000_23TeV-madgraph_2.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-3500_23TeV-madgraph_2.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-4000_23TeV-madgraph_2.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-4500_23TeV-madgraph_2.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-600_23TeV-madgraph_2.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-800_23TeV-madgraph_2.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1000_33TeV-madgraph_3.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1200_33TeV-madgraph_3.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1400_33TeV-madgraph_3.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1600_33TeV-madgraph_3.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1800_33TeV-madgraph_3.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-2000_33TeV-madgraph_3.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-2500_33TeV-madgraph_3.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-3000_33TeV-madgraph_3.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-3500_33TeV-madgraph_3.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-4000_33TeV-madgraph_3.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-4500_33TeV-madgraph_3.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-600_33TeV-madgraph_3.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-800_33TeV-madgraph_3.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1000_43TeV-madgraph_4.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1200_43TeV-madgraph_4.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1400_43TeV-madgraph_4.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1600_43TeV-madgraph_4.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-1800_43TeV-madgraph_4.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-2000_43TeV-madgraph_4.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-2500_43TeV-madgraph_4.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-3000_43TeV-madgraph_4.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-3500_43TeV-madgraph_4.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-4000_43TeV-madgraph_4.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-4500_43TeV-madgraph_4.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-600_43TeV-madgraph_4.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToWW_narrow_M-800_43TeV-madgraph_4.root",
                    ]),
            #("q/g",["dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_1.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_2.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_3.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_4.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_5.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_6.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_7.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_8.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_9.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_10.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_11.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_12.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_13.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_14.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_15.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_16.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_17.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_18.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_19.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_20.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_21.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_22.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_23.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_24.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_25.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_26.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_27.root",
            #        ]),
            ]
     ptbins=[(250,350),
            (350,450),
            (550,650),
	    #(600,800),
            (800,1200),
	    (1400,1600),
            (1600,2400),
            #(2500,3500),
            #(3500,4500),
           ]
    elif "Z" in prefix:    
     samples=[("Z",["dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-1000_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-1200_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-1400_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-1600_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-1800_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-2000_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-2500_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-3000_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-3500_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-4000_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-4500_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-600_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravToZZToZlepZhad_narrow_M-800_13TeV-madgraph_1.root",
                    ]),
            #("q/g",["dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_1.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_2.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_3.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_4.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_5.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_6.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_7.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_8.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_9.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_10.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_11.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_12.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_13.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_14.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_15.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_16.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_17.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_18.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_19.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_20.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_21.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_22.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_23.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_24.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_25.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_26.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_27.root",
            #        ]),
            ]
     ptbins=[(250,350),
            (350,450),
            (550,650),
	    #(600,800),
            (800,1200),
	    (1400,1600),
            (1600,2400),
            #(2500,3500),
            #(3500,4500),
           ]
    elif "H" in prefix:    
     samples=[("H",["dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-1000_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-1200_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-1400_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-1600_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-1800_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-2000_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-2500_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-3000_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-3500_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-4000_13TeV-madgraph_1.root",
                    "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_BulkGravTohhTohbbhbb_narrow_M-4500_13TeV-madgraph_1.root",
                    ]),
            #("q/g",["dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_1.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_2.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_3.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_4.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_5.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_6.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_7.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_8.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_9.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_10.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_11.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_12.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_13.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_14.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_15.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_16.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_17.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_18.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_19.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_20.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_21.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_22.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_23.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_24.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_25.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_26.root",
            #        "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Spring15/WtaggingWithPuppi/EXOVVTree_QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_27.root",
            #        ]),
            ]
     ptbins=[#(250,350),
            #(350,450),
            (550,650),
	    #(600,800),
            (800,1200),
	    (1400,1600),
            (1600,2400),
            #(2500,3500),
            #(3500,4500),
           ]

    variables=["prunedMass","prunedMassL23Corrected","prunedMassL123Corrected",
    "softdropMass","softdropMassL23Corrected","softdropMassL123Corrected",
    "prunedMassTau21","softdropMassTau21",
    "tau21","puppiTau21",
    "trimmedMass","trimmedMassL23Corrected","trimmedMassL123Corrected",
    "puppiPrunedMass","puppiPrunedMassL23Corrected","puppiPrunedMassL123Corrected",
    "puppiSoftdropMass","puppiSoftdropMassL23Corrected","puppiSoftdropMassL123Corrected",
    "D2","JEC",
    "genPrunedMass","genSoftDropMass",
    #"alpacaPrunedMass","alpacaSoftdropMass","alpacaTau21",
    #"puppiAlpacaPrunedMass","puppiAlpacaSoftdropMass","puppiAlpacaTau21",
    ]
    
    #samples=[("W",["flatTuple.root"])]
    
    color=[1,2,4,6,7,8,9,11,40,41]
    
    plots=[]
    for sample,filenames in samples:
      for ptbin in ptbins:
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#prunedMass',';pruned mass;N',60,0,300)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#prunedMassL23Corrected',';L23-corrected pruned mass;N',60,0,300)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#prunedMassL123Corrected',';L123-corrected pruned mass;N',60,0,300)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#softdropMass',';softdrop mass;N',60,0,300)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#softdropMassL23Corrected',';L23-corrected softdrop mass;N',60,0,300)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#softdropMassL123Corrected',';L123-corrected softdrop mass;N',60,0,300)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#prunedMassTau21',';L23-corrected pruned mass;N',60,0,300)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#softdropMassTau21',';L23-corrected softdrop mass;N',60,0,300)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#tau21',';#tau_{2}/#tau_{1};N',50,0,2)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#puppiTau21',';Puppi #tau_{2}/#tau_{1};N',50,0,2)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#trimmedMass',';trimmed mass;N',60,0,300)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#trimmedMassL23Corrected',';L23-corrected trimmed mass;N',60,0,300)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#trimmedMassL123Corrected',';L123-corrected trimmed mass;N',60,0,300)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#puppiPrunedMass',';Puppi pruned mass;N',60,0,300)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#puppiPrunedMassL23Corrected',';Puppi+L23-corrected pruned mass;N',60,0,300)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#puppiPrunedMassL123Corrected',';Puppi+L123-corrected pruned mass;N',60,0,300)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#puppiSoftdropMass',';Puppi softdrop mass;N',60,0,300)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#puppiSoftdropMassL23Corrected',';Puppi+L23-corrected softdrop mass;N',60,0,300)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#puppiSoftdropMassL123Corrected',';Puppi+L123-corrected softdrop mass;N',60,0,300)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#D2',';D2;N',50,0,10)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#JEC',';JEC;N',50,1.0,1.2)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#genPrunedMass',';generator pruned mass;N',60,0,300)]
        plots += [TH1F(prefix+sample+str(ptbin[0])+'#genSoftdropMass',';generator softdrop mass;N',60,0,300)]
        #plots += [TH1F(prefix+sample+str(ptbin[0])+'#alpacaPrunedMass',';Alpaca pruned mass;N',60,0,300)]
        #plots += [TH1F(prefix+sample+str(ptbin[0])+'#alpacaSoftdropMass',';Alpaca softdrop mass;N',60,0,300)]
        #plots += [TH1F(prefix+sample+str(ptbin[0])+'#alpacaTau21',';Alpaca #tau_{2}/#tau_{1};N',50,0,2)]
        #plots += [TH1F(prefix+sample+str(ptbin[0])+'#puppiAlpacaPrunedMass',';Puppi+Alpaca pruned mass;N',60,0,300)]
        #plots += [TH1F(prefix+sample+str(ptbin[0])+'#puppiAlpacaSoftdropMass',';Puppi+Alpaca softdrop mass;N',60,0,300)]
        #plots += [TH1F(prefix+sample+str(ptbin[0])+'#puppiAlpacaTau21',';Puppi+Alpaca #tau_{2}/#tau_{1};N',50,0,2)]
    for plot in plots:
        plot.Sumw2()
    for k in range(len(samples)):
      sample=samples[k][0]
      filenames=samples[k][1]
      for f in filenames[:]:
     	  print f
     	  fil=TFile.Open(f)
	  files+=[fil]
	  try:
     	    events=fil.Get("ntuplizer/tree")
            nevents=events.GetEntries()
	  except:
	    print "error opening"
	    continue
          event_count=0
     	  print sample,nevents
     	  for event in events:
     	    event_count+=1
	    #if event_count>1000: break
     	    if event_count%10000==1: print "event",event_count
     	    if len(event.jetAK8_pt)<1 or not event.jetAK8_IDTight[0]: continue
            for j in range(len(ptbins)):
	     offsetindex=len(variables)*j+len(variables)*len(ptbins)*k
	     if event.jetAK8_pt[0]>ptbins[j][0] and event.jetAK8_pt[0]<ptbins[j][1]:
              plots[0+offsetindex].Fill(event.jetAK8_pruned_mass[0])
              plots[1+offsetindex].Fill(event.jetAK8_pruned_massCorr[0])
              plots[2+offsetindex].Fill(event.jetAK8_pruned_mass[0]*event.jetAK8_jec[0])
              plots[3+offsetindex].Fill(event.jetAK8_softdrop_mass[0])
              plots[4+offsetindex].Fill(event.jetAK8_softdrop_massCorr[0])
              plots[5+offsetindex].Fill(event.jetAK8_softdrop_mass[0]*event.jetAK8_jec[0])
	      if event.jetAK8_tau1[0]>0 and event.jetAK8_tau2[0]/event.jetAK8_tau1[0]<0.5:
                plots[6+offsetindex].Fill(event.jetAK8_pruned_massCorr[0])
                plots[7+offsetindex].Fill(event.jetAK8_softdrop_massCorr[0])
              if event.jetAK8_tau1[0]>0 and event.jetAK8_pruned_massCorr[0]>65 and event.jetAK8_pruned_massCorr[0]<95:
	        plots[8+offsetindex].Fill(event.jetAK8_tau2[0]/event.jetAK8_tau1[0])
              if event.jetAK8_puppi_tau1[0]>0 and event.jetAK8_pruned_massCorr[0]>65 and event.jetAK8_pruned_massCorr[0]<95:
	        plots[9+offsetindex].Fill(event.jetAK8_puppi_tau2[0]/event.jetAK8_puppi_tau1[0])
              plots[10+offsetindex].Fill(event.jetAK10_trimmed_mass[0])
              plots[11+offsetindex].Fill(event.jetAK10_trimmed_massCorr[0])
              plots[12+offsetindex].Fill(event.jetAK10_trimmed_mass[0]*event.jetAK8_jec[0])
              #plots[13+offsetindex].Fill(event.jetAK8_puppi_pruned_mass[0])
              #plots[14+offsetindex].Fill(event.jetAK8_puppi_pruned_massCorr[0])
              #plots[15+offsetindex].Fill(event.jetAK8_puppi_pruned_mass[0]*event.jetAK8_jec[0])
              plots[16+offsetindex].Fill(event.jetAK8_puppi_softdrop_mass[0])
              plots[17+offsetindex].Fill(event.jetAK8_puppi_softdrop_massCorr[0])
              plots[18+offsetindex].Fill(event.jetAK8_puppi_softdrop_mass[0]*event.jetAK8_jec[0])
	      if event.jetAK10_ecf2[0]>0 and event.jetAK10_trimmed_massCorr[0]>65 and event.jetAK10_trimmed_massCorr[0]<95:
                plots[19+offsetindex].Fill(event.jetAK10_ecf3[0]*pow(event.jetAK10_ecf1[0],3)/pow(event.jetAK10_ecf2[0],3))
              plots[20+offsetindex].Fill(event.jetAK8_pruned_jec[0])
              #plots[21+offsetindex].Fill(event.genJetAK8_prunedmass[0])
              #plots[22+offsetindex].Fill(event.genJetAK8_softdropmass[0])
#      for f in filenames[:]:
#     	  print f.replace("Puppi","Alpaca")
#     	  fil=TFile.Open(f.replace("Puppi","Alpaca"))
#	  files+=[fil]
#	  try:
#     	    events=fil.Get("ntuplizer/tree")
#            nevents=events.GetEntries()
#	  except:
#	    print "error opening"
#	    continue
#          event_count=0
#     	  print sample,nevents
#     	  for event in events:
#     	    event_count+=1
#	    if event_count>1000: break
#     	    if event_count%10000==1: print "event",event_count
#     	    if len(event.jetAK8_pt)<1 or not event.jetAK8_IDTight[0]: continue
#            for j in range(len(ptbins)):
#	     offsetindex=len(variables)*j+len(variables)*len(ptbins)*k
#	     if event.jetAK8_pt[0]>ptbins[j][0] and event.jetAK8_pt[0]<ptbins[j][1]:
#              #plots[23+offsetindex].Fill(event.jetAK8_puppi_pruned_mass[0])
#              plots[24+offsetindex].Fill(event.jetAK8_puppi_softdrop_mass[0])
#              if event.jetAK8_puppi_tau1[0]>0 and event.jetAK8_pruned_massCorr[0]>65 and event.jetAK8_pruned_massCorr[0]<95:
#	        plots[25+offsetindex].Fill(event.jetAK8_puppi_tau2[0]/event.jetAK8_puppi_tau1[0])
#      for f in filenames[:]:
#     	  print f.replace("Puppi","AlpacaPuppi")
#     	  fil=TFile.Open(f.replace("Puppi","AlpacaPuppi"))
#	  files+=[fil]
#	  try:
#     	    events=fil.Get("ntuplizer/tree")
#            nevents=events.GetEntries()
#	  except:
#	    print "error opening"
#	    continue
#          event_count=0
#     	  print sample,nevents
#     	  for event in events:
#     	    event_count+=1
#	    #if event_count>1000: break
#     	    if event_count%10000==1: print "event",event_count
#     	    if len(event.jetAK8_pt)<1 or not event.jetAK8_IDTight[0]: continue
#            for j in range(len(ptbins)):
#	     offsetindex=len(variables)*j+len(variables)*len(ptbins)*k
#	     if event.jetAK8_pt[0]>ptbins[j][0] and event.jetAK8_pt[0]<ptbins[j][1]:
#              #plots[26+offsetindex].Fill(event.jetAK8_puppi_pruned_mass[0])
#              plots[27+offsetindex].Fill(event.jetAK8_puppi_softdrop_mass[0])
#              if event.jetAK8_puppi_tau1[0]>0 and event.jetAK8_pruned_massCorr[0]>65 and event.jetAK8_pruned_massCorr[0]<95:
#	        plots[28+offsetindex].Fill(event.jetAK8_puppi_tau2[0]/event.jetAK8_puppi_tau1[0])

    print plots
    for i in range(len(variables)):
        canvas = TCanvas(variables[i],"",0,0,200,200)
        if "ass" in variables[i]:
          legend1=TLegend(0.55,0.4,0.95,0.9,"sample, p_{T}, mean, width")
        else:
	  legend1=TLegend(0.55,0.4,0.95,0.9,"sample, p_{T}, mean")
        for k in range(len(samples)):
          peaks=[]
          sample=samples[k][0]
          for j in range(len(ptbins)):
	    offsetindex=len(variables)*j+len(variables)*len(ptbins)*k
            plots[i+offsetindex].SetLineColor(color[j])
            plots[i+offsetindex].SetLineStyle(k+1)
            plots[i+offsetindex].SetLineWidth(k+1)
	    if j==0 and k==0:
              plots[i+offsetindex].Draw("he")
	      integral=plots[i+offsetindex].Integral()
	    else:
	      if plots[i+offsetindex].Integral()>0:
	        plots[i+offsetindex].Scale(integral/plots[i+offsetindex].Integral())
	      plots[i+offsetindex].Draw("hesame")
	    if (sample=="W" or sample=="Z" or sample=="H") and "ass" in variables[i]:
	     #N = 1
             #res = array.array("d", [0.]*N)
             #q = array.array("d", [0.5])
             #plots[i+offsetindex].GetQuantiles(N, res, q)
	     #mean=res[0]
	     maxbin=0
	     maxcontent=0
	     for b in range(plots[i+offsetindex].GetXaxis().GetNbins()):
	       if plots[i+offsetindex].GetXaxis().GetBinCenter(b+1)>50 and plots[i+offsetindex].GetBinContent(b+1)>maxcontent:
	          maxbin=b
		  maxcontent=plots[i+offsetindex].GetBinContent(b+1)
	     mean=plots[i+offsetindex].GetXaxis().GetBinCenter(maxbin)
	     g1 = TF1("g1","gaus", mean-11.,mean+11.)
	     plots[i+offsetindex].Fit(g1, "R0")
	     if i==1:
	      print g1.GetParameter(1)
	      print plots[i+offsetindex].Integral(plots[i+offsetindex].FindBin(65),plots[i+offsetindex].FindBin(100)),plots[i+offsetindex].Integral(plots[i+offsetindex].FindBin(60),plots[i+offsetindex].FindBin(105))
             hmed = int(g1.GetParameter(1)*10.)/10.
             hres = int(g1.GetParameter(2)*10.)/10.
	    else:
	     hmed = int(plots[i+offsetindex].GetMean()*1000.)/1000.
	    peaks+=[hmed]
	    entry=sample+", "+str(int(ptbins[j][0]+ptbins[j][1])/2)
	    entry+=", "+str(hmed)
	    if (sample=="W" or sample=="Z" or sample=="H") and "ass" in variables[i]:
              entry+=", "+str(hres)
            legend1.AddEntry(plots[i+offsetindex],entry,"l")
          print variables[i],peaks
	legend1.SetTextSize(0.04)
        legend1.SetFillStyle(0)
        legend1.Draw("same")

        if "tau21" in variables[i]:
          legend2=TLegend(0.55,0.9,0.9,0.95,"65<m_{pruned,corr}<95")
          legend2.SetTextSize(0.04)
          legend2.SetFillStyle(0)
          legend2.Draw("same")

        if "Tau21" in variables[i]:
          legend2=TLegend(0.55,0.9,0.9,0.95,"#tau_{2}/#tau_{1}<0.5")
          legend2.SetTextSize(0.04)
          legend2.SetFillStyle(0)
          legend2.Draw("same")

        if "D2" in variables[i]:
          legend2=TLegend(0.55,0.9,0.9,0.95,"65<m_{trimmed,corr}<95")
          legend2.SetTextSize(0.04)
          legend2.SetFillStyle(0)
          legend2.Draw("same")

        canvas.SaveAs(prefix + variables[i] + '.pdf')
        canvas.SaveAs(prefix + variables[i] + '.eps')
        canvas.SaveAs(prefix + variables[i] + '.root')
