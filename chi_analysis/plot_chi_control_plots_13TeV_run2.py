import os, sys
import array
from ROOT import gROOT,gStyle,TFile,TCanvas,TLegend

#gROOT.Reset()
gROOT.SetStyle("Plain")
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptFit(0)
gStyle.SetTitleOffset(1.2,"Y")
gStyle.SetPadLeftMargin(0.18)
gStyle.SetPadBottomMargin(0.15)
gStyle.SetPadTopMargin(0.08)
gStyle.SetPadRightMargin(0.05)
gStyle.SetMarkerSize(0.5)
gStyle.SetHistLineWidth(1)
gStyle.SetStatFontSize(0.020)
gStyle.SetTitleSize(0.06, "XYZ")
gStyle.SetLabelSize(0.05, "XYZ")
gStyle.SetNdivisions(510, "XYZ")
gStyle.SetLegendBorderSize(0)

if __name__ == '__main__':

   variables=["#chi","y_{boost}","p_{T1}","p_{T2}","y_{1}","y_{2}","#phi_{1}","#phi_{2}","METsumET","dPtsumPt","#Delta#phi","y"]
   label=["#chi","y_{boost}","p_{T1} [GeV]","p_{T2} [GeV]","y_{1}","y_{2}","#phi_{1}","#phi_{2}","missing E_{T} / #sum E_{T}","(p_{T1}-p_{T2})/(p_{T1}+p_{T2})","#Delta#phi","y"]
   #variables=["#chi","y_{boost}","p_{T1}","p_{T2}","y_{1}","y_{2}","METsumET","dPtsumPt","#Delta#phi"]
   #label=["#chi","y_{boost}","p_{T1} [GeV]","p_{T2} [GeV]","y_{1}","y_{2}","missing E_{T} / #sum E_{T}","(p_{T1}-p_{T2})/(p_{T1}+p_{T2})","#Delta#phi"]
   #variables=["#chi","y_{boost}","p_{T1}","p_{T2}","y_{1}","y_{2}","#Delta#phi"]
   #label=["#chi","y_{boost}","p_{T1} [GeV]","p_{T2} [GeV]","y_{1}","y_{2}","#Delta#phi"]

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

   colors=[1,2,3,4,6,7,8,9,10,11,12,13]
   styles=[1,2,3,4,5,6,7,8,9,11,12,13]

   chi_binnings=[]
   for mass_bin in chi_bins:
        chi_binnings+=[array.array('d')]
        for chi_bin in mass_bin:
            chi_binnings[-1].append(chi_bin)

   mass_binning=array.array('d')
   for mass_bin in [1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000]:
     mass_binning.append(mass_bin)
   
   prefix="datacard_shapelimit13TeV_run2_"
   dire="data/"
   version=""
   version="_run2"; postfix1617="_L1prefire"; postfix18="_HEM"
   #version="_run2_noHEM_noPrefire"; postfix1617=""; postfix18=""
   use_UL=False
   compare_Run2Run3=True
   use_2024=True
   compare_EOYvsUL=False
   compare_NULvsUL=False
   compare_EOYvsUL_MC=False
   compare_RECOvsGEN=False
   compare_years=False
   if use_UL:
     version="_UL"+version
   if compare_Run2Run3:
     if not use_2024:
       version="_2022_2023"+version
     else:
       version="_2024"+version
   if compare_EOYvsUL:
     version="_EOYvsUL"+version
   if compare_NULvsUL:
     version="_NULvsUL"+version
   if compare_years:
     version="_years"+version
   if compare_EOYvsUL_MC:
     version="_EOYvsUL_MC"+version
   if compare_RECOvsGEN:
     version="_RECOvsGEN"+version

   if compare_EOYvsUL or use_UL or compare_NULvsUL or compare_years:
    data=[("UL16preVFP", 1, "Data (UL 2016)"),
          ("UL16postVFP", 1, "")
         ]
    data2=[("UL17", 1, "Data (UL 2017)")
         ]
    data3=[("UL18", 1, "Data (UL 2018)")
         ]
   elif compare_Run2Run3:
    data=[("UL16preVFP", 1, "Data (Run2)"),
          ("UL16postVFP", 1, ""),
          ("UL17", 1, ""),
          ("UL18", 1, "")
         ]
    if use_2024:
      data2=[("2024-16May2025", 1, "Data (2024)")
           ]
      data3=[]#("2024", 1, "Data (2024 old)")]
    else:
      data2=[("2023-28May2024", 1, "Data (2023)")
           ]
      data3=[("2022", 1, "Data (2022)")
           ]
   elif compare_EOYvsUL_MC or compare_RECOvsGEN:
    UL16preFactor=19.52/36.33 # 19.52 is lumi of preVFP, 36.33 is lumi of all UL16
    UL16postFactor=16.81/36.33 # 16.81 is lumi of postVFP, 36.33 is lumi of all UL16
    data=[("UL16preVFP_QCDmadgraph-HT200to300",UL16preFactor*1710000./17969592, "MG+Py QCD (UL16)"),
       ("UL16preVFP_QCDmadgraph-HT300to500",UL16preFactor*347500./13586390, ""),
       ("UL16preVFP_QCDmadgraph-HT500to700",UL16preFactor*30363.051/55497082, ""),
       ("UL16preVFP_QCDmadgraph-HT700to1000",UL16preFactor*6428.869/15242034, ""),
       ("UL16preVFP_QCDmadgraph-HT1000to1500",UL16preFactor*1122.659/13559959, ""),
       ("UL16preVFP_QCDmadgraph-HT1500to2000",UL16preFactor*108.163/9661950, ""),
       ("UL16preVFP_QCDmadgraph-HT2000toInf",UL16preFactor*22.008/4827641, ""),
       ("UL16postVFP_QCDmadgraph-HT200to300",UL16postFactor*1710000./42723038, "MG+Py QCD (UL16)"),
       ("UL16postVFP_QCDmadgraph-HT300to500",UL16postFactor*347500./45502889, ""),
       ("UL16postVFP_QCDmadgraph-HT500to700",UL16postFactor*30363.051/15066884, ""),
       ("UL16postVFP_QCDmadgraph-HT700to1000",UL16postFactor*6428.869/13714842, ""),
       ("UL16postVFP_QCDmadgraph-HT1000to1500",UL16postFactor*1122.659/12416669, ""),
       ("UL16postVFP_QCDmadgraph-HT1500to2000",UL16postFactor*108.163/9244228, ""),
       ("UL16postVFP_QCDmadgraph-HT2000toInf",UL16postFactor*22.008/4843949, ""),
       ]
    data2=[("UL17_QCDmadgraph-HT200to300",1710000./42316128, "MG+Py QCD (UL17)"),
       ("UL17_QCDmadgraph-HT300to500",347500./42914024, ""),
       ("UL17_QCDmadgraph-HT500to700",30363.051/35745565, ""),
       ("UL17_QCDmadgraph-HT700to1000",6428.869/33646855, ""),
       ("UL17_QCDmadgraph-HT1000to1500",1122.659/10136610, ""),
       ("UL17_QCDmadgraph-HT1500to2000",108.163/7528926, ""),
       ("UL17_QCDmadgraph-HT2000toInf",22.008/4089387, ""),
       ]
    data3=[("UL18_QCDmadgraph-HT200to300",1710000./56298746, "MG+Py QCD (UL18)"),
       ("UL18_QCDmadgraph-HT300to500",347500./60991701, ""),
       ("UL18_QCDmadgraph-HT500to700",30363.051/48640047, ""),
       ("UL18_QCDmadgraph-HT700to1000",6428.869/47925782, ""),
       ("UL18_QCDmadgraph-HT1000to1500",1122.659/14244456, ""),
       ("UL18_QCDmadgraph-HT1500to2000",108.163/10751607, ""),
       ("UL18_QCDmadgraph-HT2000toInf",22.008/5278880, ""),
       ]
   else:
    data=[("2016", 1, "Data (2016)")
         ]
    data2=[("2017", 1, "Data (2017)")
         ]
    data3=[("2018", 1, "Data (2018)")
         ]
   if compare_EOYvsUL:
    mc=[("2016", 1, "Data (EOY 2016)")
         ]
    mc2=[("2017", 1, "Data (EOY 2017)")
         ]
    mc3=[("2018", 1, "Data (EOY 2018)")
         ]
   elif compare_NULvsUL:
    mc=[("NUL16preVFP", 1, "Data (NUL 2016)"),
          ("NUL16postVFP", 1, "")
         ]
    mc2=[("NUL17", 1, "Data (NUL 2017)")
         ]
    mc3=[("NUL18", 1, "Data (NUL 2018)")
         ]
   elif compare_years:
    mc=[("UL16preVFP", 1, "Data (Run2)"),
        ("UL16postVFP", 1, ""),
        ("UL17", 1, ""),
        ("UL18", 1, "")]
    mc2=[("UL16preVFP", 1, "Data (Run2)"),
        ("UL16postVFP", 1, ""),
        ("UL17", 1, ""),
        ("UL18", 1, "")]
    mc3=[("UL16preVFP", 1, "Data (Run2)"),
        ("UL16postVFP", 1, ""),
        ("UL17", 1, ""),
        ("UL18", 1, "")]
   elif compare_Run2Run3:
    mc=[("UL16preVFP_QCDmadgraph-HT200to300",36.31/137.62*19.52/36.33*1710000./17969592, "MG+Py QCD (Run2)"), # 19.52 is lumi of preVFP, 36.33 is lumi of all UL16
       ("UL16preVFP_QCDmadgraph-HT300to500",36.31/137.62*19.52/36.33*347500./13586390, ""),
       ("UL16preVFP_QCDmadgraph-HT500to700",36.31/137.62*19.52/36.33*30363.051/55497082, ""),
       ("UL16preVFP_QCDmadgraph-HT700to1000",36.31/137.62*19.52/36.33*6428.869/15242034, ""),
       ("UL16preVFP_QCDmadgraph-HT1000to1500",36.31/137.62*19.52/36.33*1122.659/13559959, ""),
       ("UL16preVFP_QCDmadgraph-HT1500to2000",36.31/137.62*19.52/36.33*108.163/9661950, ""),
       ("UL16preVFP_QCDmadgraph-HT2000toInf",36.31/137.62*19.52/36.33*22.008/4827641, ""),
       ("UL16postVFP_QCDmadgraph-HT200to300",36.31/137.62*16.81/36.33*1710000./42723038, ""), # 16.81 is lumi of postVFP, 36.33 is lumi of all UL16
       ("UL16postVFP_QCDmadgraph-HT300to500",36.31/137.62*16.81/36.33*347500./45502889, ""),
       ("UL16postVFP_QCDmadgraph-HT500to700",36.31/137.62*16.81/36.33*30363.051/15066884, ""),
       ("UL16postVFP_QCDmadgraph-HT700to1000",36.31/137.62*16.81/36.33*6428.869/13714842, ""),
       ("UL16postVFP_QCDmadgraph-HT1000to1500",36.31/137.62*16.81/36.33*1122.659/12416669, ""),
       ("UL16postVFP_QCDmadgraph-HT1500to2000",36.31/137.62*16.81/36.33*108.163/9244228, ""),
       ("UL16postVFP_QCDmadgraph-HT2000toInf",36.31/137.62*16.81/36.33*22.008/4843949, ""),
       ("UL17_QCDmadgraph-HT200to300",41.48/137.62*1710000./42316128, ""),
       ("UL17_QCDmadgraph-HT300to500",41.48/137.62*347500./42914024, ""),
       ("UL17_QCDmadgraph-HT500to700",41.48/137.62*30363.051/35745565, ""),
       ("UL17_QCDmadgraph-HT700to1000",41.48/137.62*6428.869/33646855, ""),
       ("UL17_QCDmadgraph-HT1000to1500",41.48/137.62*1122.659/10136610, ""),
       ("UL17_QCDmadgraph-HT1500to2000",41.48/137.62*108.163/7528926, ""),
       ("UL17_QCDmadgraph-HT2000toInf",41.48/137.62*22.008/4089387, ""),
       ("UL18_QCDmadgraph-HT200to300",59.83/137.62*1710000./56298746, ""),
       ("UL18_QCDmadgraph-HT300to500",59.83/137.62*347500./60991701, ""),
       ("UL18_QCDmadgraph-HT500to700",59.83/137.62*30363.051/48640047, ""),
       ("UL18_QCDmadgraph-HT700to1000",59.83/137.62*6428.869/47925782, ""),
       ("UL18_QCDmadgraph-HT1000to1500",59.83/137.62*1122.659/14244456, ""),
       ("UL18_QCDmadgraph-HT1500to2000",59.83/137.62*108.163/10751607, ""),
       ("UL18_QCDmadgraph-HT2000toInf",59.83/137.62*22.008/5278880, ""),
       ]
    mc2=[("2023_QCDmadgraph-28May2024-HT200to400",17.060/27.207*1.877e+06/37041533, "MG+Py QCD (2023)"),
       ("2023_QCDmadgraph-28May2024-HT400to600",17.060/27.207*9.521e+04/36680240, ""),
       ("2023_QCDmadgraph-28May2024-HT600to800",17.060/27.207*1.345e+04/34513335, ""),
       ("2023_QCDmadgraph-28May2024-HT800to1000",17.060/27.207*3.106e+03/37794527, ""),
       ("2023_QCDmadgraph-28May2024-HT1000to1200",17.060/27.207*8.938e+02/32193981, ""),
       ("2023_QCDmadgraph-28May2024-HT1200to1500",17.060/27.207*3.808e+02/40348927, ""),
       ("2023_QCDmadgraph-28May2024-HT1500to2000",17.060/27.207*1.157e+02/39490287, ""),
       ("2023_QCDmadgraph-28May2024-HT2000",17.060/27.207*2.621e+01/42353820, ""),
       ("2023_QCDmadgraphBPix-28May2024-HT200to400",9.525/27.207*2.030e+06/17128041, ""),
       ("2023_QCDmadgraphBPix-28May2024-HT400to600",9.525/27.207*9.459e+04/20454226, ""),
       ("2023_QCDmadgraphBPix-28May2024-HT600to800",9.525/27.207*1.417e+04/19700360, ""),
       ("2023_QCDmadgraphBPix-28May2024-HT800to1000",9.525/27.207*3.013e+03/17683808, ""),
       ("2023_QCDmadgraphBPix-28May2024-HT1000to1200",9.525/27.207*9.183e+02/18270279, ""),
       ("2023_QCDmadgraphBPix-28May2024-HT1200to1500",9.525/27.207*3.929e+02/18878616, ""),
       ("2023_QCDmadgraphBPix-28May2024-HT1500to2000",9.525/27.207*1.340e+02/17117867, ""),
       ("2023_QCDmadgraphBPix-28May2024-HT2000",9.525/27.207*2.676e+01/20420289, ""),
       ]
    if use_2024:
      mc3=[]
    else:
      mc3=[("2022_QCDmadgraph-HT200to400",8.077/35.18188*2.028e+06/20642715, "MG+Py QCD (2022)"),
       ("2022_QCDmadgraph-HT400to600",8.077/35.18188*1.003e+05/19602817, ""),
       ("2022_QCDmadgraph-HT600to800",8.077/35.18188*1.332e+04/19028458, ""),
       ("2022_QCDmadgraph-HT800to1000",8.077/35.18188*3.031e+03/21602068, ""),
       ("2022_QCDmadgraph-HT1000to1200",8.077/35.18188*8.579e+02/20642169, ""),
       ("2022_QCDmadgraph-HT1200to1500",8.077/35.18188*3.893e+02/21112774, ""),
       ("2022_QCDmadgraph-HT1500to2000",8.077/35.18188*1.290e+02/21191630, ""),
       ("2022_QCDmadgraph-HT2000",8.077/35.18188*2.496e+01/18132719, ""),
       ("2022_QCDmadgraphEE-HT200to400",27.0072/35.18188*1.824e+06/70275585, ""),
       ("2022_QCDmadgraphEE-HT400to600",27.0072/35.18188*9.232e+04/68707093, ""),
       ("2022_QCDmadgraphEE-HT600to800",27.0072/35.18188*1.340e+04/63298004, ""),
       ("2022_QCDmadgraphEE-HT800to1000",27.0072/35.18188*3.011e+03/66615474, ""),
       ("2022_QCDmadgraphEE-HT1000to1200",27.0072/35.18188*9.685e+02/70854616, ""),
       ("2022_QCDmadgraphEE-HT1200to1500",27.0072/35.18188*4.131e+02/70804499, ""),
       ("2022_QCDmadgraphEE-HT1500to2000",27.0072/35.18188*1.234e+02/63307771, ""),
       ("2022_QCDmadgraphEE-HT2000",27.0072/35.18188*2.612e+01/65102938, ""),
       ]
   elif use_UL and not compare_EOYvsUL_MC or compare_RECOvsGEN:
    mc=[("UL16preVFP_QCDmadgraph-HT200to300",19.52/36.33*1710000./17969592, "MG+Py QCD (UL16)"), # 19.52 is lumi of preVFP, 36.33 is lumi of all UL16
       ("UL16preVFP_QCDmadgraph-HT300to500",19.52/36.33*347500./13586390, ""),
       ("UL16preVFP_QCDmadgraph-HT500to700",19.52/36.33*30363.051/55497082, ""),
       ("UL16preVFP_QCDmadgraph-HT700to1000",19.52/36.33*6428.869/15242034, ""),
       ("UL16preVFP_QCDmadgraph-HT1000to1500",19.52/36.33*1122.659/13559959, ""),
       ("UL16preVFP_QCDmadgraph-HT1500to2000",19.52/36.33*108.163/9661950, ""),
       ("UL16preVFP_QCDmadgraph-HT2000toInf",19.52/36.33*22.008/4827641, ""),
       ("UL16postVFP_QCDmadgraph-HT200to300",16.81/36.33*1710000./42723038, "MG+Py QCD (UL16)"), # 16.81 is lumi of postVFP, 36.33 is lumi of all UL16
       ("UL16postVFP_QCDmadgraph-HT300to500",16.81/36.33*347500./45502889, ""),
       ("UL16postVFP_QCDmadgraph-HT500to700",16.81/36.33*30363.051/15066884, ""),
       ("UL16postVFP_QCDmadgraph-HT700to1000",16.81/36.33*6428.869/13714842, ""),
       ("UL16postVFP_QCDmadgraph-HT1000to1500",16.81/36.33*1122.659/12416669, ""),
       ("UL16postVFP_QCDmadgraph-HT1500to2000",16.81/36.33*108.163/9244228, ""),
       ("UL16postVFP_QCDmadgraph-HT2000toInf",16.81/36.33*22.008/4843949, ""),
       ]
    mc2=[("UL17_QCDmadgraph-HT200to300",1710000./42316128, "MG+Py QCD (UL17)"),
       ("UL17_QCDmadgraph-HT300to500",347500./42914024, ""),
       ("UL17_QCDmadgraph-HT500to700",30363.051/35745565, ""),
       ("UL17_QCDmadgraph-HT700to1000",6428.869/33646855, ""),
       ("UL17_QCDmadgraph-HT1000to1500",1122.659/10136610, ""),
       ("UL17_QCDmadgraph-HT1500to2000",108.163/7528926, ""),
       ("UL17_QCDmadgraph-HT2000toInf",22.008/4089387, ""),
       ]
    mc3=[("UL18_QCDmadgraph-HT200to300",1710000./56298746, "MG+Py QCD (UL18)"),
       ("UL18_QCDmadgraph-HT300to500",347500./60991701, ""),
       ("UL18_QCDmadgraph-HT500to700",30363.051/48640047, ""),
       ("UL18_QCDmadgraph-HT700to1000",6428.869/47925782, ""),
       ("UL18_QCDmadgraph-HT1000to1500",1122.659/14244456, ""),
       ("UL18_QCDmadgraph-HT1500to2000",108.163/10751607, ""),
       ("UL18_QCDmadgraph-HT2000toInf",22.008/5278880, ""),
       ]
   else:
    mc=[("2016_QCDmadgraph-HT200to300",1712000./56709875, "MG+Py QCD (EOY16)"),
       ("2016_QCDmadgraph-HT300to500",347700./53096517, ""),
       ("2016_QCDmadgraph-HT500to700",32100./52906552, ""),
       ("2016_QCDmadgraph-HT700to1000",6831./36741540, ""),
       ("2016_QCDmadgraph-HT1000to1500",1207./15210939, ""),
       ("2016_QCDmadgraph-HT1500to2000",119.9/11839357, ""),
       ("2016_QCDmadgraph-HT2000toInf",25.24/5947849, ""),
       ]
    mc2=[("2017_QCDmadgraph-HT200to300",1545000./58990434, "MG+Py QCD (EOY17)"),
       ("2017_QCDmadgraph-HT300to500",323300./58748739, ""),
       ("2017_QCDmadgraph-HT500to700",30000./54366431, ""),
       ("2017_QCDmadgraph-HT700to1000",6324./46924322, ""),
       ("2017_QCDmadgraph-HT1000to1500",1090./16495598, ""),
       ("2017_QCDmadgraph-HT1500to2000",101./11196479, ""),
       ("2017_QCDmadgraph-HT2000toInf",20.43/5362513, ""),
       ]
    mc3=[("2018_QCDmadgraph-HT200to300",1461000./54289442, "MG+Py QCD (EOY18)"),
       ("2018_QCDmadgraph-HT300to500",311900./54512704, ""),
       ("2018_QCDmadgraph-HT500to700",29070./53919811, ""),
       ("2018_QCDmadgraph-HT700to1000",5962./48158738, ""),
       ("2018_QCDmadgraph-HT1000to1500",1005./14945819, ""),
       ("2018_QCDmadgraph-HT1500to2000",101.8/10707847, ""),
       ("2018_QCDmadgraph-HT2000toInf",20.54/5329144, ""),
       ]
   f_data=[]
   f_data2=[]
   f_data3=[]
   for name,xsec,l in data:
      if ("UL16" in name or "UL17" in name) and not "QCD" in name: postfix=postfix1617
      elif "UL18" in name and not "QCD" in name: postfix=postfix18
      else: postfix=""
      f_data+=[TFile.Open(dire+prefix+name+postfix+"_chi.root")]
      print(prefix+name+postfix+"_chi.root")
   for name,xsec,l in data2:
      if ("UL16" in name or "UL17" in name) and not "QCD" in name: postfix=postfix1617
      elif "UL18" in name and not "QCD" in name: postfix=postfix18
      else: postfix=""
      f_data2+=[TFile.Open(dire+prefix+name+postfix+"_chi.root")]
      print(prefix+name+postfix+"_chi.root")
   for name,xsec,l in data3:
      if ("UL16" in name or "UL17" in name) and not "QCD" in name: postfix=postfix1617
      elif "UL18" in name and not "QCD" in name: postfix=postfix18
      else: postfix=""
      f_data3+=[TFile.Open(dire+prefix+name+postfix+"_chi.root")]
      print(prefix+name+postfix+"_chi.root")
   if compare_Run2Run3:
     if not use_2024:
       lumi=[35.18,137.62,27.20]
     else:
       lumi=[109.08,137.62,0]
   else:
     lumi=[36.33,41.53,59.74]
   f_mc=[]
   f_mc2=[]
   f_mc3=[]
   for name,xsec,l in mc:
      if ("UL16" in name or "UL17" in name) and not "QCD" in name: postfix=postfix1617
      elif "UL18" in name and not "QCD" in name: postfix=postfix18
      else: postfix=""
      f_mc+=[TFile.Open(dire+prefix+name+postfix+("-GEN" if compare_RECOvsGEN else "")+"_chi.root")]
      print(prefix+name+postfix+"_chi.root")
   for name,xsec,l in mc2:
      if ("UL16" in name or "UL17" in name) and not "QCD" in name: postfix=postfix1617
      elif "UL18" in name and not "QCD" in name: postfix=postfix18
      else: postfix=""
      f_mc2+=[TFile.Open(dire+prefix+name+postfix+("-GEN" if compare_RECOvsGEN else "")+"_chi.root")]
      print(prefix+name+postfix+"_chi.root")
   for name,xsec,l in mc3:
      if ("UL16" in name or "UL17" in name) and not "QCD" in name: postfix=postfix1617
      elif "UL18" in name and not "QCD" in name: postfix=postfix18
      else: postfix=""
      f_mc3+=[TFile.Open(dire+prefix+name+postfix+("-GEN" if compare_RECOvsGEN else "")+"_chi.root")]
      print(prefix+name+postfix+"_chi.root")

   if compare_EOYvsUL_MC:
     minmass=1900
   else:
     minmass=2400
  
   for var in ["mass"]:
     canvas = TCanvas("","",0,0,200,300)
     canvas.Divide(1,3,0,0,0)
     canvas.GetPad(1).SetPad(0.0,0.40,1.0,1.0)
     canvas.GetPad(1).SetLeftMargin(0.15)
     canvas.GetPad(1).SetRightMargin(0.08)
     canvas.GetPad(1).SetTopMargin(0.08)
     canvas.GetPad(1).SetBottomMargin(0.05)
     canvas.GetPad(2).SetPad(0.0,0.25,1.0,0.40)
     canvas.GetPad(2).SetLeftMargin(0.15)
     canvas.GetPad(2).SetRightMargin(0.08)
     canvas.GetPad(2).SetTopMargin(0.08)
     canvas.GetPad(2).SetBottomMargin(0.05)
     canvas.GetPad(3).SetPad(0.0,0.0,1.0,0.25)
     canvas.GetPad(3).SetLeftMargin(0.15)
     canvas.GetPad(3).SetRightMargin(0.08)
     canvas.GetPad(3).SetTopMargin(0.08)
     canvas.GetPad(3).SetBottomMargin(0.45)
     canvas.cd(1)
     canvas.GetPad(1).SetLogy(True)

     legend=TLegend(0.5,0.5,0.95,0.9,"1<=#Chi<16 , y_{boost}<1.11")
     
     hist=f_data[0].Get(prefix+data[0][0].split("-HT")[0]+var)
     hist.Scale(data[0][1])
     for i in range(1,len(data)):
         hist.Add(f_data[i].Get(prefix+data[i][0].split("-HT")[0]+var),data[i][1])
     hist=hist.Rebin(4)#len(mass_binning)-1,hist.GetName()+"_rebin1",mass_binning)
     normfactor=hist.Integral(hist.FindBin(minmass),hist.FindBin(10000))
     if compare_EOYvsUL or compare_NULvsUL:
       normfactor=lumi[0]*1000.
     if compare_EOYvsUL_MC or compare_RECOvsGEN:
       normfactor=1.
     hist.Scale(1./normfactor)
     hist.SetLineColor(1)
     hist.SetMarkerStyle(24)
     hist.SetMarkerColor(1)
     hist.SetMarkerSize(0.2)
     hist.SetTitle("")
     hist.GetXaxis().SetLabelColor(0)
     if compare_EOYvsUL or compare_NULvsUL or compare_EOYvsUL_MC or compare_RECOvsGEN:
       hist.GetYaxis().SetTitle("Cross section [pb]")
     else:
       hist.GetYaxis().SetTitle("Normalized distribution")
     hist.GetXaxis().SetRangeUser(minmass,10000)
     if compare_EOYvsUL_MC or compare_RECOvsGEN:
       hist.GetYaxis().SetRangeUser(0.5/1e6,hist.GetMaximum()*1.5)
     else:
       hist.GetYaxis().SetRangeUser(0.5/normfactor,hist.GetMaximum()*1.5)
     hist.GetXaxis().SetTitleOffset(1.1)
     hist.GetYaxis().SetTitleOffset(1.1)
     hist.GetXaxis().SetLabelSize(0.05)
     hist.GetYaxis().SetLabelSize(0.05)
     hist.GetXaxis().SetTitleSize(0.06)
     hist.GetYaxis().SetTitleSize(0.06)
     hist.SetStats(False)
     hist.Draw("pe")
     hist.GetXaxis().SetRangeUser(minmass,10000)
     legend.AddEntry(hist,data[0][2]+(" RECO" if compare_RECOvsGEN else ""),"lpe")

     hist2=f_data2[0].Get(prefix+data2[0][0].split("-HT")[0]+var)
     hist2.Scale(data2[0][1])
     for i in range(1,len(data2)):
        hist2.Add(f_data2[i].Get(prefix+data2[i][0].split("-HT")[0]+var),data2[i][1])
     hist2=hist2.Rebin(4)#len(mass_binning)-1,hist2.GetName()+"_rebin1",mass_binning)
     normfactor2=hist2.Integral(hist2.FindBin(minmass),hist2.FindBin(10000))
     if compare_EOYvsUL or compare_NULvsUL:
       normfactor2=lumi[1]*1000.
     if compare_EOYvsUL_MC or compare_RECOvsGEN:
       normfactor2=1.
     hist2.Scale(1./normfactor2)
     hist2.SetLineColor(2)
     hist2.SetMarkerStyle(25)
     hist2.SetMarkerColor(2)
     hist2.SetMarkerSize(0.2)
     hist2.SetStats(False)
     hist2.Draw("pesame")
     legend.AddEntry(hist2,data2[0][2]+(" RECO" if compare_RECOvsGEN else ""),"lpe")
     
     if len(f_data3)>0:
       hist3=f_data3[0].Get(prefix+data3[0][0].split("-HT")[0]+var)
       hist3.Scale(data3[0][1])
       for i in range(1,len(data3)):
          hist3.Add(f_data3[i].Get(prefix+data3[i][0].split("-HT")[0]+var),data3[i][1])
       hist3=hist3.Rebin(4)#len(mass_binning)-1,hist3.GetName()+"_rebin1",mass_binning)
       normfactor3=hist3.Integral(hist3.FindBin(minmass),hist3.FindBin(10000))
       if compare_EOYvsUL or compare_NULvsUL:
         normfactor3=lumi[2]*1000.
       if compare_EOYvsUL_MC or compare_RECOvsGEN:
         normfactor3=1.
       hist3.Scale(1./normfactor3)
       hist3.SetLineColor(4)
       hist3.SetMarkerStyle(26)
       hist3.SetMarkerColor(4)
       hist3.SetMarkerSize(0.2)
       hist3.SetStats(False)
       hist3.Draw("pesame")
       legend.AddEntry(hist3,data3[0][2]+(" RECO" if compare_RECOvsGEN else ""),"lpe")
     
     hist_mc=f_mc[0].Get(prefix+mc[0][0].split("-HT")[0]+var)
     hist_mc.Scale(mc[0][1])
     for i in range(1,len(mc)):
         hist_mc.Add(f_mc[i].Get(prefix+mc[i][0].split("-HT")[0]+var),mc[i][1])
     hist_mc=hist_mc.Rebin(4)#len(mass_binning)-1,hist_mc.GetName()+"_rebin1",mass_binning)
     if compare_EOYvsUL or compare_NULvsUL:
       hist_mc.Scale(1./normfactor)
     elif compare_EOYvsUL_MC or compare_RECOvsGEN or compare_years:
       hist_mc.Scale(1.)
     else:
       hist_mc.Scale(hist.Integral(hist.FindBin(minmass),hist.GetNbinsX())/hist_mc.Integral(hist_mc.FindBin(minmass),hist_mc.GetNbinsX()))
     hist_mc.SetLineColor(1)
     hist_mc.SetStats(False)
     hist_mc.Draw("histsame")
     legend.AddEntry(hist_mc,mc[0][2]+(" GEN" if compare_RECOvsGEN else ""),"l")

     hist_mc2=f_mc2[0].Get(prefix+mc2[0][0].split("-HT")[0]+var)
     hist_mc2.Scale(mc2[0][1])
     for i in range(1,len(mc2)):
        #print(f_mc2[i],prefix+mc2[0][0].split("-HT")[0]+var, i,mc2[i][1])
        hist_mc2.Add(f_mc2[i].Get(prefix+mc2[i][0].split("-HT")[0]+var),mc2[i][1])
     hist_mc2=hist_mc2.Rebin(4)#len(mass_binning)-1,hist_mc2.GetName()+"_rebin1",mass_binning)
     if compare_EOYvsUL or compare_NULvsUL:
       hist_mc2.Scale(1./normfactor2)
     elif compare_EOYvsUL_MC or compare_RECOvsGEN:
       hist_mc2.Scale(1.)
     else:
       hist_mc2.Scale(hist2.Integral(hist2.FindBin(minmass),hist2.GetNbinsX())/hist_mc2.Integral(hist_mc2.FindBin(minmass),hist_mc2.GetNbinsX()))
     hist_mc2.SetLineColor(2)
     hist_mc2.SetStats(False)
     hist_mc2.Draw("histsame")
     legend.AddEntry(hist_mc2,mc2[0][2]+(" GEN" if compare_RECOvsGEN else ""),"l")

     if len(f_mc3)>0:
       hist_mc3=f_mc3[0].Get(prefix+mc3[0][0].split("-HT")[0]+var)
       hist_mc3.Scale(mc3[0][1])
       for i in range(1,len(mc3)):
          hist_mc3.Add(f_mc3[i].Get(prefix+mc3[i][0].split("-HT")[0]+var),mc3[i][1])
       hist_mc3=hist_mc3.Rebin(4)#len(mass_binning)-1,hist_mc3.GetName()+"_rebin1",mass_binning)
       if compare_EOYvsUL or compare_NULvsUL:
         hist_mc3.Scale(1./normfactor3)
       elif compare_EOYvsUL_MC or compare_RECOvsGEN:
         hist_mc3.Scale(1.)
       else:
         hist_mc3.Scale(hist3.Integral(hist3.FindBin(minmass),hist3.GetNbinsX())/hist_mc3.Integral(hist_mc3.FindBin(minmass),hist_mc3.GetNbinsX()))
       hist_mc3.SetLineColor(4)
       hist_mc3.SetStats(False)
       hist_mc3.Draw("histsame")
       legend.AddEntry(hist_mc3,mc3[0][2]+(" GEN" if compare_RECOvsGEN else ""),"l")

     hist.Draw("pesame")
     hist2.Draw("pesame")
     if len(f_data3)>0: hist3.Draw("pesame")
     hist.Draw("axissame")
     
     legend.SetTextSize(0.04)
     legend.SetFillStyle(0)
     legend.Draw("same")

     if compare_Run2Run3:
       histrebinned=hist2.Clone(hist.GetName()+"rebinned")
       for b in range(hist.GetNbinsX()):
           histrebinned.SetBinContent(b+1,hist.GetBinContent(b+1))
           histrebinned.SetBinError(b+1,hist.GetBinError(b+1))
       histrebinned.SetLineColor(hist.GetLineColor())
       histrebinned.SetLineStyle(hist.GetLineStyle())
       histrebinned.SetMarkerColor(hist.GetMarkerColor())
       histrebinned.SetMarkerStyle(hist.GetMarkerStyle())
       histrebinned.SetMarkerSize(hist.GetMarkerSize())
       hist=histrebinned
       hist.GetXaxis().SetRangeUser(minmass,10000)
       histbackup=hist
       hist=hist_mc # use MC as reference for first ratio plot
       hist.GetXaxis().SetRangeUser(minmass,10000)
       histbackup2=hist2
       hist2=hist_mc
       if len(f_mc3)>0:
         histbackup3=hist3
         hist3=hist_mc

     canvas.cd(2)
     ratio=hist.Clone(hist.GetName()+"ratio")
     ratio.Divide(hist,hist)
     for b in range(hist.GetNbinsX()):
       if hist.GetBinContent(b+1)>0:
         ratio.SetBinError(b+1,hist.GetBinError(b+1)/hist.GetBinContent(b+1))
     ratio.SetTitle("")
     ratio.GetYaxis().SetTitle("Run / year" if compare_years else ("Sim / "+mc[0][2] if compare_Run2Run3 else ("EOY / UL" if compare_EOYvsUL or compare_EOYvsUL_MC else ("new / old" if compare_NULvsUL else ("RECO / GEN" if compare_RECOvsGEN else "Sim / Data")))))
     ratio.GetYaxis().SetTitleSize(0.18/(len(ratio.GetYaxis().GetTitle())/10))
     ratio.GetYaxis().SetTitleOffset(0.3)
     ratio.SetMarkerSize(0.1)
     ratio.GetYaxis().SetLabelSize(0.2)
     ratio.GetYaxis().SetRangeUser(0,2)
     if compare_EOYvsUL or compare_NULvsUL or compare_EOYvsUL_MC or compare_RECOvsGEN or compare_years:
      ratio.GetYaxis().SetRangeUser(0.8,1.2)
     ratio.GetXaxis().SetNdivisions(506)
     ratio.GetYaxis().SetNdivisions(503)
     ratio.GetXaxis().SetLabelColor(0)
     ratio.GetXaxis().SetTitleSize(0.2)
     ratio.GetXaxis().SetTitleOffset(1.1)
     ratio.GetXaxis().SetLabelSize(0.18)
     ratio.Draw("histe")
     ratio_mc=hist_mc.Clone(hist_mc.GetName()+"ratio")
     ratio_mc.Divide(hist_mc,hist)
     for b in range(hist.GetNbinsX()):
       if hist.GetBinContent(b+1)>0:
         ratio_mc.SetBinError(b+1,hist.GetBinError(b+1)/hist.GetBinContent(b+1))
     ratio_mc.Draw("histsame")
     ratio_mc2=hist_mc2.Clone(hist_mc2.GetName()+"ratio")
     ratio_mc2.Divide(hist_mc2,hist2)
     for b in range(hist2.GetNbinsX()):
       if hist2.GetBinContent(b+1)>0:
         ratio_mc2.SetBinError(b+1,hist2.GetBinError(b+1)/hist2.GetBinContent(b+1))
     ratio_mc2.Draw("histsame")
     if len(f_mc3)>0:
       ratio_mc3=hist_mc3.Clone(hist_mc3.GetName()+"ratio")
       ratio_mc3.Divide(hist_mc3,hist3)
       for b in range(hist3.GetNbinsX()):
         if hist3.GetBinContent(b+1)>0:
           ratio_mc3.SetBinError(b+1,hist3.GetBinError(b+1)/hist3.GetBinContent(b+1))
       ratio_mc3.Draw("histsame")
       ratio.Draw("axissame")

     if compare_Run2Run3:
       hist=histbackup
       hist2=histbackup2
       if len(f_mc3)>0:
         hist3=histbackup3

     canvas.cd(3)
     ddratio=hist.Clone(hist.GetName()+"ddratio")
     ddratio.Divide(hist,hist)
     for b in range(hist.GetNbinsX()):
       if hist.GetBinContent(b+1)>0:
         ddratio.SetBinError(b+1,hist.GetBinError(b+1)/hist.GetBinContent(b+1))
     ddratio.SetTitle("")
     if compare_EOYvsUL_MC or compare_RECOvsGEN:
       ddratio.GetYaxis().SetTitle("Sim / "+data[0][2]+"")
     else:
       ddratio.GetYaxis().SetTitle("Data / "+data[0][2]+"")
   
     ddratio.GetYaxis().SetTitleSize(0.11/(len(ddratio.GetYaxis().GetTitle())/10))
     ddratio.GetYaxis().SetTitleOffset(0.5)
     ddratio.SetMarkerSize(0.1)
     ddratio.GetYaxis().SetLabelSize(0.13)
     ddratio.GetYaxis().SetRangeUser(0.8,1.2)
     if compare_Run2Run3:
       ddratio.GetYaxis().SetRangeUser(0,2)
     ddratio.GetXaxis().SetNdivisions(506)
     ddratio.GetYaxis().SetNdivisions(503)
     ddratio.GetXaxis().SetLabelColor(1)
     ddratio.GetXaxis().SetTitle("dijet mass [GeV]")
     ddratio.GetXaxis().SetTitleSize(0.14)
     ddratio.GetXaxis().SetTitleOffset(1.1)
     ddratio.GetXaxis().SetLabelSize(0.12)
     ddratio.Draw("histe")
     ddratio2=hist2.Clone(hist2.GetName()+"ddratio")
     ddratio2.Divide(hist2,hist)
     for b in range(hist.GetNbinsX()):
       if hist.GetBinContent(b+1)>0:
         ddratio2.SetBinError(b+1,hist.GetBinError(b+1)/hist.GetBinContent(b+1))
     ddratio2.Draw("histsame")
     if len(f_data3)>0:
       ddratio3=hist3.Clone(hist3.GetName()+"ddratio")
       ddratio3.Divide(hist3,hist)
       for b in range(hist.GetNbinsX()):
         if hist.GetBinContent(b+1)>0:
           ddratio3.SetBinError(b+1,hist.GetBinError(b+1)/hist.GetBinContent(b+1))
       ddratio3.Draw("histsame")
     ratio.Draw("axissame")

     canvas.cd(1)
     hist.GetYaxis().SetTitleOffset(1.2)
     
     canvas.SaveAs("plots/chi_control_plots_mass"+version+".root")
     canvas.SaveAs("plots/chi_control_plots_mass"+version+".pdf")

   for var in variables:
    log=(var=="p_{T1}" or var=="p_{T2}" or var=="METsumET" or var=="#Delta#phi")
    legends=[]
    for mass in range(len(massbins)):
        print(var, str(massbins[mass][0])+"_"+str(massbins[mass][1]))
        canvas = TCanvas("","",0,0,200,300)
        canvas.Divide(1,3,0,0,0)
        canvas.GetPad(1).SetPad(0.0,0.40,1.0,1.0)
        canvas.GetPad(1).SetLeftMargin(0.15)
        canvas.GetPad(1).SetRightMargin(0.08)
        canvas.GetPad(1).SetTopMargin(0.08)
        canvas.GetPad(1).SetBottomMargin(0.05)
        canvas.GetPad(2).SetPad(0.0,0.25,1.0,0.40)
        canvas.GetPad(2).SetLeftMargin(0.15)
        canvas.GetPad(2).SetRightMargin(0.08)
        canvas.GetPad(2).SetTopMargin(0.08)
        canvas.GetPad(2).SetBottomMargin(0.05)
        canvas.GetPad(3).SetPad(0.0,0.0,1.0,0.25)
        canvas.GetPad(3).SetLeftMargin(0.15)
        canvas.GetPad(3).SetRightMargin(0.08)
        canvas.GetPad(3).SetTopMargin(0.08)
        canvas.GetPad(3).SetBottomMargin(0.45)
        canvas.cd(1)
        canvas.GetPad(1).SetLogy(log)
        legend=TLegend(0.45,0.6,0.95,0.90,(str(massbins[mass][0])+"<m_{jj}<"+str(massbins[mass][1])+" GeV").replace("7000<m_{jj}<13000","m_{jj}>7000"))
        legends+=[legend]
    
        name=prefix+data[0][0].split("-HT")[0]+var+str(massbins[mass][0])+"_"+str(13600 if "2024" in data[0][0] and massbins[mass][1]==13000 else massbins[mass][1])
        if var=="#chi": name+="_rebin1"
        if var=="y":
          hist=f_data[0].Get(name.replace("y"+str(massbins[mass][0]),"y_{1}"+str(massbins[mass][0])))
          hist.Add(f_data[0].Get(name.replace("y"+str(massbins[mass][0]),"y_{2}"+str(massbins[mass][0]))))
        else:
          hist=f_data[0].Get(name)
        for i in range(1,len(data)):
            namei=prefix+data[i][0].split("-HT")[0]+var+str(massbins[mass][0])+"_"+str(13600 if "2024" in data[i][0] and massbins[mass][1]==13000 else massbins[mass][1])
            if var=="#chi": namei+="_rebin1"
            if var=="y":
              hist.Add(f_data[i].Get(namei.replace("y"+str(massbins[mass][0]),"y_{1}"+str(massbins[mass][0]))),data[i][1]/data[0][1])
              hist.Add(f_data[i].Get(namei.replace("y"+str(massbins[mass][0]),"y_{2}"+str(massbins[mass][0]))),data[i][1]/data[0][1])
            else:
              hist.Add(f_data[i].Get(namei),data[i][1]/data[0][1])
        print("mass bin",mass,"data 2016 integral",hist.Integral())
        integral=hist.Integral()
        if var=="#chi":
            hist=hist.Rebin(len(chi_binnings[mass])-1,hist.GetName()+"_rebin1",chi_binnings[mass])
            for b in range(hist.GetNbinsX()):
               hist.SetBinContent(b+1,hist.GetBinContent(b+1)/hist.GetBinWidth(b+1))
               hist.SetBinError(b+1,hist.GetBinError(b+1)/hist.GetBinWidth(b+1))
        if var.startswith("#phi"):
          hist=hist.Rebin(2)
        hist.SetLineColor(1)
        hist.SetMarkerStyle(24)
        hist.SetMarkerSize(0.2)
        miny=0
        if hist.Integral()>0:
            miny=log*0.1/hist.GetMaximum()
            hist.Scale(1./integral)
        hist.SetTitle("")
        #hist.GetXaxis().SetTitle(label[variables.index(var)])
        hist.GetXaxis().SetLabelColor(0)
        hist.GetYaxis().SetTitle("Normalized distribution")
        hist.GetYaxis().SetRangeUser(miny,hist.GetMaximum()*(1.5+log*10))
        #hist.GetXaxis().SetTitleOffset(1.1)
        hist.GetYaxis().SetTitleOffset(1.1)
        hist.GetXaxis().SetLabelSize(0.05)
        hist.GetYaxis().SetLabelSize(0.05)
        hist.GetXaxis().SetTitleSize(0.06)
        hist.GetYaxis().SetTitleSize(0.06)
        hist.SetStats(False)
        hist.Draw("pe")
        legend.AddEntry(hist,data[0][2]+(" RECO" if compare_RECOvsGEN else ""),"lpe")

        name=prefix+data2[0][0].split("-HT")[0]+var+str(massbins[mass][0])+"_"+str(13600 if "2024" in data2[0][0] and massbins[mass][1]==13000 else massbins[mass][1])
        if var=="#chi": name+="_rebin1"
        if var=="y":
          hist2=f_data2[0].Get(name.replace("y"+str(massbins[mass][0]),"y_{1}"+str(massbins[mass][0])))
          hist2.Add(f_data2[0].Get(name.replace("y"+str(massbins[mass][0]),"y_{2}"+str(massbins[mass][0]))))
        else:
          hist2=f_data2[0].Get(name)
        for i in range(1,len(data2)):
            namei=prefix+data2[i][0].split("-HT")[0]+var+str(massbins[mass][0])+"_"+str(13600 if "2024" in data2[i][0] and massbins[mass][1]==13000 else massbins[mass][1])
            if var=="#chi": namei+="_rebin1"
            if var=="y":
              hist2.Add(f_data2[i].Get(namei.replace("y"+str(massbins[mass][0]),"y_{1}"+str(massbins[mass][0]))),data2[i][1]/data2[0][1])
              hist2.Add(f_data2[i].Get(namei.replace("y"+str(massbins[mass][0]),"y_{2}"+str(massbins[mass][0]))),data2[i][1]/data2[0][1])
            else:
              hist2.Add(f_data2[i].Get(namei),data2[i][1]/data2[0][1])
        print("mass bin",mass,"data 2017 integral",hist2.Integral())
        integral=hist2.Integral()
        if var=="#chi":
            hist2=hist2.Rebin(len(chi_binnings[mass])-1,hist2.GetName()+"_rebin1",chi_binnings[mass])
            for b in range(hist2.GetNbinsX()):
               hist2.SetBinContent(b+1,hist2.GetBinContent(b+1)/hist2.GetBinWidth(b+1))
               hist2.SetBinError(b+1,hist2.GetBinError(b+1)/hist2.GetBinWidth(b+1))
        if var.startswith("#phi"):
          hist2=hist2.Rebin(2)
        hist2.SetLineColor(2)
        hist2.SetMarkerStyle(25)
        hist2.SetMarkerColor(2)
        hist2.SetMarkerSize(0.2)
        miny=0
        if hist2.Integral()>0:
            miny=log*0.1/integral
            hist2.Scale(1./integral)
        hist2.SetTitle("")
        hist2.SetStats(False)
        hist2.Draw("pesame")
        legend.AddEntry(hist2,data2[0][2]+(" RECO" if compare_RECOvsGEN else ""),"lpe")

        if len(f_data3)>0:
          name=prefix+data3[0][0].split("-HT")[0]+var+str(massbins[mass][0])+"_"+str(13600 if "2024" in data3[0][0] and massbins[mass][1]==13000 else massbins[mass][1])
          if var=="#chi": name+="_rebin1"
          if var=="y":
            hist3=f_data3[0].Get(name.replace("y"+str(massbins[mass][0]),"y_{1}"+str(massbins[mass][0])))
            hist3.Add(f_data3[0].Get(name.replace("y"+str(massbins[mass][0]),"y_{2}"+str(massbins[mass][0]))))
          else:
            hist3=f_data3[0].Get(name)
          for i in range(1,len(data3)):
              namei=prefix+data3[i][0].split("-HT")[0]+var+str(massbins[mass][0])+"_"+str(13600 if "2024-" in data3[i][0] and massbins[mass][1]==13000 else massbins[mass][1])
              if var=="#chi": namei+="_rebin1"
              if var=="y":
                hist3.Add(f_data3[i].Get(namei.replace("y"+str(massbins[mass][0]),"y_{1}"+str(massbins[mass][0]))),data3[i][1]/data3[0][1])
                hist3.Add(f_data3[i].Get(namei.replace("y"+str(massbins[mass][0]),"y_{2}"+str(massbins[mass][0]))),data3[i][1]/data3[0][1])
              else:
                hist3.Add(f_data3[i].Get(namei),data3[i][1]/data3[0][1])
          print("mass bin",mass,"data 2018 integral",hist3.Integral())
          integral=hist3.Integral()
          if var=="#chi":
              hist3=hist3.Rebin(len(chi_binnings[mass])-1,hist3.GetName()+"_rebin1",chi_binnings[mass])
              for b in range(hist3.GetNbinsX()):
                 hist3.SetBinContent(b+1,hist3.GetBinContent(b+1)/hist3.GetBinWidth(b+1))
                 hist3.SetBinError(b+1,hist3.GetBinError(b+1)/hist3.GetBinWidth(b+1))
          if var.startswith("#phi"):
            hist3=hist3.Rebin(2)
          hist3.SetLineColor(4)
          hist3.SetMarkerStyle(25)
          hist3.SetMarkerColor(4)
          hist3.SetMarkerSize(0.2)
          miny=0
          if hist3.Integral()>0:
              miny=log*0.1/integral
              hist3.Scale(1./integral)
          hist3.SetTitle("")
          hist3.SetStats(False)
          hist3.Draw("pesame")
          legend.AddEntry(hist3,data3[0][2]+(" RECO" if compare_RECOvsGEN else ""),"lpe")

        name=prefix+mc[0][0].split("-HT")[0]+var+str(massbins[mass][0])+"_"+str(massbins[mass][1])
        if var=="#chi": name+="_rebin1"
        if var=="y":
          hist_mc=f_mc[0].Get(name.replace("y"+str(massbins[mass][0]),"y_{1}"+str(massbins[mass][0])))
          hist_mc.Add(f_mc[0].Get(name.replace("y"+str(massbins[mass][0]),"y_{2}"+str(massbins[mass][0]))))
        else:
          hist_mc=f_mc[0].Get(name)
        for i in range(1,len(mc)):
            namei=prefix+mc[i][0].split("-HT")[0]+var+str(massbins[mass][0])+"_"+str(massbins[mass][1])
            if var=="#chi": namei+="_rebin1"
            if var=="y":
              hist_mc.Add(f_mc[i].Get(namei.replace("y"+str(massbins[mass][0]),"y_{1}"+str(massbins[mass][0]))),mc[i][1]/mc[0][1])
              hist_mc.Add(f_mc[i].Get(namei.replace("y"+str(massbins[mass][0]),"y_{2}"+str(massbins[mass][0]))),mc[i][1]/mc[0][1])
            else:
              hist_mc.Add(f_mc[i].Get(namei),mc[i][1]/mc[0][1])
        if var=="#chi":
            hist_mc=hist_mc.Rebin(len(chi_binnings[mass])-1,hist_mc.GetName()+"_rebin1",chi_binnings[mass])
            for b in range(hist_mc.GetNbinsX()):
               hist_mc.SetBinContent(b+1,hist_mc.GetBinContent(b+1)/hist_mc.GetBinWidth(b+1))
               hist_mc.SetBinError(b+1,hist_mc.GetBinError(b+1)/hist_mc.GetBinWidth(b+1))
        if var.startswith("#phi"):
          hist_mc=hist_mc.Rebin(2)
        if hist_mc.Integral()>0:
            hist_mc.Scale(hist.Integral()/hist_mc.Integral())
        hist_mc.SetLineColor(1)
        hist_mc.SetStats(False)
        hist_mc.Draw("histsame")
        legend.AddEntry(hist_mc,mc[0][2]+(" GEN" if compare_RECOvsGEN else ""),"l")

        name=prefix+mc2[0][0].split("-HT")[0]+var+str(massbins[mass][0])+"_"+str(massbins[mass][1])
        if var=="#chi": name+="_rebin1"
        if var=="y":
          hist_mc2=f_mc2[0].Get(name.replace("y"+str(massbins[mass][0]),"y_{1}"+str(massbins[mass][0])))
          hist_mc2.Add(f_mc2[0].Get(name.replace("y"+str(massbins[mass][0]),"y_{2}"+str(massbins[mass][0]))))
        else:
          hist_mc2=f_mc2[0].Get(name)
        for i in range(1,len(mc2)):
            namei=prefix+mc2[i][0].split("-HT")[0]+var+str(massbins[mass][0])+"_"+str(massbins[mass][1])
            if var=="#chi": namei+="_rebin1"
            if var=="y":
              hist_mc2.Add(f_mc2[i].Get(namei.replace("y"+str(massbins[mass][0]),"y_{1}"+str(massbins[mass][0]))),mc2[i][1]/mc2[0][1])
              hist_mc2.Add(f_mc2[i].Get(namei.replace("y"+str(massbins[mass][0]),"y_{2}"+str(massbins[mass][0]))),mc2[i][1]/mc2[0][1])
            else:
              hist_mc2.Add(f_mc2[i].Get(namei),mc2[i][1]/mc2[0][1])
        if var=="#chi":
            hist_mc2=hist_mc2.Rebin(len(chi_binnings[mass])-1,hist_mc2.GetName()+"_rebin1",chi_binnings[mass])
            for b in range(hist_mc2.GetNbinsX()):
               hist_mc2.SetBinContent(b+1,hist_mc2.GetBinContent(b+1)/hist_mc2.GetBinWidth(b+1))
               hist_mc2.SetBinError(b+1,hist_mc2.GetBinError(b+1)/hist_mc2.GetBinWidth(b+1))
        if var.startswith("#phi"):
          hist_mc2=hist_mc2.Rebin(2)
        if hist_mc2.Integral()>0:
            hist_mc2.Scale(hist2.Integral()/hist_mc2.Integral())
        hist_mc2.SetLineColor(2)
        hist_mc2.SetStats(False)
        hist_mc2.Draw("histsame")
        legend.AddEntry(hist_mc2,mc2[0][2]+(" GEN" if compare_RECOvsGEN else ""),"l")

        if len(f_mc3)>0:
          name=prefix+mc3[0][0].split("-HT")[0]+var+str(massbins[mass][0])+"_"+str(massbins[mass][1])
          if var=="#chi": name+="_rebin1"
          if var=="y":
            hist_mc3=f_mc3[0].Get(name.replace("y"+str(massbins[mass][0]),"y_{1}"+str(massbins[mass][0])))
            hist_mc3.Add(f_mc3[0].Get(name.replace("y"+str(massbins[mass][0]),"y_{2}"+str(massbins[mass][0]))))
          else:
            hist_mc3=f_mc3[0].Get(name)
          for i in range(1,len(mc3)):
              namei=prefix+mc3[i][0].split("-HT")[0]+var+str(massbins[mass][0])+"_"+str(massbins[mass][1])
              if var=="#chi": namei+="_rebin1"
              if var=="y":
                hist_mc3.Add(f_mc3[i].Get(namei.replace("y"+str(massbins[mass][0]),"y_{1}"+str(massbins[mass][0]))),mc3[i][1]/mc3[0][1])
                hist_mc3.Add(f_mc3[i].Get(namei.replace("y"+str(massbins[mass][0]),"y_{2}"+str(massbins[mass][0]))),mc3[i][1]/mc3[0][1])
              else:
                hist_mc3.Add(f_mc3[i].Get(namei),mc3[i][1]/mc3[0][1])
          if var=="#chi":
              hist_mc3=hist_mc3.Rebin(len(chi_binnings[mass])-1,hist_mc3.GetName()+"_rebin1",chi_binnings[mass])
              for b in range(hist_mc3.GetNbinsX()):
                 hist_mc3.SetBinContent(b+1,hist_mc3.GetBinContent(b+1)/hist_mc3.GetBinWidth(b+1))
                 hist_mc3.SetBinError(b+1,hist_mc3.GetBinError(b+1)/hist_mc3.GetBinWidth(b+1))
          if var.startswith("#phi"):
             hist_mc3=hist_mc3.Rebin(2)
          if hist_mc3.Integral()>0:
              hist_mc3.Scale(hist3.Integral()/hist_mc3.Integral())
          hist_mc3.SetLineColor(4)
          hist_mc3.SetStats(False)
          hist_mc3.Draw("histsame")
          legend.AddEntry(hist_mc3,mc3[0][2]+(" GEN" if compare_RECOvsGEN else ""),"l")

        hist.Draw("lesame")
        hist2.Draw("lesame")
        if len(f_data3)>0: hist3.Draw("lesame")
        hist.Draw("axissame")

        legend.SetTextSize(0.04)
        legend.SetFillStyle(0)
        legend.Draw("same")

        if compare_Run2Run3:
          histbackup=hist
          hist=hist_mc # use MC as reference for first ratio plot
          histbackup2=hist2
          hist2=hist_mc
          if len(f_mc3)>0:
            histbackup3=hist3
            hist3=hist_mc

        canvas.cd(2)
        ratio=hist.Clone(hist.GetName()+"ratio")
        ratio.Divide(hist,hist)
        for b in range(hist.GetNbinsX()):
          if hist.GetBinContent(b+1)>0:
            ratio.SetBinError(b+1,hist.GetBinError(b+1)/hist.GetBinContent(b+1))
        ratio.SetTitle("")
        ratio.GetYaxis().SetTitle("Run2 / year" if compare_years else ("Sim / "+mc[0][2] if compare_Run2Run3 else ("EOY / UL" if compare_EOYvsUL or compare_EOYvsUL_MC else ("new / old" if compare_NULvsUL else ("RECO / GEN" if compare_RECOvsGEN else "Sim / Data")))))
        ratio.GetYaxis().SetTitleSize(0.18/(len(ratio.GetYaxis().GetTitle())/10))
        ratio.GetYaxis().SetTitleOffset(0.3)
        ratio.SetMarkerSize(0.1)
        ratio.GetYaxis().SetLabelSize(0.2)
        ratio.GetYaxis().SetRangeUser(0,2)
        if var=="#chi" and mass<=7:
          ratio.GetYaxis().SetRangeUser(0.5,1.5)
        if var=="#chi" and mass<=5:
          ratio.GetYaxis().SetRangeUser(0.8,1.2)
        if compare_EOYvsUL or compare_NULvsUL or compare_years or compare_Run2Run3:
         if var in ["#chi","y_{boost}","p_{T1}","p_{T2}","y_{1}","y_{2}","y","#phi_{1}","#phi_{2}"] and mass<=7:
          ratio.GetYaxis().SetRangeUser(0.5,1.5)
         if var=="#chi" and mass<=5:
          ratio.GetYaxis().SetRangeUser(0.9,1.1)
         if var=="#chi" and mass<=4:
          ratio.GetYaxis().SetRangeUser(0.95,1.05)
        ratio.GetXaxis().SetNdivisions(506)
        ratio.GetYaxis().SetNdivisions(503)
        ratio.GetXaxis().SetLabelColor(0)
        ratio.GetXaxis().SetTitleSize(0.2)
        ratio.GetXaxis().SetTitleOffset(1.1)
        ratio.GetXaxis().SetLabelSize(0.18)
        ratio.Draw("histe")
        ratio_mc=hist_mc.Clone(hist_mc.GetName()+"ratio")
        ratio_mc.Divide(hist_mc,hist)
        for b in range(hist.GetNbinsX()):
          if hist.GetBinContent(b+1)>0:
            ratio_mc.SetBinError(b+1,hist.GetBinError(b+1)/hist.GetBinContent(b+1))
        ratio_mc.Draw("histsame")
        ratio_mc2=hist_mc2.Clone(hist_mc2.GetName()+"ratio")
        ratio_mc2.Divide(hist_mc2,hist2)
        for b in range(hist2.GetNbinsX()):
          if hist2.GetBinContent(b+1)>0:
            ratio_mc2.SetBinError(b+1,hist2.GetBinError(b+1)/hist2.GetBinContent(b+1))
        ratio_mc2.Draw("histsame")
        if len(f_mc3)>0: 
          ratio_mc3=hist_mc3.Clone(hist_mc3.GetName()+"ratio")
          ratio_mc3.Divide(hist_mc3,hist3)
          for b in range(hist3.GetNbinsX()):
            if hist3.GetBinContent(b+1)>0:
              ratio_mc3.SetBinError(b+1,hist3.GetBinError(b+1)/hist3.GetBinContent(b+1))
          ratio_mc3.Draw("histsame")
          ratio.Draw("axissame")
        if compare_Run2Run3:
          hist=histbackup
          hist2=histbackup2
          if len(f_mc3)>0:
            hist3=histbackup3

        canvas.cd(3)
        ddratio=hist.Clone(hist.GetName()+"ddratio")
        ddratio.Divide(hist,hist)
        for b in range(hist.GetNbinsX()):
          if hist.GetBinContent(b+1)>0:
            ddratio.SetBinError(b+1,hist.GetBinError(b+1)/hist.GetBinContent(b+1))
        ddratio.SetTitle("")
        if compare_EOYvsUL_MC or compare_RECOvsGEN:
          ddratio.GetYaxis().SetTitle("Sim / "+data[0][2]+"")
        else:
          ddratio.GetYaxis().SetTitle("Data / "+data[0][2]+"")
        ddratio.GetYaxis().SetTitleSize(0.11/(len(ddratio.GetYaxis().GetTitle())/10))
        ddratio.GetYaxis().SetTitleOffset(0.5)
        ddratio.SetMarkerSize(0.1)
        ddratio.GetYaxis().SetLabelSize(0.13)
        ddratio.GetYaxis().SetRangeUser(0,2)
        if var in ["#chi","y_{boost}","p_{T1}","p_{T2}","y_{1}","y_{2}","y","#phi_{1}","#phi_{2}"] and mass<=7:
          ddratio.GetYaxis().SetRangeUser(0.5,1.5)
        if var=="#chi" and mass<=5:
          ddratio.GetYaxis().SetRangeUser(0.9,1.1)
        if var=="#chi" and mass<=4:
          ddratio.GetYaxis().SetRangeUser(0.95,1.05)
        ddratio.GetXaxis().SetNdivisions(506)
        ddratio.GetYaxis().SetNdivisions(503)
        ddratio.GetXaxis().SetLabelColor(1)
        ddratio.GetXaxis().SetTitle(label[variables.index(var)])
        ddratio.GetXaxis().SetTitleSize(0.14)
        ddratio.GetXaxis().SetTitleOffset(1.1)
        ddratio.GetXaxis().SetLabelSize(0.12)
        ddratio.Draw("histe")
        ddratio2=hist2.Clone(hist2.GetName()+"ddratio")
        ddratio2.Divide(hist2,hist)
        for b in range(hist.GetNbinsX()):
          if hist.GetBinContent(b+1)>0:
            ddratio2.SetBinError(b+1,hist.GetBinError(b+1)/hist.GetBinContent(b+1))
        ddratio2.Draw("histsame")
        if len(f_data3)>0: 
          ddratio3=hist3.Clone(hist3.GetName()+"ddratio")
          ddratio3.Divide(hist3,hist)
          for b in range(hist.GetNbinsX()):
            if hist.GetBinContent(b+1)>0:
              ddratio3.SetBinError(b+1,hist.GetBinError(b+1)/hist.GetBinContent(b+1))
          ddratio3.Draw("histsame")
        ratio.Draw("axissame")

        canvas.cd(1)
        hist.GetYaxis().SetTitleOffset(1.2)

        canvas.SaveAs("plots/chi_control_plots_"+var.replace("_","").replace("{","").replace("}","").replace("#","").replace("+","p").replace("-","m")+"_"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+version+".root")
        canvas.SaveAs("plots/chi_control_plots_"+var.replace("_","").replace("{","").replace("}","").replace("#","").replace("+","p").replace("-","m")+"_"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+version+".pdf")
