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
    for b in range(h1.GetXaxis().GetNbins()):
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)/h1.GetBinWidth(b+1))
        h1.SetBinError(b+1,h1.GetBinError(b+1)/h1.GetBinWidth(b+1))
    fit=TF1(h1.GetName()+"smooth",f,1,16)
    h1.Fit(fit,"RNQWW")
    for chi_bin in range(h1.GetXaxis().GetNbins()):
      h1.SetBinContent(chi_bin+1,fit.Eval(h1.GetBinCenter(chi_bin+1)))
    for b in range(h1.GetXaxis().GetNbins()):
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)*h1.GetBinWidth(b+1))
        h1.SetBinError(b+1,h1.GetBinError(b+1)*h1.GetBinWidth(b+1))
    return h1

def smoothChi(h1):
    for b in range(h1.GetXaxis().GetNbins()):
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)/h1.GetBinWidth(b+1))
        h1.SetBinError(b+1,h1.GetBinError(b+1)/h1.GetBinWidth(b+1))
    fit=TF1(h1.GetName()+"smooth","pol2",3.5,16)
    h1.Fit(fit,"RNQ")
    for chi_bin in range(h1.FindBin(4.5),h1.GetXaxis().GetNbins()):
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
    use_UL=True
    responsePostfix="" # "", "herwigpp", "Test", "Train", "GS", "SysUp", "SysDown", "RECO"
    run="2" # "2" for Run2 or "3" for 13.6 TeV projection
    use_NNPDF3=True
    use_CP2=True # for ADD
    
    if use_UL:
      years=["UL16preVFP","UL16postVFP","UL17","UL18"]
    else:
      years=["2016","2017","2018"]

    if useNNLO:
      if use_NNPDF3:
        pdfset="nn31nnlo"
      else:
        pdfset="ct14nnlo"
    else:
      pdfset="ct14nlo"
    if useM2:
      muScale="m2"
      muAltScale="pt12"
    else:
      muScale="pt12"
      muAltScale="m2"

    if "ct14" in pdfset:
      PDFmembers=56
    else:
      PDFmembers=100

    if run=="2":
      if use_NNPDF3:
        prefixs=["versions/run2ULNNLO_m2_NNPDF3/datacard_shapelimit13TeV"]
      elif useM2:
        prefixs=["versions/run2ULNNLO_m2/datacard_shapelimit13TeV"]
      elif useNNLO:
        prefixs=["versions/run2ULNNLO_pt12/datacard_shapelimit13TeV"]
      else:
        prefixs=["versions/run2ULNLO_pt12/datacard_shapelimit13TeV"]
    elif run=="3":
      prefixs=["versions/run3NNLO_m2/datacard_shapelimit13TeV"]
    else:
      whatprefix
    input_prefix1="datacards/datacard_shapelimit13TeV" # GEN signals
    input_prefix2="versions/2016NLO/datacard_shapelimit13TeV" # other signals
    prefixs[0]+=responsePostfix
 
    # negligible jes sources removed
    jessources=["AbsoluteScale",
                "AbsoluteMPFBias",
                "SinglePionHCAL",
                "FlavorQCD"]+\
                ["TimePtEta"+y for y in years]+\
                ["RelativePtBB",
                "RelativeBal",
                "RelativeFSR"]+\
                ["RelativeSample"+y for y in years]+\
                ["RelativeStatEC"+y for y in years]+\
                ["RelativeJEREC1"+y for y in years]+\
                ["PileUpDataMC",
                "PileUpPtRef",
                "PileUpPtBB",
                #"PileUpPtEC1",
                ]+["PileUpPtEC1"+y for y in years]+\
                ["SumInQuadrature",
                ]
    print(len(jessources)-1,"jes sources")
    jersources=["JER"+y for y in years]+\
               ["SumInQuadrature",
                ] #["JER1"+y for y in years]+["JER2"+y for y in years]+

    if run=="3":
      years=["2024"]

    colors=[1,2,3,4,6,7,8,9,12,28,34,38,40,41,42,43,44,45,46,47,48,49] #11,20
 
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
    mass_binning=array.array('d')
    for mass_bin in massbins:
        mass_binning.append(mass_bin[0])
    mass_binning.append(13000)
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

    samples+=[("QCD",[]),
            ]
    samples+=[("QCDAntiCIplusLL12000",[]),
             ]
    samples+=[("QCDCIplusLL12000",[]),
             ]
    samplesSkip=[("QCDCIplusLL8000",[]),
             ("QCDCIplusLL9000",[]),
             ("QCDCIplusLL10000",[]),
             ("QCDCIplusLL11000",[]),
             ("QCDCIplusLL12000",[]),
             ("QCDCIplusLL13000",[]),
             ("QCDCIplusLL14000",[]),
             ("QCDCIplusLL16000",[]),
             ("QCDCIplusLL18000",[]),
             ]
    samplesSkip+=[("QCDCIminusLL8000",[]),
             ("QCDCIminusLL9000",[]),
             ("QCDCIminusLL10000",[]),
             ("QCDCIminusLL11000",[]),
             ("QCDCIminusLL12000",[]),
             ("QCDCIminusLL13000",[]),
             ("QCDCIminusLL14000",[]),
             ("QCDCIminusLL16000",[]),
             ("QCDCIminusLL18000",[]),
             ]
    if run=="3":
       minMasses=[1500,1900,2400,2800,3300,3800,4300,5200,6000] # for mass bins 1.9, 2.4, 3.0, 3.6, 4.2, 4.8, 5.4, 6.0, 7.0
       maxMasses=[1900,2400,2800,3300,3800,4300,5200,6000,13600] # for mass bins 1.9, 2.4, 3.0, 3.6, 4.2, 4.8, 5.4, 6.0, 7.0
       lambdaTes=[9000,10000,11000,12000,13000,14000,15000,16000,17000,18000,19000,20000,21000,22000,25000,30000]
       couplings=[(0,0,0,1),]
       for lambdaT in lambdaTes:
         samples1+=[("QCDADD"+str(lambdaT),[("pythia8_add_m"+str(minMass)+"_"+str(maxMasses[minMasses.index(minMass)])+"_"+str(lambdaT)+"_0_0_0_1_13p6TeV_Nov2022",1) for minMass in minMasses])]
    elif use_NNPDF3:
       minMasses=[1500,1900,2400,2800,3300,3800,4300,5200,6000] # for mass bins 1.9, 2.4, 3.0, 3.6, 4.2, 4.8, 5.4, 6.0, 7.0
       maxMasses=[1900,2400,2800,3300,3800,4300,5200,6000,13000] # for mass bins 1.9, 2.4, 3.0, 3.6, 4.2, 4.8, 5.4, 6.0, 7.0
       lambdaTes=[9000,10000,11000,12000,13000,14000,15000,16000,17000,18000,19000,20000,21000,22000,25000,30000]
       couplings=[(0,0,0,1),]
       for lambdaT in lambdaTes:
         samples1+=[("QCDADD"+str(lambdaT),[("pythia8_add_m"+str(minMass)+"_"+str(maxMasses[minMasses.index(minMass)])+"_"+str(lambdaT)+"_0_0_0_1_13TeV_Feb2025"+("CP2" if use_CP2 else "CP5"),1) for minMass in minMasses])]
    else:
      samples1+=[#("QCDADD6000",[]),
             #("QCDADD7000",[]),
             #("QCDADD8000",[]),
             ("QCDADD9000",[]),
             ("QCDADD10000",[]),
             ("QCDADD11000",[]),
             ("QCDADD12000",[]),
             ("QCDADD13000",[]),
             ("QCDADD14000",[]),
             ("QCDADD15000",[]),
             ("QCDADD16000",[]),
             ("QCDADD17000",[]),
             ("QCDADD18000",[]),
             ("QCDADD19000",[]),
             ("QCDADD20000",[]),
             #("QCDADD21000",[]),
             #("QCDADD22000",[]),
             ]

    #for m in range(5,31):
    #for m in [13]:
    #for m in [10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30]:
    for m in [10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,35,40,45,50,55,60,65,70]:
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
      for weight in ['gdmv_0_gdma_1p0_gv_0_ga_0p01', 'gdmv_0_gdma_1p0_gv_0_ga_0p05', 'gdmv_0_gdma_1p0_gv_0_ga_0p1', 'gdmv_0_gdma_1p0_gv_0_ga_0p2', 'gdmv_0_gdma_1p0_gv_0_ga_0p25', 'gdmv_0_gdma_1p0_gv_0_ga_0p3', 'gdmv_0_gdma_1p0_gv_0_ga_0p5', 'gdmv_0_gdma_1p0_gv_0_ga_0p75', 'gdmv_0_gdma_1p0_gv_0_ga_1', 'gdmv_0_gdma_1p0_gv_0_ga_1p5', 'gdmv_0_gdma_1p0_gv_0_ga_2p0', 'gdmv_0_gdma_1p0_gv_0_ga_2p5', 'gdmv_0_gdma_1p0_gv_0_ga_3p0']:
      #for weight in ['gdmv_0_gdma_1p0_gv_0_ga_0p1', 'gdmv_0_gdma_1p0_gv_0_ga_0p2', 'gdmv_0_gdma_1p0_gv_0_ga_0p25', 'gdmv_0_gdma_1p0_gv_0_ga_0p3', 'gdmv_0_gdma_1p0_gv_0_ga_0p5', 'gdmv_0_gdma_1p0_gv_0_ga_0p75', 'gdmv_0_gdma_1p0_gv_0_ga_1', 'gdmv_0_gdma_1p0_gv_0_ga_1p5']:
      #for weight in ['gdmv_0_gdma_1p0_gv_0_ga_1p5']:
         samples4+=[("DMAxial_Dijet_LO_Mphi_"+str(mass)+"_"+str(mDM)+"_1p0_1p0_Mar5_"+weight,[("DMAxial_Dijet_LO_Mphi_"+str(mass)+"_"+str(mDM)+"_1p0_1p0_Mar5_"+weight,0)]),
             ]

    if use_NNPDF3:
       xsecs_qbhadd6={7500.0: 0.01255937, 8000.0: 0.003966422, 8500.0: 0.001199453, 9000.0: 0.0003439181, 9500.0: 9.214433e-05, 10000.0: 2.255567e-05, 11000.0: 8.795348e-07, 12000.0: 1.020805e-08}
       for signalMass in [7500,8000,8500,9000,9500,10000,11000,12000]:
         samples3+=[("QBH_ADD6_"+str(signalMass),[('qbh_ADD6_'+str(+signalMass)+"_Feb2025CP5",1.*xsecs_qbhadd6[signalMass]),])] # Factor 3 to reproduce old results. To be fixed.
       xsecs_qbhrs1={4500.0: 0.04722447, 5000.0: 0.0164621, 5500.0: 0.005754282, 6000.0: 0.001997539, 6500.0: 0.0006832241, 7000.0: 0.0002286843, 7500.0: 7.442768e-05, 8000.0: 2.339223e-05}
       for signalMass in [4500,5000,5500,6000,6500,7000,7500,8000]:
         samples3+=[("QBH_RS1_"+str(signalMass),[('qbh_RS1_'+str(signalMass)+"_Feb2025CP5",1.*xsecs_qbhrs1[signalMass]),])] # Factor 3 to reproduce old results. To be fixed.
    else:
       for m in [[7500,0.01678352],[8000,0.004871688],[8500,0.001292072],[9000,0.0003054339],[9500,0.00006221544],[10000,0.00001040396],[10500,0.000001327101],[11000,0.0000001145733]]:
         samples3+=[("QBH_"+str(m[0])+"_6",[("QBH_"+str(m[0])+"_6",m[1])]),]

       for m in [[4500,0.05148],[5000,0.01829],[5500,0.006472],[6000,0.002250],[6500,0.0007599],[7000,0.0002461]]:
         samples3+=[("QBH_"+str(m[0])+"_RS1",[("QBH_"+str(m[0])+"_RS1",m[1])]),]

    for weight in ['fa1000','fa1500','fa2000','fa2500','fa3000','fa3500','fa4000','fa4500','fa5000','fa50000']:
         samples5+=[("alp_QCD_"+weight,[("alp_QCD_"+weight,0)]),]

    for weight in ["CG0p1","CG0p05","CG0p04","CG0p03","CG0p025","CG0p02","CG0p015","CG0p01","CG0p0075","CG0p005","CG0p0025","CG0p0"]:
         samples6+=[("tripleG_QCD_"+weight,[("tripleG_QCD_"+weight,0)]),]

    samples7=[]
    xsecs_qbh={'QBH_MD4000_MBH5000_n4': 1.42769, 'QBH_MD4000_MBH5000_n6': 3.15044, 'QBH_MD4000_MBH5000_n2': 0.316956, 'QBH_MD2000_MBH3000_n2': 57.4928, 'QBH_MD2000_MBH3000_n4': 232.017, 'QBH_MD2000_MBH3000_n6': 496.996, 'QBH_MD6000_MBH7000_n2': 0.00238888, 'QBH_MD6000_MBH7000_n4': 0.0112931, 'QBH_MD6000_MBH7000_n6': 0.0253378, 'QBH_MD7000_MBH8000_n2': 0.000175315, 'QBH_MD7000_MBH8000_n6': 0.00189162, 'QBH_MD7000_MBH8000_n4': 0.000834913, 'QBH_MD8000_MBH9000_n4': 4.61159e-05, 'QBH_MD8000_MBH9000_n6': 0.000104591, 'QBH_MD8000_MBH9000_n2': 9.4628e-06, 'QBH_MD3000_MBH4000_n2': 3.74809, 'QBH_MD3000_MBH4000_n6': 35.319, 'QBH_MD3000_MBH4000_n4': 16.1407, 'QBH_MD5000_MBH6000_n2': 0.0283503, 'QBH_MD5000_MBH6000_n6': 0.287511, 'QBH_MD5000_MBH6000_n4': 0.129882, 'QBH_MD9000_MBH10000_n4': 1.23232e-06, 'QBH_MD9000_MBH10000_n6': 2.79188e-06, 'QBH_MD9000_MBH10000_n2': 2.53752e-07}
    for md in [2000,3000,4000,5000,6000,7000,8000,9000]:
      for n in [2,4,6]:
        samples7+=[("QBH_MD"+str(md)+"_MBH"+str(md+1000)+"_n"+str(n),[("UL18_QBH_MD"+str(md)+"_MBH"+str(md+1000)+"_n"+str(n),xsecs_qbh["QBH_MD"+str(md)+"_MBH"+str(md+1000)+"_n"+str(n)])])]

    # all samples
    samples=samples+samples1+samples2+samples3+samples4+samples5+samples6+samples7
    # for add
    #samples=samples1
    # for ci
    #samples=samples2
    # for qbh
    #samples=samples3
    # for alp+tripleG
    #samples=samples5
    #samples=samples6
    # for QBH check
    #samples=samples7
    # for postfit plots
    #samples=[("DMAxial_Dijet_LO_Mphi_7000_4000_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_1",[("DMAxial_Dijet_LO_Mphi_7000_4000_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_1",0)]), ]
    #samples=[("DMAxial_Dijet_LO_Mphi_2000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_1",[("DMAxial_Dijet_LO_Mphi_2000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_1",0)]), ]
    #samples=[("DMAxial_Dijet_LO_Mphi_5000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_1",[("DMAxial_Dijet_LO_Mphi_5000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_1",0)]), ]
    #samples=[("cs_ct14nnlo_14000_V-A-",[])]
    
    if len(sys.argv)>1:
      if int(sys.argv[1])>=len(samples):
        samples=[]
      else:
        samples=[samples[int(sys.argv[1])]]

    print(len(samples))

    dataevents={}
    datahist={}
    data=None
    dataplot={}
    qcdnorm={}
    qcd13TeVnorm={}
    nloqcdnormraw={}
    nloqcdnorm={}
    cinorm={}
    for prefix in prefixs: 
     # signal cards
     for i in range(len(samples)):
      if samples[i][0]=="QCDCIplusLL8000":
        sample=prefix + '_GENnp-0-run'+run+'_chi.root'
      elif samples[i][0]=="QCDCIplusLL9000":
        sample=prefix + '_GENnp-1-run'+run+'_chi.root'
      elif samples[i][0]=="QCDCIplusLL10000":
        sample=prefix + '_GENnp-2-run'+run+'_chi.root'
      elif samples[i][0]=="QCDCIplusLL11000":
        sample=prefix + '_GENnp-3-run'+run+'_chi.root'
      elif samples[i][0]=="QCDCIplusLL12000":
        sample=prefix + '_GENnp-4-run'+run+'_chi.root'
      elif samples[i][0]=="QCDCIplusLL13000":
        sample=prefix + '_GENnp-5-run'+run+'_chi.root'
      elif samples[i][0]=="QCDCIplusLL14000":
        sample=prefix + '_GENnp-6-run'+run+'_chi.root'
      elif samples[i][0]=="QCDCIplusLL16000":
        sample=prefix + '_GENnp-7-run'+run+'_chi.root'
      elif samples[i][0]=="QCDCIplusLL18000":
        sample=prefix + '_GENnp-8-run'+run+'_chi.root'
      elif samples[i][0]=="QCDCIminusLL8000":
        sample=prefix + '_GENnp-9-run'+run+'_chi.root'
      elif samples[i][0]=="QCDCIminusLL9000":
        sample=prefix + '_GENnp-10-run'+run+'_chi.root'
      elif samples[i][0]=="QCDCIminusLL10000":
        sample=prefix + '_GENnp-11-run'+run+'_chi.root'
      elif samples[i][0]=="QCDCIminusLL11000":
        sample=prefix + '_GENnp-12-run'+run+'_chi.root'
      elif samples[i][0]=="QCDCIminusLL12000":
        sample=prefix + '_GENnp-13-run'+run+'_chi.root'
      elif samples[i][0]=="QCDCIminusLL13000":
        sample=prefix + '_GENnp-14-run'+run+'_chi.root'
      elif samples[i][0]=="QCDCIminusLL14000":
        sample=prefix + '_GENnp-15-run'+run+'_chi.root'
      elif samples[i][0]=="QCDCIminusLL16000":
        sample=prefix + '_GENnp-16-run'+run+'_chi.root'
      elif samples[i][0]=="QCDCIminusLL18000":
        sample=prefix + '_GENnp-17-run'+run+'_chi.root'
      elif "ADD" in samples[i][0] and run=="3":
        sample=prefix + "_GENnp-run"+run+"-" + samples[i][0] + '_chi.root'
      elif ("ADD" in samples[i][0] or "QBH" in samples[i][0]) and use_NNPDF3:
        sample=prefix + "_GENnp-"+("CP2" if use_CP2 else "CP5")+"-run"+run+"-" + samples[i][0] + '_chi.root'
      elif samples[i][0]=="QCDADD6000":
        sample=prefix + '_GENnp-18-run'+run+'_chi.root'
      elif samples[i][0]=="QCDADD7000":
        sample=prefix + '_GENnp-19-run'+run+'_chi.root'
      elif samples[i][0]=="QCDADD8000":
        sample=prefix + '_GENnp-20-run'+run+'_chi.root'
      elif samples[i][0]=="QCDADD9000":
        sample=prefix + '_GENnp-21-run'+run+'_chi.root'
      elif samples[i][0]=="QCDADD10000":
        sample=prefix + '_GENnp-22-run'+run+'_chi.root'
      elif samples[i][0]=="QCDADD11000":
        sample=prefix + '_GENnp-23-run'+run+'_chi.root'
      elif samples[i][0]=="QCDADD12000":
        sample=prefix + '_GENnp-24-run'+run+'_chi.root'
      elif samples[i][0]=="QCDADD13000":
        sample=prefix + '_GENnp-25-run'+run+'_chi.root'
      elif samples[i][0]=="QCDADD14000":
        sample=prefix + '_GENnp-26-run'+run+'_chi.root'
      elif samples[i][0]=="QCDADD15000":
        sample=prefix + '_GENnp-27-run'+run+'_chi.root'
      elif samples[i][0]=="QCDADD16000":
        sample=prefix + '_GENnp-28-run'+run+'_chi.root'
      elif samples[i][0]=="QCDADD17000":
        sample=prefix + '_GENnp-29-run'+run+'_chi.root'
      elif samples[i][0]=="QCDADD18000":
        sample=prefix + '_GENnp-30-run'+run+'_chi.root'
      elif samples[i][0]=="QCDADD19000":
        sample=prefix + '_GENnp-31-run'+run+'_chi.root'
      elif samples[i][0]=="QCDADD20000":
        sample=prefix + '_GENnp-32-run'+run+'_chi.root'
      elif samples[i][0]=="QCDADD21000":
        sample=prefix + '_GENnp-33-run'+run+'_chi.root'
      elif samples[i][0]=="QCDADD22000":
        sample=prefix + '_GENnp-34-run'+run+'_chi.root'
      elif samples[i][0]=="QCDAntiCIplusLL12000":
        sample=prefix + '_GENnp-antici-run'+run+'_chi.root'
      elif samples[i][0]=="QCD":
        sample=prefix + '_GEN-QCD-run'+run+'_chi.root'
      elif "QBH_MD" in samples[i][0]:
        sample=prefix + "_run2_UL18_" + samples[i][0] + '-GEN_chi.root'
      elif "alp" in samples[i][0] or "tripleG" in samples[i][0] or "DM" in samples[i][0] or "ll" in samples[i][0] or "cs" in samples[i][0] or "wide" in samples[i][0] or "QBH" in samples[i][0]:
        sample=prefix + "_" + samples[i][0] + '-run'+run+'_chi.root'
      #if "ADD" in samples[i][0]:
      #  sample=prefix + '_GENaddv3_chi2016.root'
      #elif "CIplus" in samples[i][0]:
      #  sample=prefix + '_GENciv3_chi2016.root'
      #elif "CIminus" in samples[i][0]:
      #  sample=prefix + '_GENciminusv3_chi2016.root'
      else:
        sample=prefix + '_GEN-QCD-run'+run+'_chi.root'
      print("output file", sample)
      if "GEN" in sample:
        insignalsample=(sample.replace(prefix,input_prefix1).replace("_1-run3","_1-run2")) ###### FIXME when run3 samples exist
      else:
        insignalsample=(input_prefix2 + "_cs_ct14nnlo_30000_V-A+-run2_chi.root" if "cs" in sample else sample.replace(prefix,input_prefix2).replace("_1-run3","_1-run2")) ###### FIXME when run3 samples exist
      print("input signal file", insignalsample)
  
      insignalfile=TFile(insignalsample,'READ')
      closefiles=[insignalfile]
      out=TFile(sample,'RECREATE')
      closefiles+=[out]

      # RECO sample
      if responsePostfix=="RECO":
        recosample=insignalsample.replace("-GEN","")
        print("open reco file", recosample)
        freco=TFile(recosample,'READ')
        closefiles+=[freco]
      
      # LO QCD file
      sample2='datacards/datacard_shapelimit13TeV_GEN-QCD-run'+run+'_chi.root'
      print(sample2)
      in2=TFile(sample2,'READ')

      # LO QCD 13 TeV CP5 file to compute 13.6 TeV / 13 TeV k-factor
      sample13TeV='datacards/datacard_shapelimit13TeV_GEN-QCD-CP5-run2_chi.root'
      print(sample13TeV)
      in13TeV=TFile(sample13TeV,'READ')

     # Madgraph QCD file
      sampleMadgraph=input_prefix2 + '_alp_QCD_fa50000-run2_chi.root'
      print(sampleMadgraph)
      inMadgraph=TFile(sampleMadgraph,'READ')

      # data file
      #insample='datacards/chiHist_dataReReco_v3_PFHT900.root' #2016
      #insample='datacards/datacard_shapelimit13TeV_25nsData13combi_chi.root' # buggy data
      if run=="3":
        insample="data/datacard_shapelimit13TeV_run2_2024_chi.root" # ultra legacy reco
      elif use_UL:
        insample="datacards/datacard_shapelimit13TeV_run2_UL16preVFP_L1prefire_chi.root" # ultra legacy reco
      else:
        insample="datacards/datacard_shapelimit13TeV_run2_2016_L1prefire_chi.root" # Aug rereco
      print(insample)
      infile=TFile(insample,'READ')
      insample16postVFP="datacards/datacard_shapelimit13TeV_run2_UL16postVFP_L1prefire_chi.root"
      print(insample16postVFP)
      infile16postVFP=TFile(insample16postVFP,'READ')
      #insample17='datacards/uhh2.AnalysisModuleRunner.DATA.Run2017_RunBCDEF_17Nov2017-v1.root' #2017
      if use_UL:
        insample17="datacards/datacard_shapelimit13TeV_run2_UL17_L1prefire_chi.root"
      else:
        insample17="datacards/datacard_shapelimit13TeV_run2_2017_L1prefire_chi.root"
      print(insample17)
      infile17=TFile(insample17,'READ')
      #insample18='datacards/uhh2.AnalysisModuleRunner.DATA.Run2018_RunABCD_RunII_102X_v1.root' #2018
      if use_UL:
        insample18="datacards/datacard_shapelimit13TeV_run2_UL18_HEM_chi.root"
      else:
        insample18="datacards/datacard_shapelimit13TeV_run2_2018_HEM_chi.root"
      print(insample18)
      infile18=TFile(insample18,'READ')

      # unfolded data file
      #unfoldsample='datacards/Unfolded_chiNtuple_dataReReco_v3_Coarse_PFHT900_fromCB_AK4SF_pythia8_Pt_170toInf.root'
      unfoldsample='datacards/Unfolded_chiNtuple_dataReReco_v3_Coarse_PFHT900_fromCB_AK4SF_pythia8_Pt_170toInf_MatrixInvert.root'
      print(unfoldsample)
      unfoldfile=TFile(unfoldsample,'READ')

      # (N)NLO correction
      #filename1nu2="fastnlo/RunII/fnl5662j_v23_fix_CT14nlo_allmu_ak4.root"
      if useNNLO:
        #filename1nu2="fastnlo/NNLO/2jet.NNLO.fnl5662j_mjj_chi_ct14nnlo_cppread_mu_"+muScale+".root"
        filename1nu2="fastnlo/NNLO/2jet.NNLO.fnl5662j_mjj_chi_norm_v25_"+pdfset+"_cppread_mu_"+muScale+".root"
      else:
        filename1nu2="fastnlo/NNLO/2jet.NNLO.fnl5662j_mjj_chi_"+pdfset+"_cppread_mu_pt12.root"
      print(filename1nu2)
      nlofile2 = TFile.Open(filename1nu2)
      closefiles+=[nlofile2]

      # (N)NLO uncertainties
      filename1nusys="fastnlo/NNLO/newcifnl5662j_cs_"+pdfset+("_"+muScale if useNNLO else "")+"_70000_V-A-.root"
      print(filename1nusys)
      nlofilesys = TFile.Open(filename1nusys)
      closefiles+=[nlofilesys]
      if useNNLO:
        filename1altscale="fastnlo/NNLO/newcifnl5662j_cs_"+pdfset+("_"+muAltScale if useNNLO else "")+"_70000_V-A-.root"
      else:
        filename1altscale="fastnlo/NNLO/newcifnl5662j_cs_ct14nnlo_pt12_70000_V-A-.root" # NNLO is alternative to NLO
      print(filename1altscale)
      nlofilealtscale = TFile.Open(filename1altscale)
      closefiles+=[nlofilealtscale]
      
       # DM uncertainties
      filename1dmpdf="datacards/chi_dm_pdf_plots6000_13TeV_run2.root"
      print(filename1dmpdf)
      dmpdffile = TFile.Open(filename1dmpdf)
      closefiles+=[dmpdffile]

      # EWK correction
      filename1ewk="fastnlo/RunII/DijetAngularCMS13_ewk.root"
      print(filename1ewk)
      ewkfile = TFile.Open(filename1ewk)
      closefiles+=[ewkfile]

      # JES uncertainty QCD
      jesfiles=[]
      for n in range(5):
        filename1jes="plots/chi_systematic_plotschi_QCDmadgraphJES_"+str(n)+"_13TeV"+("_UL" if use_UL else "")+"_run2_PileUpPtEC1uncorrelated.root"
        print(filename1jes)
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
            print("JES source not found", source)
            error

      # JES uncertainty CI
      jescifiles=[]
      for n in range(5):
        filename1jesci="plots/chi_systematic_plotschi_QCDmadgraphJES_"+str(n)+"_13TeV"+("_UL" if use_UL else "")+"_run2_PileUpPtEC1uncorrelated.root"
        print(filename1jesci)
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
            print("JES source not found", source)
            error

      # JER uncertainty QCD
      jerfiles=[]
      for n in range(2):
        filename1jer="plots/chi_systematic_plotschi_QCDJER_"+str(n)+"_13TeV"+("_UL" if use_UL else "")+"_run2.root"
        print(filename1jer)
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
            print("JER source not found", source)
            error

      # JER uncertainty CI
      jercifiles=[]
      for n in range(2):
        filename1jerci="plots/chi_systematic_plotschi_QCDJER_"+str(n)+"_13TeV"+("_UL" if use_UL else "")+"_run2.root"
        print(filename1jerci)
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
            print("JER source not found", source)
            error

      # prefire uncertainty
      prefirefiles=[]
      filename1prefire="plots/chi_systematic_plotschi_run2prefire_13TeV_run2.root"
      print(filename1prefire)
      prefirefiles += [TFile.Open(filename1prefire)]
      closefiles+=[prefirefiles[-1]]
      prefirehists={}
      for j in range(len(massbins)):
          prefirehists[str(j)]=[]
          for prefirefile in prefirefiles:
            prefirepad=prefirefile.Get("prefire")
            prefire=prefirepad.GetListOfPrimitives()[j+useUnfoldedData*1]
            prefirehists[str(j)]+=[a for a in prefire.GetListOfPrimitives() if "TF1" in str(a) or "pol" in str(a)]
          if len(prefirehists[str(j)])!=2:
            print("prefire source not found",prefirehists)
            error
      
      # trigger efficiency
      trigger_files=[]
      trigger_files+=[TFile("plots/chi_trigger_plots_2016_SingleMuonHLT_PFHT900orHLT_PFHT800orHLT_PFJet450orHLT_PFJet500orHLT_CaloJet500_NoJetID_run2.root")]
      trigger_files+=[TFile("plots/chi_trigger_plots_2017_SingleMuonHLT_PFHT1050orHLT_PFJet500orHLT_PFJet550orHLT_CaloJet500_NoJetIDorHLT_CaloJet550_NoJetID_run2.root")]
      trigger_files+=[TFile("plots/chi_trigger_plots_2018_SingleMuonHLT_PFHT1050orHLT_PFJet500orHLT_PFJet550orHLT_CaloJet500_NoJetIDorHLT_CaloJet550_NoJetID_run2.root")]
      print(trigger_files)
      mass_spectrum=in2.Get("QCDmass")
      trigger_canvases=[]
      for trigger_file in trigger_files:
        trigger_canvases+=[trigger_file.Get("trigger")]
      trigger_efficiency=[]
      for tc in range(len(trigger_canvases)):
       trigger_efficiency+=[{}]
       for j in range(len(massbins)):
        trigger_efficiency[tc][j]={}
       for ft in trigger_canvases[tc].GetListOfPrimitives():
        if ("TF1" in str(ft) or "erf" in str(ft)) and "chi-" in str(ft):
          #ft.SetRange(massbins[0][0],massbins[-1][1])
          #ft.Update()
          #print ft.GetXmin(),ft.GetXmax()
          chi=int(str(ft).split("-")[-1].split('"')[0])
          mass_spectrum_eff=mass_spectrum.Clone(mass_spectrum.GetName()+"eff"+str(tc)+"_"+str(j)+"_"+str(chi))
          mass_spectrum_eff.Multiply(ft)
          for j in range(len(massbins)):
            trigger_efficiency[tc][j][int(chi)]=mass_spectrum_eff.Integral(mass_spectrum_eff.FindBin(massbins[j][0]),mass_spectrum_eff.FindBin(min(massbins[j][1],ft.GetXmax())))
            if ft.GetXmax()<massbins[j][1]:
              trigger_efficiency[tc][j][int(chi)]+=mass_spectrum.Integral(mass_spectrum.FindBin(ft.GetXmax()),mass_spectrum.FindBin(massbins[j][1]))
            #print j, chi, trigger_efficiency[tc][j][chi]
            #print massbins[j][0], massbins[j][1], ft.Eval(massbins[j][0]),ft.Eval(min(massbins[j][1],ft.GetXmax()))
            #if mass_spectrum.GetBinContent(mass_spectrum.FindBin(massbins[j][0]))>0 and mass_spectrum.GetBinContent(mass_spectrum.FindBin(massbins[j][1])):
            #  print mass_spectrum_eff.GetBinContent(mass_spectrum_eff.FindBin(massbins[j][0]))/mass_spectrum.GetBinContent(mass_spectrum.FindBin(massbins[j][0])),mass_spectrum_eff.GetBinContent(mass_spectrum_eff.FindBin(min(massbins[j][1],ft.GetXmax())))/mass_spectrum.GetBinContent(mass_spectrum.FindBin(min(massbins[j][1],ft.GetXmax())))
            refmass=mass_spectrum.Integral(mass_spectrum.FindBin(massbins[j][0]),mass_spectrum.FindBin(min(massbins[j][1],ft.GetXmax())))
            if ft.GetXmax()<massbins[j][1]:
              refmass+=mass_spectrum.Integral(mass_spectrum.FindBin(ft.GetXmax()),mass_spectrum.FindBin(massbins[j][1]))
            trigger_efficiency[tc][j][int(chi)]/=refmass
            #print j, chi, trigger_efficiency[tc][j][chi]

      canvas = TCanvas("","",0,0,800,600)
      canvas.Divide(4,3)
      plots=[]
      legends=[]
      
      datahist2d_chiBins2=TH2F("dijet_mass2_chi2","M_{jj} vs #chi -- DATA",len(mass_binning)-1, mass_binning, len(chi_binnings[0])-1, chi_binnings[0])
      datahist2d_chiBins3=TH2F("dijet_mass2_chi3","M_{jj} vs #chi -- DATA",len(mass_binning)-1, mass_binning, len(chi_binnings[-1])-1, chi_binnings[-1])

      for j in range(len(massbins)):
        col=1

        # QCD (empty background, not used in limit)
        histname='QCD#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_backup"
        print(histname)
        if only6000:
          qcd=in2.Get(histname.replace("6000_13000","6000_7000"))
        else:
          qcd=in2.Get(histname)
        out.cd()
        qcd.Write(histname.replace("_backup",""))
        qcd=qcd.Rebin(len(chi_binnings[j])-1,histname,chi_binnings[j])
        qcd.Scale(1e10*1e9) #mb -> pb and factor 1e-10 for backup
        if j in qcdnorm.keys():
           qcd.Scale(qcdnorm[j]/qcd.Integral())
        else:
           qcdnorm[j]=qcd.Integral()
        
        # Calculate ratio of 13.6 TeV and 13 TeV QCD LO cross section
        qcd13TeV=in13TeV.Get(histname)
        qcd13TeV=qcd13TeV.Rebin(len(chi_binnings[j])-1,histname,chi_binnings[j])
        qcd13TeV.Scale(1e10*1e9) #mb -> pb and factor 1e-10 for backup
        if j in qcd13TeVnorm.keys():
           qcd13TeV.Scale(qcd13TeVnorm[j]/qcd13TeV.Integral())
        else:
           qcd13TeVnorm[j]=qcd13TeV.Integral()
        factor13p6=qcd.Integral()/qcd13TeV.Integral()
        print("k-factor LO 13.6 TeV / 13 TeV", factor13p6)

        # NLO
        nloqcd=None
        for k in mass_bins_nlo_list[j]:
         #histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])
         histname='qcd_chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-1.0"
         print(histname)
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
        #if not useUnfoldedData:
        nloqcdbackup=nloqcd.Clone(nloqcd.GetName()+"_backup")
        nloqcd=smoothChi(nloqcd) # SMOOTH NNLO PREDICTION (FIX ME)
        if run=="3": # APPLY A SCALE FACTOR TO EXTRAPOLATE DATA AND PREDICTION TO 13.6 TEV
           nloqcdbackup.Scale(factor13p6)
        print("theory integral (pb):", nloqcdbackup.Integral())
        print("k-factor NNLO / LO", nloqcdbackup.Integral()/qcd.Integral())

        # NLO normalized
        nloqcdnormd=None
        for k in mass_bins_nlo_list[j]:
         histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])
         print(histname)
         hnlo = TH1F(nlofilesys.Get(histname))
         #hnlo.Scale(float(mass_bins_nlo3[k+1]-mass_bins_nlo3[k]))
         hnlo=rebin(hnlo,len(chi_binnings[j])-1,chi_binnings[j])
         if nloqcdnormd:
            nloqcdnormd.Add(hnlo)
         else:
            nloqcdnormd=hnlo

        # EWK corrections
        histname='chi-'+str(massbins[j]).strip("()").replace(',',"-").replace(' ',"").replace("1200-1500","1900-2400").replace("1500-1900","1900-2400").replace("6000-7000","6000-6600").replace("6000-13000","6000-6600").replace("7000-13000","6600-13000") # FIX compute lower and higher mass bins
        print(histname)
        ewk=ewkfile.Get(histname)
        for b in range(nloqcd.GetXaxis().GetNbins()):
            low_bin=ewk.FindBin(nloqcd.GetXaxis().GetBinLowEdge(b+1))
            up_bin=ewk.FindBin(nloqcd.GetXaxis().GetBinUpEdge(b+1))
            correction=ewk.Integral(low_bin,up_bin-1)/(up_bin-low_bin)
            #print "EWK", correction
            if not "EWK" in samples[i][0]:
               nloqcd.SetBinContent(b+1,nloqcd.GetBinContent(b+1)*correction)
        nloqcdnormraw[j]=nloqcdbackup.Integral()
        nloqcdnorm[j]=nloqcdbackup.Integral()
        nloqcdnorm[j]*=(1.-59.83/137.6*(1.57-0.87)/(2.*pi)) # SCALE CROSS SECTION TO ACCOUNT FOR HEM VETO
        nloqcdnorm[j]*=138000. # SCALE CROSS SECTION BY LUMI
        nloqcd.Scale(1./nloqcd.Integral())
        ewk.SetName("ewk-"+histname)

        # DATA
        histname="dijet_"+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_chi"
        print(histname)
        if useLensData:
          if "13000" in str(massbins[j]):
            histname2="dijet_m_chi_2__projY_"+str(massbins[j]).strip("()").replace(',',"-").replace(' ',"")
          else:
            histname2="dijet_m_chi_2__projY_"+str(massbins[j]).strip("()").replace(',',"-").replace(' ',"")
          print(histname2)
          data = TH1D(unfoldfile.Get(histname2))
          data.SetName(histname)
          data=data.Rebin(len(chi_binnings[j])-1,data.GetName()+"_rebin1",chi_binnings[j])
        elif useUnfoldedData:
          histname2="dijet_mass_"+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_chi_unfolded"
          print(histname2)
          data = TH1F(unfoldfile.Get(histname2))
          data.SetName(histname)
          data=data.Rebin(len(chi_binnings[j])-1,data.GetName()+"_rebin1",chi_binnings[j])
        else:
          #histname2="dijet_"+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"").replace("1200_1500","1900_2400").replace("1500_1900","1900_2400").replace("6000_7000","6000_13000").replace("7000_13000","6000_13000")+"_chi"
          histname2="datacard_shapelimit13TeV_run2_"+years[0]+"#chi"+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print(histname2)
          data = TH1F(infile.Get(histname2))
          data.SetName(histname)
          for b in range(data.GetXaxis().GetNbins()):
            datahist2d_chiBins2.Fill(massbins[j][0]+0.1,data.GetXaxis().GetBinCenter(b+1),data.GetBinContent(b+1))
            datahist2d_chiBins3.Fill(massbins[j][0]+0.1,data.GetXaxis().GetBinCenter(b+1),data.GetBinContent(b+1))
          data=data.Rebin(len(chi_binnings[j])-1,data.GetName()+"_rebin1",chi_binnings[j])
          if run=="2":
           if use_UL:
            histname16postVFP="datacard_shapelimit13TeV_run2_"+years[1]+"#chi"+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
            print(histname16postVFP)
            data16postVFP = TH1F(infile16postVFP.Get(histname16postVFP))
            for b in range(data.GetXaxis().GetNbins()):
              datahist2d_chiBins2.Fill(massbins[j][0]+0.1,data.GetXaxis().GetBinCenter(b+1),data16postVFP.GetBinContent(b+1))
              datahist2d_chiBins3.Fill(massbins[j][0]+0.1,data.GetXaxis().GetBinCenter(b+1),data16postVFP.GetBinContent(b+1))
            data16postVFP=data16postVFP.Rebin(len(chi_binnings[j])-1,data.GetName()+"_rebin1",chi_binnings[j])
            data.Add(data16postVFP)
           if j==0: data.Scale(40.77) # trigger prescales
           if j==1: data.Scale(40.77) # trigger prescales
           if j==2: data.Scale(14.44) # trigger prescales
           #histname17="Dijet/chi_"+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"").replace("1200_1500","1900_2400").replace("1500_1900","1900_2400").replace("6000_7000","6000_6600").replace("7000_13000","6600_13000")
           histname17="datacard_shapelimit13TeV_run2_"+years[-2]+"#chi"+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
           print(histname17)
           data17 = TH1F(infile17.Get(histname17))
           for b in range(data.GetXaxis().GetNbins()):
            datahist2d_chiBins2.Fill(massbins[j][0]+0.1,data.GetXaxis().GetBinCenter(b+1),data17.GetBinContent(b+1))
            datahist2d_chiBins3.Fill(massbins[j][0]+0.1,data.GetXaxis().GetBinCenter(b+1),data17.GetBinContent(b+1))
           data17=data17.Rebin(len(chi_binnings[j])-1,data.GetName()+"_rebin1",chi_binnings[j])
           if j==0: data17.Scale(60.17) # trigger prescales
           if j==1: data17.Scale(60.17) # trigger prescales
           if j==2: data17.Scale(23.76) # trigger prescales
           data.Add(data17)
           #histname18="Dijet/chi_"+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")
           histname18="datacard_shapelimit13TeV_run2_"+years[-1]+"#chi"+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
           print(histname18)
           data18 = TH1F(infile18.Get(histname18))
           for b in range(data.GetXaxis().GetNbins()):
            datahist2d_chiBins2.Fill(massbins[j][0]+0.1,data.GetXaxis().GetBinCenter(b+1),data18.GetBinContent(b+1))
            datahist2d_chiBins3.Fill(massbins[j][0]+0.1,data.GetXaxis().GetBinCenter(b+1),data18.GetBinContent(b+1))
           data18=data18.Rebin(len(chi_binnings[j])-1,data.GetName()+"_rebin1",chi_binnings[j])
           if j==0: data18.Scale(83.37) # trigger prescales
           if j==1: data18.Scale(83.37) # trigger prescales
           if j==2: data18.Scale(29.54) # trigger prescales
           data.Add(data18)
        #if run=="3": # APPLY A SCALE FACTOR TO EXTRAPOLATE DATA AND PREDICTION TO 13.6 TEV
        #   lumiratio=35./138. # LUMI ASSUMPTION FOR RUN3
        #   data.Scale(factor13p6*lumiratio)
        #   dn=data.Integral()
        #   data=nloqcd.Clone(histname)
        #   data.Scale(dn/nloqcd.Integral())
        dataevents[j]=data.Integral()
        print("data events:", dataevents[j], ", in pb:", dataevents[j]/138000.)
        out.cd()
        if not injectSignal:
         histname='data_obs#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
         data.Write(histname)

        # Madgraph reference sample
        histname='alp_QCD_fa50000#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_backup"
        print(histname)
        if only6000:
          qcdMadgraph=inMadgraph.Get(histname.replace("6000_13000","6000_7000"))
        else:
          qcdMadgraph=inMadgraph.Get(histname)

        # CI (=LO CI+NLO QCD)
        histname=samples[i][0].replace("Anti","")+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_backup"
        print(histname)
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
          filenamecinlo="fastnlo/NNLO/"+("newci" if massbins[j][0]>=3600 else "")+"fnl5662j_"+samples[i][0].replace("QCD","").replace("nnlo","nnlo_"+muScale)+".root" # calcTheoryUncert.py from Jingyu already gives normalized signals
          if massbins[j][0]<3600: # No CI calculation available
            filenamecinlo=filenamecinlo.replace("35000","30000").replace("40000","30000").replace("45000","30000").replace("50000","30000").replace("55000","30000").replace("60000","30000").replace("65000","30000").replace("70000","30000")
          print(filenamecinlo)
          cinlofile = TFile.Open(filenamecinlo)
          closefiles+=[cinlofile]
          filenamecinloaltscale="fastnlo/NNLO/"+("newci" if massbins[j][0]>=3600 else "")+"fnl5662j_"+samples[i][0].replace("QCD","").replace("nnlo","nnlo_"+muAltScale)+".root" # calcTheoryUncert.py from Jingyu already gives normalized signals
          if massbins[j][0]<3600:# No CI calculation available
            filenamecinloaltscale=filenamecinloaltscale.replace("35000","30000").replace("40000","30000").replace("45000","30000").replace("50000","30000").replace("55000","30000").replace("60000","30000").replace("65000","30000").replace("70000","30000")
          print(filenamecinloaltscale)
          cinlofilealtscale = TFile.Open(filenamecinloaltscale)
          closefiles+=[cinlofilealtscale]
          ci=None
          for k in mass_bins_nlo_list[j]:
            histname2='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])
            print(histname2)
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
          print(filenamecinlo)
          cinlofile = TFile.Open(filenamecinlo)
          closefiles+=[cinlofile]
          histname2="chi-"+str(massbinsci[0])+"-"+str(massbinsci[1])
          print(histname2)
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
            cibackup=insignalfile.Get(histname.replace("6000_13000","6000_7000"))
          else:
            cibackup=insignalfile.Get(histname)
          histname=histname.replace("_backup","")
          ci=cibackup.Clone(histname)
          #ci=smooth(ci,"pol3") # SMOOTH DM PREDICTION (FIX ME)
          ci=smoothChi(ci) # SMOOTH DM PREDICTION (FIX ME)
          ci=ci.Rebin(len(chi_binnings[j])-1,ci.GetName(),chi_binnings[j])
          ci.Scale(1./nloqcdbackup.Integral())
          print(histname,"signal fraction in first bin", ci.GetBinContent(1)/nloqcd.GetBinContent(1))
          ci.Add(nloqcd)
        elif "alp" in samples[i][0]:
          if only6000:
            cibackup=insignalfile.Get(histname.replace("6000_13000","6000_7000"))
          else:
            cibackup=insignalfile.Get(histname)
          histname=histname.replace("_backup","")
          ci=cibackup.Clone(histname)
          print("ALP:",ci.Integral(), "QCDMG:",qcdMadgraph.Integral(), "QCDPY:",qcd.Integral(), "QCDNLO:",nloqcdbackup.Integral())
          ci.Add(qcdMadgraph,-1)
          #ci=smooth(ci,"pol3") # SMOOTH ALP PREDICTION (FIX ME)
          ci=smoothChi(ci) # SMOOTH ALP PREDICTION (FIX ME)
          ci.Scale(nloqcdbackup.Integral()/qcdMadgraph.Integral()) # CORRECT FOR MISSING FACTOR LO->NNLO k-factor of ~10 in Madgraph QCD cross sections
          ci=ci.Rebin(len(chi_binnings[j])-1,ci.GetName(),chi_binnings[j])
          ci.Scale(1./nloqcdbackup.Integral())
          print(histname,"signal fraction in first bin", ci.GetBinContent(1)/nloqcd.GetBinContent(1))
          ci.Add(nloqcd)
        elif "QBH_MD" in samples[i][0]:
            histname=samples[i][0]+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
            histnamein='datacard_shapelimit13TeV_run2_UL18_QBH#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
            print(histnamein)
            cibackup=insignalfile.Get(histnamein)
            ci=cibackup.Clone(histname)
            ci=ci.Rebin(len(chi_binnings[j])-1,ci.GetName(),chi_binnings[j])
            ci.Scale(samples[i][1][0][1])
            ci.Scale(1./nloqcdbackup.Integral())
            ci.Add(nloqcd)
        elif "QBH" in samples[i][0] and use_NNPDF3:
            histname=samples[i][0]+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
            cibackup=insignalfile.Get(histname)
            histname=histname.replace("_backup","")
            ci=cibackup.Clone(histname)
            ci=ci.Rebin(len(chi_binnings[j])-1,ci.GetName(),chi_binnings[j])
            ci.Scale(samples[i][1][0][1])
            ci.Scale(1./nloqcdbackup.Integral())
            ci.Add(nloqcd)
        elif "QBH" in samples[i][0]:
            histname=samples[i][0]+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
            if j<5:
                ci=nloqcd.Clone(histname)
            else:
                if samples[i][0].split("_")[2]=='6':
                    histnamein='QCDADD'+samples[i][0].split("_")[2]+samples[i][0].split("_")[0]+samples[i][0].split("_")[1]+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")
                else:
                    histnamein='QCD'+samples[i][0].split("_")[2]+samples[i][0].split("_")[0]+samples[i][0].split("_")[1]+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")
                print(histnamein)
                cibackup=insignalfile.Get(histnamein.replace("6000_7000","6000_13000").replace("7000_13000","6000_13000")) # FIX compute QBH for highest mass bins
                ci=cibackup.Clone(histname)
                ci=ci.Rebin(len(chi_binnings[j])-1,ci.GetName(),chi_binnings[j])
                ci.Scale(samples[i][1][0][1]/1000000)
                ci.Scale(1./nloqcdbackup.Integral())
                ci.Add(nloqcd)
        elif "wide" in samples[i][0]:
          cibackup=insignalfile.Get(histname)
          try:
              histname=cibackup.GetName().replace("_backup","")
          except:
            print("problem reading", histname)
            break
          ci=cibackup.Clone(histname)
          ci=ci.Rebin(len(chi_binnings[j])-1,ci.GetName(),chi_binnings[j])
          ci.Scale(10./nloqcdbackup.Integral()) # make in units if 10pb
          ci.Add(nloqcd)
        elif "CIplus" in samples[i][0]:
          print("CREATE FAKE SIGNAL")
          if True:#j<=3: # fake signal for signficances
            histnamealt=samples[i][0].replace("Anti","")+'#chi'+str((4800,5400)).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_backup"
          else:
            histnamealt=histname
          cibackup=insignalfile.Get(histnamealt).Clone(histname)
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
        elif samples[i][0]=="QCD":
          if only6000:
            cibackup=insignalfile.Get(histname.replace("6000_13000","6000_7000"))
          else:  
            cibackup=insignalfile.Get(histname)
          histname=histname.replace("_backup","")
          ci=cibackup.Clone(histname)
          ci=ci.Rebin(len(chi_binnings[j])-1,ci.GetName(),chi_binnings[j])
          ci.Scale(1e9) #mb -> pb
          ci.Scale(1./ci.Integral())
        else:
          if only6000:
            cibackup=insignalfile.Get(histname.replace("6000_13000","6000_7000"))
          else:  
            cibackup=insignalfile.Get(histname)
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
        #cibackup.Write()
        ci.Clone(histname+"_nosmear").Write(histname+"_nosmear")
        cinonorm=ci.Clone(histname+"_nonorm_nosmear")
        cinonorm.Scale(nloqcdnorm[j]/dataevents[j])
        cinonorm.Write(histname+"_nonorm_nosmear")
        ci.Write(histname)

        # ALT (=NLO QCD)
        histname=samples[i][0].replace("Anti","")+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
        print(histname)
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
        alt.SetLineColor(colors[col])
        alt.SetTitle("")
        alt.GetYaxis().SetTitle("dN/d#chi")
        alt.Clone(histname+"_nosmear").Write(histname+"_nosmear")
        altnonorm=alt.Clone(histname+"_nonorm_nosmear")
        altnonorm.Scale(nloqcdnorm[j]/dataevents[j])
        altnonorm.Write(histname+"_nonorm_nosmear")
        alt.Write(histname)
        col+=1
        
        if injectSignal:
         data.Add(ci,1)
         data.Add(alt,-1)
         histname='data_obs#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
         data.Write(histname)

        # model, JERtail, sim uncertainty
        uncertainties=["model","JERtail","sim"]
        slopes={}
        slopes[1200]=[0.002,0.01*0.5,0.010] # FIX compute
        slopes[1500]=[0.002,0.01*0.5,0.010] # FIX compute
        slopes[1900]=[0.002,0.006*0.5,0.005] # Updated 27.June 2023 from v4
        slopes[2400]=[0.014,0.007*0.5,0.005] # model updated 28.4.2024
        slopes[3000]=[0.006,0.010*0.5,0.011]
        slopes[3600]=[0.012,0.014*0.5,0.007]
        slopes[4200]=[0.010,0.017*0.5,0.006]
        slopes[4800]=[0.004,0.005*0.5,0.010]
        slopes[5400]=[0.016,0.045*0.5,0.005]
        slopes[6000]=[0.010,0.072*0.5,0.005]
        slopes[7000]=[0.056,0.105*0.5,0.004]
        modelup={}
        modeldown={}
        cimodelup={}
        cimodeldown={}
        uncertaintyvariants=[[str(mn[0]) for mn in massbins]+[""],[str(mn[0]) for mn in massbins]+[""],[str(mn[0]) for mn in massbins]+[""]]
        for unc in uncertainties:
         print(unc)
         for modeln in uncertaintyvariants[uncertainties.index(unc)]:
               histname=samples[i][0]+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
               clone=ci.Clone(histname)
               clone=clone.Rebin(len(chi_binnings[j])-1,clone.GetName(),chi_binnings[j])
               cimodelup[unc]=clone.Clone(histname+"_"+unc+modeln+"Up")
               cimodeldown[unc]=clone.Clone(histname+"_"+unc+modeln+"Down")
               if modeln=="" or str(massbins[j][0])==modeln:
                 for b in range(clone.GetNbinsX()):
                   cimodelup[unc].SetBinContent(b+1,clone.GetBinContent(b+1)*(1.+(clone.GetBinCenter(b+1)-8.5)/7.5*slopes[massbins[j][0]][uncertainties.index(unc)]))
                   cimodeldown[unc].SetBinContent(b+1,clone.GetBinContent(b+1)*(1.-(clone.GetBinCenter(b+1)-8.5)/7.5*slopes[massbins[j][0]][uncertainties.index(unc)]))
               #cimodelup[unc].Scale(dataevents[j]/cimodelup[unc].Integral())
               #cimodeldown[unc].Scale(dataevents[j]/cimodeldown[unc].Integral())
               cimodelup[unc].SetLineColor(colors[col])
               cimodelup[unc].SetLineStyle(2)
               cimodeldown[unc].SetLineColor(colors[col])
               cimodeldown[unc].SetLineStyle(3)
               out.cd()
               cimodelup[unc].Write()
               cimodeldown[unc].Write()

               histname=samples[i][0]+'_ALT#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
               clone=alt.Clone(histname)
               clone=clone.Rebin(len(chi_binnings[j])-1,clone.GetName(),chi_binnings[j])
               modelup[unc]=clone.Clone(histname+"_"+unc+modeln+"Up")
               modeldown[unc]=clone.Clone(histname+"_"+unc+modeln+"Down")
               if modeln=="" or str(massbins[j][0])==modeln:
                 for b in range(clone.GetNbinsX()):
                   modelup[unc].SetBinContent(b+1,clone.GetBinContent(b+1)*(1.+(clone.GetBinCenter(b+1)-8.5)/7.5*slopes[massbins[j][0]][uncertainties.index(unc)]))
                   modeldown[unc].SetBinContent(b+1,clone.GetBinContent(b+1)*(1.-(clone.GetBinCenter(b+1)-8.5)/7.5*slopes[massbins[j][0]][uncertainties.index(unc)]))
               #modelup[unc].Scale(dataevents[j]/modelup[unc].Integral())
               #modeldown[unc].Scale(dataevents[j]/modeldown[unc].Integral())
               modelup[unc].SetLineColor(colors[col])
               modelup[unc].SetLineStyle(2)
               modeldown[unc].SetLineColor(colors[col])
               modeldown[unc].SetLineStyle(3)
               out.cd()
               modelup[unc].Write()
               modeldown[unc].Write()
               if modeln=="": col+=1

        # jes uncertainty
        print("jes")
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
          cijesup.SetLineColor(colors[col])
          cijesup.SetLineStyle(2)
          cijesdown.SetLineColor(colors[col])
          cijesdown.SetLineStyle(3)
          out.cd()
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
          jesup.SetLineColor(colors[col])
          jesup.SetLineStyle(2)
          jesdown.SetLineColor(colors[col])
          jesdown.SetLineStyle(3)
          out.cd()
          jesup.Write()
          jesdown.Write()
          if jesname=="": col+=1

        # jer uncertainty
        print("jer")
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
          cijerup.SetLineColor(colors[col])
          cijerup.SetLineStyle(2)
          cijerdown.SetLineColor(colors[col])
          cijerdown.SetLineStyle(3)
          out.cd()
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
          jerup.SetLineColor(colors[col])
          jerup.SetLineStyle(2)
          jerdown.SetLineColor(colors[col])
          jerdown.SetLineStyle(3)
          out.cd()
          jerup.Write()
          jerdown.Write()
          if jername=="": col+=1

        # prefire uncertainty
        print("prefire")
        histname=samples[i][0]+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
        clone=ci.Clone(histname)
        clone=clone.Rebin(len(chi_binnings[j])-1,clone.GetName(),chi_binnings[j])
        ciprefireup=clone.Clone(histname+"_prefireUp")
        ciprefiredown=clone.Clone(histname+"_prefireDown")
        for b in range(clone.GetNbinsX()):
            ciprefireup.SetBinContent(b+1,clone.GetBinContent(b+1)*((prefirehists[str(j)][0].Eval(clone.GetXaxis().GetBinCenter(b+1))-1.)*0.5+1.))
            ciprefiredown.SetBinContent(b+1,clone.GetBinContent(b+1)*((prefirehists[str(j)][1].Eval(clone.GetXaxis().GetBinCenter(b+1))-1.)*0.5+1.))
        ciprefireup.SetLineColor(colors[col])
        ciprefireup.SetLineStyle(2)
        ciprefiredown.SetLineColor(colors[col])
        ciprefiredown.SetLineStyle(3)
        out.cd()
        ciprefireup.Write()
        ciprefiredown.Write()

        histname=samples[i][0]+'_ALT#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
        clone=alt.Clone(histname)
        clone=clone.Rebin(len(chi_binnings[j])-1,clone.GetName(),chi_binnings[j])
        prefireup=clone.Clone(histname+"_prefireUp")
        prefiredown=clone.Clone(histname+"_prefireDown")
        for b in range(clone.GetNbinsX()):
            prefireup.SetBinContent(b+1,clone.GetBinContent(b+1)*((prefirehists[str(j)][0].Eval(clone.GetXaxis().GetBinCenter(b+1))-1.)*0.5+1.))
            prefiredown.SetBinContent(b+1,clone.GetBinContent(b+1)*((prefirehists[str(j)][1].Eval(clone.GetXaxis().GetBinCenter(b+1))-1.)*0.5+1.))
        prefireup.SetLineColor(colors[col])
        prefireup.SetLineStyle(2)
        prefiredown.SetLineColor(colors[col])
        prefiredown.SetLineStyle(3)
        out.cd()
        prefireup.Write()
        prefiredown.Write()
        col+=1

        # trigger uncertainty
        print("trigger")
        histname=samples[i][0]+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
        clone=ci.Clone(histname)
        clone=clone.Rebin(len(chi_binnings[j])-1,clone.GetName(),chi_binnings[j])
        citriggerup=clone.Clone(histname+"_triggerUp")
        citriggerdown=clone.Clone(histname+"_triggerDown")
        for b in range(clone.GetNbinsX()):
            efficiency=trigger_efficiency[0][j][chi_bins[j][b]]*36.33/137.64+trigger_efficiency[1][j][chi_bins[j][b]]*41.48/137.64+trigger_efficiency[2][j][chi_bins[j][b]]*59.83/137.64
            #print j,b,efficiency
            citriggerup.SetBinContent(b+1,clone.GetBinContent(b+1)*(1./efficiency))
            citriggerdown.SetBinContent(b+1,clone.GetBinContent(b+1)*(efficiency))
            if not useUnfoldedData:
              data.SetBinContent(b+1,data.GetBinContent(b+1)*(1./efficiency)) # Correct for trigger efficiency in data
        citriggerup.Scale(clone.Integral()/citriggerup.Integral())
        citriggerdown.Scale(clone.Integral()/citriggerdown.Integral())
        citriggerup.SetLineColor(colors[col])
        citriggerup.SetLineStyle(2)
        citriggerdown.SetLineColor(colors[col])
        citriggerdown.SetLineStyle(3)
        out.cd()
        citriggerup.Write()
        citriggerdown.Write()

        histname=samples[i][0]+'_ALT#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
        clone=alt.Clone(histname)
        clone=clone.Rebin(len(chi_binnings[j])-1,clone.GetName(),chi_binnings[j])
        triggerup=clone.Clone(histname+"_triggerUp")
        triggerdown=clone.Clone(histname+"_triggerDown")
        for b in range(clone.GetNbinsX()):
            efficiency=trigger_efficiency[0][j][chi_bins[j][b]]*36.33/137.64+trigger_efficiency[1][j][chi_bins[j][b]]*41.48/137.64+trigger_efficiency[2][j][chi_bins[j][b]]*59.83/137.64
            triggerup.SetBinContent(b+1,clone.GetBinContent(b+1)*(1./efficiency))
            triggerdown.SetBinContent(b+1,clone.GetBinContent(b+1)*(efficiency))
        triggerup.Scale(clone.Integral()/triggerup.Integral())
        triggerdown.Scale(clone.Integral()/triggerdown.Integral())
        triggerup.SetLineColor(colors[col])
        triggerup.SetLineStyle(2)
        triggerdown.SetLineColor(colors[col])
        triggerdown.SetLineStyle(3)
        out.cd()
        triggerup.Write()
        triggerdown.Write()
        col+=1

        # NLO PDFup/down
        nloPDFupqcd=None
        for k in mass_bins_nlo_list[j]:
         histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"PDFUp"
         print(histname)
         hnloPDFup = TH1F(nlofilesys.Get(histname))
         #hnloPDFup.Scale(float(mass_bins_nlo3[k+1]-mass_bins_nlo3[k]))
         hnloPDFup=rebin(hnloPDFup,len(chi_binnings[j])-1,chi_binnings[j])
         if nloPDFupqcd:
            nloPDFupqcd.Add(hnloPDFup)
         else:
            nloPDFupqcd=hnloPDFup
        nloPDFupqcd.Add(nloqcdnormd,-1)
        nloPDFupqcd.Scale(1./nloqcdnormd.Integral())
        nloPDFupqcd=smooth(nloPDFupqcd,"pol3") # SMOOTH NNLO PREDICTION (FIX ME)

        nloPDFdownqcd=None
        for k in mass_bins_nlo_list[j]:
         histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"PDFDown"
         print(histname)
         hnloPDFdown = TH1F(nlofilesys.Get(histname))
         #hnloPDFdown.Scale(float(mass_bins_nlo3[k+1]-mass_bins_nlo3[k]))
         hnloPDFdown=rebin(hnloPDFdown,len(chi_binnings[j])-1,chi_binnings[j])
         if nloPDFdownqcd:
            nloPDFdownqcd.Add(hnloPDFdown)
         else:
            nloPDFdownqcd=hnloPDFdown
        nloPDFdownqcd.Add(nloqcdnormd,-1)
        nloPDFdownqcd.Scale(1./nloqcdnormd.Integral())
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
        pdfup.SetLineColor(colors[col])
        pdfup.SetLineStyle(2)
        pdfdown.SetLineColor(colors[col])
        pdfdown.SetLineStyle(3)
        out.cd()
        pdfup.Write()
        pdfdown.Write()

        if "lo" in samples[i][0] or "cteq66" in samples[i][0] or "cteq6ll" in samples[i][0]:
           nloPDFupci=None
           for k in mass_bins_nlo_list[j]:
            histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"PDFUp"
            print(histname)
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
            print(histname)
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
           dmpdfplot=dmpdfcanvas.GetListOfPrimitives()[j+useUnfoldedData*1]
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
        cipdfup.SetLineColor(colors[col])
        cipdfup.SetLineStyle(2)
        cipdfdown.SetLineColor(colors[col])
        cipdfdown.SetLineStyle(3)
        out.cd()
        cipdfup.Write()
        cipdfdown.Write()
        col+=1

        # NLO Scaleup/down
        for scaleVariation in ["MuR","MuF","Alt",""]:
          nloScaleupqcd=None
          for k in mass_bins_nlo_list[j]:
           histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"Scale"+scaleVariation+"Up"
           print(histname)
           if scaleVariation=="":
             if useNNLO:
               hnloScaleup = TH1F(nlofilesys.Get(histname))
             else:
               hnloScaleup = TH1F(nlofilesys.Get(histname.replace("Scale","scale")).Clone(histname))
           elif scaleVariation=="Alt":
             hnloScaleup = TH1F(nlofilealtscale.Get("CIJET_fnl5662j_cs_001_"+pdfset+"_"+muAltScale+"_0_"+str(PDFmembers)+"_70000_V-A-_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-1.0_addmu")).Clone(histname)
           elif scaleVariation=="MuR":
             hnloScaleup = TH1F(nlofilesys.Get("CIJET_fnl5662j_cs_001_"+pdfset+"_"+muScale+"_0_"+str(PDFmembers)+"_70000_V-A-_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-2.0-1.0_addmu")).Clone(histname)
           elif scaleVariation=="MuF":
             hnloScaleup = TH1F(nlofilesys.Get("CIJET_fnl5662j_cs_001_"+pdfset+"_"+muScale+"_0_"+str(PDFmembers)+"_70000_V-A-_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-2.0_addmu")).Clone(histname)
           #hnloScaleup.Scale(float(mass_bins_nlo3[k+1]-mass_bins_nlo3[k]))
           hnloScaleup=rebin(hnloScaleup,len(chi_binnings[j])-1,chi_binnings[j])
           if nloScaleupqcd:
              nloScaleupqcd.Add(hnloScaleup)
           else:
              nloScaleupqcd=hnloScaleup
          nloScaleupqcd.Add(nloqcdnormd,-1)
          nloScaleupqcd.Scale(1./nloqcdnormd.Integral())
          nloScaleupqcd=smooth(nloScaleupqcd,"pol3") # SMOOTH NNLO PREDICTION (FIX ME)

          nloScaledownqcd=None
          for k in mass_bins_nlo_list[j]:
           histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"Scale"+scaleVariation+"Down"
           print(histname)
           if scaleVariation=="":
             if useNNLO:
               hnloScaledown = TH1F(nlofilesys.Get(histname))
             else:
               hnloScaledown = TH1F(nlofilesys.Get(histname.replace("Scale","scale")).Clone(histname))
           elif scaleVariation=="Alt":
             hnloScaledown = TH1F(nlofilesys.Get("CIJET_fnl5662j_cs_001_"+pdfset+"_"+muScale+"_0_"+str(PDFmembers)+"_70000_V-A-_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-1.0_addmu")).Clone(histname)
           elif scaleVariation=="MuR":
             hnloScaledown = TH1F(nlofilesys.Get("CIJET_fnl5662j_cs_001_"+pdfset+"_"+muScale+"_0_"+str(PDFmembers)+"_70000_V-A-_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-0.5-1.0_addmu")).Clone(histname)
           elif scaleVariation=="MuF":
             hnloScaledown = TH1F(nlofilesys.Get("CIJET_fnl5662j_cs_001_"+pdfset+"_"+muScale+"_0_"+str(PDFmembers)+"_70000_V-A-_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-0.5_addmu")).Clone(histname)
           #hnloScaledown.Scale(float(mass_bins_nlo3[k+1]-mass_bins_nlo3[k]))
           hnloScaledown=rebin(hnloScaledown,len(chi_binnings[j])-1,chi_binnings[j])
           if nloScaledownqcd:
              nloScaledownqcd.Add(hnloScaledown)
           else:
              nloScaledownqcd=hnloScaledown
          nloScaledownqcd.Add(nloqcdnormd,-1)
          nloScaledownqcd.Scale(1./nloqcdnormd.Integral())
          nloScaledownqcd=smooth(nloScaledownqcd,"pol3") # SMOOTH NNLO PREDICTION (FIX ME)

          scaleup=alt.Clone(alt.GetName()+"_scale"+scaleVariation+"Up")
          scaledown=alt.Clone(alt.GetName()+"_scale"+scaleVariation+"Down")
          scaleup.Add(nloScaleupqcd,dataevents[j])
          scaledown.Add(nloScaledownqcd,dataevents[j])
          for b in range(scaleup.GetXaxis().GetNbins()):
              scaleup.SetBinError(b+1,0)
              scaledown.SetBinError(b+1,0)
              if scaleVariation=="" and (b+1)==scaleup.FindBin(scaleup.GetMean()):
                  tmp=scaleup.GetBinContent(b+1)
                  scaleup.SetBinContent(b+1,(scaleup.GetBinContent(b+1)+scaledown.GetBinContent(b+1))/2.)
                  scaledown.SetBinContent(b+1,(tmp+scaledown.GetBinContent(b+1))/2.)
              if scaleVariation=="" and (b+1)>scaleup.FindBin(scaleup.GetMean()):
                  tmp=scaleup.GetBinContent(b+1)
                  scaleup.SetBinContent(b+1,scaledown.GetBinContent(b+1))
                  scaledown.SetBinContent(b+1,tmp)
          #scaleup.Scale(dataevents[j]/scaleup.Integral())
          #scaledown.Scale(dataevents[j]/scaledown.Integral())
          scaleup.SetLineColor(colors[col])
          scaleup.SetLineStyle(2)
          scaledown.SetLineColor(colors[col])
          scaledown.SetLineStyle(3)
          if scaleVariation=="":
           for mn in massbins:
            if massbins.index(mn)==j:
              scaleupmn=scaleup.Clone(scaleup.GetName().replace("scale","scale"+str(mn[0])))
              scaledownmn=scaledown.Clone(scaledown.GetName().replace("scale","scale"+str(mn[0])))
            else:
              scaleupmn=alt.Clone(scaleup.GetName().replace("scale","scale"+str(mn[0])))
              scaledownmn=alt.Clone(scaledown.GetName().replace("scale","scale"+str(mn[0])))
            scaleupmn.Write()
            scaledownmn.Write()
          out.cd()
          scaleup.Write()
          scaledown.Write()
          if scaleVariation=="Alt":
            scalealt=scaleup
          
          if "lo" in samples[i][0] or "cteq66" in samples[i][0] or "cteq6ll" in samples[i][0]:
             signalmassname="_".join(samples[i][0].split("_")[2:4])
             if massbins[j][0]<3600: ### No CI calculation available
               signalmassname=signalmassname.replace("35000","30000").replace("40000","30000").replace("45000","30000").replace("50000","30000").replace("55000","30000").replace("60000","30000").replace("65000","30000").replace("70000","30000")
             nloScaleupci=None
             for k in mass_bins_nlo_list[j]:
              histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"Scale"+scaleVariation+"Up"
              print(histname)
              if scaleVariation=="":
                if useNNLO:
                  hnloScaleup = TH1F(cinlofile.Get(histname))
                else:
                  hnloScaleup = TH1F(cinlofile.Get(histname.replace("Scale","scale")).Clone(histname))
              elif scaleVariation=="Alt":
                hnloScaleup = TH1F(cinlofilealtscale.Get("CIJET_fnl5662j_cs_001_"+("ct14nlo_0_56" if massbins[j][0]<3600 else pdfset+"_"+muAltScale+"_0_"+str(PDFmembers))+"_"+signalmassname+"_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-1.0_addmu")).Clone(histname)
              elif scaleVariation=="MuR":
                hnloScaleup = TH1F(cinlofile.Get("CIJET_fnl5662j_cs_001_"+("ct14nlo_0_56" if massbins[j][0]<3600 else pdfset+"_"+muScale+"_0_"+str(PDFmembers))+"_"+signalmassname+"_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-2.0-1.0_addmu")).Clone(histname)
              elif scaleVariation=="MuF":
                hnloScaleup = TH1F(cinlofile.Get("CIJET_fnl5662j_cs_001_"+("ct14nlo_0_56" if massbins[j][0]<3600 else pdfset+"_"+muScale+"_0_"+str(PDFmembers))+"_"+signalmassname+"_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-2.0_addmu")).Clone(histname)
              hnloScaleup=hnloScaleup.Rebin(len(chi_binnings[j])-1,histname.replace("_backup",""),chi_binnings[j])
              if nloScaleupci:
                  nloScaleupci.Add(hnloScaleup)
              else:
                  nloScaleupci=hnloScaleup
             nloScaleupci.Add(cibackup,-1)
             nloScaleupci.Scale(1./cibackup.Integral())

             nloScaledownci=None
             for k in mass_bins_nlo_list[j]:
              histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"Scale"+scaleVariation+"Down"
              print(histname)
              if scaleVariation=="":
                if useNNLO:
                  hnloScaledown = TH1F(cinlofile.Get(histname))
                else:
                  hnloScaledown = TH1F(cinlofile.Get(histname.replace("Scale","scale")).Clone(histname))
              elif scaleVariation=="Alt":
                hnloScaledown = TH1F(cinlofile.Get("CIJET_fnl5662j_cs_001_"+("ct14nlo_0_56" if massbins[j][0]<3600 else pdfset+"_"+muScale+"_0_"+str(PDFmembers))+"_"+signalmassname+"_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-1.0_addmu")).Clone(histname)
              elif scaleVariation=="MuR":
                hnloScaledown = TH1F(cinlofile.Get("CIJET_fnl5662j_cs_001_"+("ct14nlo_0_56" if massbins[j][0]<3600 else pdfset+"_"+muScale+"_0_"+str(PDFmembers))+"_"+signalmassname+"_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-0.5-1.0_addmu")).Clone(histname)
              elif scaleVariation=="MuF":
                hnloScaledown = TH1F(cinlofile.Get("CIJET_fnl5662j_cs_001_"+("ct14nlo_0_56" if massbins[j][0]<3600 else pdfset+"_"+muScale+"_0_"+str(PDFmembers))+"_"+signalmassname+"_mu_qcd_chi-"+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-0.5_addmu")).Clone(histname)
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
              if scaleVariation=="" and (b+1)==ciscaleup.FindBin(ciscaleup.GetMean()):
                 tmp=ciscaleup.GetBinContent(b+1)
                 ciscaleup.SetBinContent(b+1,(ciscaleup.GetBinContent(b+1)+ciscaledown.GetBinContent(b+1))/2.)
                 ciscaledown.SetBinContent(b+1,(tmp+ciscaledown.GetBinContent(b+1))/2.)
              if scaleVariation=="" and (b+1)>ciscaleup.FindBin(ciscaleup.GetMean()):
                 tmp=ciscaleup.GetBinContent(b+1)
                 ciscaleup.SetBinContent(b+1,ciscaledown.GetBinContent(b+1))
                 ciscaledown.SetBinContent(b+1,tmp)
          #ciscaleup.Scale(dataevents[j]/ciscaleup.Integral())
          #ciscaledown.Scale(dataevents[j]/ciscaledown.Integral())
          ciscaleup.SetLineColor(colors[col])
          ciscaleup.SetLineStyle(2)
          ciscaledown.SetLineColor(colors[col])
          ciscaledown.SetLineStyle(3)
          out.cd()
          ciscaleup.Write()
          ciscaledown.Write()
          if scaleVariation=="":
           col+=1
           for mn in massbins:
            if massbins.index(mn)==j:
              ciscaleupmn=scaleup.Clone(ciscaleup.GetName().replace("scale","scale"+str(mn[0])))
              ciscaledownmn=scaledown.Clone(ciscaledown.GetName().replace("scale","scale"+str(mn[0])))
            else:
              ciscaleupmn=ci.Clone(ciscaleup.GetName().replace("scale","scale"+str(mn[0])))
              ciscaledownmn=ci.Clone(ciscaledown.GetName().replace("scale","scale"+str(mn[0])))
            ciscaleupmn.Write()
            ciscaledownmn.Write()
          if scaleVariation=="Alt":
            col+=1
            ciscalealt=ciscaleup
        
        # theory stat uncertainties
        theorystatup={}
        theorystatdown={}
        citheorystatup={}
        citheorystatdown={}
        nloStatupqcd=None
        for k in mass_bins_nlo_list[j]:
         histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"StatUp"
         print(histname)
         if useNNLO:
           hnloStatup = TH1F(nlofilesys.Get(histname))
         else:
           # take statistical uncertainty from NNLO (since not available for NLO)
           hnloStatup = TH1F(nlofilealtscale.Get(histname).Clone(histname))
           #hnloStatup.Scale(1.+slopes[massbins[j][0]][0]) # take statistical uncertainty from Pythia/Herwig response matrix
         hnloStatup=rebin(hnloStatup,len(chi_binnings[j])-1,chi_binnings[j])
         if nloStatupqcd:
            nloStatupqcd.Add(hnloStatup)
         else:
            nloStatupqcd=hnloStatup
         if useNNLO:
           hnloref = TH1F(nlofilesys.Get(histname.replace("StatUp","")))
         else:
           hnloref = TH1F(nlofilealtscale.Get(histname.replace("StatUp","")))
         hnloref=rebin(hnloref,len(chi_binnings[j])-1,chi_binnings[j])
         nloStatupqcd.Add(hnloref,-1)
        nloStatupqcd.Scale(dataevents[j]/nloqcdnormd.Integral())
        cihistname=samples[i][0]+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
        ciclone=ci.Clone(cihistname)
        ciclone=ciclone.Rebin(len(chi_binnings[j])-1,ciclone.GetName(),chi_binnings[j])
        citheorystatup[""]=ciclone.Clone(cihistname+"_"+"statUp")
        citheorystatdown[""]=ciclone.Clone(cihistname+"_"+"statDown")
        citheorystatup[""].SetLineColor(colors[col])
        citheorystatup[""].SetLineStyle(2)
        citheorystatdown[""].SetLineColor(colors[col])
        citheorystatdown[""].SetLineStyle(3)
        histname=samples[i][0]+'_ALT#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
        clone=alt.Clone(histname)
        clone=clone.Rebin(len(chi_binnings[j])-1,clone.GetName(),chi_binnings[j])
        theorystatup[""]=clone.Clone(histname+"_"+"statUp")
        theorystatdown[""]=clone.Clone(histname+"_"+"statDown")
        theorystatup[""].SetLineColor(colors[col])
        theorystatup[""].SetLineStyle(2)
        theorystatdown[""].SetLineColor(colors[col])
        theorystatdown[""].SetLineStyle(3)
        for massbin in range(len(massbins)):
         for chibin in range(len(chi_binnings[massbin])-1):
          unc="stat"+str(massbin)+"_"+str(chibin)
          citheorystatup[unc]=ciclone.Clone(cihistname+"_"+unc+"Up")
          citheorystatdown[unc]=ciclone.Clone(cihistname+"_"+unc+"Down")
          if j==massbin:
            theorystatvalue=nloStatupqcd.GetBinContent(chibin+1)
            citheorystatup[""].SetBinContent(chibin+1,ciclone.GetBinContent(chibin+1)+theorystatvalue)
            citheorystatdown[""].SetBinContent(chibin+1,ciclone.GetBinContent(chibin+1)-theorystatvalue)
            citheorystatup[unc].SetBinContent(chibin+1,ciclone.GetBinContent(chibin+1)+theorystatvalue)
            citheorystatdown[unc].SetBinContent(chibin+1,ciclone.GetBinContent(chibin+1)-theorystatvalue)
          out.cd()
          citheorystatup[unc].Write()
          citheorystatdown[unc].Write()

          theorystatup[unc]=clone.Clone(histname+"_"+unc+"Up")
          theorystatdown[unc]=clone.Clone(histname+"_"+unc+"Down")
          if j==massbin:
             theorystatup[""].SetBinContent(chibin+1,clone.GetBinContent(chibin+1)+theorystatvalue)
             theorystatdown[""].SetBinContent(chibin+1,clone.GetBinContent(chibin+1)-theorystatvalue)
             theorystatup[unc].SetBinContent(chibin+1,clone.GetBinContent(chibin+1)+theorystatvalue)
             theorystatdown[unc].SetBinContent(chibin+1,clone.GetBinContent(chibin+1)-theorystatvalue)
          out.cd()
          theorystatup[unc].Write()
          theorystatdown[unc].Write()
        out.cd()
        citheorystatup[""].Write()
        citheorystatdown[""].Write()
        theorystatup[""].Write()
        theorystatdown[""].Write()
        col+=1

        # DATA BLINDED
        #data=alt.Clone("data_blinded")
        #for b in range(data.GetXaxis().GetNbins()):
        #    data.SetBinError(b+1,sqrt(data.GetBinContent(b+1)))
        #out.cd()
        #data.Write('data_obs#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1")
      
        # FAKE SIGNAL
        #ci=alt.Clone("fake_signal")
        #out.cd()
        #ci.Write(samples[i][0]+'#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1")
      
        # PLOTS
        #if j<3:
        #   continue
        canvas.cd(j+2)#j-2
        legend1=TLegend(0.2,0.6,0.9,0.95,(str(massbins[j][0])+"<m_{jj}<"+str(massbins[j][1])+" GeV").replace("7000<m_{jj}<13000","m_{jj}>7000"))
        legends+=[legend1]
        legend1.AddEntry(data,"data","lpe")
        alt=cloneNormalize(alt)
        plots+=[alt]
        alt.Draw("he")
        legend1.AddEntry(alt,"QCD","l")
        
        # background uncertainties
        for unc in uncertainties:
          modelup[unc]=cloneNormalize(modelup[unc])
          plots+=[modelup[unc]]
          #modelup[unc].Draw("hesame")
          legend1.AddEntry(modelup[unc],unc,"l")
          modeldown[unc]=cloneNormalize(modeldown[unc])
          plots+=[modeldown[unc]]
          #modeldown[unc].Draw("hesame")
        jesup=cloneNormalize(jesup)
        plots+=[jesup]
        #jesup.Draw("hesame")
        legend1.AddEntry(jesup,"JES","l")
        jesdown=cloneNormalize(jesdown)
        plots+=[jesdown]
        #jesdown.Draw("hesame")
        jerup=cloneNormalize(jerup)
        plots+=[jerup]
        #jerup.Draw("hesame")
        legend1.AddEntry(jerup,"JER","l")
        jerdown=cloneNormalize(jerdown)
        plots+=[jerdown]
        #jerdown.Draw("hesame")
        prefireup=cloneNormalize(prefireup)
        plots+=[prefireup]
        #prefireup.Draw("hesame")
        legend1.AddEntry(prefireup,"prefire","l")
        prefiredown=cloneNormalize(prefiredown)
        plots+=[prefiredown]
        #prefiredown.Draw("hesame")
        triggerup=cloneNormalize(triggerup)
        plots+=[triggerup]
        #triggerup.Draw("hesame")
        legend1.AddEntry(triggerup,"trigger","l")
        triggerdown=cloneNormalize(triggerdown)
        plots+=[triggerdown]
        #triggerdown.Draw("hesame")
        pdfup=cloneNormalize(pdfup)
        plots+=[pdfup]
        #pdfup.Draw("hesame")
        legend1.AddEntry(pdfup,"PDF","l")
        pdfdown=cloneNormalize(pdfdown)
        plots+=[pdfdown]
        #pdfdown.Draw("hesame")
        scalealt=cloneNormalize(scalealt)
        plots+=[scalealt]
        #scalealt.Draw("hesame")
        legend1.AddEntry(scalealt,"scale alt","l")
        scaleup=cloneNormalize(scaleup)
        plots+=[scaleup]
        #scaleup.Draw("hesame")
        legend1.AddEntry(scaleup,"scale","l")
        scaledown=cloneNormalize(scaledown)
        plots+=[scaledown]
        #scaledown.Draw("hesame")
        theorystatup[""]=cloneNormalize(theorystatup[""])
        plots+=[theorystatup[""]]
        #theorystatup[""].Draw("hesame")
        legend1.AddEntry(theorystatup[""],"stat","l")
        theorystatdown[""]=cloneNormalize(theorystatdown[""])
        plots+=[theorystatdown[""]]
        #theorystatdown[""].Draw("hesame")
        
        if samples[i][0]!="QCD":
         ci=cloneNormalize(ci)
         plots+=[ci]
         ci.Draw("hesame")
         legend1.AddEntry(ci,"New Physics","l")
        
         # signal uncertainties
         for unc in uncertainties:
           cimodelup[unc]=cloneNormalize(cimodelup[unc])
           plots+=[cimodelup[unc]]
           cimodelup[unc].Draw("hesame")
           cimodeldown[unc]=cloneNormalize(cimodeldown[unc])
           plots+=[cimodeldown[unc]]
           cimodeldown[unc].Draw("hesame")
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
         ciprefireup=cloneNormalize(ciprefireup)
         plots+=[ciprefireup]
         ciprefireup.Draw("hesame")
         ciprefiredown=cloneNormalize(ciprefiredown)
         plots+=[ciprefiredown]
         ciprefiredown.Draw("hesame")
         citriggerup=cloneNormalize(citriggerup)
         plots+=[citriggerup]
         citriggerup.Draw("hesame")
         citriggerdown=cloneNormalize(citriggerdown)
         plots+=[citriggerdown]
         citriggerdown.Draw("hesame")
         cipdfup=cloneNormalize(cipdfup)
         plots+=[cipdfup]
         cipdfup.Draw("hesame")
         cipdfdown=cloneNormalize(cipdfdown)
         plots+=[cipdfdown]
         cipdfdown.Draw("hesame")
         ciscalealt=cloneNormalize(ciscalealt)
         plots+=[ciscalealt]
         ciscalealt.Draw("hesame")
         ciscaleup=cloneNormalize(ciscaleup)
         plots+=[ciscaleup]
         ciscaleup.Draw("hesame")
         ciscaledown=cloneNormalize(ciscaledown)
         plots+=[ciscaledown]
         ciscaledown.Draw("hesame")
         citheorystatup[""]=cloneNormalize(citheorystatup[""])
         plots+=[citheorystatup[""]]
         citheorystatup[""].Draw("hesame")
         citheorystatdown[""]=cloneNormalize(citheorystatdown[""])
         plots+=[citheorystatdown[""]]
         citheorystatdown[""].Draw("hesame")
 
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

        #miny=0#min(cloneNormalize(datahist[j]).GetMinimum(),alt.GetMinimum())
        #maxy=max(cloneNormalize(datahist[j]).GetMaximum(),alt.GetMaximum())
        #alt.GetYaxis().SetRangeUser(miny*0.8,(maxy-miny)*2.+miny)
        miny=alt.GetMinimum()
        maxy=alt.GetMaximum()
        alt.GetYaxis().SetRangeUser(miny*0.8,(maxy-miny)*2.+miny)
        
        legend1.SetTextSize(0.04)
        legend1.SetFillStyle(0)
        legend1.Draw("same")

      canvas.SaveAs(prefix + "_"+samples[i][0].replace("QCD","") + '_sys_run'+run+'.pdf')
      #canvas.SaveAs(prefix + "_"+samples[i][0].replace("QCD","") + '_sys_run2.eps')

      out.cd()
      datahist2d_chiBins2.Write()
      datahist2d_chiBins3.Write()

      #out.Close()
      #out=TFile(sample,'UPDATE')
      
      if not useUnfoldedData:
        canvas = TCanvas("","",0,0,800,600)
        canvas.Divide(4,3)
        print("Applying response matrix")
        #matrix=TFile("datacards/responseMatrices_20190930.root",'READ')
        #matrix1=matrix.Get("TMatrixT<double>;2") #low mass chi binning
        #matrix2=matrix.Get("TMatrixT<double>;3") #highest mass chi binning
        #matrixMassBins=[1000,1200,1500,1900,2400,3000,3600,4200,4800,5400,6000,7000,13000]
        #print (len(matrixMassBins)-1)*(len(chi_bins[0])-1)
        #matrix1.Print()
        #print (len(matrixMassBins)-1)*(len(chi_bins[-1])-1)
        #matrix2.Print()
        if "herwigp" in responsePostfix:
          responses=TFile("datacards/Response_"+responsePostfix+"_Pt_170toInf_CB2_AK4SF_March24.root",'READ')
        elif "madgraphMLM" in responsePostfix:
          responses=TFile("datacards/Response_"+responsePostfix+"_HT_200toInf_CB2_AK4SF_March24.root",'READ')
        elif "T" in responsePostfix:
          responses=TFile("datacards/Response_pythia8_Pt_170toInf_CB2_AK4SF_March24_"+responsePostfix+".root",'READ')
        elif "Sys" in responsePostfix:
          responses=TFile("datacards/Response_pythia8_Pt_170toInf_GS_AK4SF_"+responsePostfix+"_March24.root",'READ')
        elif "GS" in responsePostfix:
          responses=TFile("datacards/Response_pythia8_Pt_170toInf_GS_AK4SF_March24.root",'READ')
        elif responsePostfix=="" or responsePostfix=="RECO":
          responses=TFile("datacards/Response_pythia8_Pt_170toInf_CB2_AK4SF_March24.root",'READ')

        # modify histograms
        histogram_names=[
          (samples[i][0]+'_ALT#chi',"_rebin1"),
          (samples[i][0]+'#chi',"_rebin1"),
        ]
        for pre,post in histogram_names:
          althists=[]
          althistsclones=[]
          histsunfold={}
          for j in range(len(massbins)):
            histname=pre+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+post
            #print histname
            althists+=[out.Get(histname)]
            althistsclones+=[althists[-1].Clone(althists[-1].GetName()+"original")]
            althistsclones[-1].Scale(nloqcdnorm[j]/dataevents[j]) # Use (N)NLO norm rather than data (since event numbers in mjj<2.4 TeV are not accurate)
            althists[-1].Scale(0)

          # integral response matrix files from Jingyu
          response_integrals={}
          for j1 in range(len(massbins)):
            for b1 in range(althists[j1].GetNbinsX()):
              response_integrals[str(j1)+"_"+str(b1)]=0
              for j2 in range(len(massbins)):
                for b2 in range(althists[j2].GetNbinsX()):
                  response=0
                  response_name="QCD_ALT#chi"+str(massbins[j2]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_bin_"+str(j1)+"_"
                  if j2==(len(massbins)-1):
                    # highest reco mass has coarse gen-chi-binning in Jingyu's files
                    bin1=althists[j2].GetXaxis().FindBin(althists[j1].GetBinCenter(b1+1))-1
                    response=responses.Get(response_name+str(bin1)).GetBinContent(b2+1)/althists[j2].GetXaxis().GetBinWidth(bin1+1)*althists[j1].GetXaxis().GetBinWidth(b1+1)
                  elif j1==(len(massbins)-1):
                    # highest gen mass bin has fine gen-chi-binning in Jingyu's files
                    for bin1 in range(althists[0].GetNbinsX()):
                      if althists[0].GetXaxis().GetBinUpEdge(bin1+1)>althists[j1].GetXaxis().GetBinLowEdge(b1+1) and\
                         althists[0].GetXaxis().GetBinLowEdge(bin1+1)<althists[j1].GetXaxis().GetBinUpEdge(b1+1):
                        response+=responses.Get(response_name+str(bin1)).GetBinContent(b2+1)
                  else:
                    response=responses.Get(response_name+str(b1)).GetBinContent(b2+1)
                  if abs(j1-j2)>4: response=0 # REMOVE FAR OFF-DIAGONAL
                  response_integrals[str(j1)+"_"+str(b1)]+=response

          if samples[i][0]=="QCD": # for unfolding
            for j in range(len(massbins)):
              for j1 in range(len(massbins)-3):
                for b1 in range(althists[j1].GetNbinsX()):
                  unfoldbin="_bin_"+str(j1)+"_"+str(b1)+"_"
                  histsunfold[str(j)+unfoldbin]=althists[j].Clone(althists[j].GetName()+unfoldbin)
          
          countgen=0
          for j1 in range(len(massbins)):
            for b1 in range(althists[j1].GetNbinsX()):
              for j2 in range(len(massbins)):
                for b2 in range(althists[j2].GetNbinsX()):
                  # map massbins on matrixMassBins
                  #j1m=j1+1
                  #j2m=j2+1
                  #response=0
                  #if j1==(len(massbins)-1) and j2==(len(massbins)-1):
                  #  # both in highest mass bin
                  #  response=matrix2[j2m+b2*(len(matrixMassBins)-1)][j1m+b1*(len(matrixMassBins)-1)]
                  #elif j1!=(len(massbins)-1) and j2!=(len(massbins)-1):
                  #  # both in lower mass bins
                  #  response=matrix1[j2m+b2*(len(matrixMassBins)-1)][j1m+b1*(len(matrixMassBins)-1)]
                  #elif j1==(len(massbins)-1):
                  #  # one in highest, one in lower mass bin
                  #  for bin1 in range(althists[0].GetNbinsX()):
                  #    if althists[0].GetXaxis().GetBinUpEdge(bin1+1)>althists[j1].GetXaxis().GetBinLowEdge(b1+1) and\
                  #       althists[0].GetXaxis().GetBinLowEdge(bin1+1)<althists[j1].GetXaxis().GetBinUpEdge(b1+1):
                  #      response+=matrix1[j2m+b2*(len(matrixMassBins)-1)][j1m+bin1*(len(matrixMassBins)-1)]/\
                  #          (althists[j1].GetXaxis().GetBinUpEdge(b1+1)-althists[j1].GetXaxis().GetBinLowEdge(b1+1))*\
                  #          (althists[0].GetXaxis().GetBinUpEdge(bin1+1)-althists[0].GetXaxis().GetBinLowEdge(bin1+1))
                  #elif j2==(len(massbins)-1):
                  #  # one in highest, one in lower mass bin
                  #  for bin2 in range(althists[0].GetNbinsX()):
                  #    if althists[0].GetXaxis().GetBinUpEdge(bin2+1)>althists[j2].GetXaxis().GetBinLowEdge(b2+1) and\
                  #       althists[0].GetXaxis().GetBinLowEdge(bin2+1)<althists[j2].GetXaxis().GetBinUpEdge(b2+1):
                  #      response+=matrix1[j2m+bin2*(len(matrixMassBins)-1)][j1m+b1*(len(matrixMassBins)-1)]
                  response=0
                  response_name="QCD_ALT#chi"+str(massbins[j2]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_bin_"+str(j1)+"_"
                  if j2==(len(massbins)-1):
                    # highest reco mass has coarse gen-chi-binning in Jingyu's files
                    bin1=althists[j2].GetXaxis().FindBin(althists[j1].GetBinCenter(b1+1))-1
                    response=responses.Get(response_name+str(bin1)).GetBinContent(b2+1)/althists[j2].GetXaxis().GetBinWidth(bin1+1)*althists[j1].GetXaxis().GetBinWidth(b1+1)
                  elif j1==(len(massbins)-1):
                    # highest gen mass bin has fine gen-chi-binning in Jingyu's files
                    for bin1 in range(althists[0].GetNbinsX()):
                      if althists[0].GetXaxis().GetBinUpEdge(bin1+1)>althists[j1].GetXaxis().GetBinLowEdge(b1+1) and\
                         althists[0].GetXaxis().GetBinLowEdge(bin1+1)<althists[j1].GetXaxis().GetBinUpEdge(b1+1):
                        response+=responses.Get(response_name+str(bin1)).GetBinContent(b2+1)
                  else:
                    response=responses.Get(response_name+str(b1)).GetBinContent(b2+1)
                  if abs(j1-j2)>4: response=0 # REMOVE FAR OFF-DIAGONAL
                  response/=response_integrals[str(j1)+"_"+str(b1)]
                  #print "gen",j1,b1,"reco",j2,b2,"old-response",response,"new-response",response_name,response
                  althists[j2].Fill(althists[j2].GetBinCenter(b2+1),althistsclones[j1].GetBinContent(b1+1)*response)
                  if samples[i][0]=="QCD": # for unfolding
                    unfoldbin="_bin_"+str(max(0,j1-3))+"_"+str(b1)+"_" # 2400 is underflow bin
                    histsunfold[str(j2)+unfoldbin].Fill(althists[j2].GetBinCenter(b2+1),max(0,althistsclones[j1].GetBinContent(b1+1)*response))
          out.cd()
          for j in range(len(massbins)):
            if samples[i][0]=="QCD": # for unfolding
              for j1 in range(len(massbins)-3):
                for b1 in range(althists[j1].GetNbinsX()):
                   unfoldbin="_bin_"+str(j1)+"_"+str(b1)+"_"
                   histsunfold[str(j)+unfoldbin].Scale(dataevents[j]/althists[j].Integral())
                   for b in range(althists[j].GetXaxis().GetNbins()):
                     histsunfold[str(j)+unfoldbin].SetBinError(b+1,0)
            histnonorm=althists[j].Clone(althists[j].GetName()+"_nonorm").Write(althists[j].GetName()+"_nonorm")

            if responsePostfix=="RECO" and pre==samples[i][0]+'#chi':
              print("Using RECO QBH instead of smeared QBH prediction")
              althists[j]=out.Get(althists[j].GetName().replace(samples[i][0]+'#chi',samples[i][0]+'_ALT#chi')+"_nonorm").Clone(althists[j].GetName())
              histnamein='datacard_shapelimit13TeV_run2_UL18_QBH#chi'+str(massbins[j]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
              print(histnamein)
              reco=freco.Get(histnamein)
              reco=reco.Rebin(len(chi_binnings[j])-1,reco.GetName(),chi_binnings[j])
              reco.Scale(samples[i][1][0][1]*nloqcdnorm[j]/nloqcdnormraw[j])
              print(reco.Integral(), althists[j].Integral())
              althists[j].Add(reco)

            #print althists[j].GetName()+"_nonorm"
            print(dataevents[j],althistsclones[j].Integral(),althists[j].Integral())
            althists[j].Scale(dataevents[j]/althists[j].Integral())
            althistsclones[j].Scale(dataevents[j]/althistsclones[j].Integral())
            for b in range(althists[j].GetXaxis().GetNbins()):
              althists[j].SetBinError(b+1,0)
            #print althists[j].GetBinContent(1),althistsclones[j].GetBinContent(1)
              
          for hist in althists:
            hist.Write()
          if samples[i][0]=="QCD": # for unfolding
            for j in range(len(massbins)):
              for j1 in range(len(massbins)-3):
                for b1 in range(althists[j1].GetNbinsX()):
                  unfoldbin="_bin_"+str(j1)+"_"+str(b1)+"_"
                  histsunfold[str(j)+unfoldbin].Write()
          for j in range(len(massbins)):
            canvas.cd(j+2)#j-2
            alt=cloneNormalize(althists[j])
            if "ALT" in pre:
              alt.Draw("he")
              miny=0#min(cloneNormalize(datahist[j]).GetMinimum(),alt.GetMinimum())
              maxy=max(cloneNormalize(datahist[j]).GetMaximum(),alt.GetMaximum())
              alt.GetYaxis().SetRangeUser(miny*0.8,(maxy-miny)*2.+miny)
            else:
              if samples[i][0]!="QCD":
               alt.Draw("hesame")
            plots+=[alt]
          syshists={}
          syss=["model","JERtail","sim","pdf"]
          syss+=["model"+str(mn[0]) for mn in massbins]
          skipInSum=["model"+str(mn[0]) for mn in massbins]
          syss+=["JERtail"+str(mn[0]) for mn in massbins]
          skipInSum=["JERtail"+str(mn[0]) for mn in massbins]
          syss+=["sim"+str(mn[0]) for mn in massbins]
          skipInSum=["sim"+str(mn[0]) for mn in massbins]
          syss+=["jes"]
          for n in range(len(jessources)-1):
            syss+=["jes"+str(n+1)]
            skipInSum+=["jes"+str(n+1)]
          syss+=["jer"]
          for n in range(len(jersources)-1):
            syss+=["jer"+str(n+1)]
            skipInSum+=["jer"+str(n+1)]
          syss+=["prefire"]
          syss+=["trigger"]
          syss+=["scale"]
          syss+=["scaleAlt","scaleMuR","scaleMuF"]
          skipInSum+=["scaleAlt","scaleMuR","scaleMuF"]
          syss+=["scale"+str(mn[0]) for mn in massbins]
          skipInSum+=["scale"+str(mn[0]) for mn in massbins]
          syss+=["stat"]
          for mn in range(len(massbins)):
            for cn in range(len(chi_binnings[mn])-1):
              syss+=["stat"+str(mn)+"_"+str(cn)]
              skipInSum+=["stat"+str(mn)+"_"+str(cn)]
          for j in range(len(massbins)):
            print(massbins[j])
            for sys in syss:
              #print sys
              for shift in ["Up","Down"]:
                histname=althists[j].GetName()+"_"+sys+shift
                #print histname
                sysHist=out.Get(histname)#.Clone()
                #print sysHist
                #print dir(sysHist)
                sysNorm=sysHist.Integral()

                for b in range(sysHist.GetNbinsX()):
                  if althistsclones[j].GetBinContent(b+1)*althists[j].GetBinContent(b+1)>0:
                    sysHist.SetBinContent(b+1,sysHist.GetBinContent(b+1)/althistsclones[j].GetBinContent(b+1)*althists[j].GetBinContent(b+1))
                if sysHist.Integral()>0:
                  sysHist.Scale(sysNorm/althistsclones[j].Integral()*althists[j].Integral()/sysHist.Integral())
                #if j==3:
                #  print "norm", sys, shift, althistsclones[j].Integral(), althists[j].Integral(), sysHist.Integral()
                #if "model2400" in sys and j==3:
                #  print "shape", sys, shift, sysHist.GetBinContent(1), sysHist.GetBinContent(sysHist.GetNbinsX())
                #sysHist.Scale(dataevents[j]/sysHist.Integral())
                sysHist.Write()
                syshists[sys+shift]=sysHist

                if samples[i][0]=="QCD": # for unfolding
                  integral=0
                  for j1 in range(len(massbins)-3):
                    for b1 in range(althists[j1].GetNbinsX()):
                      unfoldbin="_bin_"+str(j1)+"_"+str(b1)+"_"
                      sysHistunfoldbin=histsunfold[str(j)+unfoldbin].Clone(histname+unfoldbin)
                      for b in range(sysHist.GetNbinsX()):
                        if althistsclones[j].GetBinContent(b+1)*histsunfold[str(j)+unfoldbin].GetBinContent(b+1)>0:
                          sysHistunfoldbin.SetBinContent(b+1,sysHist.GetBinContent(b+1)/althists[j].GetBinContent(b+1)*histsunfold[str(j)+unfoldbin].GetBinContent(b+1))
                      #if "model2400" in sys and j==3:
                      #   print "bin shape", sys, shift, sysNorm, j1, b1, histsunfold[str(j)+unfoldbin].Integral(), sysHistunfoldbin.Integral(), sysHistunfoldbin.GetBinContent(1), sysHistunfoldbin.GetBinContent(sysHistunfoldbin.GetNbinsX())
                      sysHistunfoldbin.Write()
                      integral+=sysHistunfoldbin.Integral()
                  #if j==3:
                  #  print "integral",integral

                canvas.cd(j+2)#j-2
                alt=cloneNormalize(sysHist)
                if (samples[i][0]!="QCD" or "ALT" in pre) and not sys in skipInSum:
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
              print("stat",sqrt(unc2)/datahist[j].GetBinContent(b+1))
              for sys in syss:
               if not sys in skipInSum:
                unc2+=max(pow(syshists[sys+"Up"].GetBinContent(b+1)-althists[j].GetBinContent(b+1),2),pow(syshists[sys+"Down"].GetBinContent(b+1)-althists[j].GetBinContent(b+1),2))
              print("stat+sys",sqrt(unc2)/datahist[j].GetBinContent(b+1))
              chi2+=pow(datahist[j].GetBinContent(b+1)-althists[j].GetBinContent(b+1),2)/unc2
            pvalue=stats.distributions.chi2.sf(chi2, datahist[j].GetNbinsX())
            sign=-stats.norm.ppf(pvalue)
            print("sign",sign,"chi2/ndof",chi2/datahist[j].GetNbinsX())
        canvas.SaveAs(prefix + "_"+samples[i][0].replace("QCD","") + '_sys_run'+run+'_smear.pdf')
        #canvas.SaveAs(prefix + "_"+samples[i][0].replace("QCD","") + '_sys_run2_smear.eps')

      for closefile in closefiles:
          closefile.Close()
