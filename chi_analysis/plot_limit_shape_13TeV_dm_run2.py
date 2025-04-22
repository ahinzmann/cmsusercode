from ROOT import *
import ROOT,sys,os,math
from math import *
import numpy as np
import array
import ctypes

#gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
gROOT.SetBatch(True)
gROOT.ForceStyle()
gStyle.SetOptStat(0)
gStyle.SetOptFit(0)
gStyle.SetTitleOffset(0.9,"XY")
gStyle.SetPadLeftMargin(0.12)
gStyle.SetPadBottomMargin(0.12)
gStyle.SetPadTopMargin(0.075)
gStyle.SetPadRightMargin(0.12)
gStyle.SetMarkerSize(1.5)
gStyle.SetHistLineWidth(1)
gStyle.SetStatFontSize(0.020)
gStyle.SetNdivisions(513, "Y")
gStyle.SetNdivisions(510, "X")
gStyle.SetLegendBorderSize(0)

def medWidth(gq):
  return 6*gq**2/(4*3.141592653)+1/(12*3.141592653)

if __name__=="__main__":
  
 #for style in ["DMVector","DMAxial"]:
 for style in ["DMAxial"]:

  testStat="LHC"
  asym="a" # asymptotic CLS
  #testStat="LEP"
  #asym=""
  version="_v6b"

  isGen=False
  isCB=False
  
  isInjection=False
  injectiontext="Injection1p0"
  
  if isGen:
    prefix="limitsGen"
  elif isCB:
    prefix="limitsDetCB"
  else:
    prefix="limitsDet"
    
  signalCounter={}
  if style=="DMVector":
    counter=100
  else:
    counter=1100

  gs=["0p01","0p05","0p1","0p2","0p25","0p3","0p5","0p75","1","1p5","2p0","2p5","3p0"]
  gsplot=["0p1","0p2","0p25","0p3","0p5","0p75","1","1p5"]
  #mdms=["1","3000"]
  mdms=["1"]
  #signalMasses=[1000,1500,1750,2000,2250,2500,3000,3500,4000,4500,5000,6000,7000]
  #signalMasses=[2000,2250,2500,3000,3500,4000,4500,5000,6000,7000]
  #signalMasses=[2000,2250,2500,3000,3500,4000,4500,5000,6000]
  signalMasses=[4000,4500,5000,6000,7000]
 
  for mdm in mdms:
    for g in gs:
      signalCounter[style+"_mdm"+mdm+"_g"+g]=counter
      counter+=1

  #for mdm in mdms:
  for mdm in ["1"]:
    g_q_out=TFile(prefix+testStat+asym+"_"+style+"_mdm"+mdm+version+'.root',"RECREATE")
    
    g_q=TGraph()
    g_q_exp=TGraph()
    g_q_exp1m=TGraph()
    g_q_exp1p=TGraph()
    g_q_exp2m=TGraph()
    g_q_exp2p=TGraph()
    g_q_band=TGraphAsymmErrors()
    g_q_band_2sigma=TGraphAsymmErrors()
    for signalMass in signalMasses:
      print("--------------------")
      signal=style+"_Mphi"+str(signalMass)+"_mdm"+mdm

      if isInjection:
        signal=signal.replace("Mphi","Mphi"+injectiontext)

      limits=[]
      for g in gsplot:
        if testStat!="LEP":
          try:
            if isInjection:
              f=open(prefix+testStat+asym+str(signalCounter[style+"_mdm"+mdm+"_g"+g])+"_"+style+"_Dijet_LO_Mphi"+injectiontext+"_exp_"+str(signalMass)+"_run2"+version+".txt")
            else:
              f=open(prefix+testStat+asym+str(signalCounter[style+"_mdm"+mdm+"_g"+g])+"_"+style+"_Dijet_LO_Mphi_exp_"+str(signalMass)+"_run2"+version+".txt")
          except:
            print("can't open",prefix+testStat+asym+str(signalCounter[style+"_mdm"+mdm+"_g"+g])+"_"+style+"_Dijet_LO_Mphi_exp_"+str(signalMass)+"_run2"+version+".txt")
            continue
        else:
          try:
            if isInjection:
              f=open(prefix+testStat+asym+str(signalCounter[style+"_mdm"+mdm+"_g"+g])+"_"+style+"_Dijet_LO_Mphi"+injectiontext+"_"+str(signalMass)+"_run2"+version+".txt")
            else:
              f=open(prefix+testStat+asym+str(signalCounter[style+"_mdm"+mdm+"_g"+g])+"_"+style+"_Dijet_LO_Mphi_"+str(signalMass)+"_run2"+version+".txt")
          except:
            print("can't open",prefix+testStat+asym+str(signalCounter[style+"_mdm"+mdm+"_g"+g])+"_"+style+"_Dijet_LO_Mphi_"+str(signalMass)+"_run2"+version+".txt")
            continue
          
        limits+=[[]]
        for lin in f.readlines():
         for line in lin.split("\\n"):
 
          if "Observed Limit" in line and asym:
            limits[-1]=[float(g.replace('p','.')),float(line.strip().split(" ")[-1]),0]
          
          if "CLs = " in line and testStat=="LEP":
            limits[-1]=[float(g.replace('p','.')),float(line.strip().split(" ")[-3]),float(line.strip().split(" ")[-1])]
            if limits[-1][-1]==0:
              limits[-1][-2]=1e-6
          
          if "Observed CLs = " in line and testStat!="LEP":
            limits[-1]=[float(g.replace('p','.')),float(line.strip().split(" ")[-1]),0]
 
          if "Significance:" in line and asym:
            print("observed signficance (p-value): ",line.strip().split(" ")[-1].strip(")"))
              
          if "CLb      = " in line and testStat=="LEP":
            print("observed signficance (p-value): ",ROOT.Math.normal_quantile_c((1.-float(line.strip().split(" ")[-3]))/2.,1),"(",(1.-float(line.strip().split(" ")[-3])),")")

          if "Observed CLb = " in line and testStat!="LEP":
            print("observed signficance (p-value): ",ROOT.Math.normal_quantile_c((1.-float(line.strip().split(" ")[-1]))/2.,1),"(",(1.-float(line.strip().split(" ")[-1])),")")
            
        if len(limits[-1])==0:
          limits[-1]+=[float(g.replace('p','.'))]
          limits[-1]+=[1e-6]
          limits[-1]+=[0]
          
        try:
          if isInjection:
            f=open(prefix+testStat+asym+str(signalCounter[style+"_mdm"+mdm+"_g"+g])+"_"+style+"_Dijet_LO_Mphi"+injectiontext+"_exp_"+str(signalMass)+"_run2"+version+".txt")
          else:
            f=open(prefix+testStat+asym+str(signalCounter[style+"_mdm"+mdm+"_g"+g])+"_"+style+"_Dijet_LO_Mphi_exp_"+str(signalMass)+"_run2"+version+".txt")
        except:
          print("can't open",prefix+testStat+asym+str(signalCounter[style+"_mdm"+mdm+"_g"+g])+"_"+style+"_Dijet_LO_Mphi_exp_"+str(signalMass)+"_run2"+version+".txt")
          del limits[-1]
          continue
        #for line in f.readlines():
        for lin in f.readlines():
         for line in lin.split("\\n"):
          if "Expected" in line and asym:
              limits[-1]+=[float(line.strip().split(" ")[-1])]
          if "Expected CLs" in line:
            try:
              limits[-1]+=[float(line.strip().split(" ")[-1])]
            except:
              print("didn't find one point")
            
        for i in range(len(limits[-1]),8):
          limits[-1]+=[0]

        for i in range(len(limits[-1])):
          if limits[-1][i]==0 and i>=3:
            limits[-1][i]=1e-6

      if asym: # reorder expected limits to exp,-1,+1,-2,+2
        for g in range(len(limits)):
          limits[g]=[limits[g][0],limits[g][1],limits[g][2],limits[g][5],limits[g][6],limits[g][4],limits[g][7],limits[g][3]]

      print(limits)

      canvas = TCanvas("","",0,0,300,300)
      #canvas.GetPad(0).SetLogy()
      canvas.GetPad(0).SetLogx()
      mg=TMultiGraph()

      min_x=0.1
      max_x=11

      g=TGraph()
      g_exp_=TGraph()
      g_exp1m=TGraph()
      g_exp1p=TGraph()
      g_exp2m=TGraph()
      g_exp2p=TGraph()
      max_mass=0
      max_mass_exp=0
      max_mass_exp1m=0
      max_mass_exp1p=0
      for lim in limits:
        mass,limit,error,exp,exp1m,exp1p,exp2m,exp2p =lim
        if limit>0:
          g.SetPoint(g.GetN(),mass,log10(limit))
          max_mass=mass
        if exp>0:
          g_exp_.SetPoint(g_exp_.GetN(),mass,log10(exp))
          max_mass_exp=mass
        if exp1m>0:
          g_exp1m.SetPoint(g_exp1m.GetN(),mass,log10(exp1m))
          max_mass_exp1m=mass
        if exp1p>0:
          g_exp1p.SetPoint(g_exp1p.GetN(),mass,log10(exp1p))
          max_mass_exp1p=mass
        if exp2m>0:
          g_exp2m.SetPoint(g_exp2m.GetN(),mass,log10(exp2m))
          max_mass_exp2m=mass
        if exp2p>0:
          g_exp2p.SetPoint(g_exp2p.GetN(),mass,log10(exp2p))
          max_mass_exp2p=mass

      g.SetMarkerStyle(24)
      g.SetMarkerSize(0.5)
      g.SetLineColor(1)
      g.SetLineWidth(3)
      mg.Add(g)
      g_exp_.SetMarkerStyle(24)
      g_exp_.SetMarkerSize(0.5)
      g_exp_.SetLineColor(2)
      g_exp_.SetLineWidth(3)
      mg.Add(g_exp_)
      g_exp1m.SetMarkerStyle(24)
      g_exp1m.SetMarkerSize(0.5)
      g_exp1m.SetLineColor(3)
      g_exp1m.SetLineWidth(3)
      mg.Add(g_exp1m)
      g_exp1p.SetMarkerStyle(24)
      g_exp1p.SetMarkerSize(0.5)
      g_exp1p.SetLineColor(3)
      g_exp1p.SetLineWidth(3)
      mg.Add(g_exp1p)

      g_exp2m.SetMarkerStyle(24)
      g_exp2m.SetMarkerSize(0.5)
      g_exp2m.SetLineColor(5)
      g_exp2m.SetLineWidth(3)
      mg.Add(g_exp2m)
      g_exp2p.SetMarkerStyle(24)
      g_exp2p.SetMarkerSize(0.5)
      g_exp2p.SetLineColor(5)
      g_exp2p.SetLineWidth(3)
      mg.Add(g_exp2p)
  
      mg.Draw("apl")
      mg.SetTitle("")
      mg.GetXaxis().SetTitle("g_{q}")
      if asym:
        mg.GetYaxis().SetTitle("log_{10}(signal strength)")
        miny=-1
        maxy=1
      else:
        mg.GetYaxis().SetTitle("log_{10}(CL_{S})")
        miny=-6
        maxy=0
      mg.GetYaxis().SetRangeUser(miny,maxy)
      mg.GetXaxis().SetLimits(min_x, max_x)

      min_x=0.01
  
      if asym:
        cut=1
      else:
        cut=0.05

      l=TLine(min_x,log10(cut),max_x,log10(cut))
      l.SetLineColor(2)
      l.SetLineStyle(2)
      l.Draw("same")
      
      if asym:
        l1=TLatex((max_x-min_x)*0.5+min_x,log10(cut)*1,"signal strength = 1")
      else:
        l1=TLatex((max_x-min_x)*0.5+min_x,log10(cut)*1,"CL_{S}=0.05")
      l1.Draw("same")

      limit=0
      exp=0
      exp1m=0
      exp1p=0
      exp2m=0
      exp2p=0
      nseg=10000
      masses=np.linspace(max_x,min_x,num=nseg)      
      if (signalMass==7000 or signalMass==6000):
        for mass in masses:  
          #if limit==0 and g.Eval(mass,0)>log10(cut) and mass<=max_mass:
          if limit==0 and g.Eval(mass,0)>log10(cut):
            limit=mass
          #if exp==0 and g_exp_.Eval(mass,0)>log10(cut) and mass<=max_mass_exp:
          if exp==0 and g_exp_.Eval(mass,0)>log10(cut):
            exp=mass
            #print "exp:",i,mass
          #if exp1m==0 and g_exp1m.Eval(mass,0)>log10(cut) and mass<=max_mass_exp1m:
          if exp1m==0 and g_exp1m.Eval(mass,0)>log10(cut):
            #print "exp1m:",i,mass
            exp1m=mass
          #if exp1p==0 and g_exp1p.Eval(mass,0)>log10(cut) and mass<=max_mass_exp1p:
          if exp1p==0 and g_exp1p.Eval(mass,0)>log10(cut):
            exp1p=mass
            #print "exp1p:",i,mass
          #if exp2m==0 and g_exp2m.Eval(mass,0)>log10(cut) and mass<=max_mass_exp2m:
          if exp2m==0 and g_exp2m.Eval(mass,0)>log10(cut):
            #print "exp2m:",i,mass
            exp2m=mass
          #if exp2p==0 and g_exp2p.Eval(mass,0)>log10(cut) and mass<=max_mass_exp2p:
          if exp2p==0 and g_exp2p.Eval(mass,0)>log10(cut):  
            exp2p=mass
            #print "exp2p:",i,mass        
      else:
        for mass in masses:  
          if limit==0 and g.Eval(mass,0)>log10(cut) and mass<=max_mass:
          #if limit==0 and g.Eval(mass,0)>log10(cut):
            limit=mass
          if exp==0 and g_exp_.Eval(mass,0)>log10(cut) and mass<=max_mass_exp:
            #if exp==0 and g_exp_.Eval(mass,0)>log10(cut):
            exp=mass
            #print "exp:",i,mass
          if exp1m==0 and g_exp1m.Eval(mass,0)>log10(cut) and mass<=max_mass_exp1m:
            #if exp1m==0 and g_exp1m.Eval(mass,0)>log10(cut):
            #print "exp1m:",i,mass
            exp1m=mass
          if exp1p==0 and g_exp1p.Eval(mass,0)>log10(cut) and mass<=max_mass_exp1p:
          #if exp1p==0 and g_exp1p.Eval(mass,0)>log10(cut):
            exp1p=mass
            #print "exp1p:",i,mass
          if exp2m==0 and g_exp2m.Eval(mass,0)>log10(cut) and mass<=max_mass_exp2m:
          #if exp2m==0 and g_exp2m.Eval(mass,0)>log10(cut):
            #print "exp2m:",i,mass
            exp2m=mass
          if exp2p==0 and g_exp2p.Eval(mass,0)>log10(cut) and mass<=max_mass_exp2p:
          #if exp2p==0 and g_exp2p.Eval(mass,0)>log10(cut):  
            exp2p=mass
            #print "exp2p:",i,mass


      print("mass: %.0f" % (signalMass), "g_q limit: %.2f" % (limit), "& %.2f" % (exp), "$\pm$ %.2f" % (max(exp-exp1p,exp1m-exp)))

      #print "limit: %.6f," % (limit), "%.6f," % (exp), "%.6f, %.6f, 0, 0" % ((-max(exp-exp1p,exp1m-exp)+exp),(exp+max(exp-exp1p,exp1m-exp)))
      print("mass: %.0f" % (signalMass), "g_q limit: %.6f," % (limit), "%.6f," % (exp), "%.6f, %.6f, %.6f, %.6f" % (exp-min(exp1p,exp1m),max(exp1p,exp1m)-exp, exp-min(exp2p,exp2m), max(exp2p,exp2m)-exp))

      signalMass=signalMass/1000.
 
      g_q.SetPoint(g_q.GetN(),signalMass,limit)
      g_q_exp.SetPoint(g_q_exp.GetN(),signalMass,exp)
      g_q_exp1m.SetPoint(g_q_exp1m.GetN(),signalMass,exp1m)
      g_q_exp1p.SetPoint(g_q_exp1p.GetN(),signalMass,exp1p)
      g_q_exp2m.SetPoint(g_q_exp2m.GetN(),signalMass,exp2m)
      g_q_exp2p.SetPoint(g_q_exp2p.GetN(),signalMass,exp2p)
      g_q_band.SetPoint(g_q_band.GetN(),signalMass,exp)
      g_q_band_2sigma.SetPoint(g_q_band_2sigma.GetN(),signalMass,exp)

      #g_q_band.SetPointError(g_q_band.GetN()-1,0,0,max(exp-exp1p,exp1m-exp),max(exp-exp1p,exp1m-exp))
      g_q_band.SetPointError(g_q_band.GetN()-1,0,0,exp-min(exp1p,exp1m),max(exp1p,exp1m)-exp)
      #print(signalMass,exp,exp1p,exp1m,exp-min(exp1p,exp1m),max(exp1p,exp1m)-exp)
      #print(signalMass,exp,exp2p,exp2m,exp-min(exp2p,exp2m),max(exp2p,exp2m)-exp)
      g_q_band_2sigma.SetPointError(g_q_band_2sigma.GetN()-1,0,0,exp-min(exp2p,exp2m),max(exp2p,exp2m)-exp)
      #print "g_q_band_2sigma.GetErrorYlow(g_q_band_2sigma.GetN()-1),g_q_band.GetErrorYhigh(g_q_band_2sigma.GetN()-1):",g_q_band_2sigma.GetErrorYlow(g_q_band_2sigma.GetN()-1),g_q_band_2sigma.GetErrorYhigh(g_q_band_2sigma.GetN()-1)
  
      l2=TLine(limit,miny,limit,log10(cut))
      l2.SetLineColor(1)
      l2.SetLineStyle(2)
      l2.Draw("same")
      
      l2a=TLine(exp,miny,exp,log10(cut))
      l2a.SetLineColor(2)
      l2a.SetLineStyle(2)
      l2a.Draw("same")
  
      l2b=TLine(exp1m,miny,exp1m,log10(cut))
      l2b.SetLineColor(3)
      l2b.SetLineStyle(2)
      l2b.Draw("same")
  
      l2c=TLine(exp1p,miny,exp1p,log10(cut))
      l2c.SetLineColor(3)
      l2c.SetLineStyle(2)
      l2c.Draw("same")
  
      l2d=TLine(exp2m,miny,exp2m,log10(cut))
      l2d.SetLineColor(5)
      l2d.SetLineStyle(2)
      l2d.Draw("same")
  
      l2e=TLine(exp2p,miny,exp2p,log10(cut))
      l2e.SetLineColor(5)
      l2e.SetLineStyle(2)
      l2e.Draw("same")
  
      canvas.SaveAs(prefix+testStat+asym+signal+version+'.pdf')

    # from https://gitlab.cern.ch/cms-exo-mci/exo-dmsummaryplots/-/blob/master/GQSummary/data/EXO16046_obs.dat?ref_type=heads
    old_limits=[(2000.0,0.56026026026), (2010.0,0.554554554555), (2020.0,0.548848848849), (2030.0,0.545045045045), (2040.0,0.539339339339), (2050.0,0.535535535536), (2060.0,0.52982982983), (2070.0,0.524124124124), (2080.0,0.52032032032), (2090.0,0.514614614615), (2100.0,0.510810810811), (2110.0,0.505105105105), (2120.0,0.499399399399), (2130.0,0.495595595596), (2140.0,0.48988988989), (2150.0,0.486086086086), (2160.0,0.48038038038), (2170.0,0.474674674675), (2180.0,0.470870870871), (2190.0,0.465165165165), (2200.0,0.461361361361), (2210.0,0.455655655656), (2220.0,0.44994994995), (2230.0,0.446146146146), (2240.0,0.44044044044), (2250.0,0.434734734735), (2260.0,0.432832832833), (2270.0,0.429029029029), (2280.0,0.427127127127), (2290.0,0.423323323323), (2300.0,0.41951951952), (2310.0,0.417617617618), (2320.0,0.413813813814), (2330.0,0.411911911912), (2340.0,0.408108108108), (2350.0,0.406206206206), (2360.0,0.402402402402), (2370.0,0.398598598599), (2380.0,0.396696696697), (2390.0,0.392892892893), (2400.0,0.390990990991), (2410.0,0.387187187187), (2420.0,0.383383383383), (2430.0,0.381481481481), (2440.0,0.377677677678), (2450.0,0.375775775776), (2460.0,0.371971971972), (2470.0,0.368168168168), (2480.0,0.366266266266), (2490.0,0.362462462462), (2500.0,0.360560560561), (2510.0,0.362462462462), (2520.0,0.362462462462), (2530.0,0.364364364364), (2540.0,0.366266266266), (2550.0,0.368168168168), (2560.0,0.37007007007), (2570.0,0.371971971972), (2580.0,0.373873873874), (2590.0,0.375775775776), (2600.0,0.377677677678), (2610.0,0.37957957958), (2620.0,0.381481481481), (2630.0,0.383383383383), (2640.0,0.385285285285), (2650.0,0.385285285285), (2660.0,0.387187187187), (2670.0,0.389089089089), (2680.0,0.390990990991), (2690.0,0.392892892893), (2700.0,0.394794794795), (2710.0,0.396696696697), (2720.0,0.398598598599), (2730.0,0.400500500501), (2740.0,0.402402402402), (2750.0,0.404304304304), (2760.0,0.406206206206), (2770.0,0.408108108108), (2780.0,0.408108108108), (2790.0,0.41001001001), (2800.0,0.411911911912), (2810.0,0.413813813814), (2820.0,0.415715715716), (2830.0,0.417617617618), (2840.0,0.41951951952), (2850.0,0.421421421421), (2860.0,0.423323323323), (2870.0,0.425225225225), (2880.0,0.427127127127), (2890.0,0.429029029029), (2900.0,0.430930930931), (2910.0,0.430930930931), (2920.0,0.432832832833), (2930.0,0.434734734735), (2940.0,0.436636636637), (2950.0,0.438538538539), (2960.0,0.44044044044), (2970.0,0.442342342342), (2980.0,0.444244244244), (2990.0,0.446146146146), (3000.0,0.448048048048), (3010.0,0.446146146146), (3020.0,0.446146146146), (3030.0,0.446146146146), (3040.0,0.444244244244), (3050.0,0.444244244244), (3060.0,0.444244244244), (3070.0,0.444244244244), (3080.0,0.442342342342), (3090.0,0.442342342342), (3100.0,0.442342342342), (3110.0,0.44044044044), (3120.0,0.44044044044), (3130.0,0.44044044044), (3140.0,0.438538538539), (3150.0,0.438538538539), (3160.0,0.438538538539), (3170.0,0.436636636637), (3180.0,0.436636636637), (3190.0,0.436636636637), (3200.0,0.434734734735), (3210.0,0.434734734735), (3220.0,0.434734734735), (3230.0,0.434734734735), (3240.0,0.432832832833), (3250.0,0.432832832833), (3260.0,0.432832832833), (3270.0,0.430930930931), (3280.0,0.430930930931), (3290.0,0.430930930931), (3300.0,0.429029029029), (3310.0,0.429029029029), (3320.0,0.429029029029), (3330.0,0.427127127127), (3340.0,0.427127127127), (3350.0,0.427127127127), (3360.0,0.427127127127), (3370.0,0.425225225225), (3380.0,0.425225225225), (3390.0,0.425225225225), (3400.0,0.423323323323), (3410.0,0.423323323323), (3420.0,0.423323323323), (3430.0,0.421421421421), (3440.0,0.421421421421), (3450.0,0.421421421421), (3460.0,0.41951951952), (3470.0,0.41951951952), (3480.0,0.41951951952), (3490.0,0.41951951952), (3500.0,0.417617617618), (3510.0,0.417617617618), (3520.0,0.41951951952), (3530.0,0.41951951952), (3540.0,0.41951951952), (3550.0,0.41951951952), (3560.0,0.41951951952), (3570.0,0.41951951952), (3580.0,0.41951951952), (3590.0,0.421421421421), (3600.0,0.421421421421), (3610.0,0.421421421421), (3620.0,0.421421421421), (3630.0,0.421421421421), (3640.0,0.421421421421), (3650.0,0.421421421421), (3660.0,0.423323323323), (3670.0,0.423323323323), (3680.0,0.423323323323), (3690.0,0.423323323323), (3700.0,0.423323323323), (3710.0,0.423323323323), (3720.0,0.423323323323), (3730.0,0.425225225225), (3740.0,0.425225225225), (3750.0,0.425225225225), (3760.0,0.425225225225), (3770.0,0.425225225225), (3780.0,0.425225225225), (3790.0,0.427127127127), (3800.0,0.427127127127), (3810.0,0.427127127127), (3820.0,0.427127127127), (3830.0,0.427127127127), (3840.0,0.427127127127), (3850.0,0.427127127127), (3860.0,0.429029029029), (3870.0,0.429029029029), (3880.0,0.429029029029), (3890.0,0.429029029029), (3900.0,0.429029029029), (3910.0,0.429029029029), (3920.0,0.429029029029), (3930.0,0.430930930931), (3940.0,0.430930930931), (3950.0,0.430930930931), (3960.0,0.430930930931), (3970.0,0.430930930931), (3980.0,0.430930930931), (3990.0,0.432832832833), (4000.0,0.432832832833), (4010.0,0.438538538539), (4020.0,0.446146146146), (4030.0,0.453753753754), (4040.0,0.461361361361), (4050.0,0.467067067067), (4060.0,0.474674674675), (4070.0,0.482282282282), (4080.0,0.487987987988), (4090.0,0.495595595596), (4100.0,0.503203203203), (4110.0,0.508908908909), (4120.0,0.516516516517), (4130.0,0.524124124124), (4140.0,0.531731731732), (4150.0,0.537437437437), (4160.0,0.545045045045), (4170.0,0.552652652653), (4180.0,0.558358358358), (4190.0,0.565965965966), (4200.0,0.573573573574), (4210.0,0.579279279279), (4220.0,0.586886886887), (4230.0,0.594494494494), (4240.0,0.6002002002), (4250.0,0.607807807808), (4260.0,0.615415415415), (4270.0,0.621121121121), (4280.0,0.628728728729), (4290.0,0.634434434434), (4300.0,0.642042042042), (4310.0,0.64964964965), (4320.0,0.655355355355), (4330.0,0.662962962963), (4340.0,0.670570570571), (4350.0,0.676276276276), (4360.0,0.683883883884), (4370.0,0.691491491491), (4380.0,0.697197197197), (4390.0,0.704804804805), (4400.0,0.710510510511), (4410.0,0.718118118118), (4420.0,0.725725725726), (4430.0,0.731431431431), (4440.0,0.739039039039), (4450.0,0.744744744745), (4460.0,0.752352352352), (4470.0,0.75995995996), (4480.0,0.765665665666), (4490.0,0.773273273273), (4500.0,0.780880880881), (4510.0,0.792292292292), (4520.0,0.805605605606), (4530.0,0.818918918919), (4540.0,0.832232232232), (4550.0,0.845545545546), (4560.0,0.858858858859), (4570.0,0.872172172172), (4580.0,0.885485485485), (4590.0,0.898798798799), (4600.0,0.912112112112), (4610.0,0.925425425425), (4620.0,0.938738738739), (4630.0,0.95015015015), (4640.0,0.963463463463), (4650.0,0.976776776777), (4660.0,0.99009009009), (4670.0,1.0034034034), (4680.0,1.01671671672), (4690.0,1.03003003003), (4700.0,1.04334334334), (4710.0,1.05475475475), (4720.0,1.06806806807), (4730.0,1.08138138138), (4740.0,1.09469469469), (4750.0,1.10800800801), (4760.0,1.12132132132), (4770.0,1.13463463463), (4780.0,1.14604604605), (4790.0,1.15935935936), (4800.0,1.17267267267), (4810.0,1.18598598599), (4820.0,1.1992992993), (4830.0,1.21261261261), (4840.0,1.22402402402), (4850.0,1.23733733734), (4860.0,1.25065065065), (4870.0,1.26396396396), (4880.0,1.27727727728), (4890.0,1.28868868869), (4900.0,1.302002002), (4910.0,1.31531531532), (4920.0,1.32862862863), (4930.0,1.34194194194), (4940.0,1.35525525526), (4950.0,1.36666666667), (4960.0,1.37997997998), (4970.0,1.39329329329), (4980.0,1.40660660661), (4990.0,1.41991991992), (5000.0,1.43133133133), (5010.0,1.44464464464), (5020.0,1.45605605606), (5030.0,1.46936936937), (5040.0,1.48268268268), (5050.0,1.49409409409), (5060.0,1.50740740741), (5070.0,1.51881881882), (5080.0,1.53213213213), (5090.0,1.54354354354), (5100.0,1.55685685686), (5110.0,1.56826826827), (5120.0,1.58158158158), (5130.0,1.59489489489), (5140.0,1.60630630631), (5150.0,1.61961961962), (5160.0,1.63103103103), (5170.0,1.64434434434), (5180.0,1.65575575576), (5190.0,1.66906906907), (5200.0,1.68048048048), (5210.0,1.69379379379), (5220.0,1.70520520521), (5230.0,1.71851851852), (5240.0,1.72992992993), (5250.0,1.74324324324), (5260.0,1.75655655656), (5270.0,1.76796796797), (5280.0,1.78128128128), (5290.0,1.79269269269), (5300.0,1.80600600601), (5310.0,1.81741741742), (5320.0,1.83073073073), (5330.0,1.84214214214), (5340.0,1.85545545546), (5350.0,1.86686686687), (5360.0,1.88018018018), (5370.0,1.89159159159), (5380.0,1.9049049049), (5390.0,1.91631631632), (5400.0,1.92962962963), (5410.0,1.94104104104), (5420.0,1.95435435435), (5430.0,1.96576576577), (5440.0,1.97907907908), (5450.0,1.99239239239)]
    # from https://gitlab.cern.ch/cms-exo-mci/exo-dmsummaryplots/-/blob/master/GQSummary/data/cl90/EXO19012_obs.dat?ref_type=heads
    exo19012limit=[(1800,0.101873), (1900,0.13199), (2000,0.136569), (2100,0.114894), (2200,0.0729184), (2300,0.0586678), (2400,0.0678256), (2500,0.103493), (2600,0.145027), (2700,0.166081), (2800,0.183731), (2900,0.199127), (3000,0.21702), (3100,0.211135), (3200,0.19009), (3300,0.181343), (3400,0.19753), (3500,0.234712), (3600,0.263661), (3700,0.303623), (3800,0.345381), (3900,0.359206), (4000,0.371823), (4100,0.355785), (4200,0.344932), (4300,0.349693), (4400,0.372859), (4500,0.451645), (4600,0.528722), (4700,0.581656), (4800,0.650748)]

    bins=array.array('d',[3.75,4.25,4.75,5.25,6.75,7.25])
    h_band_low=TH1F("low","low",len(signalMasses),0,7)
    h_band_low.SetBins(len(signalMasses),bins)
    h_band_high=TH1F("high","high",len(signalMasses),0,7)
    h_band_high.SetBins(len(signalMasses),bins)
    h_band_low2=TH1F("low2","low2",len(signalMasses),0,7)
    h_band_low2.SetBins(len(signalMasses),bins)
    h_band_high2=TH1F("high2","high2",len(signalMasses),0,7)
    h_band_high2.SetBins(len(signalMasses),bins)
    for m in range(len(signalMasses)):
      d1, d2 = ctypes.c_double(0.), ctypes.c_double(0.)
      g_q_band.GetPoint(m,d1,d2)
      h_band_low.SetBinContent(m+1,d2.value-g_q_band.GetErrorYlow(m))
      h_band_high.SetBinContent(m+1,d2.value+g_q_band.GetErrorYhigh(m))
      h_band_low2.SetBinContent(m+1,d2.value-g_q_band_2sigma.GetErrorYlow(m))
      h_band_high2.SetBinContent(m+1,d2.value+g_q_band_2sigma.GetErrorYhigh(m))
    
    g_q_old=TGraph()
    for mass, limit in old_limits:
      g_q_old.SetPoint(g_q_old.GetN(),mass/1000.,limit)

    g_q_resonance=TGraph()
    for mass, limit in exo19012limit:
      g_q_resonance.SetPoint(g_q_resonance.GetN(),mass/1000.,limit)

    canvas = TCanvas("","",0,0,1800,1550)
    mg=TMultiGraph()

    min_x_new=signalMasses[0]/1000
    #min_x_new=2000

    ymin=0.0
    ymax=1.42788
    #ymax=3
    xs=np.linspace(4,7,num=20000)
    for x in xs:
      if g_q_exp.Eval(x)>=1.42788:
        max_x_new=x
        break
    max_x_new=signalMasses[-1]/1000

    limit=0
    exp=0
    exp1m=0
    exp1p=0
    exp2m=0
    exp2p=0
    cut=1
    nseg=10000
    masses=np.linspace(0.8*min_x_new,max_x_new*1.2,num=nseg)      
    for mass in masses:
      if limit==0 and g_q.Eval(mass,0)>cut:
        limit=mass
      if exp==0 and g_q_exp.Eval(mass,0)>cut:
        exp=mass
      if exp1m==0 and g_q_exp1m.Eval(mass,0)>cut:
        exp1m=mass
      if exp1p==0 and g_q_exp1p.Eval(mass,0)>cut:
        exp1p=mass
      if exp2m==0 and g_q_exp2m.Eval(mass,0)>cut:
        exp2m=mass
      if exp2p==0 and g_q_exp2p.Eval(mass,0)>cut:  
        exp2p=mass
    print("mass limit: %.1f" % (limit), "& %.1f" % (exp), "$\pm$ %.1f" % (max(exp-exp1m,exp1p-exp)))

    print("mass limit: %.2f," % (limit), "%.2f," % (exp), "%.2f, %.2f, %.2f, %.2f" % (exp-min(exp1p,exp1m),max(exp1p,exp1m)-exp, exp-min(exp2p,exp2m), max(exp2p,exp2m)-exp))
      
    g_q_band_2sigma.SetFillStyle(1001)
    g_q_band_2sigma.SetFillColor(kOrange)
    g_q_band_2sigma.SetLineColor(1)
    g_q_band_2sigma.SetLineStyle(1)
    g_q_band_2sigma.SetLineWidth(0)
    #mg.Add(g_q_band_2sigma,"3")

    g_q_band.SetFillStyle(1001)
    g_q_band.SetFillColor(kGreen-3)
    g_q_band.SetLineColor(1)
    g_q_band.SetLineStyle(1)
    g_q_band.SetLineWidth(0)
    #mg.Add(g_q_band,"3")

    g_q_exp.SetLineColor(1)
    g_q_exp.SetLineWidth(3)
    g_q_exp.SetLineStyle(kDashed)

    #mg.Add(g_q_old,"l")
    g_q_old.SetLineColor(4)
    g_q_old.SetLineWidth(2)
    g_q_old.SetLineStyle(kDotted)

    #mg.Add(g_q_resonance,"l")
    g_q_resonance.SetLineColor(2)
    g_q_resonance.SetLineWidth(2)
    g_q_resonance.SetLineStyle(kDotted)

    #mg.Add(g_q_exp,"l")
    g_q.SetMarkerStyle(24)
    g_q.SetMarkerSize(0)
    g_q.SetLineColor(1)
    g_q.SetLineWidth(3)
    #mg.Add(g_q,"pl")

    g_q.SetName("gq_obs")
    g_q_exp.SetName("gq_exp")
    g_q.Write()
    g_q_exp.Write()
    g_q_band.SetName("gq_exp_1sigma")
    g_q_band.Write()
    g_q_band_2sigma.SetName("g_q_exp_2sigma")
    g_q_band_2sigma.Write()
    mg.SetName("chi")
    mg.Write()
    
    mg=g_q
    mg.Draw("apl")
    mg.SetTitle("")
    mg.GetXaxis().SetTitle("M_{Med} [TeV]")
    mg.GetYaxis().SetTitle("g_{q}")
    mg.GetYaxis().SetLabelSize(0.04)
    mg.GetYaxis().SetTitleSize(0.05)
    mg.GetYaxis().SetTitleOffset(0.95)
    mg.GetXaxis().SetLimits(min_x_new,max_x_new)
    #mg.GetXaxis().SetRangeUser(min_x_new,max_x_new)
    mg.GetYaxis().SetRangeUser(ymin,ymax)
    mg.GetYaxis().SetNdivisions(510)
    mg.GetXaxis().SetNdivisions(510)
    mg.GetXaxis().SetLabelSize(0.04)
    mg.GetXaxis().SetLabelOffset(0.015)
    mg.GetXaxis().SetTitleSize(0.05)
    mg.GetXaxis().SetTitleOffset(1.1)

    h_band_high2.SetFillStyle(1001)
    h_band_high2.SetFillColor(kOrange)
    h_band_high2.SetLineColor(kOrange)
    h_band_high2.SetLineStyle(1)
    h_band_high2.SetLineWidth(0)
    h_band_high2.Draw("lf2same")
    h_band_high.SetFillStyle(1001)
    h_band_high.SetFillColor(kGreen-3)
    h_band_high.SetLineColor(kGreen-3)
    h_band_high.SetLineStyle(1)
    h_band_high.SetLineWidth(0)
    h_band_high.Draw("lf2same")
    h_band_low.SetFillStyle(1001)
    h_band_low.SetFillColor(kOrange)
    h_band_low.SetLineColor(kOrange)
    h_band_low.SetLineStyle(1)
    h_band_low.SetLineWidth(0)
    h_band_low.Draw("lf2same")
    h_band_low2.SetFillStyle(1001)
    h_band_low2.SetFillColor(10)
    h_band_low2.SetLineColor(1)
    h_band_low2.SetLineStyle(1)
    h_band_low2.SetLineWidth(0)
    h_band_low2.Draw("lf2same")

    g_q_old.Draw("plsame")
    g_q_resonance.Draw("plsame")
    g_q_exp.Draw("plsame")
    g_q.Draw("plsame")

    y1=TGaxis(min_x_new, ymin, min_x_new, ymax, ymin,ymax,510,"")
    y1.SetLabelSize(0)
    y1.SetTitleSize(0)
    y1.Draw()

    x2=TGaxis(min_x_new, ymax, max_x_new, ymax, min_x_new,max_x_new,510,"-")
    x2.SetLabelSize(0)
    x2.SetTitleSize(0)
    x2.Draw()

    minwidth=medWidth(ymin)
    myFunc=TF1("myFunc","pow((x-1/(12*3.141592653))*(4*3.141592653)/6,0.5)",minwidth,1)
    y2=TGaxis(max_x_new, ymin, max_x_new, ymax,"myFunc",510,"+L")
    y2.SetTitle("#Gamma/M_{Med}")
    y2.SetLabelSize(0.04)
    y2.SetTitleSize(0.05)
    y2.SetTitleOffset(1.)
    y2.SetTitleFont(42)
    y2.SetLabelFont(42)
    y2.SetLabelOffset(0.005)
    y2.Draw()

    l0p5=TLine(min_x_new,0.5,max_x_new,0.5)
    l0p5.SetLineColor(kGray+1)
    l0p5.SetLineStyle(kDashed)
    #l0p5.Draw("same")
      
    l0p5T=TLatex((max_x_new-min_x_new)*0.35+min_x_new,0.5-0.05,"g_{q}=0.5")
    l0p5T.SetTextSize(0.03)
    l0p5T.SetTextColor(kGray+1)
    #l0p5T.Draw("same")

    l1p0=TLine(min_x_new,1,max_x_new,1)
    l1p0.SetLineColor(kGray+3)
    l1p0.SetLineStyle(7)
    l1p0.Draw("same")
      
    l1p0T=TLatex((max_x_new-min_x_new)*0.5+min_x_new,1.0-0.07,"g_{q}=1.0")
    l1p0T.SetTextSize(0.04)
    l1p0T.SetTextFont(42)
    l1p0T.SetTextColor(kGray+3)
    l1p0T.Draw("same")
    
    if style=="DMAxial":
      #lt=TLatex(signalMasses[0]+100,1.27,"#splitline{Vector/Axial-Vector Mediator}{m_{DM} = "+mdm+" GeV, g_{DM} = 1.0}")
      lt=TLatex((signalMasses[0]+100)/1000.,0.14,"#splitline{#bf{Vector/Axial-Vector Mediator}}{#bf{m_{DM} = "+mdm+" GeV, g_{DM} = 1.0}}")
      #lt=TLatex(signalMasses[0]+100,1.31,"m_{DM} = "+mdm+" GeV, g_{DM} = 1.0")
    else:
      lt=TLatex((signalMasses[0]+100)/1000.,1.26,"#splitline{#bf{Vector Mediator}}{#bf{m_{DM} = "+mdm+" GeV, g_{DM} = 1.0}}")
    lt.SetTextSize(0.04)
    lt.SetTextFont(42)
    lt.Draw("same")
    
    #l=TLegend(0.13,0.5,0.43,0.72,"95% CL upper limits")
    l=TLegend(0.15,0.52,0.4,0.79,"95% CL upper limits")
    l.SetFillColor(0)
    l.SetTextFont(42)
    l.SetFillStyle(0)
    l.SetTextSize(0.04)
    l.SetShadowColor(0)
    l.AddEntry(g_q,"Observed","LP")
    l.AddEntry(g_q_exp,"Expected","LP")
    l.AddEntry(g_q_old,"2016 Observed","L")
    l.AddEntry(g_q_resonance,"Resonance Search Observed","L")
    l.AddEntry(g_q_band,"Expected #pm 1 s.d.","F")
    l.AddEntry(g_q_band_2sigma,"Expected #pm 2 s.d.","F")
    l.Draw()
    

    # CMS
    #leg2=TLatex(min_x_new,ymax+0.03,"#bf{CMS} #it{Preliminary}")
    cmsPos=min_x_new+220/1000.
    leg2=TLatex(cmsPos,ymax-0.17,"#bf{CMS}")
    leg2.SetTextFont(42)
    leg2.SetTextSize(0.06)
    # lumi
    #lumiPos=max_x_new-1700/1000.
    lumiPos=max_x_new-1000/1000.
    leg3=TLatex(lumiPos,ymax+0.03,"138 fb^{-1} (13 TeV)")
    leg3.SetTextFont(42)
    leg3.SetTextSize(0.045)
    leg2.Draw("same")
    leg3.Draw("same")
    
    #canvas.SaveAs('limits'+testStat+asym+"_"+style+"_mdm"+mdm+version+'_noSys.pdf')
    if isInjection:
      canvas.SaveAs(prefix+testStat+asym+"_"+style+"_mdm"+mdm+version+injectiontext+'_run2.pdf')
    else:
      canvas.SaveAs(prefix+testStat+asym+"_"+style+"_mdm"+mdm+version+'_run2.pdf')
