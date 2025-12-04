from ROOT import TF1,gROOT,gStyle,TGraphAsymmErrors,TFile,TH1F,TCanvas,kCyan,kRed,kAzure,kBlue,kMagenta,kTeal,kOrange,TPad,TLegend,TLatex,kGreen,kGray
import ROOT
import array, math
import os,sys
from math import *

def getFitResults(fitFile, treename):
    tree=fitFile.Get(treename)
    fitParameters=[]
    fitConstraints=[]
    for name in uncertaintynames:
      if tree.floatParsFinal().find(name):
        print(name, tree.floatParsFinal().find(name).getVal(), tree.floatParsFinal().find(name).getError())
        fitParameters.append(tree.floatParsFinal().find(name).getVal())
        fitConstraints.append(tree.floatParsFinal().find(name).getError())
      else:
        print(name, "not found")
        fitParameters.append(0)
        fitConstraints.append(1)
    return fitParameters, fitConstraints

def applyFitResults(fitParameters,fitConstraints,uncertainties,hist,hdata):
    hdataG=TGraphAsymmErrors(hdata.Clone(hdata.GetName()+"G"))
    alpha=1.-0.6827
    nevents=0
    for b in range(hdataG.GetN()):
        if unfoldedData:
            N=hdata.GetBinContent(b+1)
        elif hdata.GetBinContent(b+1)>0:
            N=1./pow(hdata.GetBinError(b+1)/hdata.GetBinContent(b+1),2)
        else:
            N=0
        print(N)
        nevents+=N
        L=0
        if N>0:
            L=ROOT.Math.gamma_quantile(alpha/2.,N,1.)
        U=ROOT.Math.gamma_quantile_c(alpha/2.,N+1,1.)
        if N>0:
          hdataG.SetPointEYlow(b,(N-L)/N*hdata.GetBinContent(b+1))
          hdataG.SetPointEYhigh(b,(U-N)/N*hdata.GetBinContent(b+1))
    print("data events:", nevents)

    hdataGsysstat=hdataG.Clone(hdataG.GetName()+"_sysstat")
    
    histbackup=hist.Clone(hist.GetName()+"_backup")
    for nn in range(len(fitParameters)):
        for b in range(hist.GetXaxis().GetNbins()):
            if fitParameters[nn]>0:
                hist.SetBinContent(b+1,hist.GetBinContent(b+1)*(1+fitParameters[nn]*(uncertainties[nn][0].GetBinContent(b+1)-uncertainties[nn][2].GetBinContent(b+1))/uncertainties[nn][2].GetBinContent(b+1)))
                histbackup.SetBinContent(b+1,histbackup.GetBinContent(b+1)*(1+fitParameters[nn]*(uncertainties[nn][0].GetBinContent(b+1)-uncertainties[nn][2].GetBinContent(b+1))/uncertainties[nn][2].GetBinContent(b+1)))
            else:
                hist.SetBinContent(b+1,hist.GetBinContent(b+1)*(1-fitParameters[nn]*(uncertainties[nn][1].GetBinContent(b+1)-uncertainties[nn][2].GetBinContent(b+1))/uncertainties[nn][2].GetBinContent(b+1)))
                histbackup.SetBinContent(b+1,histbackup.GetBinContent(b+1)*(1-fitParameters[nn]*(uncertainties[nn][1].GetBinContent(b+1)-uncertainties[nn][2].GetBinContent(b+1))/uncertainties[nn][2].GetBinContent(b+1)))
        	    
    h2new=histbackup.Clone("down"+str(massbins[massbin]))
    h3new=histbackup.Clone("up"+str(massbins[massbin]))
    for b in range(histbackup.GetXaxis().GetNbins()):
        theory_sumdown=0
        theory_sumup=0
        exp_sumdown=0
        exp_sumup=0
        nn=0
        for up,down,central in uncertainties:
            addup=fitConstraints[nn]*pow(max(0,up.GetBinContent(b+1)-central.GetBinContent(b+1),down.GetBinContent(b+1)-central.GetBinContent(b+1)),2)/pow(central.GetBinContent(b+1),2)
            adddown=fitConstraints[nn]*pow(max(0,central.GetBinContent(b+1)-up.GetBinContent(b+1),central.GetBinContent(b+1)-down.GetBinContent(b+1)),2)/pow(central.GetBinContent(b+1),2)
            if uncertaintynames[uncertainties.index([up,down,central])] in ["jer","prefire"] or "jes" in uncertaintynames[uncertainties.index([up,down,central])] or "model" in uncertaintynames[uncertainties.index([up,down,central])] or "JERtail" in uncertaintynames[uncertainties.index([up,down,central])] or "sim" in uncertaintynames[uncertainties.index([up,down,central])]:
                exp_sumup+=addup
                exp_sumdown+=adddown
                #print uncertaintynames[uncertainties.index([up,down,central])]
            else:
                theory_sumup+=addup
                theory_sumdown+=adddown
            nn+=1
        theory_sumdown=sqrt(theory_sumdown)
        theory_sumup=sqrt(theory_sumup)
        exp_sumdown=sqrt(exp_sumdown)
        exp_sumup=sqrt(exp_sumup)
        total_sumup=sqrt(pow(theory_sumup,2)+pow(exp_sumup,2))
        total_sumdown=sqrt(pow(theory_sumdown,2)+pow(exp_sumdown,2))

        hdataGsysstat.SetPointEXlow(b,0)
        hdataGsysstat.SetPointEXhigh(b,0)
        hdataGsysstat.SetPointEYlow(b,sqrt(pow(exp_sumdown*hdataG.GetY()[b],2)+pow(hdataG.GetErrorYlow(b),2)))
        hdataGsysstat.SetPointEYhigh(b,sqrt(pow(exp_sumup*hdataG.GetY()[b],2)+pow(hdataG.GetErrorYhigh(b),2)))
        
        h2new.SetBinContent(b+1,histbackup.GetBinContent(b+1)-total_sumdown*histbackup.GetBinContent(b+1))
        h3new.SetBinContent(b+1,histbackup.GetBinContent(b+1)+total_sumup*histbackup.GetBinContent(b+1))
        
    return h2new,h3new,hdataG,hdataGsysstat

def getUncertainties(fsys,basename,pname,masstext):
    uncertainties=[]
    for u in uncertaintynames:
        hup_orig=fsys.Get(basename+"_"+u+"Up")
        up=hup_orig.Clone(pname+'#chi'+masstext+"_rebin1_"+u+"Up")
        divBinWidth(up)

        hdown_orig=fsys.Get(basename+"_"+u+"Down")
        down=hdown_orig.Clone(pname+'#chi'+masstext+"_rebin1_"+u+"Down")
        divBinWidth(down)

        hcentral_orig=fsys.Get(basename)
        central=hcentral_orig.Clone(pname+'#chi'+masstext+"_rebin1")
        divBinWidth(central)
        
        uncertainties+=[[up,down,central]]

    return uncertainties
    
def divBinWidth(hist):
    #hist.SetLineColor(linecolor)
    #hist.SetLineStyle(linestyle)
    #hist.SetLineWidth(linewidth)
    #hist.Scale(1./hist.Integral())
    for b in range(hist.GetNbinsX()):
        hist.SetBinContent(b+1,hist.GetBinContent(b+1)/hist.GetBinWidth(b+1))
        hist.SetBinError(b+1,hist.GetBinError(b+1)/hist.GetBinWidth(b+1))
    return

def rebin(h1,nbins,binning):
    for b in range(h1.GetXaxis().GetNbins()):
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)*h1.GetBinWidth(b+1))
        h1.SetBinError(b+1,h1.GetBinError(b+1)*h1.GetBinWidth(b+1))
    h1=h1.Rebin(nbins,h1.GetName()+"_rebin",binning)
    h1.Scale(1./h1.Integral())
    for b in range(h1.GetXaxis().GetNbins()):
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)/h1.GetBinWidth(b+1))
        h1.SetBinError(b+1,h1.GetBinError(b+1)/h1.GetBinWidth(b+1))
    return h1

def rebin2(h1,nbins,binning):
    h1=h1.Rebin(nbins,h1.GetName()+"_rebin",binning)
    h1.Scale(1./h1.Integral())
    for b in range(h1.GetXaxis().GetNbins()):
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)/h1.GetBinWidth(b+1))
        h1.SetBinError(b+1,h1.GetBinError(b+1)/h1.GetBinWidth(b+1))
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

def setUpDMHists(hist,linecolor,linestyle,linewidth):
    hist.SetLineColor(linecolor)
    hist.SetLineStyle(linestyle)
    hist.SetLineWidth(linewidth)
    hist.Scale(1./hist.Integral())
    for b in range(hist.GetNbinsX()):
        hist.SetBinContent(b+1,hist.GetBinContent(b+1)/hist.GetBinWidth(b+1))
    return 

def divideAsymErrors(g1,h1,doEX):
    if not g1.GetN()==h1.GetNbinsX():
        print("divideAsymErrors Function Fails!!!", g1.GetN(),h1.GetNbinsX())

    g=TGraphAsymmErrors()
    for b in range(h1.GetNbinsX()):
        den=h1.GetBinContent(b+1)
        if doEX:
            bwidth=h1.GetBinWidth(b+1)/2
        else:
            bwidth=0
        if b>=g1.GetN() or g1.GetX()[b]!=h1.GetBinCenter(b+1): continue
        g.SetPoint(g.GetN(),g1.GetX()[b],g1.GetY()[b]/den)
        g.SetPointError(g.GetN()-1,bwidth,bwidth,g1.GetEYlow()[b]/den,g1.GetEYhigh()[b]/den)
        
    return g

def setupAsymErrors(g):
    g.SetMarkerStyle(20)
    g.SetMarkerSize(0.1)
    g.SetMarkerColor(1)
    g.SetLineColor(1)
    g.SetLineWidth(1)
    return

if __name__=="__main__":

  useNNLO=True # choice for QCD
  useM2=True # choice of mu-scale for QCD
  use_NNPDF3=True
  normalize=True
  
  if use_NNPDF3:
    pdfset="nn31nnlo"
  else:
    pdfset="ct14nnlo"
  if not useNNLO:
    pdfset="ct14nlo"
  if useM2:
    muScale="m2"
    muAltScale="pt12"
  else:
    muScale="pt12"
    muAltScale="m2"

  unfoldedData=True
  oldMeasurements=False
  compareRun3=False
  oldTheory=True
  signalsBSM=False
  signalsDM=False
  compareScales=True
  comparePDFs=True
  compareMu30=False
  compareMadgraph=False
  compareNoEWK=False
  comparePostfit=False
  
  lumi=138000*(1.-59.83/137.6*(1.57-0.87)/(2.*pi)) # SCALE CROSS SECTION TO ACCOUNT FOR HEM VETO

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
        (7000,13000),
  ]

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
  chi_binnings=[]
  for mass_bin in chi_bins:
      chi_binnings+=[array.array('d')]
      for chi_bin in mass_bin:
        chi_binnings[-1].append(chi_bin)

  #gROOT.Reset()
  gROOT.SetBatch(True)
  #gROOT.SetStyle("Plain")
  gStyle.SetOptStat(0)
  gStyle.SetOptFit(0)
  gStyle.SetTitleOffset(1.2,"Y")
  gStyle.SetPadLeftMargin(0.18)
  gStyle.SetPadBottomMargin(0.2)
  gStyle.SetPadTopMargin(0.055)
  gStyle.SetPadRightMargin(0.05)
  gStyle.SetMarkerSize(1.5)
  gStyle.SetHistLineWidth(1)
  gStyle.SetStatFontSize(0.025)
  gStyle.SetTitleSize(0.075, "XYZ")
  gStyle.SetLabelSize(0.0625, "XYZ")
  gStyle.SetNdivisions(510, "XYZ")
  gStyle.SetLegendBorderSize(0)
  gStyle.SetPadTickX(1)
  gStyle.SetPadTickY(1)
  gStyle.SetEndErrorSize(8)

  print("load CMS style")
  import CMS_lumi
  CMS_lumi.writeExtraText=False
  CMS_lumi.relPosX=0.13
  CMS_lumi.cmsTextSize=1.0
  #gROOT.LoadMacro("CMS_lumi.C");
  iPeriod = 5;        #// 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV 
  iPos = 33;
  #// second parameter in example_plot is iPos, which drives the position of the CMS logo in the plot
  #// iPos=11 : top-left, left-aligned
  #// iPos=33 : top-right, right-aligned
  #// iPos=22 : center, centered
  #// mode generally : 
  #//   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)
  print("make plots")

  for massbin in range(len(massbins)):
      
    massbintext = str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")
    
    prefix="datacard_shapelimit13TeV"

    c = TCanvas("combined", "combined", 0, 0, 650, 800)
    c.Divide(1,1)
    new_hists=[]
    if True:
        fdir='versions/run2UL'+('NNLO' if useNNLO else 'NLO')+'_'+muScale+("_NNPDF3" if use_NNPDF3 else "")+'/'
            
        if unfoldedData:  
            if useNNLO:
              #filename1nu2="fastnlo/NNLO/2jet.NNLO.fnl5662j_mjj_chi_ct14nnlo_cppread_mu_"+muScale+".root"
              filename1nu2="fastnlo/NNLO/2jet.NNLO.fnl5662j_mjj_chi_norm_v25_"+pdfset+"_cppread_mu_"+muScale+".root"
            else:
              filename1nu2="fastnlo/NNLO/2jet.NNLO.fnl5662j_mjj_chi_ct14nlo_cppread_mu_pt12.root"
              #filename1nu2="fastnlo/NNLO/2jet.NNLO.fnl5662j_mjj_chi_norm_v25_ct14nlo_cppread_mu_pt12.root"
            print(filename1nu2)
            nlofile2 = TFile.Open(filename1nu2)
            new_hists+=[nlofile2]
            hNloQcd=None
            for k in mass_bins_nlo_list[massbin]:
             #histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])
             histname='qcd_chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-1.0"
             print(histname)
             hnlo = TH1F(nlofile2.Get(histname))
             #hnlo.Scale(float(mass_bins_nlo3[k+1]-mass_bins_nlo3[k]))
             hnlo=hnlo.Rebin(len(chi_binnings[massbin])-1,hnlo.GetName()+"_rebin1",chi_binnings[massbin])
             #hnlo=rebin(hnlo,len(chi_binnings[j])-1,chi_binnings[j])
             if hNloQcd:
                hNloQcd.Add(hnlo)
             else:
                hNloQcd=hnlo
            hNloQcd=smoothChi(hNloQcd) # SMOOTH NNLO PREDICTION (FIX ME)
            hNloQcd.SetLineColor(5)
            hNloQcd.SetLineStyle(3)
            hNloQcd.SetLineWidth(2)

            hNloQcdNoEWK=hNloQcd.Clone(hNloQcd.GetName()+"noEWK")
            hNloQcdNoEWK.SetLineColor(8)
            hNloQcdNoEWK.SetLineStyle(4)
            hNloQcdNoEWK.SetLineWidth(2)
            if normalize:
              hNloQcdNoEWK.Scale(1./hNloQcdNoEWK.Integral())
            else:
              hNloQcdNoEWK.Scale(lumi)
            for b in range(hNloQcdNoEWK.GetXaxis().GetNbins()):
                hNloQcdNoEWK.SetBinContent(b+1,hNloQcdNoEWK.GetBinContent(b+1)/hNloQcdNoEWK.GetBinWidth(b+1))
    
    
            filename="fastnlo/RunII/DijetAngularCMS13_ewk.root"
            print(filename)
            fEWK = TFile.Open(filename)
            new_hists+=[fEWK]
            histname='chi-'+str(massbins[massbin]).strip("()").replace(',',"-").replace(' ',"").replace("1200-1500","1900-2400").replace("1500-1900","1900-2400").replace("6000-7000","6000-6600").replace("6000-13000","6000-6600").replace("7000-13000","6600-13000")
            print(histname)
            hEWK=fEWK.Get(histname)
            print("EWK hist: ")
            print(hEWK)
            for b in range(hNloQcd.GetXaxis().GetNbins()):
                low_bin=hEWK.FindBin(hNloQcd.GetXaxis().GetBinLowEdge(b+1))
                up_bin=hEWK.FindBin(hNloQcd.GetXaxis().GetBinUpEdge(b+1))
                correction=hEWK.Integral(low_bin,up_bin-1)/(up_bin-low_bin)
                print("correction: ")
                print(correction)
                hNloQcd.SetBinContent(b+1,hNloQcd.GetBinContent(b+1)*correction)
            if normalize:
              hNloQcd.Scale(1./hNloQcd.Integral())
            else:
              hNloQcd.Scale(lumi)
            for b in range(hNloQcd.GetXaxis().GetNbins()):
                hNloQcd.SetBinContent(b+1,hNloQcd.GetBinContent(b+1)/hNloQcd.GetBinWidth(b+1))

            if compareMu30:
              filename1nu30="fastnlo/NNLO/2jet.NNLO.fnl5662j_mjj_chi_norm_v25_"+pdfset+"_cppread_mu30_m2.root"
              print(filename1nu30)
              nlofile30 = TFile.Open(filename1nu30)
              new_hists+=[nlofile30]
              hnloQcd30=None
              hnloQcd30up=None
              hnloQcd30down=None
              for k in mass_bins_nlo_list[massbin]:
               histname='qcd_chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-30-central"
               print(histname)
               hnlo30 = TH1F(nlofile30.Get(histname))
               hnlo30=hnlo30.Rebin(len(chi_binnings[massbin])-1,hnlo30.GetName()+"_rebin1",chi_binnings[massbin])
               if hnloQcd30:
                  hnloQcd30.Add(hnlo30)
               else:
                  hnloQcd30=hnlo30
               histname='qcd_chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-30-up"
               print(histname)
               hnlo30up = TH1F(nlofile30.Get(histname))
               hnlo30up=hnlo30up.Rebin(len(chi_binnings[massbin])-1,hnlo30up.GetName()+"_rebin1",chi_binnings[massbin])
               if hnloQcd30up:
                  hnloQcd30up.Add(hnlo30up)
               else:
                  hnloQcd30up=hnlo30up
               histname='qcd_chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-30-down"
               print(histname)
               hnlo30down = TH1F(nlofile30.Get(histname))
               hnlo30down=hnlo30down.Rebin(len(chi_binnings[massbin])-1,hnlo30down.GetName()+"_rebin1",chi_binnings[massbin])
               if hnloQcd30down:
                  hnloQcd30down.Add(hnlo30down)
               else:
                  hnloQcd30down=hnlo30down
              g30=TGraphAsymmErrors()
              for b in range(hnloQcd30.GetNbinsX()):
                g30.SetPoint(g30.GetN(),hNloQcd.GetXaxis().GetBinCenter(b+1),hNloQcd.GetBinContent(b+1))
                print(hnloQcd30.GetBinContent(b+1),hnloQcd30down.GetBinContent(b+1),hnloQcd30up.GetBinContent(b+1))
                g30.SetPointError(g30.GetN()-1,0,0,(hnloQcd30.GetBinContent(b+1)-hnloQcd30down.GetBinContent(b+1)),(hnloQcd30up.GetBinContent(b+1)-hnloQcd30.GetBinContent(b+1)))
              g30.SetLineColor(6)
              g30.SetLineStyle(1)
              g30.SetLineWidth(2)

            if oldTheory:
             filename1nu2="fastnlo/NNLO/2jet.NNLO.fnl5662j_mjj_chi_ct14nlo_cppread_mu_pt12.root"
             #filename1nu2="fastnlo/NNLO/2jet.NNLO.fnl5662j_mjj_chi_norm_v25_ct14nlo_cppread_mu_pt12.root"
             print(filename1nu2)
             nlofile2 = TFile.Open(filename1nu2)
             new_hists+=[nlofile2]
             hNloQcdOld=None
             for k in mass_bins_nlo_list[massbin]:
              #histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])
              histname='qcd_chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-1.0"
              print(histname)
              hnlo = TH1F(nlofile2.Get(histname))
              #hnlo.Scale(float(mass_bins_nlo3[k+1]-mass_bins_nlo3[k]))
              hnlo=hnlo.Rebin(len(chi_binnings[massbin])-1,hnlo.GetName()+"_rebin1",chi_binnings[massbin])
              #hnlo=rebin(hnlo,len(chi_binnings[j])-1,chi_binnings[j])
              if hNloQcdOld:
                 hNloQcdOld.Add(hnlo)
              else:
                 hNloQcdOld=hnlo
             hNloQcdOld=smoothChi(hNloQcdOld) # SMOOTH NNLO PREDICTION (FIX ME)
             hNloQcdOld.SetLineColor(2)
             hNloQcdOld.SetLineStyle(4)
             hNloQcdOld.SetLineWidth(2)
    
             filename="fastnlo/RunII/DijetAngularCMS13_ewk.root"
             print(filename)
             fEWK = TFile.Open(filename)
             new_hists+=[fEWK]
             histname='chi-'+str(massbins[massbin]).strip("()").replace(',',"-").replace(' ',"").replace("1200-1500","1900-2400").replace("1500-1900","1900-2400").replace("6000-7000","6000-6600").replace("6000-13000","6000-6600").replace("7000-13000","6600-13000")
             print(histname)
             hEWK=fEWK.Get(histname)
             print("EWK hist: ")
             print(hEWK)
             for b in range(hNloQcdOld.GetXaxis().GetNbins()):
                 low_bin=hEWK.FindBin(hNloQcdOld.GetXaxis().GetBinLowEdge(b+1))
                 up_bin=hEWK.FindBin(hNloQcdOld.GetXaxis().GetBinUpEdge(b+1))
                 correction=hEWK.Integral(low_bin,up_bin-1)/(up_bin-low_bin)
                 print("correction: ")
                 print(correction)
                 hNloQcdOld.SetBinContent(b+1,hNloQcdOld.GetBinContent(b+1)*correction)
             if normalize:
               hNloQcdOld.Scale(1./hNloQcdOld.Integral())
             else:
               hNloQcdOld.Scale(lumi)
             for b in range(hNloQcdOld.GetXaxis().GetNbins()):
                hNloQcdOld.SetBinContent(b+1,hNloQcdOld.GetBinContent(b+1)/hNloQcdOld.GetBinWidth(b+1))

            if compareScales:
             if useNNLO:
               #filename1nu2="fastnlo/NNLO/2jet.NNLO.fnl5662j_mjj_chi_ct14nnlo_cppread_mu_"+muAltScale+".root"
               filename1nu2="fastnlo/NNLO/2jet.NNLO.fnl5662j_mjj_chi_norm_v25_"+pdfset+"_cppread_mu_"+muAltScale+".root"
             else:
               filename1nu2="fastnlo/NNLO/2jet.NNLO.fnl5662j_mjj_chi_ct14nlo_cppread_mu_pt12.root"
               #filename1nu2="fastnlo/NNLO/2jet.NNLO.fnl5662j_mjj_chi_norm_v25_ct14nlo_cppread_mu_pt12.root"
             print(filename1nu2)
             nlofile2 = TFile.Open(filename1nu2)
             new_hists+=[nlofile2]
             hNloQcdAlt=None
             for k in mass_bins_nlo_list[massbin]:
              #histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])
              histname='qcd_chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-1.0"
              print(histname)
              hnlo = TH1F(nlofile2.Get(histname))
              #hnlo.Scale(float(mass_bins_nlo3[k+1]-mass_bins_nlo3[k]))
              hnlo=hnlo.Rebin(len(chi_binnings[massbin])-1,hnlo.GetName()+"_rebin1",chi_binnings[massbin])
              #hnlo=rebin(hnlo,len(chi_binnings[j])-1,chi_binnings[j])
              if hNloQcdAlt:
                 hNloQcdAlt.Add(hnlo)
              else:
                 hNloQcdAlt=hnlo
             hNloQcdAlt=smoothChi(hNloQcdAlt) # SMOOTH NNLO PREDICTION (FIX ME)
             hNloQcdAlt.SetLineColor(4)
             hNloQcdAlt.SetLineStyle(2)
             hNloQcdAlt.SetLineWidth(2)
             
             filename="fastnlo/RunII/DijetAngularCMS13_ewk.root"
             print(filename)
             fEWK = TFile.Open(filename)
             new_hists+=[fEWK]
             histname='chi-'+str(massbins[massbin]).strip("()").replace(',',"-").replace(' ',"").replace("1200-1500","1900-2400").replace("1500-1900","1900-2400").replace("6000-7000","6000-6600").replace("6000-13000","6000-6600").replace("7000-13000","6600-13000")
             print(histname)
             hEWK=fEWK.Get(histname)
             print("EWK hist: ")
             print(hEWK)
             for b in range(hNloQcdAlt.GetXaxis().GetNbins()):
                 low_bin=hEWK.FindBin(hNloQcdAlt.GetXaxis().GetBinLowEdge(b+1))
                 up_bin=hEWK.FindBin(hNloQcdAlt.GetXaxis().GetBinUpEdge(b+1))
                 correction=hEWK.Integral(low_bin,up_bin-1)/(up_bin-low_bin)
                 print("correction: ")
                 print(correction)
                 hNloQcdAlt.SetBinContent(b+1,hNloQcdAlt.GetBinContent(b+1)*correction)
             if normalize:
               hNloQcdAlt.Scale(1./hNloQcdAlt.Integral())
             else:
               hNloQcdAlt.Scale(lumi)
             for b in range(hNloQcdAlt.GetXaxis().GetNbins()):
                 hNloQcdAlt.SetBinContent(b+1,hNloQcdAlt.GetBinContent(b+1)/hNloQcdAlt.GetBinWidth(b+1))

            if comparePDFs:
             filename1nu2="fastnlo/NNLO/2jet.NNLO.fnl5662j_mjj_chi_norm_v25_ct14nnlo_cppread_mu_"+muScale+".root"
             print(filename1nu2)
             nlofile3 = TFile.Open(filename1nu2)
             new_hists+=[nlofile3]
             hNloQcdPdf=None
             for k in mass_bins_nlo_list[massbin]:
              #histname='chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])
              histname='qcd_chi-'+str(mass_bins_nlo3[k])+"-"+str(mass_bins_nlo3[k+1])+"scale-1.0-1.0"
              print(histname)
              hnlo = TH1F(nlofile3.Get(histname))
              #hnlo.Scale(float(mass_bins_nlo3[k+1]-mass_bins_nlo3[k]))
              hnlo=hnlo.Rebin(len(chi_binnings[massbin])-1,hnlo.GetName()+"_rebin1",chi_binnings[massbin])
              #hnlo=rebin(hnlo,len(chi_binnings[j])-1,chi_binnings[j])
              if hNloQcdPdf:
                 hNloQcdPdf.Add(hnlo)
              else:
                 hNloQcdPdf=hnlo
             hNloQcdPdf=smoothChi(hNloQcdPdf) # SMOOTH NNLO PREDICTION (FIX ME)
             hNloQcdPdf.SetLineColor(6)
             hNloQcdPdf.SetLineStyle(5)
             hNloQcdPdf.SetLineWidth(2)
    
             filename="fastnlo/RunII/DijetAngularCMS13_ewk.root"
             print(filename)
             fEWK = TFile.Open(filename)
             new_hists+=[fEWK]
             histname='chi-'+str(massbins[massbin]).strip("()").replace(',',"-").replace(' ',"").replace("1200-1500","1900-2400").replace("1500-1900","1900-2400").replace("6000-7000","6000-6600").replace("6000-13000","6000-6600").replace("7000-13000","6600-13000")
             print(histname)
             hEWK=fEWK.Get(histname)
             print("EWK hist: ")
             print(hEWK)
             for b in range(hNloQcdPdf.GetXaxis().GetNbins()):
                 low_bin=hEWK.FindBin(hNloQcdPdf.GetXaxis().GetBinLowEdge(b+1))
                 up_bin=hEWK.FindBin(hNloQcdPdf.GetXaxis().GetBinUpEdge(b+1))
                 correction=hEWK.Integral(low_bin,up_bin-1)/(up_bin-low_bin)
                 print("correction: ")
                 print(correction)
                 hNloQcdPdf.SetBinContent(b+1,hNloQcdPdf.GetBinContent(b+1)*correction)
             if normalize:
               hNloQcdPdf.Scale(1./hNloQcdPdf.Integral())
             else:
               hNloQcdPdf.Scale(lumi)
             for b in range(hNloQcdPdf.GetXaxis().GetNbins()):
                 hNloQcdPdf.SetBinContent(b+1,hNloQcdPdf.GetBinContent(b+1)/hNloQcdPdf.GetBinWidth(b+1))

        else:
            filename=fdir+'datacard_shapelimit13TeV_GEN-QCD-run2_chi.root'
            print(filename)
            f = TFile.Open(filename)
            new_hists+=[f]
            histname='QCD_ALT#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
            print(filename)
            f = TFile.Open(filename)
            new_hists+=[f]
            print(histname)
            hNloQcd=f.Get(histname)
            hNloQcd.SetLineColor(5)
            hNloQcd.SetLineStyle(5)
            hNloQcd.SetLineWidth(2)
            if normalize:
              hNloQcd.Scale(1./hNloQcd.Integral())
            else:
              hNloQcd.Scale(lumi)
            for b in range(hNloQcd.GetNbinsX()):
                hNloQcd.SetBinContent(b+1,hNloQcd.GetBinContent(b+1)/hNloQcd.GetBinWidth(b+1))

            
        hNloQcdbackup=hNloQcd.Clone(hNloQcd.GetName()+"backup")    

        if oldMeasurements:
         measurements=["13 TeV, 35.9/fb",
                       "13 TeV, 2.6/fb",
                       #"8 TeV, 19.7/fb",
        	       #"7 TeV, 2.2/fb",
        	       #"7 TeV, 36/pb",
        	       ]
         filenames=["hepdata/HEPData-ins1663452-v1-root.root",
                    "hepdata/HEPData-ins1519995-v2-root.root",
                    #"hepdata/HEPData-ins1327224-v1-root.root",
                    #"hepdata/HEPData-ins1090423-v1-root.root",
                    #"hepdata/HEPData-ins889175-v1-root.root",
        	    ]
         bins={}
         bins[0]=[6000,5400,4800,4200,3600,3000,2400]
         bins[1]=[4800,4200,3600,3000,2400,1900]
         bins[2]=[4200,3600,3000,2400,1900]
         bins[3]=[3000,2400,1900,1500,1200,1000,800,600,400]
         #bins[4]=[2200,1800,1400,1100,850,650,500,350,250]
         hDatas=[]
         for filename in filenames:	             
          print(filename)
          fData2 = TFile.Open(filename)
          new_hists+=[fData2]
          fnum=filenames.index(filename)
          if massbins[massbin][0] in bins[fnum]:
           index=bins[fnum].index(massbins[massbin][0])
           histname="Table "+str(index+1)+"/Graph1D_y1"
           print(histname)
           hData2=fData2.Get(histname)
           new_hists+=[hData2]
           hData2.SetMarkerStyle(fnum+24)
           hData2.SetMarkerSize(0.8)
           hData2.SetMarkerColor(fnum+6)
           hData2.SetLineColor(fnum+6)
           hData2.SetLineWidth(1)
          else:
           hData2=TGraphAsymmErrors() 
          hDatas+=[hData2]

        if compareRun3:
         measurements=["13.6 TeV, 2023, 26.6/fb",
                       "13.6 TeV, 2022, 35.2/fb",
                       ]
         filenames=["data/datacard_shapelimit13TeV_run2_2023_chi.root",
                    "data/datacard_shapelimit13TeV_run2_2022_chi.root",
                    ]
         bins={}
         bins[0]=[7000,6000,5400,4800,4200,3600,3000]
         bins[1]=[7000,6000,5400,4800,4200,3600,3000]
         hDatas=[]
         for filename in filenames:	             
          print(filename)
          fData2 = TFile.Open(filename)
          new_hists+=[fData2]
          fnum=filenames.index(filename)
          if massbins[massbin][0] in bins[fnum]:
           histname=filename.replace("data/","").replace("_chi.root","")+"#chi"+massbintext+"_rebin1"
           print(histname)
           hData2=fData2.Get(histname)
           hData2.SetMarkerStyle(fnum+24)
           hData2.SetMarkerSize(0.8)
           hData2.SetMarkerColor(fnum+6)
           hData2.SetLineColor(fnum+6)
           hData2.SetLineWidth(1)
           if normalize:
             hData2.Scale(1./hData2.Integral())
           hData2h=rebin2(hData2,len(chi_binnings[massbin])-1,chi_binnings[massbin])
           hData2=TGraphAsymmErrors(hData2h.Clone(histname+"G"))
           hDatas+=[hData2]

        if compareMadgraph:
         for mass1,mass2 in massbins:	             
          filename="plots/chi_model_plots_chi_"+str(mass1)+"_"+str(mass2)+"_UL_run2_GEN"+("_nonorm" if not normalize else "")+".root"
          print(filename)
          fMC2 = TFile.Open(filename)
          new_hists+=[fMC2]
          cv=fMC2.Get("cv")
          h=cv.GetListOfPrimitives()[0].GetListOfPrimitives()[1]
          print(h.GetName())
          hMC2=h.Clone(h.GetName()+"Clone")
          hMC2.SetMarkerStyle(24)
          hMC2.SetMarkerSize(0.8)
          hMC2.SetMarkerColor(6)
          hMC2.SetLineColor(6)
          hMC2.SetLineWidth(1)
          if normalize:
            hMC2.Scale(1./hMC2.Integral())
          print(hMC2.Integral())
          #hMC2=rebin2(hMC2,len(chi_binnings[massbin])-1,chi_binnings[massbin])
          #hMC2=TGraphAsymmErrors(hMC2h.Clone(histname+"G"))

        #filename="datacard_shapelimit13TeV_QCD_run2_chi.root"
        #print filename
        #f = TFile.Open(filename)
        #new_hists+=[f]
        #histname='QCD#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
        #print histname
        #h0=f.Get(histname)
        h0=hNloQcd
        h0.SetLineColor(1)
        h0.SetLineStyle(3)
        h0.SetLineWidth(4)

        if unfoldedData:
            #filename="datacards/Unfolded_chiNtuple_dataReReco_v3_Coarse_PFHT900_fromCB_AK4SF_pythia8_Pt_170toInf_MatrixInvert.root"
            #histname="dijet_mass_"+massbintext.replace("1200_1500","2400_3000").replace("1500_1900","2400_3000").replace("1900_2400","2400_3000").replace("6000_7000","5400_6000").replace("7000_13000","6000_13000")+"_chi_unfolded"
            #filename=fdir+'datacard_shapelimit13TeV_GEN-QCD-run2_chi.root'
            #histname="data_obs#chi"+massbintext+"_rebin1" # TODO UNFOLD
            filename='datacards/datacard_shapelimit13TeV_unfold_withUncertainties_run2.root'
            histname="QCD_ALT#chi"+massbintext.replace("1200_1500","2400_3000").replace("1500_1900","2400_3000").replace("1900_2400","2400_3000")+"_rebin1_nosmearpostfit"
            filenamestat='datacards/datacard_shapelimit13TeV_unfold_run2.root'
        else:
            filename=fdir+'datacard_shapelimit13TeV_GEN-QCD-run2_chi.root'
            histname="data_obs#chi"+massbintext+"_rebin1"
        print(filename)
        fData = TFile.Open(filename)
        new_hists+=[fData]
        print(histname)
        h14=fData.Get(histname)
        if unfoldedData:
          print(filenamestat)
          fDataStat = TFile.Open(filenamestat)
          new_hists+=[fDataStat]
          h14stat=fDataStat.Get(histname)
          #h14stat=rebin2(h14stat,len(chi_binnings[massbin])-1,chi_binnings[massbin])
          h14stat=h14stat.Rebin(len(chi_binnings[massbin])-1,h14stat.GetName()+"_rebin",chi_binnings[massbin])
          if normalize:
            h14stat.Scale(1./h14stat.Integral())
          for b in range(h14stat.GetXaxis().GetNbins()):
               h14stat.SetBinContent(b+1,h14stat.GetBinContent(b+1)/h14stat.GetBinWidth(b+1))
               h14stat.SetBinError(b+1,h14stat.GetBinError(b+1)/h14stat.GetBinWidth(b+1))
        
        #if not unfoldedData:
        #  for b in range(h14.GetXaxis().GetNbins()):
        #    h14.SetBinContent(b+1,h14.GetBinContent(b+1)*h14.GetBinWidth(b+1))
        #    h14.SetBinError(b+1,h14.GetBinError(b+1)*h14.GetBinWidth(b+1))
        #origh14=h14.Rebin(len(chi_binnings[massbin])-1,h14.GetName()+"rebinorig",chi_binnings[massbin])
        #h14=rebin2(h14,len(chi_binnings[massbin])-1,chi_binnings[massbin])
        h14=h14.Rebin(len(chi_binnings[massbin])-1,h14.GetName()+"_rebin",chi_binnings[massbin])
        if normalize:
          h14.Scale(1./h14.Integral())
        for b in range(h14.GetXaxis().GetNbins()):
             h14.SetBinContent(b+1,h14.GetBinContent(b+1)/h14.GetBinWidth(b+1))
             h14.SetBinError(b+1,h14.GetBinError(b+1)/h14.GetBinWidth(b+1))

        #if unfoldedData:
        #  hQCDraw=fData.Get("QCD_ALT#chi"+massbintext+"_rebin1_nosmear") # TEMPORARY BIN-BY-BIN UNFOLDING
        #  hQCDsmear=fData.Get("QCD_ALT#chi"+massbintext+"_rebin1")
        #  hQCDraw=hQCDraw.Rebin(len(chi_binnings[massbin])-1,hQCDraw.GetName()+"rebin",chi_binnings[massbin])
        #  hQCDsmear=hQCDsmear.Rebin(len(chi_binnings[massbin])-1,hQCDsmear.GetName()+"rebin",chi_binnings[massbin])
        #  for b in range(h14.GetXaxis().GetNbins()):
        #    print hQCDsmear.GetBinContent(b+1)/hQCDraw.GetBinContent(b+1)
        #    h14.SetBinContent(b+1,h14.GetBinContent(b+1)/hQCDsmear.GetBinContent(b+1)*hQCDraw.GetBinContent(b+1))
        #    h14.SetBinError(b+1,h14.GetBinError(b+1)/hQCDsmear.GetBinContent(b+1)*hQCDraw.GetBinContent(b+1))
        #    origh14.SetBinContent(b+1,origh14.GetBinContent(b+1)/hQCDsmear.GetBinContent(b+1)*hQCDraw.GetBinContent(b+1))
        #    origh14.SetBinError(b+1,origh14.GetBinError(b+1)/hQCDsmear.GetBinContent(b+1)*hQCDraw.GetBinContent(b+1))
        #  h14.Scale(hQCDsmear.Integral()/hQCDraw.Integral())
        #  origh14.Scale(hQCDsmear.Integral()/hQCDraw.Integral())
      
        h14G=TGraphAsymmErrors(h14.Clone(histname+"G"))
        new_hists+=[h14G]
        setupAsymErrors(h14G)
        alpha=1.-0.6827
        nevents=0
        for b in range(h14G.GetN()):
            if h14.GetBinContent(b+1)>0:
                 N=1./pow(h14.GetBinError(b+1)/h14.GetBinContent(b+1),2)
            else:
                 N=0
            #print N
            nevents+=N
            L=0
            if N>0:
                L=ROOT.Math.gamma_quantile(alpha/2.,N,1.)
            U=ROOT.Math.gamma_quantile_c(alpha/2.,N+1,1.)
            if N>0 and not unfoldedData:
              h14G.SetPointEYlow(b,(N-L)/N*h14.GetBinContent(b+1))
              h14G.SetPointEYhigh(b,(U-N)/N*h14.GetBinContent(b+1))
            #print N, sqrt(N)/N, h14.GetBinError(b+1)/h14.GetBinContent(b+1), (N-L)/N, (U-N)/N
        print("data events:", nevents)

        h14Gsys=h14G.Clone(histname+"sys")
        new_hists+=[h14Gsys]
        h14Gsysstat=h14G.Clone(histname+"sysstat")
        new_hists+=[h14Gsysstat]

        filename=fdir+'datacard_shapelimit13TeV_GEN-QCD-run2_chi.root'
        print(filename)
        fsys = TFile.Open(filename)
        new_hists+=[fsys]
        uncertaintynames=["jes","jer","JERtail","prefire","model","sim","scale","pdf"]
        if useNNLO and not compareMu30 and not compareScales: uncertaintynames+=["stat","scaleAlt"]
        uncertainties=[]
        for u in uncertaintynames:
            histname1='QCD_ALT#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_"+u+"Up"
            print(histname1)
            histname2='QCD_ALT#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_"+u+"Down"
            print(histname2)
            histname3='QCD_ALT#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
            print(histname3)
            up=fsys.Get(histname1)
            down=fsys.Get(histname2)
            central=fsys.Get(histname3)
            uncertainties+=[[up,down,central]]
        h2new=h14.Clone("down"+str(massbins[massbin]))
        h3new=h14.Clone("up"+str(massbins[massbin]))
        chi2=0
        for b in range(h14.GetXaxis().GetNbins()):
            if b==0:# or b==h14G.GetXaxis().GetNbins()-1:
                print(massbins[massbin],b,("total" if unfoldedData else "stat"),h14G.GetErrorYlow(b)/h14G.GetY()[b],h14G.GetErrorYhigh(b)/h14G.GetY()[b])
            exp_sumdown=0
            exp_sumup=0
            theory_sumdown=0
            theory_sumup=0
            for up,down,central in uncertainties:
                if b==0:# or b==h14G.GetXaxis().GetNbins()-1:
                    print(massbins[massbin],b,uncertaintynames[uncertainties.index([up,down,central])],abs(up.GetBinContent(b+1)-central.GetBinContent(b+1))/central.GetBinContent(b+1),abs(down.GetBinContent(b+1)-central.GetBinContent(b+1))/central.GetBinContent(b+1))
                addup=pow(max(0,up.GetBinContent(b+1)-central.GetBinContent(b+1),down.GetBinContent(b+1)-central.GetBinContent(b+1)),2)/pow(central.GetBinContent(b+1),2)
                adddown=pow(max(0,central.GetBinContent(b+1)-up.GetBinContent(b+1),central.GetBinContent(b+1)-down.GetBinContent(b+1)),2)/pow(central.GetBinContent(b+1),2)
                if uncertaintynames[uncertainties.index([up,down,central])] in ["jes","jer","JERtail","prefire","model","sim"]:
                    exp_sumup+=addup
                    exp_sumdown+=adddown
                else:
                    theory_sumup+=addup
                    theory_sumdown+=adddown
            exp_sumdown=sqrt(exp_sumdown)
            exp_sumup=sqrt(exp_sumup)
            theory_sumdown=sqrt(theory_sumdown)
            theory_sumup=sqrt(theory_sumup)

            if h14G.GetY()[b]>0:
              chi2+=pow(hNloQcdbackup.GetBinContent(b+1)-h14G.GetY()[b],2) / \
                   (pow(max(exp_sumdown,exp_sumup)*hNloQcdbackup.GetBinContent(b+1),2)+pow(max(theory_sumdown,theory_sumup)*hNloQcdbackup.GetBinContent(b+1),2)+pow(max(h14G.GetErrorYlow(b),h14G.GetErrorYhigh(b)),2))
              print(chi2)

            h14Gsys.SetPointEXlow(b,0)
            h14Gsys.SetPointEXhigh(b,0)
            h14Gsysstat.SetPointEXlow(b,0)
            h14Gsysstat.SetPointEXhigh(b,0)
            if unfoldedData: # unfolded distribution already has total error
              h14Gsys.SetPointEYlow(b,sqrt(max(0,pow(h14G.GetErrorYlow(b),2)-pow(h14stat.GetBinError(b+1),2))))
              h14Gsys.SetPointEYhigh(b,sqrt(max(0,pow(h14G.GetErrorYhigh(b),2)-pow(h14stat.GetBinError(b+1),2))))
              h14Gsysstat.SetPointEYlow(b,h14G.GetErrorYlow(b))
              h14Gsysstat.SetPointEYhigh(b,h14G.GetErrorYhigh(b))
              h2new.SetBinContent(b+1,hNloQcdbackup.GetBinContent(b+1)-theory_sumdown*hNloQcdbackup.GetBinContent(b+1))
              h3new.SetBinContent(b+1,hNloQcdbackup.GetBinContent(b+1)+theory_sumup*hNloQcdbackup.GetBinContent(b+1))
              stat_up=h14stat.GetBinError(b+1)
              stat_down=h14stat.GetBinError(b+1)
            else:
              ### draw exp. sys. on data
              #h14Gsys.SetPointEYlow(b,exp_sumdown*h14G.GetY()[b])
              #h14Gsys.SetPointEYhigh(b,exp_sumup*h14G.GetY()[b])
              #h14Gsysstat.SetPointEYlow(b,sqrt(pow(exp_sumdown*h14G.GetY()[b],2)+pow(h14G.GetErrorYlow(b),2)))
              #h14Gsysstat.SetPointEYhigh(b,sqrt(pow(exp_sumup*h14G.GetY()[b],2)+pow(h14G.GetErrorYhigh(b),2)))
              #h2new.SetBinContent(b+1,hNloQcdbackup.GetBinContent(b+1)-theory_sumdown*hNloQcdbackup.GetBinContent(b+1))
              #h3new.SetBinContent(b+1,hNloQcdbackup.GetBinContent(b+1)+theory_sumup*hNloQcdbackup.GetBinContent(b+1))
              ### draw exp. sys. on background
              h14Gsys.SetPointEYlow(b,0)
              h14Gsys.SetPointEYhigh(b,0)
              h14Gsysstat.SetPointEYlow(b,h14G.GetErrorYlow(b))
              h14Gsysstat.SetPointEYhigh(b,h14G.GetErrorYhigh(b))
              h2new.SetBinContent(b+1,hNloQcdbackup.GetBinContent(b+1)-sqrt(pow(exp_sumdown,2)+pow(theory_sumdown,2))*hNloQcdbackup.GetBinContent(b+1))
              h3new.SetBinContent(b+1,hNloQcdbackup.GetBinContent(b+1)+sqrt(pow(exp_sumup,2)+pow(theory_sumup,2))*hNloQcdbackup.GetBinContent(b+1))
              stat_up=h14G.GetErrorYhigh(b)
              stat_down=h14G.GetErrorYlow(b)
            h14G.SetPointEYlow(b,0)
            h14G.SetPointEYhigh(b,0)
            print("{0:.1f} TO {1:.1f}; {2:.4f} +{3:.4f},-{4:.4f} (DSYS=+{5:.4f},-{6:.4f})".format(h2new.GetXaxis().GetBinLowEdge(b+1),h2new.GetXaxis().GetBinUpEdge(b+1),h14G.GetY()[b],sqrt(pow(stat_up,2)),sqrt(pow(stat_down,2)),sqrt(pow(exp_sumup*h14G.GetY()[b],2)),sqrt(pow(exp_sumdown*h14G.GetY()[b],2))))
        print("chi2/ndof",chi2/hNloQcdbackup.GetXaxis().GetNbins())
        new_hists+=[h2new]
        new_hists+=[h3new]
        h2new.SetLineStyle(1)
        h3new.SetLineStyle(1)
        bandcolor=(kGray if unfoldedData else kCyan)
        h2new.SetLineColor(bandcolor)
        h2new.SetFillColor(10)
        h3new.SetLineColor(bandcolor)
        h3new.SetFillColor(bandcolor)

        #LOAD POSTFIT QCD
        if comparePostfit:
                # QCD systematics
                pname='QCD'
                histnameprefix="DMAxial_Dijet_LO_Mphi_7000_4000_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_1"
                filenameprefix=fdir+prefix+"_"+histnameprefix
                fpostfit=TFile.Open(filenameprefix+"-run2_chi.root")
                masstext=str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")
                basename=histnameprefix+"_ALT#chi"+masstext+"_rebin1"
                uncertainties=getUncertainties(fpostfit,basename,pname,masstext)
                # Shift theory predictions according to fitted nuisance parameters
                fitFile=TFile("limitsDetLHCa1108_DMAxial_Dijet_LO_Mphi_v6b/fitDiagnostics"+histnameprefix+".root")
                bfitParameters,bfitConstraints=getFitResults(fitFile, 'fit_b')
                hbPostfit=hNloQcd.Clone("QCDPostfit"+str(massbins[massbin]))
                applyFitResults(bfitParameters,bfitConstraints,uncertainties,hbPostfit,h14)
                hbPostfit.SetLineWidth(2)
                hbPostfit.SetLineStyle(2)

        if massbin>=0:
         if True: #FIX
          filename=fdir+"datacard_shapelimit13TeV_cs_"+pdfset+"_30000_LL+-run2_chi.root"
          histname='cs_'+pdfset+'_30000_LL+#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print(filename)
          f = TFile.Open(filename)
          new_hists+=[f]
          print(histname)
          hci=f.Get(histname)
          hci.SetLineColor(kRed)
          hci.SetLineWidth(3)
          hci.SetLineStyle(8)
          if normalize:
            hci.Scale(1./hci.Integral())
          else:
            hci.Scale(lumi)
          for b in range(hci.GetNbinsX()):
               hci.SetBinContent(b+1,hci.GetBinContent(b+1)/hci.GetBinWidth(b+1))

          filename=fdir+"datacard_shapelimit13TeV_cs_"+pdfset+"_30000_LL--run2_chi.root"
          histname='cs_'+pdfset+'_30000_LL-#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print(filename)
          f = TFile.Open(filename)
          new_hists+=[f]
          print(histname)
          hcib=f.Get(histname)
          hcib.SetLineColor(kAzure+3)
          hcib.SetLineStyle(6)
          hcib.SetLineWidth(3)
          if normalize:
            hcib.Scale(1./hcib.Integral())
          else:
            hcib.Scale(lumi)
          for b in range(hcib.GetNbinsX()):
               hcib.SetBinContent(b+1,hcib.GetBinContent(b+1)/hcib.GetBinWidth(b+1))

          filename=fdir+"datacard_shapelimit13TeV_cs_"+pdfset+"_30000_VV+-run2_chi.root"
          histname='cs_'+pdfset+'_30000_VV+#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print(filename)
          f = TFile.Open(filename)
          new_hists+=[f]
          print(histname)
          hcic=f.Get(histname)
          hcic.SetLineColor(36)
          hcic.SetLineStyle(5)
          hcic.SetLineWidth(3)
          if normalize:
            hcic.Scale(1./hcic.Integral())
          else:
            hcic.Scale(lumi)
          for b in range(hcic.GetNbinsX()):
               hcic.SetBinContent(b+1,hcic.GetBinContent(b+1)/hcic.GetBinWidth(b+1))

          filename=fdir+"datacard_shapelimit13TeV_cs_"+pdfset+"_30000_VV--run2_chi.root"
          histname='cs_'+pdfset+'_30000_VV-#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print(filename)
          f = TFile.Open(filename)
          new_hists+=[f]
          print(histname)
          hcid=f.Get(histname)
          hcid.SetLineColor(28)
          hcid.SetLineStyle(5)
          hcid.SetLineWidth(3)
          if normalize:
            hcid.Scale(1./hcid.Integral())
          else:
            hcid.Scale(lumi)
          for b in range(hcid.GetNbinsX()):
               hcid.SetBinContent(b+1,hcid.GetBinContent(b+1)/hcid.GetBinWidth(b+1))

          filename=fdir+"datacard_shapelimit13TeV_cs_"+pdfset+"_30000_AA+-run2_chi.root"
          histname='cs_'+pdfset+'_30000_AA+#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print(filename)
          f = TFile.Open(filename)
          new_hists+=[f]
          print(histname)
          hcie=f.Get(histname)
          hcie.SetLineColor(8)
          hcie.SetLineStyle(5)
          hcie.SetLineWidth(3)
          if normalize:
            hcie.Scale(1./hcie.Integral())
          else:
            hcie.Scale(lumi)
          for b in range(hcie.GetNbinsX()):
               hcie.SetBinContent(b+1,hcie.GetBinContent(b+1)/hcie.GetBinWidth(b+1))

          filename=fdir+"datacard_shapelimit13TeV_cs_"+pdfset+"_30000_AA--run2_chi.root"
          histname='cs_'+pdfset+'_30000_AA-#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print(filename)
          f = TFile.Open(filename)
          new_hists+=[f]
          print(histname)
          hcif=f.Get(histname)
          hcif.SetLineColor(9)
          hcif.SetLineStyle(5)
          hcif.SetLineWidth(3)
          if normalize:
            hcif.Scale(1./hcif.Integral())
          else:
            hcif.Scale(lumi)
          for b in range(hcif.GetNbinsX()):
               hcif.SetBinContent(b+1,hcif.GetBinContent(b+1)/hcif.GetBinWidth(b+1))

          filename=fdir+"datacard_shapelimit13TeV_cs_"+pdfset+"_30000_V-A+-run2_chi.root"
          histname='cs_'+pdfset+'_30000_V-A+#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print(filename)
          f = TFile.Open(filename)
          new_hists+=[f]
          print(histname)
          hcig=f.Get(histname)
          hcig.SetLineColor(kAzure+3)
          hcig.SetLineStyle(5)
          hcig.SetLineWidth(3)
          if normalize:
            hcig.Scale(1./hcig.Integral())
          else:
            hcig.Scale(lumi)
          for b in range(hcig.GetNbinsX()):
               hcig.SetBinContent(b+1,hcig.GetBinContent(b+1)/hcig.GetBinWidth(b+1))

          filename=fdir+"datacard_shapelimit13TeV_cs_"+pdfset+"_30000_V-A--run2_chi.root"
          histname='cs_'+pdfset+'_30000_V-A-#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print(filename)
          f = TFile.Open(filename)
          new_hists+=[f]
          print(histname)
          hcih=f.Get(histname)
          hcih.SetLineColor(kAzure+3)
          hcih.SetLineStyle(5)
          hcih.SetLineWidth(3)
          if normalize:
            hcih.Scale(1./hcih.Integral())
          else:
            hcih.Scale(lumi)
          for b in range(hcih.GetNbinsX()):
               hcih.SetBinContent(b+1,hcih.GetBinContent(b+1)/hcih.GetBinWidth(b+1))
         if True:
          filename=fdir+"datacard_shapelimit13TeV_GENnp-22-CP2-run2_chi.root"
          print(filename)
          f = TFile.Open(filename)
          new_hists+=[f]
          histname='QCDADD13000#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print(histname)
          hgrw=f.Get(histname)
          hgrw.SetLineColor(kAzure+7)
          hgrw.SetLineStyle(7)
          hgrw.SetLineWidth(4)
          if normalize:
            hgrw.Scale(1./hgrw.Integral())
          else:
            hgrw.Scale(lumi)
          for b in range(hgrw.GetNbinsX()):
               hgrw.SetBinContent(b+1,hgrw.GetBinContent(b+1)/hgrw.GetBinWidth(b+1))
            
          coupling="0p5"
          filename=fdir+'datacard_shapelimit13TeV_DMAxial_Dijet_LO_Mphi_2000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_'+coupling+'-run2_chi.root'
          histname='DMAxial_Dijet_LO_Mphi_2000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_'+coupling+'#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print(filename)
          f7 = TFile.Open(filename)
          new_hists+=[f7]
          print(histname)
          hdm=f7.Get(histname)
          setUpDMHists(hdm,kBlue,4,2)

          filename=fdir+'datacard_shapelimit13TeV_DMAxial_Dijet_LO_Mphi_3000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_'+coupling+'-run2_chi.root'
          histname='DMAxial_Dijet_LO_Mphi_3000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_'+coupling+'#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print(filename)
          f7a = TFile.Open(filename)
          new_hists+=[f7a]
          print(histname)
          hdma=f7a.Get(histname)
          setUpDMHists(hdma,kMagenta,7,2)

          filename=fdir+'datacard_shapelimit13TeV_DMAxial_Dijet_LO_Mphi_4000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_'+coupling+'-run2_chi.root'
          histname='DMAxial_Dijet_LO_Mphi_4000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_'+coupling+'#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print(filename)
          f7b = TFile.Open(filename)
          new_hists+=[f7b]
          print(histname)
          hdmb=f7b.Get(histname)
          setUpDMHists(hdmb,kTeal+4,10,2)

          filename=fdir+'datacard_shapelimit13TeV_DMAxial_Dijet_LO_Mphi_5000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_'+coupling+'-run2_chi.root'
          histname='DMAxial_Dijet_LO_Mphi_5000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_'+coupling+'#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print(filename)
          f7c = TFile.Open(filename)
          new_hists+=[f7c]
          print(histname)
          hdmc=f7c.Get(histname)
          setUpDMHists(hdmc,kOrange+7,1,2)

          filename=fdir+'datacard_shapelimit13TeV_DMAxial_Dijet_LO_Mphi_6000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_'+coupling+'-run2_chi.root'
          histname='DMAxial_Dijet_LO_Mphi_6000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_'+coupling+'#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print(filename)
          f7d = TFile.Open(filename)
          new_hists+=[f7d]
          print(histname)
          hdmd=f7d.Get(histname)
          setUpDMHists(hdmd,kAzure+8,1,2)

        if massbin>3:
         if True:
          filename=fdir+"datacard_shapelimit13TeV_GENnp-CP5-run2-QBH_ADD6_9000_chi.root"
          print(filename)
          f = TFile.Open(filename)
          new_hists+=[f]
          histname='QBH_ADD6_9000#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print(histname)
          hqbh=f.Get(histname)
          hqbh=hqbh.Rebin(len(chi_binnings[massbin])-1,hqbh.GetName()+"_rebin",chi_binnings[massbin])
          hqbh.SetLineColor(kGreen+3)
          hqbh.SetLineStyle(5)
          hqbh.SetLineWidth(4)
          if normalize:
            hqbh.Scale(1./hqbh.Integral())
          else:
            hqbh.Scale(lumi)
          for b in range(hqbh.GetNbinsX()):
              hqbh.SetBinContent(b+1,hqbh.GetBinContent(b+1)/hqbh.GetBinWidth(b+1))

        if True:
          filename=fdir+"datacard_shapelimit13TeV_alp_QCD_fa2500-run2_chi.root"
          histname='alp_QCD_fa2500#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print(filename)
          f = TFile.Open(filename)
          new_hists+=[f]
          print(histname)
          halp=f.Get(histname)
          halp.SetLineColor(kBlue)
          halp.SetLineWidth(3)
          halp.SetLineStyle(9)
          if normalize:
            halp.Scale(1./halp.Integral())
          else:
            halp.Scale(lumi)
          for b in range(halp.GetNbinsX()):
               halp.SetBinContent(b+1,halp.GetBinContent(b+1)/halp.GetBinWidth(b+1))

          filename=fdir+"datacard_shapelimit13TeV_tripleG_QCD_CG0p01-run2_chi.root"
          histname='tripleG_QCD_CG0p01#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print(filename)
          f = TFile.Open(filename)
          new_hists+=[f]
          print(histname)
          htripleG=f.Get(histname)
          htripleG.SetLineColor(kRed)
          htripleG.SetLineWidth(3)
          htripleG.SetLineStyle(10)
          if normalize:
            htripleG.Scale(1./htripleG.Integral())
          else:
            htripleG.Scale(lumi)
          for b in range(htripleG.GetNbinsX()):
               htripleG.SetBinContent(b+1,htripleG.GetBinContent(b+1)/htripleG.GetBinWidth(b+1))

#        if massbin>4:
#            h0.GetYaxis().SetRangeUser(0.02,0.22)
#        elif massbin>3:
#            h0.GetYaxis().SetRangeUser(0.045,0.14)
#        elif massbin>1:
#            h0.GetYaxis().SetRangeUser(0.045,0.12)
#        else:
#            h0.GetYaxis().SetRangeUser(0.055,0.105)

        if massbin>9:
            h0.GetYaxis().SetRangeUser(0.0,0.8*h14.Integral())
        elif massbin>8:
            h0.GetYaxis().SetRangeUser(0.03*h14.Integral(),0.26*h14.Integral())
        elif massbin>7:
            h0.GetYaxis().SetRangeUser(0.04*h14.Integral(),0.18*h14.Integral())
        elif massbin>6:
            h0.GetYaxis().SetRangeUser(0.055*h14.Integral(),0.15*h14.Integral())
        elif massbin>4:
            h0.GetYaxis().SetRangeUser(0.065*h14.Integral(),0.12*h14.Integral())
        else:
            h0.GetYaxis().SetRangeUser(0.065*h14.Integral(),0.12*h14.Integral())
            
        h0.GetXaxis().SetTitle("#chi_{dijet}")
        if unfoldedData:
          h0.GetYaxis().SetTitle("1/#sigma_{dijet} d#sigma_{dijet}/d#chi_{dijet}")
        else:
          h0.GetYaxis().SetTitle("1/N_{dijet} dN_{dijet}/d#chi_{dijet}")
        h0.GetYaxis().SetTitleOffset(1.4)
        h0.GetXaxis().SetTitleOffset(0.8)
        h0.GetYaxis().SetTitleSize(0.0625)
        h0.GetYaxis().SetLabelSize(0.05)
        h0.GetXaxis().SetTitleSize(0.0625)
        h0.GetXaxis().SetLabelSize(0.05)
        h0.GetXaxis().SetTickLength(0.03)
        h0.GetYaxis().SetNdivisions(505)
        
        c.cd()

        pad1=TPad("","",0,0,1,0.3)
        pad1.SetTopMargin(0.08)
        pad1.SetBottomMargin(0.30)
        pad1.Draw()
        pad1.cd()
        
        h0Div=h0.Clone(h0.GetName()+"_ratio")
        h3newDiv=h3new.Clone(h3new.GetName()+"_ratio")
        h2newDiv=h2new.Clone(h2new.GetName()+"_ratio")
        if oldTheory:
          hNloQcdOldDiv=hNloQcdOld.Clone(hNloQcdOld.GetName()+"_ratio")
        if compareScales:
          hNloQcdAltDiv=hNloQcdAlt.Clone(hNloQcdAlt.GetName()+"_ratio")
        if comparePDFs:
          hNloQcdPdfDiv=hNloQcdPdf.Clone(hNloQcdAlt.GetName()+"_ratio")

        h14GDiv=h14G.Clone(h14G.GetName()+"_ratio")
        h14GsysDiv=h14Gsys.Clone(h14Gsys.GetName()+"_ratio")
        h14GsysstatDiv=h14Gsysstat.Clone(h14Gsysstat.GetName()+"_ratio")

        if compareMadgraph:
          hMC2Div=hMC2.Clone(hMC2.GetName()+"_ratio")
        if unfoldedData:
          hNloQcdNoEWKDiv=hNloQcdNoEWK.Clone(hNloQcdNoEWK.GetName()+"_ratio")

        hDataDivs=[]
        h0Div.Divide(h0)
        h3newDiv.Divide(h0)
        h2newDiv.Divide(h0)
        if oldTheory:
          hNloQcdOldDiv.Divide(h0)
        if compareMu30:
          g30Div=divideAsymErrors(g30,h0,False)
          g30Div.SetMarkerStyle(g30.GetMarkerStyle())
          g30Div.SetMarkerSize(g30.GetMarkerSize())
          g30Div.SetMarkerColor(g30.GetMarkerColor())
          g30Div.SetLineColor(g30.GetLineColor())
          g30Div.SetLineWidth(g30.GetLineWidth())
          hDataDivs+=[g30Div]
        if compareScales:
          hNloQcdAltDiv.Divide(h0)
        if comparePDFs:
          hNloQcdPdfDiv.Divide(h0)
        if oldMeasurements or compareRun3:
         for hData2 in hDatas:
          hData2Div=divideAsymErrors(hData2,h0,True)
          hData2Div.SetMarkerStyle(hData2.GetMarkerStyle())
          hData2Div.SetMarkerSize(hData2.GetMarkerSize())
          hData2Div.SetMarkerColor(hData2.GetMarkerColor())
          hData2Div.SetLineColor(hData2.GetLineColor())
          hData2Div.SetLineWidth(hData2.GetLineWidth())
          hDataDivs+=[hData2Div]
        if compareMadgraph:
          hMC2Div.Divide(h0)
        if compareNoEWK:
          hNloQcdNoEWKDiv.Divide(h0)

        #print h14G.GetN(),h14G.Eval(1),h14G.Eval(1.5), h14.GetBinContent(1)
        #print h14G.GetX()[0], h14G.GetY()[0]
        h14GDiv=divideAsymErrors(h14G,h0,True)
        setupAsymErrors(h14GDiv)
        h14GsysDiv=divideAsymErrors(h14Gsys,h0,False)
        setupAsymErrors(h14GsysDiv)
        h14GsysstatDiv=divideAsymErrors(h14Gsysstat,h0,True)
        setupAsymErrors(h14GsysstatDiv)

#        if massbin>4:
#            h0Div.GetYaxis().SetRangeUser(0.,2.)
#        elif massbin>3:
#            h0Div.GetYaxis().SetRangeUser(0.4,1.6)
#        elif massbin>1:
#            h0Div.GetYaxis().SetRangeUser(0.7,1.3)
#        else:
#            h0Div.GetYaxis().SetRangeUser(0.8,1.2)

        if massbin>9:
            h0Div.GetYaxis().SetRangeUser(0.,3.5)
        elif massbin>8:
            h0Div.GetYaxis().SetRangeUser(0.,2.)
        elif massbin>6:
            h0Div.GetYaxis().SetRangeUser(0.7,1.3)
        else:
            h0Div.GetYaxis().SetRangeUser(0.85,1.15)
        h0Div.GetYaxis().SetTitle("#frac{Data}{NNLO QCD+EW}")
        h0Div.GetXaxis().SetTitle("#chi_{dijet}")
        h0Div.GetYaxis().SetTitleOffset(0.65)
        h0Div.GetYaxis().SetTitleSize(0.125)
        h0Div.GetYaxis().SetLabelSize(0.118)
        h0Div.GetXaxis().SetTitleSize(0.15)
        h0Div.GetXaxis().SetLabelSize(0.118)
        h0Div.GetYaxis().SetTickLength(0.04)
        h0Div.GetXaxis().SetTickLength(0.08)
        
        h0Div.Draw("axis")
        h3newDiv.Draw("histsame")
        h2newDiv.Draw("histsame")
        if oldTheory:
         hNloQcdOldDiv.Draw("histsame")
        if compareMu30:
         g30Div.Draw("esame")
        if compareScales:
         hNloQcdAltDiv.Draw("histsame")
        if comparePDFs:
         hNloQcdPdfDiv.Draw("histsame")
        if oldMeasurements or compareRun3:
         for hData2Div in hDataDivs:
          hData2Div.Draw("pzesame")
        if compareMadgraph:
         hMC2Div.Draw("histsame")
        if compareNoEWK:
         hNloQcdNoEWKDiv.Draw("histsame")
        h0Div.Draw("histsame")
        #h14GDiv.Draw("pzesame")
        h14GsysDiv.Draw("||same")
        h14GsysstatDiv.Draw("zesame")
        h0Div.Draw("axissame")
        h14GDiv.SetMarkerSize(0.8)
        h14GDiv.Draw("pzesame")

        c.cd()

        pad2=TPad("","",0,0.3,1,1)
        pad2.SetBottomMargin(0.001)
        pad2.SetTopMargin(0.07)
        pad2.Draw()
        pad2.cd()
        
        if massbin==len(massbins)-1:
            h0.Draw("axis")
        else:
            h0.Draw("axissame")
        h3new.Draw("histsame")
        h2new.Draw("histsame")
        if oldTheory:
         hNloQcdOld.Draw("histsame")
        if compareMu30:
         g30.Draw("esame")
        if compareScales:
         hNloQcdAlt.Draw("histsame")
        if comparePDFs:
         hNloQcdPdf.Draw("histsame")
        if oldMeasurements or compareRun3:
         for hData2 in hDatas:
          hData2.Draw("pzesame")
        if compareMadgraph:
          hMC2.Draw("histsame")
        if compareNoEWK:
          hNloQcdNoEWK.Draw("histsame")
        h0.Draw("histsame")
        if comparePostfit:
          hbPostfit.Draw("histsame")
        if massbin>=5 and signalsBSM:
            #hci.Draw("histsame")
            hcib.Draw("histsame")
            hcic.Draw("histsame")
            #hcid.Draw("histsame")
            #hcie.Draw("histsame")
            #hcif.Draw("histsame")
            #hcig.Draw("histsame")
            #hcih.Draw("histsame")
        if massbin>5 and signalsBSM:
            hgrw.Draw("histsame")
        if massbin>8 and signalsBSM:
            hqbh.Draw("histsame")
        if massbin == 3 and signalsDM:
            #hdm.Draw("histsame")
            #hdma.Draw("histsame")
            hdmb.Draw("histsame")
            hdmc.Draw("histsame")
        if massbin == 4 and signalsDM:
            #hdma.Draw("histsame")
            hdmb.Draw("histsame")
            hdmc.Draw("histsame")
            hdmd.Draw("histsame")
        if massbin == 5 and signalsDM:
            hdmb.Draw("histsame")
            hdmc.Draw("histsame")
            hdmd.Draw("histsame")
        if massbin == 6 and signalsDM:
            hdmc.Draw("histsame")
            hdmd.Draw("histsame")
        if massbin == 7 and signalsDM:
            hdmc.Draw("histsame")
            hdmd.Draw("histsame")
        if massbin == 8 and signalsDM:
            hdmc.Draw("histsame")
            hdmd.Draw("histsame")
        if massbin == 9 and signalsDM:
            hdmd.Draw("histsame")
        if massbin>3 and signalsBSM:
            halp.Draw("histsame")
            htripleG.Draw("histsame")
        #h14G.Draw("pzesame")
        h14Gsys.Draw("||same")
        h14Gsysstat.Draw("zesame")
        h14G.SetMarkerSize(0.8)
        h14G.Draw("pzesame")
        h0.Draw("axissame")
        #h14G.Print()
        #h14Gsys.Print()
        #h14Gsysstat.Print()

        if massbin>=9:
            ylabel=0.35
        else:
            ylabel=0.10
        
        if massbin==0: title="1.2 < #font[52]{M}_{#font[42]{jj}} < 1.5 TeV"
        if massbin==1: title="1.5 < #font[52]{M}_{#font[42]{jj}} < 1.9 TeV"
        if massbin==2: title="1.9 < #font[52]{M}_{#font[42]{jj}} < 2.4 TeV"
        if massbin==3: title="2.4 < #font[52]{M}_{#font[42]{jj}} < 3.0 TeV"
        if massbin==4: title="3.0 < #font[52]{M}_{#font[42]{jj}} < 3.6 TeV"
        if massbin==5: title="3.6 < #font[52]{M}_{#font[42]{jj}} < 4.2 TeV"
        if massbin==6: title="4.2 < #font[52]{M}_{#font[42]{jj}} < 4.8 TeV"
        if massbin==7: title="4.8 < #font[52]{M}_{#font[42]{jj}} < 5.4 TeV"
        if massbin==8: title="5.4 < #font[52]{M}_{#font[42]{jj}} < 6.0 TeV"
        if massbin==9: title="6.0 < #font[52]{M}_{#font[42]{jj}} < 7.0 TeV"
        if massbin==10: title=" #font[52]{M}_{#font[42]{jj}} > 7.0 TeV"

        #title+=" TeV"
        #if offsets[massbin]==0: titleo=""
        #elif offsets[massbin]<0: titleo="("+str(offsets[massbin])+")"
        #else: titleo="(+"+str(offsets[massbin])+")"

        l=TLegend(0.61,ylabel,1.0,ylabel+0.005,title)
        l.SetTextSize(0.044)
        l.SetFillStyle(0)
        l.Draw("same")
        new_hists+=[l]

        #lo=TLegend(0.82,ylabel,1.4,ylabel+0.01,titleo)
        #lo.SetTextSize(0.033)
        #lo.SetFillStyle(0)
        #lo.Draw("same")
        #new_hists+=[lo]

    l5=TLegend(0.74,0.92,1.0,0.92,"CMS")
    l5.SetTextSize(0.044)
    l5.SetFillStyle(0)
    #l5.Draw("same")
    new_hists+=[l5]
     
    l5=TLegend(0.7,0.88,1.0,0.88,"#sqrt{s} = 13 TeV")
    l5.SetTextSize(0.044)
    l5.SetFillStyle(0)
    #l5.Draw("same")
    new_hists+=[l5]
     
    l=TLegend(0.7,0.82,1.0,0.82,"L = 138 fb^{-1}")
    l.SetTextSize(0.044)
    l.SetFillStyle(0)
    #l.Draw("same")
    new_hists+=[l]
    
    banner=TLatex(0.33,0.96,"CMS,   L = 138 fb^{-1},   #sqrt{s} = 13 TeV")
    banner.SetNDC()
    banner.SetTextSize(0.044)
    #banner.Draw()

    h3newnew=h3new.Clone()
    h3newnew.SetLineColor(1)
    h3newnew.SetLineStyle(3)
    h3newnew.SetLineWidth(4)

    if massbin>=7:
        l2=TLegend(0.25,0.45,0.5,0.91,"")
    elif massbin>=4:
        l2=TLegend(0.25,0.55,0.5,0.91,"")
    else:
        l2=TLegend(0.25,0.55,0.5,0.91,"")
    l2.SetTextSize(0.044)
    l2.SetMargin(0.33)
    if oldMeasurements or compareRun3:
      l2.AddEntry(h14G,"13 TeV 138/fb","ple")
      if oldMeasurements or compareRun3:
       for hData2 in hDatas:
        if hData2.GetN()>0:
         l2.AddEntry(hData2,measurements[hDatas.index(hData2)],"ple")
    else:
      l2.AddEntry(h14G,"Data","ple")
    if compareScales:
      if compareMu30:
        l2.AddEntry(h3newnew,"NNLO QCD+EW (#mu=m_{jj}, 6)","fl")
        l2.AddEntry(g30,"NNLO QCD+EW (#mu=m_{jj}, 30)","le")
      else:
        l2.AddEntry(h3newnew,"NNLO QCD+EW (#mu=m_{jj}, NNPDF3.1)","fl")
      l2.AddEntry(hNloQcdAlt,"NNLO QCD+EW (#mu=<p_{T}>, NNPDF3.1)","fl")
    else:
      if useNNLO:
        l2.AddEntry(h3newnew,"NNLO QCD+EW","fl")
      else:
        l2.AddEntry(h3newnew,"NLO QCD+EW","fl")
    if comparePDFs:
      l2.AddEntry(hNloQcdPdf,"NNLO QCD+EW (#mu=m_{jj}, CT14)","fl")
    if oldTheory:
      l2.AddEntry(hNloQcdOld,"NLO QCD+EW (#mu=<p_{T}>, CT14)","fl")
    if compareMadgraph:
      l2.AddEntry(hMC2,"LO+MLM QCD (Madgraph)","l")
    if compareNoEWK:
      l2.AddEntry(hNloQcdNoEWK,"NNLO QCD (#mu=m_{jj}, NNPDF3.1)","l")
    if comparePostfit:
      l2.AddEntry(hbPostfit,"Postfit QCD","l")
    if massbin>=5 and signalsBSM:
        #l2.AddEntry(hci,"#Lambda_{LL}^{#font[122]{+}} (CI) = 30 TeV","l")
        l2.AddEntry(hcib,"#Lambda_{LL}^{#font[122]{-}} (CI) = 30 TeV","l")
        l2.AddEntry(hcic,"#Lambda_{VV}^{#font[122]{+}} (CI) = 30 TeV","l")
        #l2.AddEntry(hcid,"#Lambda_{VV}^{#font[122]{-}} (CI) = 30 TeV","l")
        #l2.AddEntry(hcie,"#Lambda_{AA}^{#font[122]{+}} (CI) = 30 TeV","l")
        #l2.AddEntry(hcif,"#Lambda_{AA}^{#font[122]{-}} (CI) = 30 TeV","l")
        #l2.AddEntry(hcig,"#Lambda_{V-A} (CI) = 30 TeV","l")
        #l2.AddEntry(hcih,"#Lambda_{V-A}^{#font[122]{-}} (CI) = 30 TeV","l")
    if massbin > 5 and signalsBSM:
        l2.AddEntry(hgrw,"#Lambda_{T} (GRW) = 13 TeV","l")
    if massbin > 6 and signalsBSM:
        l2.AddEntry(hqbh,"M_{QBH} (n_{ED} = 6 ADD) = 9 TeV","l")
    if massbin == 3 and signalsDM:
        #l2.AddEntry(hdm,"M_{Med} = 2 TeV (g_{q} = 0.5)","l")
        #l2.AddEntry(hdma,"M_{Med} = 3 TeV (g_{q} = 0.5)","l")
        l2.AddEntry(hdmb,"M_{Med} = 4 TeV (g_{q} = 0.5)","l")
        l2.AddEntry(hdmc,"M_{Med} = 5 TeV (g_{q} = 0.5)","l")
    if massbin == 4 and signalsDM:
        #l2.AddEntry(hdma,"M_{Med} = 3 TeV (g_{q} = 0.5)","l")
        l2.AddEntry(hdmb,"M_{Med} = 4 TeV (g_{q} = 0.5)","l")
        l2.AddEntry(hdmc,"M_{Med} = 5 TeV (g_{q} = 0.5)","l")
        l2.AddEntry(hdmd,"M_{Med} = 6 TeV (g_{q} = 0.5)","l")
    if massbin == 5 and signalsDM:
        l2.AddEntry(hdmb,"M_{Med} = 4 TeV (g_{q} = 0.5)","l")
        l2.AddEntry(hdmc,"M_{Med} = 5 TeV (g_{q} = 0.5)","l")
        l2.AddEntry(hdmd,"M_{Med} = 6 TeV (g_{q} = 0.5)","l")
    if massbin == 6 and signalsDM:
        l2.AddEntry(hdmc,"M_{Med} = 5 TeV (g_{q} = 0.5)","l")
        l2.AddEntry(hdmd,"M_{Med} = 6 TeV (g_{q} = 0.5)","l")
    if massbin == 7 and signalsDM:
        l2.AddEntry(hdmc,"M_{Med} = 5 TeV (g_{q} = 0.5)","l")
        l2.AddEntry(hdmd,"M_{Med} = 6 TeV (g_{q} = 0.5)","l")
    if massbin == 8 and signalsDM:
        l2.AddEntry(hdmc,"M_{Med} = 5 TeV (g_{q} = 0.5)","l")
        l2.AddEntry(hdmd,"M_{Med} = 6 TeV (g_{q} = 0.5)","l")
    if massbin == 9 and signalsDM:
        l2.AddEntry(hdmd,"M_{Med} = 6 TeV (g_{q} = 0.5)","l")
    if massbin > 3 and signalsBSM:
        l2.AddEntry(halp,"#font[52]{f}_{a} / #font[52]{c}_{g} (ALP) = 2.5 TeV","l")
        l2.AddEntry(htripleG,"#Lambda / #sqrt{#font[52]{C}_{G}} (SMEFT) = 10 TeV","l")
    l2.SetFillStyle(0)
    l2.Draw("same")

#    l2b=TLegend(0.3,0.54,0.6,0.93,"")
#    l2b.SetTextSize(0.035)
#    l2b.AddEntry(h14G," ","")
#    l2b.AddEntry(h0," ","l")
#    #l2b.AddEntry(hqbh," ","")
#    l2b.AddEntry(hci," ","")
#    l2b.AddEntry(hcib," ","")
#    l2b.AddEntry(hcic," ","")
#    l2b.AddEntry(hcid," ","")
#    #l2b.AddEntry(hcie," ","")
#    #l2b.AddEntry(hcif," ","")
#    l2b.AddEntry(hcig," ","")
#    #l2b.AddEntry(hcih," ","")
#    l2b.AddEntry(hgrw," ","")
#    l2b.SetFillStyle(0)
#    #l2b.Draw("same")
    
    #tl1=TLine(16,0,16,0.24)
    #tl1.SetLineColor(1)
    #tl1.SetLineStyle(1)
    #tl1.SetLineWidth(1)
    #tl1.Draw("same")

    #// writing the lumi information and the CMS "logo"
    CMS_lumi.CMS_lumi( c, iPeriod, iPos );

    postfix=""
    if not unfoldedData:
      postfix="_detector"
    if oldMeasurements:
      postfix+="_compare"
    if compareScales:
      postfix+="_scales"
    if comparePDFs:
      postfix+="_pdfs"
    if signalsDM:
      postfix+="_dm"
    if compareRun3:
      postfix+="_run3"
    if compareMu30:
      postfix+="_mu30"
    if compareMadgraph:
      postfix+="_madgraph"
    if not normalize:
      postfix+="_nonorm"
    c.SaveAs(fdir+prefix + "_combined_theory"+str(massbin)+postfix+"_run2.pdf")
    #c.SaveAs(fdir+prefix + "_combined_theory"+str(massbin)+postfix+"_run2.eps")
