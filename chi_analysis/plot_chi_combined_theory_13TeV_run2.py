from ROOT import *
import ROOT
import array, math
import os,sys
from math import *

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
        print "divideAsymErrors Function Fails!!!"
        sys.exit()

    g=TGraphAsymmErrors()
    for b in range(h1.GetNbinsX()):
        den=h1.GetBinContent(b+1)
        if doEX:
            bwidth=h1.GetBinWidth(b+1)/2
        else:
            bwidth=0
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

  unfoldedData=False

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

  for massbin in range(len(massbins)):
      
    massbintext = str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")
    
    #gROOT.Reset()
    gROOT.SetBatch(True)
    gROOT.SetStyle("Plain")
    gStyle.SetOptStat(0)
    gStyle.SetOptFit(0)
    gStyle.SetTitleOffset(1.2,"Y")
    gStyle.SetPadLeftMargin(0.18)
    gStyle.SetPadBottomMargin(0.11)
    gStyle.SetPadTopMargin(0.055)
    gStyle.SetPadRightMargin(0.05)
    gStyle.SetMarkerSize(1.5)
    gStyle.SetHistLineWidth(1)
    gStyle.SetStatFontSize(0.020)
    gStyle.SetTitleSize(0.06, "XYZ")
    gStyle.SetLabelSize(0.05, "XYZ")
    gStyle.SetNdivisions(510, "XYZ")
    gStyle.SetLegendBorderSize(0)
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetEndErrorSize(8)

    print "load CMS style"
    gROOT.LoadMacro("CMS_lumi.C");
    iPeriod = 5;	#// 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV 
    iPos = 33;
    #// second parameter in example_plot is iPos, which drives the position of the CMS logo in the plot
    #// iPos=11 : top-left, left-aligned
    #// iPos=33 : top-right, right-aligned
    #// iPos=22 : center, centered
    #// mode generally : 
    #//   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)

    prefix="datacard_shapelimit13TeV"

    c = TCanvas("combined", "combined", 0, 0, 650, 800)
    c.Divide(1,1)
    new_hists=[]
    if True:
        if unfoldedData:
            fdir='invertMatrixOct20/'
        else:
            fdir='versions/run2NNLOMar25/'
            
        if unfoldedData:  
            filename="fastnlo/RunII/fnl5662j_v23_fix_CT14nlo_allmu_ak4.root"
            print filename
            fNloQcd = TFile.Open(filename)
            new_hists+=[fNloQcd]
            histname='chi-'+str(massbins[massbin]).strip("()").replace(',',"-").replace(' ',"").replace("6000-13000","6000-6600")
            print histname
            hNloQcd=fNloQcd.Get(histname)
	    print "NLO QCD hist: "
            print fNloQcd
	    hNloQcd=rebin(hNloQcd,len(chi_binnings[massbin])-1,chi_binnings[massbin])
            hNloQcd.SetLineColor(5)
            hNloQcd.SetLineStyle(3)
            hNloQcd.SetLineWidth(2)
    
            filename="fastnlo/RunII/DijetAngularCMS13_ewk.root"
            print filename
            fEWK = TFile.Open(filename)
            new_hists+=[fEWK]
            histname='chi-'+str(massbins[massbin]).strip("()").replace(',',"-").replace(' ',"").replace("6000-13000","6000-6600")
            print histname
            hEWK=fEWK.Get(histname)
	    print "EWK hist: "
            print hEWK
            for b in range(hNloQcd.GetXaxis().GetNbins()):
	        low_bin=hEWK.FindBin(hNloQcd.GetXaxis().GetBinLowEdge(b+1))
	        up_bin=hEWK.FindBin(hNloQcd.GetXaxis().GetBinUpEdge(b+1))
	        correction=hEWK.Integral(low_bin,up_bin-1)/(up_bin-low_bin)
	        print "correction: "
                print correction
                hNloQcd.SetBinContent(b+1,hNloQcd.GetBinContent(b+1)*correction*hNloQcd.GetBinWidth(b+1))
            hNloQcd.Scale(1./hNloQcd.Integral())
            for b in range(hNloQcd.GetXaxis().GetNbins()):
	        hNloQcd.SetBinContent(b+1,hNloQcd.GetBinContent(b+1)/hNloQcd.GetBinWidth(b+1))

        else:
            filename=fdir+'datacard_shapelimit13TeV_GEN-QCD-run2_chi.root'
            print filename
            f = TFile.Open(filename)
            new_hists+=[f]
            histname='QCD_ALT#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
            print filename
            f = TFile.Open(filename)
            new_hists+=[f]
            print histname
            hNloQcd=f.Get(histname)
            hNloQcd.SetLineColor(5)
            hNloQcd.SetLineStyle(3)
            hNloQcd.SetLineWidth(2)
	    hNloQcd.Scale(1./hNloQcd.Integral())
	    for b in range(hNloQcd.GetNbinsX()):
	        hNloQcd.SetBinContent(b+1,hNloQcd.GetBinContent(b+1)/hNloQcd.GetBinWidth(b+1))
            
        hNloQcdbackup=hNloQcd.Clone(hNloQcd.GetName()+"backup")    

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
            filename="datacards/Unfolded_chiNtuple_dataReReco_v3_Coarse_PFHT900_fromCB_AK4SF_pythia8_Pt_170toInf_MatrixInvert.root"
            histname="dijet_mass_"+massbintext+"_chi_unfolded"
	else:
            histname="data_obs#chi"+massbintext+"_rebin1"
        print filename
        fData = TFile.Open(filename)
        new_hists+=[fData]
        print histname
        h14=fData.Get(histname)
	#if not unfoldedData:
        #  for b in range(h14.GetXaxis().GetNbins()):
        #    h14.SetBinContent(b+1,h14.GetBinContent(b+1)*h14.GetBinWidth(b+1))
        #    h14.SetBinError(b+1,h14.GetBinError(b+1)*h14.GetBinWidth(b+1))
	origh14=h14.Rebin(len(chi_binnings[massbin])-1,h14.GetName()+"rebinorig",chi_binnings[massbin])
	h14=rebin2(h14,len(chi_binnings[massbin])-1,chi_binnings[massbin])
	    
	h14G=TGraphAsymmErrors(h14.Clone(histname+"G"))
	new_hists+=[h14G]
        setupAsymErrors(h14G)
	alpha=1.-0.6827
	nevents=0
	for b in range(h14G.GetN()):
	    if unfoldedData:
	      N=origh14.GetBinContent(b+1)
	    else:
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
	    if N>0:
              h14G.SetPointEYlow(b,(N-L)/N*h14.GetBinContent(b+1))
              h14G.SetPointEYhigh(b,(U-N)/N*h14.GetBinContent(b+1))
            #print N, sqrt(N)/N, origh14.GetBinError(b+1)/origh14.GetBinContent(b+1), h14.GetBinError(b+1)/h14.GetBinContent(b+1), (N-L)/N, (U-N)/N
        print "data events:", nevents
	
	h14Gsys=h14G.Clone(histname+"sys")
	new_hists+=[h14Gsys]
	h14Gsysstat=h14G.Clone(histname+"sysstat")
	new_hists+=[h14Gsysstat]

        #filename="datacard_shapelimit13TeV_QCD_chi.root"
        print filename
        fsys = TFile.Open(filename)
        new_hists+=[fsys]
        uncertaintynames=["jer","jes","pdf","scaleAlt","model"]
        uncertainties=[]
        for u in uncertaintynames:
            histname1='QCD_ALT#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_"+u+"Up"
            print histname1
            histname2='QCD_ALT#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1_"+u+"Down"
            print histname2
            histname3='QCD_ALT#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
            print histname3
            up=fsys.Get(histname1)
            down=fsys.Get(histname2)
            central=fsys.Get(histname3)
            uncertainties+=[[up,down,central]]
        h2new=h14.Clone("down"+str(massbins[massbin]))
        h3new=h14.Clone("up"+str(massbins[massbin]))
	chi2=0
        for b in range(h14.GetXaxis().GetNbins()):
            if b==0:# or b==h14G.GetXaxis().GetNbins()-1:
                print massbins[massbin],b,"stat",h14G.GetErrorYlow(b)/h14G.GetY()[b],h14G.GetErrorYhigh(b)/h14G.GetY()[b]
	    exp_sumdown=0
	    exp_sumup=0
            theory_sumdown=0
            theory_sumup=0
            for up,down,central in uncertainties:
	        if b==0:# or b==h14G.GetXaxis().GetNbins()-1:
                    print massbins[massbin],b,uncertaintynames[uncertainties.index([up,down,central])],abs(up.GetBinContent(b+1)-central.GetBinContent(b+1))/central.GetBinContent(b+1),abs(down.GetBinContent(b+1)-central.GetBinContent(b+1))/central.GetBinContent(b+1)
                addup=pow(max(0,up.GetBinContent(b+1)-central.GetBinContent(b+1),down.GetBinContent(b+1)-central.GetBinContent(b+1)),2)/pow(central.GetBinContent(b+1),2)
		adddown=pow(max(0,central.GetBinContent(b+1)-up.GetBinContent(b+1),central.GetBinContent(b+1)-down.GetBinContent(b+1)),2)/pow(central.GetBinContent(b+1),2)
                if uncertaintynames[uncertainties.index([up,down,central])]=="jer" or uncertaintynames[uncertainties.index([up,down,central])]=="jes" or uncertaintynames[uncertainties.index([up,down,central])]=="model":
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
	           (pow(max(exp_sumdown,exp_sumup)*h14G.GetY()[b],2)+pow(max(theory_sumdown,theory_sumup)*h14G.GetY()[b],2)+pow(max(h14G.GetErrorYlow(b),h14G.GetErrorYhigh(b)),2))

            h14Gsys.SetPointEXlow(b,0)
            h14Gsys.SetPointEXhigh(b,0)
            h14Gsys.SetPointEYlow(b,exp_sumdown*h14G.GetY()[b])
            h14Gsys.SetPointEYhigh(b,exp_sumup*h14G.GetY()[b])
            h14Gsysstat.SetPointEXlow(b,0)
            h14Gsysstat.SetPointEXhigh(b,0)
            h14Gsysstat.SetPointEYlow(b,sqrt(pow(exp_sumdown*h14G.GetY()[b],2)+pow(h14G.GetErrorYlow(b),2)))
            h14Gsysstat.SetPointEYhigh(b,sqrt(pow(exp_sumup*h14G.GetY()[b],2)+pow(h14G.GetErrorYhigh(b),2)))
	    stat_up=h14G.GetErrorYhigh(b)
            stat_down=h14G.GetErrorYlow(b)
            h14G.SetPointEYlow(b,0)
            h14G.SetPointEYhigh(b,0)
            h2new.SetBinContent(b+1,hNloQcdbackup.GetBinContent(b+1)-theory_sumdown*hNloQcdbackup.GetBinContent(b+1))
            h3new.SetBinContent(b+1,hNloQcdbackup.GetBinContent(b+1)+theory_sumup*hNloQcdbackup.GetBinContent(b+1))
	    print "{0:.1f} TO {1:.1f}; {2:.4f} +{3:.4f},-{4:.4f} (DSYS=+{5:.4f},-{6:.4f})".format(h2new.GetXaxis().GetBinLowEdge(b+1),h2new.GetXaxis().GetBinUpEdge(b+1),h14G.GetY()[b],sqrt(pow(stat_up,2)),sqrt(pow(stat_down,2)),sqrt(pow(exp_sumup*h14G.GetY()[b],2)),sqrt(pow(exp_sumdown*h14G.GetY()[b],2)))
        print "chi2/ndof",chi2/h14G.GetXaxis().GetNbins()
	new_hists+=[h2new]
        new_hists+=[h3new]
        h2new.SetLineStyle(1)
        h3new.SetLineStyle(1)
        h2new.SetLineColor(kGray)
        h2new.SetFillColor(10)
        h3new.SetLineColor(kGray)
        h3new.SetFillColor(kGray)
        
        if massbin>=0:
         if True: #FIX
	  filename=fdir+"datacard_shapelimit13TeV_cs_ct14nnlo_14000_LL+-run2_chi.root"
          histname='cs_ct14nnlo_14000_LL+#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print filename
          f = TFile.Open(filename)
          new_hists+=[f]
          print histname
          hci=f.Get(histname)
          hci.SetLineColor(kRed)
          hci.SetLineWidth(3)
          hci.SetLineStyle(8)
	  hci.Scale(1./hci.Integral())
	  for b in range(hci.GetNbinsX()):
	       hci.SetBinContent(b+1,hci.GetBinContent(b+1)/hci.GetBinWidth(b+1))

          filename=fdir+"datacard_shapelimit13TeV_cs_ct14nnlo_14000_LL--run2_chi.root"
          histname='cs_ct14nnlo_14000_LL-#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print filename
          f = TFile.Open(filename)
          new_hists+=[f]
          print histname
          hcib=f.Get(histname)
          hcib.SetLineColor(kAzure+3)
          hcib.SetLineStyle(6)
          hcib.SetLineWidth(3)
	  hcib.Scale(1./hcib.Integral())
	  for b in range(hcib.GetNbinsX()):
	       hcib.SetBinContent(b+1,hcib.GetBinContent(b+1)/hcib.GetBinWidth(b+1))

          filename=fdir+"datacard_shapelimit13TeV_cs_ct14nnlo_14000_VV+-run2_chi.root"
          histname='cs_ct14nnlo_14000_VV+#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print filename
          f = TFile.Open(filename)
          new_hists+=[f]
          print histname
          hcic=f.Get(histname)
          hcic.SetLineColor(36)
          hcic.SetLineStyle(5)
          hcic.SetLineWidth(3)
	  hcic.Scale(1./hcic.Integral())
	  for b in range(hcic.GetNbinsX()):
	       hcic.SetBinContent(b+1,hcic.GetBinContent(b+1)/hcic.GetBinWidth(b+1))

          filename=fdir+"datacard_shapelimit13TeV_cs_ct14nnlo_14000_VV--run2_chi.root"
          histname='cs_ct14nnlo_14000_VV-#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print filename
          f = TFile.Open(filename)
          new_hists+=[f]
          print histname
          hcid=f.Get(histname)
          hcid.SetLineColor(28)
          hcid.SetLineStyle(5)
          hcid.SetLineWidth(3)
	  hcid.Scale(1./hcid.Integral())
	  for b in range(hcid.GetNbinsX()):
	       hcid.SetBinContent(b+1,hcid.GetBinContent(b+1)/hcid.GetBinWidth(b+1))

          filename=fdir+"datacard_shapelimit13TeV_cs_ct14nnlo_14000_AA+-run2_chi.root"
          histname='cs_ct14nnlo_14000_AA+#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print filename
          f = TFile.Open(filename)
          new_hists+=[f]
          print histname
          hcie=f.Get(histname)
          hcie.SetLineColor(8)
          hcie.SetLineStyle(5)
          hcie.SetLineWidth(3)
	  hcie.Scale(1./hcie.Integral())
	  for b in range(hcie.GetNbinsX()):
	       hcie.SetBinContent(b+1,hcie.GetBinContent(b+1)/hcie.GetBinWidth(b+1))

          filename=fdir+"datacard_shapelimit13TeV_cs_ct14nnlo_14000_AA--run2_chi.root"
          histname='cs_ct14nnlo_14000_AA-#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print filename
          f = TFile.Open(filename)
          new_hists+=[f]
          print histname
          hcif=f.Get(histname)
          hcif.SetLineColor(9)
          hcif.SetLineStyle(5)
          hcif.SetLineWidth(3)
	  hcif.Scale(1./hcif.Integral())
	  for b in range(hcif.GetNbinsX()):
	       hcif.SetBinContent(b+1,hcif.GetBinContent(b+1)/hcif.GetBinWidth(b+1))

          filename=fdir+"datacard_shapelimit13TeV_cs_ct14nnlo_14000_V-A+-run2_chi.root"
          histname='cs_ct14nnlo_14000_V-A+#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print filename
          f = TFile.Open(filename)
          new_hists+=[f]
          print histname
          hcig=f.Get(histname)
          hcig.SetLineColor(kAzure+3)
          hcig.SetLineStyle(5)
          hcig.SetLineWidth(3)
	  hcig.Scale(1./hcig.Integral())
	  for b in range(hcig.GetNbinsX()):
	       hcig.SetBinContent(b+1,hcig.GetBinContent(b+1)/hcig.GetBinWidth(b+1))

          filename=fdir+"datacard_shapelimit13TeV_cs_ct14nnlo_14000_V-A--run2_chi.root"
          histname='cs_ct14nnlo_14000_V-A-#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print filename
          f = TFile.Open(filename)
          new_hists+=[f]
          print histname
          hcih=f.Get(histname)
          hcih.SetLineColor(kAzure+3)
          hcih.SetLineStyle(5)
          hcih.SetLineWidth(3)
	  hcih.Scale(1./hcih.Integral())
	  for b in range(hcih.GetNbinsX()):
	       hcih.SetBinContent(b+1,hcih.GetBinContent(b+1)/hcih.GetBinWidth(b+1))
         if True:
          filename=fdir+"datacard_shapelimit13TeV_GENnp-24-run2_chi.root"
          print filename
          f = TFile.Open(filename)
          new_hists+=[f]
          histname='QCDADD12000#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print histname
          hgrw=f.Get(histname)
          hgrw.SetLineColor(kAzure+7)
          hgrw.SetLineStyle(7)
          hgrw.SetLineWidth(4)
	  hgrw.Scale(1./hgrw.Integral())
	  for b in range(hgrw.GetNbinsX()):
	       hgrw.SetBinContent(b+1,hgrw.GetBinContent(b+1)/hgrw.GetBinWidth(b+1))
            
          coupling="1"
          filename=fdir+'datacard_shapelimit13TeV_DMAxial_Dijet_LO_Mphi_2000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_'+coupling+'-run2_chi.root'
          histname='DMAxial_Dijet_LO_Mphi_2000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_'+coupling+'#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print filename
          f7 = TFile.Open(filename)
          new_hists+=[f7]
          print histname
          hdm=f7.Get(histname)
          setUpDMHists(hdm,kBlue,4,2)

          filename=fdir+'datacard_shapelimit13TeV_DMAxial_Dijet_LO_Mphi_3000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_'+coupling+'-run2_chi.root'
          histname='DMAxial_Dijet_LO_Mphi_3000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_'+coupling+'#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print filename
          f7a = TFile.Open(filename)
          new_hists+=[f7a]
          print histname
          hdma=f7a.Get(histname)
          setUpDMHists(hdma,kMagenta,7,2)

          filename=fdir+'datacard_shapelimit13TeV_DMAxial_Dijet_LO_Mphi_1000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_'+coupling+'-run2_chi.root'
          histname='DMAxial_Dijet_LO_Mphi_1000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_'+coupling+'#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print filename
          f7b = TFile.Open(filename)
          new_hists+=[f7b]
          print histname
          hdmb=f7b.Get(histname)
          setUpDMHists(hdmb,kTeal+4,10,2)

          filename=fdir+'datacard_shapelimit13TeV_DMAxial_Dijet_LO_Mphi_5000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_'+coupling+'-run2_chi.root'
          histname='DMAxial_Dijet_LO_Mphi_5000_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_'+coupling+'#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print filename
          f7c = TFile.Open(filename)
          new_hists+=[f7c]
          print histname
          hdmc=f7c.Get(histname)
          setUpDMHists(hdmc,kOrange+7,1,2)

          filename=fdir+'datacard_shapelimit13TeV_DMAxial_Dijet_LO_Mphi_1500_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_'+coupling+'-run2_chi.root'
          histname='DMAxial_Dijet_LO_Mphi_1500_1_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_'+coupling+'#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print filename
          f7d = TFile.Open(filename)
          new_hists+=[f7d]
          print histname
          hdmd=f7d.Get(histname)
          setUpDMHists(hdmd,kAzure+8,1,2)

        if massbin>3:
         if False: #FIX
	  if unfoldedData:  
              filename=fdir+"datacard_shapelimit13TeV_QBH_8000_6_chi_v1.root"
          else:
              filename=fdir+"datacard_shapelimit13TeV_QBH_8000_6-run2_chi.root"
          print filename
          f = TFile.Open(filename)
          new_hists+=[f]
          if unfoldedData:
              histname='QCDADD6QBH8000#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          else:
              histname='QBH_8000_6#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print histname
          hqbh=f.Get(histname)
          hqbh=hqbh.Rebin(len(chi_binnings[massbin])-1,hqbh.GetName()+"_rebin",chi_binnings[massbin])
          hqbh.SetLineColor(kGreen+3)
          hqbh.SetLineStyle(5)
          hqbh.SetLineWidth(4)
          hqbh.Scale(1./hqbh.Integral())
          for b in range(hqbh.GetNbinsX()):
              hqbh.SetBinContent(b+1,hqbh.GetBinContent(b+1)/hqbh.GetBinWidth(b+1))

        if True: #FIX
	  filename=fdir+"datacard_shapelimit13TeV_alp_QCD_fa2500-run2_chi.root"
          histname='alp_QCD_fa2500#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print filename
          f = TFile.Open(filename)
          new_hists+=[f]
          print histname
          halp=f.Get(histname)
          halp.SetLineColor(kBlue)
          halp.SetLineWidth(3)
          halp.SetLineStyle(9)
	  halp.Scale(1./halp.Integral())
	  for b in range(halp.GetNbinsX()):
	       halp.SetBinContent(b+1,halp.GetBinContent(b+1)/halp.GetBinWidth(b+1))

	  filename=fdir+"datacard_shapelimit13TeV_tripleG_QCD_CG0p01-run2_chi.root"
          histname='tripleG_QCD_CG0p01#chi'+str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")+"_rebin1"
          print filename
          f = TFile.Open(filename)
          new_hists+=[f]
          print histname
          htripleG=f.Get(histname)
          htripleG.SetLineColor(kRed)
          htripleG.SetLineWidth(3)
          htripleG.SetLineStyle(10)
	  htripleG.Scale(1./htripleG.Integral())
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

        if massbin>7:
            h0.GetYaxis().SetRangeUser(0.02,0.22)
        elif massbin>6:
            h0.GetYaxis().SetRangeUser(0.02,0.22)
        elif massbin>4:
            h0.GetYaxis().SetRangeUser(0.045,0.12)
        else:
            h0.GetYaxis().SetRangeUser(0.045,0.12)
            
        h0.GetXaxis().SetTitle("#chi_{dijet}")
        h0.GetYaxis().SetTitle("1/#sigma_{dijet} d#sigma_{dijet}/d#chi_{dijet}")
        h0.GetYaxis().SetTitleOffset(1.4)
        h0.GetXaxis().SetTitleOffset(0.8)
        h0.GetYaxis().SetTitleSize(0.05)
        h0.GetYaxis().SetLabelSize(0.04)
        h0.GetXaxis().SetTitleSize(0.05)
        h0.GetXaxis().SetLabelSize(0.04)
        h0.GetXaxis().SetTickLength(0.03)
        h0.GetYaxis().SetNdivisions(505)
        
        c.cd()

        pad1=TPad("","",0,0,1,0.3)
        pad1.SetTopMargin(0.08)
        pad1.SetBottomMargin(0.25)
        pad1.Draw()
        pad1.cd()
        
        h0Div=h0.Clone(h0.GetName()+"_ratio")
        h3newDiv=h3new.Clone(h3new.GetName()+"_ratio")
        h2newDiv=h2new.Clone(h2new.GetName()+"_ratio")

        h14GDiv=h14G.Clone(h14G.GetName()+"_ratio")
        h14GsysDiv=h14Gsys.Clone(h14Gsys.GetName()+"_ratio")
        h14GsysstatDiv=h14Gsysstat.Clone(h14Gsysstat.GetName()+"_ratio")

        h0Div.Divide(h0)
        h3newDiv.Divide(h0)
        h2newDiv.Divide(h0)
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

        if massbin>7:
            h0Div.GetYaxis().SetRangeUser(0.,2.)
        elif massbin>6:
            h0Div.GetYaxis().SetRangeUser(0.,2.)
        elif massbin>4:
            h0Div.GetYaxis().SetRangeUser(0.7,1.3)
        else:
            h0Div.GetYaxis().SetRangeUser(0.7,1.3)
        h0Div.GetYaxis().SetTitle("#frac{Data}{NNLO QCD + EW}")
        h0Div.GetXaxis().SetTitle("#chi_{dijet}")
        h0Div.GetYaxis().SetTitleOffset(0.65)
        h0Div.GetYaxis().SetTitleSize(0.1)
        h0Div.GetYaxis().SetLabelSize(0.095)
        h0Div.GetXaxis().SetTitleSize(0.12)
        h0Div.GetXaxis().SetLabelSize(0.095)
        h0Div.GetYaxis().SetTickLength(0.04)
        h0Div.GetXaxis().SetTickLength(0.08)
        
        h0Div.Draw("axis")
        h3newDiv.Draw("histsame")
        h2newDiv.Draw("histsame")
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
        h0.Draw("histsame")
#        if massbin>=5: #FIX
#            hci.Draw("histsame") #FIX
#            hcib.Draw("histsame") #FIX
            #hcic.Draw("histsame")
            #hcid.Draw("histsame")
            #hcie.Draw("histsame")
            #hcif.Draw("histsame")
            #hcig.Draw("histsame")
            #hcih.Draw("histsame")
        if massbin>5:
            hgrw.Draw("histsame")
#        if massbin>6: #FIX
#            hqbh.Draw("histsame") #FIX
        if massbin == 3 and False:
            hdm.Draw("histsame")
            hdma.Draw("histsame")
            #hdmb.Draw("histsame")
            hdmc.Draw("histsame")
            #hdmd.Draw("histsame")
        if massbin == 4 and False:
            #hdm.Draw("histsame")
            hdma.Draw("histsame")
            #hdmb.Draw("histsame")
            hdmc.Draw("histsame")
            #hdmd.Draw("histsame")
        if massbin == 5 and False:
            #hdm.Draw("histsame")
            #hdma.Draw("histsame")
            #hdmb.Draw("histsame")
            hdmc.Draw("histsame")
            #hdmd.Draw("histsame")
        if massbin == 6 and False:
            #hdm.Draw("histsame")
            #hdma.Draw("histsame")
            #hdmb.Draw("histsame")
            hdmc.Draw("histsame")
            #hdmd.Draw("histsame")
        if massbin == 7 and False:
            #hdm.Draw("histsame")
            #hdma.Draw("histsame")
            #hdmb.Draw("histsame")
            hdmc.Draw("histsame")
            #hdmd.Draw("histsame")
        if massbin == 8 and False:
            #hdm.Draw("histsame")
            #hdma.Draw("histsame")
            #hdmb.Draw("histsame")
            hdmc.Draw("histsame")
            #hdmd.Draw("histsame")
        #if massbin == 9:
            #hdm.Draw("histsame")
            #hdma.Draw("histsame")
            #hdmb.Draw("histsame")
            #hdmc.Draw("histsame")
            #hdmd.Draw("histsame")
        if True:
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

        if massbin>=6:
            ylabel=0.35
        else:
            ylabel=0.4
        
        if massbin==0: title="1.2 < #font[72]{M_{jj}} < 1.5 TeV"
        if massbin==1: title="1.5 < #font[72]{M_{jj}} < 1.9 TeV"
        if massbin==2: title="1.9 < #font[72]{M_{jj}} < 2.4 TeV"
        if massbin==3: title="2.4 < #font[72]{M_{jj}} < 3.0 TeV"
        if massbin==4: title="3.0 < #font[72]{M_{jj}} < 3.6 TeV"
        if massbin==5: title="3.6 < #font[72]{M_{jj}} < 4.2 TeV"
        if massbin==6: title="4.2 < #font[72]{M_{jj}} < 4.8 TeV"
        if massbin==7: title="4.8 < #font[72]{M_{jj}} < 5.4 TeV"
        if massbin==8: title="5.4 < #font[72]{M_{jj}} < 6.0 TeV"
        if massbin==9: title="6.0 < #font[72]{M_{jj}} < 7.0 TeV"
        if massbin==10: title=" #font[72]{M_{jj}} > 7.0 TeV"

        #title+=" TeV"
        #if offsets[massbin]==0: titleo=""
        #elif offsets[massbin]<0: titleo="("+str(offsets[massbin])+")"
        #else: titleo="(+"+str(offsets[massbin])+")"

        l=TLegend(0.6,ylabel,1.0,ylabel+0.005,title)
        l.SetTextSize(0.035)
        l.SetFillStyle(0)
        l.Draw("same")
        new_hists+=[l]

        #lo=TLegend(0.82,ylabel,1.4,ylabel+0.01,titleo)
        #lo.SetTextSize(0.033)
        #lo.SetFillStyle(0)
        #lo.Draw("same")
        #new_hists+=[lo]

    l5=TLegend(0.74,0.92,1.0,0.92,"CMS")
    l5.SetTextSize(0.035)
    l5.SetFillStyle(0)
    #l5.Draw("same")
    new_hists+=[l5]
     
    l5=TLegend(0.7,0.88,1.0,0.88,"#sqrt{s} = 13 TeV")
    l5.SetTextSize(0.035)
    l5.SetFillStyle(0)
    #l5.Draw("same")
    new_hists+=[l5]
     
    l=TLegend(0.7,0.82,1.0,0.82,"L = 137 fb^{-1}")
    l.SetTextSize(0.035)
    l.SetFillStyle(0)
    #l.Draw("same")
    new_hists+=[l]
    
    banner=TLatex(0.33,0.96,"CMS,   L = 137 fb^{-1},   #sqrt{s} = 13 TeV")
    banner.SetNDC()
    banner.SetTextSize(0.035)
    #banner.Draw()

    h3newnew=h3new.Clone()
    h3newnew.SetLineColor(1)
    h3newnew.SetLineStyle(3)
    h3newnew.SetLineWidth(4)

    if massbin>=7:
        l2=TLegend(0.3,0.45,0.6,0.91,"")
    elif massbin>=4:
        l2=TLegend(0.3,0.55,0.6,0.91,"")
    else:
        l2=TLegend(0.3,0.55,0.6,0.91,"")
    l2.SetTextSize(0.033)
    l2.SetMargin(0.33)
    l2.AddEntry(h14G,"Data","ple")
    l2.AddEntry(h3newnew,"NNLO QCD + EW","fl")
#    if massbin>=5: #FIX
#        l2.AddEntry(hci,"#Lambda_{LL}^{#font[122]{+}} (CI) = 14 TeV","l") #FIX
#        l2.AddEntry(hcib,"#Lambda_{LL}^{#font[122]{-}} (CI) = 14 TeV","l") #FIX
        #l2.AddEntry(hcic,"#Lambda_{VV}^{#font[122]{+}} (CI) = 14 TeV","l")
        #l2.AddEntry(hcid,"#Lambda_{VV}^{#font[122]{-}} (CI) = 14 TeV","l")
        #l2.AddEntry(hcie,"#Lambda_{AA}^{#font[122]{+}} (CI) = 14 TeV","l")
        #l2.AddEntry(hcif,"#Lambda_{AA}^{#font[122]{-}} (CI) = 14 TeV","l")
        #l2.AddEntry(hcig,"#Lambda_{V-A} (CI) = 14 TeV","l")
        #l2.AddEntry(hcih,"#Lambda_{V-A}^{#font[122]{-}} (CI) = 14 TeV","l")
    if massbin > 5:
        l2.AddEntry(hgrw,"#Lambda_{T} (GRW) = 12 TeV","l")
#    if massbin > 6: #FIX
#        l2.AddEntry(hqbh,"M_{QBH} (n_{ED} = 6 ADD) = 8 TeV","l") #FIX
    if massbin == 3 and False:
        #l2.AddEntry(hdmb,"M_{Med} (DM g_{q} = 1.0) = 1 TeV","l")
        #l2.AddEntry(hdmd,"M_{Med} (DM g_{q} = 1.0) = 1.5 TeV","l")
        l2.AddEntry(hdm,"M_{Med} (DM g_{q} = 1.0) = 2 TeV","l")
        l2.AddEntry(hdma,"M_{Med} (DM g_{q} = 1.0) = 3 TeV","l")
        l2.AddEntry(hdmc,"M_{Med} (DM g_{q} = 1.0) = 5 TeV","l")
    if massbin == 4 and False:
        #l2.AddEntry(hdm,"M_{Med} (DM g_{q} = 1.0) = 2.0 TeV","l")
        l2.AddEntry(hdma,"M_{Med} (DM g_{q} = 1.0) = 3 TeV","l")
        #l2.AddEntry(hdmb,"M_{Med} (DM g_{q} = 1.0) = 4.0 TeV","l")
        l2.AddEntry(hdmc,"M_{Med} (DM g_{q} = 1.0) = 5 TeV","l")
        #l2.AddEntry(hdmd,"M_{Med} (DM g_{q} = 1.0) = 6.0 TeV","l")
    if massbin == 5 and False:
        #l2.AddEntry(hdm,"M_{Med} (DM g_{q} = 1.0) = 2.0 TeV","l")
        #l2.AddEntry(hdma,"M_{Med} (DM g_{q} = 1.0) = 3.0 TeV","l")
        #l2.AddEntry(hdmb,"M_{Med} (DM g_{q} = 1.0) = 4.0 TeV","l")
        l2.AddEntry(hdmc,"M_{Med} (DM g_{q} = 1.0) = 5 TeV","l")
        #l2.AddEntry(hdmd,"M_{Med} (DM g_{q} = 1.0) = 6.0 TeV","l")
    if massbin == 6 and False:
        #l2.AddEntry(hdm,"M_{Med} (DM g_{q} = 1.0) = 2.0 TeV","l")
        #l2.AddEntry(hdma,"M_{Med} (DM g_{q} = 1.0) = 3.0 TeV","l")
        #l2.AddEntry(hdmb,"M_{Med} (DM g_{q} = 1.0) = 4.0 TeV","l")
        l2.AddEntry(hdmc,"M_{Med} (DM g_{q} = 1.0) = 5 TeV","l")
        #l2.AddEntry(hdmd,"M_{Med} (DM g_{q} = 1.0) = 6.0 TeV","l")
    if massbin == 7 and False:
        #l2.AddEntry(hdm,"M_{Med} (DM g_{q} = 1.0) = 2.0 TeV","l")
        #l2.AddEntry(hdma,"M_{Med} (DM g_{q} = 1.0) = 3.0 TeV","l")
        #l2.AddEntry(hdmb,"M_{Med} (DM g_{q} = 1.0) = 4.0 TeV","l")
        l2.AddEntry(hdmc,"M_{Med} (DM g_{q} = 1.0) = 5 TeV","l")
        #l2.AddEntry(hdmd,"M_{Med} (DM g_{q} = 1.0) = 6.0 TeV","l")
    if massbin == 8 and False:
        #l2.AddEntry(hdm,"M_{Med} (DM g_{q} = 1.0) = 2.0 TeV","l")
        #l2.AddEntry(hdma,"M_{Med} (DM g_{q} = 1.0) = 3.0 TeV","l")
        #l2.AddEntry(hdmb,"M_{Med} (DM g_{q} = 1.0) = 4.0 TeV","l")
        l2.AddEntry(hdmc,"M_{Med} (DM g_{q} = 1.0) = 5 TeV","l")
        #l2.AddEntry(hdmd,"M_{Med} (DM g_{q} = 1.0) = 6.0 TeV","l")
    #if massbin == 9:
        #l2.AddEntry(hdm,"M_{Med} (DM g_{q} = 1.0) = 2.0 TeV","l")
        #l2.AddEntry(hdma,"M_{Med} (DM g_{q} = 1.0) = 3.0 TeV","l")
        #l2.AddEntry(hdmb,"M_{Med} (DM g_{q} = 1.0) = 4.0 TeV","l")
        #l2.AddEntry(hdmc,"M_{Med} (DM g_{q} = 1.0) = 5.0 TeV","l")
        #l2.AddEntry(hdmd,"M_{Med} (DM g_{q} = 1.0) = 6.0 TeV","l")
    if True:
        l2.AddEntry(halp,"f_{a} = 2.5 TeV","l")
        l2.AddEntry(htripleG,"C_{G}/#Lambda^{2} = 0.010 TeV","l")
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
    CMS_lumi( c, iPeriod, iPos );

    c.SaveAs(fdir+prefix + "_combined_theory"+str(massbin)+"_run2.pdf")
    #c.SaveAs(fdir+prefix + "_combined_theory"+str(massbin)+"_run2.eps")