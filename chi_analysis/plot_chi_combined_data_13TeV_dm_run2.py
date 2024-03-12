from ROOT import *
import ROOT
import array, math, sys
import os
from math import *
from scipy import stats

def rebin(h1,nbins,binning):
    for b in range(h1.GetXaxis().GetNbins()):
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)*h1.GetBinWidth(b+1))
        h1.SetBinError(b+1,h1.GetBinError(b+1)*h1.GetBinWidth(b+1))
    h1=h1.Rebin(nbins,h1.GetName()+"_rebin",binning)
    #h1.Scale(1./h1.Integral())
    for b in range(h1.GetXaxis().GetNbins()):
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)/h1.GetBinWidth(b+1))
        h1.SetBinError(b+1,h1.GetBinError(b+1)/h1.GetBinWidth(b+1))
    return h1

def rebin2(h1,nbins,binning):
    h1=h1.Rebin(nbins,h1.GetName()+"_rebin",binning)
    #h1.Scale(1./h1.Integral())
    for b in range(h1.GetXaxis().GetNbins()):
        h1.SetBinContent(b+1,h1.GetBinContent(b+1)/h1.GetBinWidth(b+1))
        h1.SetBinError(b+1,h1.GetBinError(b+1)/h1.GetBinWidth(b+1))
    return h1

def divBinWidth(hist):
    #hist.SetLineColor(linecolor)
    #hist.SetLineStyle(linestyle)
    #hist.SetLineWidth(linewidth)
    #hist.Scale(1./hist.Integral())
    for b in range(hist.GetNbinsX()):
        hist.SetBinContent(b+1,hist.GetBinContent(b+1)/hist.GetBinWidth(b+1))
        hist.SetBinError(b+1,hist.GetBinError(b+1)/hist.GetBinWidth(b+1))
    return

def getFitResults(fitFile, treename):
    tree=fitFile.Get(treename)
    fitParameters=[]
    fitConstraints=[]
    for name in uncertaintynames:
        print name, tree.floatParsFinal().find(name).getVal()
        fitParameters.append(tree.floatParsFinal().find(name).getVal())
        fitConstraints.append(tree.floatParsFinal().find(name).getError())
    return fitParameters, fitConstraints

def applyFitResults(fitParameters,fitConstraints,uncertainties,hist,hdata):
    hdataG=TGraphAsymmErrors(hdata.Clone(hdata.GetName()+"G"))
    alpha=1.-0.6827
    nevents=0
    for b in range(hdataG.GetN()):
        if unfoldedData:
            N=h14.GetBinContent(b+1)
        else:
            N=1./pow(h14.GetBinError(b+1)/h14.GetBinContent(b+1),2)
        print N
        nevents+=N
        L=0
        if N>0:
            L=ROOT.Math.gamma_quantile(alpha/2.,N,1.)
        U=ROOT.Math.gamma_quantile_c(alpha/2.,N+1,1.)
        hdataG.SetPointEYlow(b,(N-L)/N*hdata.GetBinContent(b+1))
        hdataG.SetPointEYhigh(b,(U-N)/N*hdata.GetBinContent(b+1))
    print "data events:", nevents

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
            if uncertaintynames[uncertainties.index([up,down,central])] in ["jer","JERtail","prefire","model"] or "jes" in uncertaintynames[uncertainties.index([up,down,central])] or "sim" in uncertaintynames[uncertainties.index([up,down,central])]:
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

        hdataGsysstat.SetPointEXlow(b,0)
        hdataGsysstat.SetPointEXhigh(b,0)
        hdataGsysstat.SetPointEYlow(b,sqrt(pow(exp_sumdown*hdataG.GetY()[b],2)+pow(hdataG.GetErrorYlow(b),2)))
        hdataGsysstat.SetPointEYhigh(b,sqrt(pow(exp_sumup*hdataG.GetY()[b],2)+pow(hdataG.GetErrorYhigh(b),2)))
        
        h2new.SetBinContent(b+1,histbackup.GetBinContent(b+1)-theory_sumdown*histbackup.GetBinContent(b+1))
        h3new.SetBinContent(b+1,histbackup.GetBinContent(b+1)+theory_sumup*histbackup.GetBinContent(b+1))
        
    return h2new,h3new,hdataGsysstat

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
    
def shiftWRTmu(uncertainties,mu,hDm,hbPrefit,hsPrefit):
    hDm.Add(hDm,-1)
    hDm.Add(hsPrefit,mu**2)
    hDm.Add(hbPrefit,(1-mu**2))
    return
    

if __name__=="__main__":

    useNNLO=False # choice for QCD
    useM2=False # choice of mu-scale for QCD
    runs="2" # "2" or "3" or "23"
    run=runs[-1]
  
    unfoldedData=False
    isCB=False
    version="v6b"

    print "start ROOT"
    #gROOT.Reset()
    gROOT.SetStyle("Plain")
    gROOT.SetBatch(True)
    gStyle.SetOptStat(0)
    gStyle.SetOptFit(0)
    gStyle.SetTitleOffset(1.2,"Y")
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetPadBottomMargin(0.11)
    gStyle.SetPadTopMargin(0.05)
    gStyle.SetPadRightMargin(0.05)
    gStyle.SetMarkerSize(2.5)
    gStyle.SetHistLineWidth(1)
    gStyle.SetStatFontSize(0.020)
    gStyle.SetTitleSize(0.06, "XYZ")
    gStyle.SetLabelSize(0.05, "XYZ")
    gStyle.SetLegendBorderSize(0)
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetEndErrorSize(5)

    if unfoldedData:
        SaveDir="./fitcheckGENv6b/"
    elif isCB:
        SaveDir="./fitcheckDETCBv6b/"
    else:
        SaveDir="./"
    
    if os.path.exists(SaveDir)==False:
        os.mkdir(SaveDir)

    signalMasses=[1000,1500,1750,2000,2250,2500,3000,3500,4000,4500,5000,6000,7000]
    signalMasses=[7000]

    gas=["0p01","0p05","0p1","0p2","0p25","0p3","0p5","0p75","1","1p5","2p0","2p5","3p0"]

    counter=1100
    signalExtraName={}
    GA={}
    for mdm in ["1"]:
        for ga in gas:
            signalExtraName[counter]="_"+mdm+"_1p0_1p0_Mar5_gdmv_0_gdma_1p0_gv_0_ga_"+ga
            GA[counter]=ga
            counter+=1
    print (counter-1)
    #print signalName,signalExtraName

    if run=="2":
      if useM2:
        prefix="versions/run2ULNNLO_m2/datacard_shapelimit13TeV"
      elif useNNLO:
        prefix="versions/run2ULNNLO_pt12/datacard_shapelimit13TeV"
      else:
        prefix="versions/run2ULNLO_pt12/datacard_shapelimit13TeV"
    elif run=="3":
      prefix="versions/run3ULNNLO_pt12/datacard_shapelimit13TeV"
    else:
      whatprefix
    if unfoldedData:
            prefix+="_unfolded"

    new_hists=[]
    #for j in [1102,1103,1104,1105,1106,1107,1108]:
    for j in [1108]:
        for signalMass in signalMasses:
	    if version=="v6":
              if signalMass<=2500:
            	  massbins=[(2400,3000)]
              elif signalMass<=3000:
            	  massbins=[(2400,3000),(3000,3600)]
              elif signalMass<=3500:
            	  massbins=[(2400,3000),(3000,3600),(3600,4200)]
              elif signalMass<=4000:
            	  massbins=[(2400,3000),(3000,3600),(3600,4200)]
              elif signalMass<=4500:
            	  massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800)]  
              elif signalMass<=5000:
            	  massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000)]
              elif signalMass<=6000:
            	  massbins=[(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000)]
              else:
            	  massbins=[(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
            if version=="v6b":
              if signalMass<=2500:
            	  massbins=[(2400,3000)]
              elif signalMass<=3000:
            	  massbins=[(2400,3000),(3000,3600)]
              elif signalMass<=3500:
            	  massbins=[(2400,3000),(3000,3600),(3600,4200)]
              elif signalMass<=4000:
            	  massbins=[(2400,3000),(3000,3600),(3600,4200)]
              elif signalMass<=4500:
            	  massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800)]  
              elif signalMass<=5000:
            	  massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000)]
              elif signalMass<=6000:
            	  massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000)]
              else:
            	  massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
            elif version=="v9":
  	      if signalMass<=2500:
  	    	massbins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000)]
  	      elif signalMass<=3000:
  	    	massbins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600)]
  	      elif signalMass<=3500:
  	    	massbins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200)]
  	      elif signalMass<=4000:
  	    	massbins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200)]
  	      elif signalMass<=4500:
  	    	massbins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800)]
  	      elif signalMass<=5000:
  	    	massbins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000)]
  	      elif signalMass<=6000:
  	    	massbins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000)]
  	      else:
  	    	massbins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]

            histnameprefix=("DMAxial_Dijet_LO_Mphi_"+str(signalMass)+signalExtraName[j]).replace("7000_1","7000_4000") # FIX produce 7000 mdm=1 sample
            filenameprefix=prefix+"_"+histnameprefix

            uncertaintynames=["pdf","scale","JERtail","model","jer","prefire"] # "scaleAlt", "sim"
            #uncertaintynames.append("model")
            for m in massbins:
                uncertaintynames.append("sim"+str(m[0]))
            for i in range(1,28):
                uncertaintynames.append("jes"+str(i))
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
            if useNNLO:
              for massbin in massbins:
                massindex=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)].index(massbin)
                for chibin in range(len(chi_bins[massindex])-1):
                  uncertaintynames.append("stat"+str(massindex)+"_"+str(chibin))
        
            for massbin in range(len(massbins)):
            
                masstext=str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")

                # Unfolded data

                #filename="datacards/Unfolded_chiNtuple_dataReReco_v3_Coarse_PFHT900_fromCB_AK4SF_pythia8_Pt_170toInf_MatrixInvert.root"
                #masstext=str(massbins[massbin]).strip("()").replace(',',"_").replace(' ',"")
                #if unfoldedData:
                #    histname="dijet_mass_"+masstext+"_chi_unfolded"
    	        #else:
                #    histname="dijet_mass_"+masstext+"_chi"
              
                #print filename
                #fData = TFile.Open(filename)
                #new_hists+=[fData]
                #print histname
                #h14=fData.Get(histname)
    	        #divBinWidth(h14)
                #h14.SetLineColor(4)
                #new_hists+=[h14]

                # open datacards
                if unfoldedData:
                    filename="./invertMatrixAug30/"+filenameprefix+"-run2_chi.root"
                elif isCB:
                    filename="./crystalBallSmearedAug30/"+filenameprefix+"-run2_chi.root"
                else:
                    filename=SaveDir+filenameprefix+"-run2_chi.root"

                if unfoldedData:
                    fitFile=TFile("limitsGenLHCa"+str(j)+"_DMAxial_Dijet_LO_Mphi_"+version+"/fitDiagnostics"+histnameprefix+".root")
                elif isCB:
                    fitFile=TFile("limitsDetCBLHCa"+str(j)+"_DMAxial_Dijet_LO_Mphi_"+version+"/fitDiagnostics"+histnameprefix+".root")
                else:
                    fitFile=TFile("limitsDetLHCa"+str(j)+"_DMAxial_Dijet_LO_Mphi_"+version+"/fitDiagnostics"+histnameprefix+".root")

                print "Fit Parameters from:", fitFile.GetName()

                f=TFile.Open(filename)
                new_hists+=[f]

		dataname="data_obs#chi"+masstext+"_rebin1"
                print "Data Hist Name:",dataname
		
                h14=f.Get(dataname)
		divBinWidth(h14)
                h14.SetLineColor(4)
                new_hists+=[h14]
                
		basename=histnameprefix+"_ALT#chi"+masstext+"_rebin1"
                print "Background Hist Name:",basename

                # NLO OCD plot (with EWK correction)
                
                hNloQcd_orig=f.Get(basename)
                new_hists+=[hNloQcd_orig]
                hNloQcd=hNloQcd_orig.Clone()
                new_hists+=[hNloQcd]
    
                for b in range(hNloQcd.GetXaxis().GetNbins()):
    	            hNloQcd.SetBinContent(b+1,hNloQcd.GetBinContent(b+1)/hNloQcd.GetBinWidth(b+1))
    
                hNloQcd.SetLineColor(1)
                hNloQcd.SetLineStyle(2)
                hNloQcd.SetLineWidth(1)
            
                hNloQcd.SetTitle("")
                hNloQcd.GetXaxis().SetTitle("#chi")
                hNloQcd.GetXaxis().SetTitleOffset(0.7)
                hNloQcd.GetYaxis().SetTitle("N")

                # QCD systematics

                pname='QCD'
                uncertainties=getUncertainties(f,basename,pname,masstext)
                
                new_hists+=[uncertainties]

                hbPrefit=hNloQcd_orig.Clone("QCDprefit"+str(massbins[massbin]))
                divBinWidth(hbPrefit)
                hbPrefit.SetLineWidth(1)
                hbPrefit.SetLineColor(1)
                hbPrefit.SetLineStyle(1)
                new_hists+=[hbPrefit]

                
                # Shift theory predictions according to fitted nuisance parameters
                
                bfitParameters,bfitConstraints=getFitResults(fitFile, 'fit_b')

                h2bnew,h3bnew,h14Gsysstat=applyFitResults(bfitParameters,bfitConstraints,uncertainties,hNloQcd,h14)

                # DM plot
                basename=histnameprefix+"#chi"+masstext+"_rebin1"
                print "Signal Hist Name:",basename
                
                hDm_orig=f.Get(basename)
                new_hists+=[hDm_orig]
                hDm=hDm_orig.Clone()
                new_hists+=[hDm]
    
                for b in range(hDm.GetXaxis().GetNbins()):
    	            hDm.SetBinContent(b+1,hDm.GetBinContent(b+1)/hDm.GetBinWidth(b+1))
    
                hDm.SetLineColor(2)
                hDm.SetLineStyle(5)
                hDm.SetLineWidth(2)
            
                # DM systematics

                pname='DM'
                uncertainties=getUncertainties(f,basename,pname,masstext)
                
                new_hists+=[uncertainties]

                hsPrefit=hDm_orig.Clone("DMprefit"+str(massbins[massbin]))
                divBinWidth(hsPrefit)
                hsPrefit.SetLineWidth(2)
                hsPrefit.SetLineColor(2)
                hsPrefit.SetLineStyle(1)
                new_hists+=[hsPrefit]

                # Shift theory predictions according to fitted signal strength
                tree=fitFile.Get('fit_s')
                mu=tree.floatParsFinal().find("x").getVal()
                shiftWRTmu(uncertainties,mu,hDm,hbPrefit,hsPrefit)
                
                # Shift theory prediction according to fitted nuisance parameters
                
                sfitParameters,sfitConstraints=getFitResults(fitFile, 'fit_s')

                h2snew,h3snew,h14Gsysstat_sb=applyFitResults(sfitParameters,sfitConstraints,uncertainties,hDm,h14)
        
                # Plotting
                h2bnew.SetLineStyle(1)
                h2bnew.SetLineColor(15)
                h2bnew.SetFillColor(10)
                
                h3bnew.SetLineStyle(1)
                h3bnew.SetLineColor(15)
                h3bnew.SetFillColor(15)
    
                h3newnew=h3bnew.Clone()
                h3newnew.SetLineStyle(2)
                h3newnew.SetLineColor(1)
        
                new_hists+=[h2bnew]
                new_hists+=[h3bnew]
        
                canvas=TCanvas("post-fit", "post-fit", 0, 0, 300, 350)
                canvas.cd()

   	        if massbins[massbin][0]>=7000:
   	   	    hNloQcd.SetMinimum(0)
   	   	    hNloQcd.SetMaximum(0.22*hNloQcd.Integral()*3.)
   	        elif massbins[massbin][0]>=5400:
   	   	    hNloQcd.SetMinimum(0.02*hNloQcd.Integral()*15/12)
   	   	    hNloQcd.SetMaximum(0.22*hNloQcd.Integral()*15/12)
   	        else:
   	   	    hNloQcd.SetMinimum(0.045*hNloQcd.Integral()*15/12)
   	   	    hNloQcd.SetMaximum(0.12*hNloQcd.Integral()*15/12)

                hNloQcd.Draw("axissame")
                h3bnew.Draw("histsame")
                h2bnew.Draw("histsame")
                h14.Draw("zesame")
                h14Gsysstat.Draw("zesame")
                hNloQcd.Draw("histsame")
                hbPrefit.Draw("histsame")
                hDm.Draw("histsame")
                hsPrefit.Draw("histsame")
                hNloQcd.Draw("axissame")

                l1=TLegend(0.19,0.6,0.45,0.95,masstext.replace("_","<m_{jj}<"))
                l1.SetFillStyle(0)
                l1.SetTextSize(0.035)
                l1.Draw("same")

                l3=TLegend(0.19,0.8,0.45,0.95,"M_{Med}="+str(signalMass)+", g_{q}="+GA[j].replace("p","."))
                l3.SetFillStyle(0)
                l3.SetTextSize(0.035)
                l3.Draw("same")
                
                if unfoldedData:
                    l2=TLegend(0.45,0.6,0.95,0.93,"Particle-Level")
                elif isCB:
                    l2=TLegend(0.45,0.6,0.95,0.93,"Detector-Level CB Smeared")
                else:
                    l2=TLegend(0.45,0.6,0.95,0.93,"Detector-Level")
                
                l2.SetTextSize(0.035)
                l2.AddEntry(h14,"Data","ple")
                l2.AddEntry(h3newnew,"Background post-fit","fl")
                l2.AddEntry(hbPrefit,"Background pre-fit","l")
                l2.AddEntry(hDm,"Background+Signal post-fit","fl")
                l2.AddEntry(hsPrefit,"Background+Signal pre-fit","l")
                l2.SetFillStyle(0)
                l2.Draw("same")
        
                h2bnew.SetLineStyle(1)
                h2bnew.SetLineColor(15)
                h2bnew.SetFillColor(10)
                
                h3bnew.SetLineStyle(1)
                h3bnew.SetLineColor(15)
                h3bnew.SetFillColor(15)

                print "systematic uncertainty:", (h3bnew.GetBinContent(1)-hNloQcd.GetBinContent(1))/hNloQcd.GetBinContent(1)
        
                canvas.SaveAs(SaveDir + prefix + "_combined_fit_"+histnameprefix+"_"+masstext+"_run2.pdf")

                #sys.exit()
