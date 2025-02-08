import os, sys
from ROOT import * 
import array

#gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
#gROOT.Reset()
gROOT.SetStyle("Plain")
gROOT.SetBatch(True)
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

if __name__ == '__main__':

    run="2"
    prefix="/data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/versions/run"+run+"ULNNLO_m2/datacard_shapelimit13TeV"

    massbins=[(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
    all_mass_bins=[(1200,1500),(1500,1900),(1900,2400),(2400,3000),(3000,3600),(3600,4200),(4200,4800),(4800,5400),(5400,6000),(6000,7000),(7000,13000)]
    mass_bin_offset=len(all_mass_bins)-len(massbins)

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

    canvas = TCanvas("","",0,0,600,600)
    canvas.SetRightMargin(0.15)
    gStyle.SetPalette(kLightTemperature)
    plots=[]
    legends=[]
    new_hists=[]

    name="response"
    normed=True

    filename=prefix+"_GEN-QCD-run2_chi.root"
    print filename
    f = TFile.Open(filename)

    binsj=0
    for j in range(len(massbins)):
      for chibinj in range(len(chi_bins[j+mass_bin_offset])-1):
        binsj+=1

    response=TH2F("response","response",binsj, 0, binsj-1, binsj, 0, binsj-1)
    
    binsj=0
    for j in range(len(massbins)):
      for chibinj in range(len(chi_bins[j+mass_bin_offset])-1):
        responsesum=0
        for i in range(len(massbins)):
         for chibini in range(len(chi_bins[i+mass_bin_offset])-1):
          hQCD=f.Get("QCD_ALT#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1"+"_bin_"+str(j)+"_"+str(chibinj)+"_")
          responsesum+=hQCD.GetBinContent(chibini+1)
        binsi=0
        for i in range(len(massbins)):
         for chibini in range(len(chi_bins[i+mass_bin_offset])-1):
          hQCD=f.Get("QCD_ALT#chi"+str(massbins[i][0])+"_"+str(massbins[i][1])+"_rebin1"+"_bin_"+str(j)+"_"+str(chibinj)+"_")
          if responsesum>0:
           if abs(i-j)<1 and abs(chibini-chibinj)<1:
            response.SetBinContent(binsj+1,binsi+1,hQCD.GetBinContent(chibini+1)/(responsesum if normed else 1.0))
          binsi+=1
        binsj+=1
   
    response.GetZaxis().SetRangeUser(0,1)
    response.GetXaxis().SetTitle("Gen bin")
    response.GetYaxis().SetTitle("Reco bin")
    response.Draw("colz")

    canvas.SaveAs(prefix + "_"+name+"_run"+run+".pdf")

    m=TMatrixD(response.GetXaxis().GetNbins(), response.GetYaxis().GetNbins())
    for x in range(response.GetXaxis().GetNbins()):
      for y in range(response.GetYaxis().GetNbins()):
        m[x,y]=response.GetBinContent(x+1,y+1)
    m.Print()
    
    from ROOT import TDecompSVD
    svd=TDecompSVD(m)
    #svd.GetMatrix().Print()
    print [s for s in svd.GetSig()]
    print "condition", max(svd.GetSig())/max(0,min(svd.GetSig()))
