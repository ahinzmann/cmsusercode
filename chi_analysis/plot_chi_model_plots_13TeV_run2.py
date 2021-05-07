import os, sys
import array
from ROOT import * 

#gROOT.Reset()
gROOT.SetStyle("Plain")
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

   variables=["#chi"]
   label=["#chi"]

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
   
   prefix="datacard_shapelimit13TeV_run2_"
   postfix="_run2"
   
   for scenario in ["model"]:

     data=["2016",
  	   "2017",
  	   "2018",
  	   ]
     lumifactor={}
     lumifactor["2016"]=36.33/137.6
     lumifactor["2017"]=41.53/137.6
     lumifactor["2018"]=59.74/137.6
     mclegend={}
     mc={}
     if scenario=="model":
       mclegend[1]="Madgraph+Pythia"
       mc[1]=[("2016_QCDmadgraph-HT200to300",1712000./56709875*lumifactor["2016"]),
     	   ("2016_QCDmadgraph-HT300to500",347700./53096517*lumifactor["2016"]),
     	   ("2016_QCDmadgraph-HT500to700",32100./52906552*lumifactor["2016"]),
     	   ("2016_QCDmadgraph-HT700to1000",6831./36741540*lumifactor["2016"]),
     	   ("2016_QCDmadgraph-HT1000to1500",1207./15210939*lumifactor["2016"]),
     	   ("2016_QCDmadgraph-HT1500to2000",119.9/11839357*lumifactor["2016"]),
     	   ("2016_QCDmadgraph-HT2000toInf",25.24/5947849*lumifactor["2016"]),
     	   ("2017_QCDmadgraph-HT200to300",1545000./58990434*lumifactor["2017"]),
     	   ("2017_QCDmadgraph-HT300to500",323300./58748739*lumifactor["2017"]),
     	   ("2017_QCDmadgraph-HT500to700",30000./54366431*lumifactor["2017"]),
     	   ("2017_QCDmadgraph-HT700to1000",6324./46924322*lumifactor["2017"]),
     	   ("2017_QCDmadgraph-HT1000to1500",1090./16495598*lumifactor["2017"]),
     	   ("2017_QCDmadgraph-HT1500to2000",101./11196479*lumifactor["2017"]),
     	   ("2017_QCDmadgraph-HT2000toInf",20.43/5362513*lumifactor["2017"]),
     	   ("2018_QCDmadgraph-HT200to300",1461000./54289442*lumifactor["2018"]),
     	   ("2018_QCDmadgraph-HT300to500",311900./54512704*lumifactor["2018"]),
     	   ("2018_QCDmadgraph-HT500to700",29070./53919811*lumifactor["2018"]),
     	   ("2018_QCDmadgraph-HT700to1000",5962./48158738*lumifactor["2018"]),
     	   ("2018_QCDmadgraph-HT1000to1500",1005./14945819*lumifactor["2018"]),
     	   ("2018_QCDmadgraph-HT1500to2000",101.8/10707847*lumifactor["2018"]),
     	   ("2018_QCDmadgraph-HT2000toInf",20.54/5329144*lumifactor["2018"]),
     	   ]
       mclegend[2]="Pythia"
       mc[2]=[("2018_QCDpythia",1),
     	   ]
       mclegend[3]="Herwig7"
       mc[3]=[("2018_QCDherwig",1),
     	   ]
       #mclegend[4]="MG+Py QCD (smeared)"
       #mc[2]=[("QCDmadgraph_JER_2016_1",1),
     #	    ("QCDmadgraph_JER_2016_2",1),
     #	    ("QCDmadgraph_JER_2016_3",1),
     #	    ("QCDmadgraph_JER_2016_4",1),
     #	    ("QCDmadgraph_JER_2016_5",1),
     #	    ("QCDmadgraph_JER_2016_6",1),
     #	    ("QCDmadgraph_JER_2016_7",1),
     #	   ]
       #mclegend[3]="Py QCD (smeared)"
       #mc[3]=[("QCD_JER_2016_1",1),
     	#    ("QCD_JER_2016_2",1),
     	#    ("QCD_JER_2016_3",1),
     	#    ("QCD_JER_2016_4",1),
     	#    ("QCD_JER_2016_5",1),
     	#    ("QCD_JER_2016_6",1),
     	#    ("QCD_JER_2016_7",1),
     	#   ]
       #mclegend[4]="Py QCD (particle)"
       #mc[4]=[("QCD_JES_2017_1",1),
     #	    ("QCD_JES_2017_2",1),
     #	    ("QCD_JES_2017_3",1),
     #	    ("QCD_JES_2017_4",1),
     #	    ("QCD_JES_2017_5",1),
     #	    ("QCD_JES_2017_6",1),
     #	    ("QCD_JES_2017_7",1),
     	#   ]
       #mclegend[5]="Py QCD (2018 smeared)"
       #mc[5]=[("QCD_JER_2018_1",1),
     #	    ("QCD_JER_2018_2",1),
     #	    ("QCD_JER_2018_3",1),
     #	    ("QCD_JER_2018_4",1),
     #	    ("QCD_JER_2018_5",1),
     #	    ("QCD_JER_2018_6",1),
     #	    ("QCD_JER_2018_7",1),
     #	   ]
       #mclegend[4]="Hw QCD (2016 smeared)"
       #mc[4]=[("QCDherwig_JER_2016_1",1),
     	#    ("QCDherwig_JER_2016_2",1),
     	#    ("QCDherwig_JER_2016_3",1),
     	#    ("QCDherwig_JER_2016_4",1),
     	#    ("QCDherwig_JER_2016_5",1),
     	#    ("QCDherwig_JER_2016_6",1),
     	#    ("QCDherwig_JER_2016_7",1),
     	#   ]
     f_data=[]
     for name in data:
  	f_data+=[TFile.Open(prefix+name+"_chi.root")]
     f_mc={}
     for m in range(1,len(mc.keys())+1):
       f_mc[m]=[]
     for m in [1,2,3]:
       for name,xsec in mc[m]:
  	f_mc[m]+=[TFile.Open(prefix+name+"_chi.root")]
     for m in range(4,len(mc.keys())+1):
       for name,xsec in mc[m]:
  	f_mc[m]+=[TFile.Open(prefix.replace("_run2","")+name+"_chi.root")]

     for var in ["mass"]:
       canvas = TCanvas("","",0,0,200,300)
       canvas.Divide(1,2,0,0,0)
       canvas.GetPad(1).SetPad(0.0,0.30,1.0,1.0)
       canvas.GetPad(1).SetLeftMargin(0.15)
       canvas.GetPad(1).SetRightMargin(0.08)
       canvas.GetPad(1).SetTopMargin(0.08)
       canvas.GetPad(1).SetBottomMargin(0.05)
       canvas.GetPad(2).SetPad(0.0,0.0,1.0,0.30)
       canvas.GetPad(2).SetLeftMargin(0.15)
       canvas.GetPad(2).SetRightMargin(0.08)
       canvas.GetPad(2).SetTopMargin(0.08)
       canvas.GetPad(2).SetBottomMargin(0.45)
       canvas.cd(1)
       canvas.GetPad(1).SetLogy(True)

       legend=TLegend(0.5,0.5,0.95,0.9,"1<=#Chi<16 , y_{boost}<1.11")

       hist=f_data[0].Get(prefix+data[0]+var)
       for i in range(1,len(data)):
  	   hist.Add(f_data[i].Get(prefix+data[i]+var))
       hist.SetLineColor(1)
       hist.SetMarkerStyle(24)
       hist.SetMarkerColor(1)
       hist.SetMarkerSize(0.2)
       hist.SetTitle("")
       hist.GetXaxis().SetLabelColor(0)
       hist.GetYaxis().SetTitle("N")
       hist.GetXaxis().SetRangeUser(2400,8400)
       hist.GetYaxis().SetRangeUser(0.5,hist.GetMaximum()*1.5)
       hist.GetXaxis().SetTitleOffset(1.1)
       hist.GetYaxis().SetTitleOffset(1.1)
       hist.GetXaxis().SetLabelSize(0.05)
       hist.GetYaxis().SetLabelSize(0.05)
       hist.GetXaxis().SetTitleSize(0.06)
       hist.GetYaxis().SetTitleSize(0.06)
       hist.SetStats(False)
       hist.Draw("pe")
       legend.AddEntry(hist,"Data","lpe")

       hist_mc={}
       for m in [1,2,3]:
         v=prefix+mc[m][0][0].replace("-HT200to300","").replace("-0","")+var
	 print v
         hist_mc[m]=f_mc[m][0].Get(v)
         for i in range(1,len(mc[m])):
   	     v=prefix+mc[m][0][0].replace("-HT200to300","").replace("2016",mc[m][i][0].split("_")[0])+var
	     print v
  	     hist_mc[m].Add(f_mc[m][i].Get(v),mc[m][i][1]/mc[m][0][1])
         hist_mc[m].Scale(hist.Integral(hist.FindBin(2400),hist.GetNbinsX())/hist_mc[m].Integral(hist_mc[m].FindBin(2400),hist_mc[m].GetNbinsX()))
         hist_mc[m].SetLineColor(colors[m])
         hist_mc[m].SetStats(False)
         hist_mc[m].Draw("histesame")
         legend.AddEntry(hist_mc[m],mclegend[m],"l")

       for m in range(4,len(mc.keys())+1):
         v=mc[m][0][0].replace("_JER_2016_1","").replace("_JER_2017_1","").replace("_JER_2018_1","").replace("_JES_2016_1","").replace("_JES_2017_1","").replace("_JES_2018_1","")+var
         hist_mc[m]=f_mc[m][0].Get(v)
         for i in range(1,len(mc[m])):
  	   hist_mc[m].Add(f_mc[m][i].Get(mc[m][0][0].replace("_JER_2016_1","").replace("_JER_2017_1","").replace("_JER_2018_1","").replace("_JES_2016_1","").replace("_JES_2017_1","").replace("_JES_2018_1","")+var),mc[m][i][1]/mc[m][0][1])
         hist_mc[m].Scale(hist.Integral(hist.FindBin(2400),hist.GetNbinsX())/hist_mc[m].Integral(hist_mc[m].FindBin(2400),hist_mc[m].GetNbinsX()))
         hist_mc[m].SetLineColor(colors[m])
         hist_mc[m].SetStats(False)
         hist_mc[m].Draw("histesame")
         legend.AddEntry(hist_mc[m],mclegend[m],"l")

       hist.Draw("pesame")
       hist.Draw("axissame")
       
       legend.SetTextSize(0.04)
       legend.SetFillStyle(0)
       legend.Draw("same")

       canvas.cd(2)
       ratio=hist.Clone(hist.GetName()+"ratio")
       ratio.Divide(hist,hist)
       for b in range(hist.GetNbinsX()):
  	 if hist.GetBinContent(b+1)>0:
  	   ratio.SetBinError(b+1,hist.GetBinError(b+1)/hist.GetBinContent(b+1))
       ratio.SetTitle("")
       ratio.GetYaxis().SetTitle("Sim / Data")
       ratio.GetYaxis().SetTitleSize(0.11)
       ratio.GetYaxis().SetTitleOffset(0.5)
       ratio.SetMarkerSize(0.1)
       ratio.GetYaxis().SetLabelSize(0.13)
       ratio.GetYaxis().SetRangeUser(0,2)
       ratio.GetXaxis().SetNdivisions(506)
       ratio.GetYaxis().SetNdivisions(503)
       ratio.GetXaxis().SetLabelColor(1)
       ratio.GetXaxis().SetTitle("dijet mass [GeV]")
       ratio.GetXaxis().SetTitleSize(0.14)
       ratio.GetXaxis().SetTitleOffset(1.1)
       ratio.GetXaxis().SetLabelSize(0.12)
       ratio.Draw("histe")
       ratio_mc={}
       for m in range(1,len(mc.keys())+1):
         ratio_mc[m]=hist_mc[m].Clone(hist_mc[m].GetName()+"ratio")
         ratio_mc[m].Divide(hist_mc[m],hist)
         for b in range(hist.GetNbinsX()):
  	   if hist.GetBinContent(b+1)>0:
  	     ratio_mc[m].SetBinError(b+1,hist.GetBinError(b+1)/hist.GetBinContent(b+1))
         ratio_mc[m].Draw("histesame")
       ratio.Draw("axissame")

       canvas.cd(1)
       hist.GetYaxis().SetTitleOffset(1.2)
       
       canvas.SaveAs("chi_"+scenario+"_plots_mass"+postfix+".root")
       canvas.SaveAs("chi_"+scenario+"_plots_mass"+postfix+".pdf")

     for var in variables:
      log=(var=="p_{T1}" or var=="p_{T2}" or var=="METsumET" or var=="#Delta#phi")
      legends=[]
      for mass in range(len(massbins)):
  	  print var, str(massbins[mass][0])+"_"+str(massbins[mass][1])
  	  canvas = TCanvas("","",0,0,200,300)
  	  canvas.Divide(1,2,0,0,0)
  	  canvas.GetPad(1).SetPad(0.0,0.30,1.0,1.0)
  	  canvas.GetPad(1).SetLeftMargin(0.15)
  	  canvas.GetPad(1).SetRightMargin(0.08)
  	  canvas.GetPad(1).SetTopMargin(0.08)
  	  canvas.GetPad(1).SetBottomMargin(0.05)
  	  canvas.GetPad(2).SetPad(0.0,0.0,1.0,0.3)
  	  canvas.GetPad(2).SetLeftMargin(0.15)
  	  canvas.GetPad(2).SetRightMargin(0.08)
  	  canvas.GetPad(2).SetTopMargin(0.08)
  	  canvas.GetPad(2).SetBottomMargin(0.45)
  	  canvas.cd(1)
  	  canvas.GetPad(1).SetLogy(log)
  	  legend=TLegend(0.45,0.6,0.95,0.90,(str(massbins[mass][0])+"<m_{jj}<"+str(massbins[mass][1])+" GeV").replace("7000<m_{jj}<13000","m_{jj}>7000"))
  	  legends+=[legend]
      
  	  name=prefix+data[0]+var+str(massbins[mass][0])+"_"+str(massbins[mass][1])
  	  if var=="#chi": name+="_rebin1"
  	  hist=f_data[0].Get(name)
  	  for i in range(1,len(data)):
  	      hist.Add(f_data[i].Get(name.replace("2016",data[i])))
  	  if var=="#chi":
  	      hist=hist.Rebin(len(chi_binnings[mass])-1,hist.GetName()+"_rebin1",chi_binnings[mass])
  	      for b in range(hist.GetNbinsX()):
  		 hist.SetBinContent(b+1,hist.GetBinContent(b+1)/hist.GetBinWidth(b+1))
  		 hist.SetBinError(b+1,hist.GetBinError(b+1)/hist.GetBinWidth(b+1))
  	  hist.SetLineColor(1)
  	  hist.SetMarkerStyle(24)
  	  hist.SetMarkerSize(0.2)
  	  miny=0
  	  if hist.Integral()>0:
  	      miny=log*0.1/hist.GetMaximum()
  	      hist.Scale(1./hist.Integral())
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
  	  legend.AddEntry(hist,"Data","lpe")

  	  print "mass bin",mass,"data integral",hist.Integral()

          for m in [1,2,3]:
  	    name=prefix+mc[m][0][0].replace("-HT200to300","").replace("-0","")+var+str(massbins[mass][0])+"_"+str(massbins[mass][1])
  	    if var=="#chi": name+="_rebin1"
  	    hist_mc[m]=f_mc[m][0].Get(name)
  	    for i in range(1,len(mc[m])):
  	   	hist_mc[m].Add(f_mc[m][i].Get(name.replace("2016",mc[m][i][0].split("_")[0])),mc[m][i][1]/mc[m][0][1])
  	    if var=="#chi":
  	   	hist_mc[m]=hist_mc[m].Rebin(len(chi_binnings[mass])-1,hist_mc[m].GetName()+"_rebin1",chi_binnings[mass])
  	   	for b in range(hist_mc[m].GetNbinsX()):
  	 	   hist_mc[m].SetBinContent(b+1,hist_mc[m].GetBinContent(b+1)/hist_mc[m].GetBinWidth(b+1))
  	 	   hist_mc[m].SetBinError(b+1,hist_mc[m].GetBinError(b+1)/hist_mc[m].GetBinWidth(b+1))
  	    if hist_mc[m].Integral()>0:
  	   	hist_mc[m].Scale(hist.Integral()/hist_mc[m].Integral())
  	    hist_mc[m].SetLineColor(colors[m])
  	    hist_mc[m].SetStats(False)
  	    hist_mc[m].Draw("histesame")
  	    legend.AddEntry(hist_mc[m],mclegend[m],"l")

          for m in range(4,len(mc.keys())+1):
  	    name=mc[m][0][0].replace("_JER_2016_1","").replace("_JER_2017_1","").replace("_JER_2018_1","").replace("_JES_2016_1","").replace("_JES_2017_1","").replace("_JES_2018_1","")+var+str(massbins[mass][0])+"_"+str(massbins[mass][1])
  	    if var=="#chi": name+="_rebin1"
  	    hist_mc[m]=f_mc[m][0].Get(name)
  	    for i in range(1,len(mc[m])):
  	      hist_mc[m].Add(f_mc[m][i].Get(name),mc[m][i][1]/mc[m][0][1])
  	    if var=="#chi":
  	      hist_mc[m]=hist_mc[m].Rebin(len(chi_binnings[mass])-1,hist_mc[m].GetName()+"_rebin1",chi_binnings[mass])
  	      for b in range(hist_mc[m].GetNbinsX()):
  		 hist_mc[m].SetBinContent(b+1,hist_mc[m].GetBinContent(b+1)/hist_mc[m].GetBinWidth(b+1))
  		 hist_mc[m].SetBinError(b+1,hist_mc[m].GetBinError(b+1)/hist_mc[m].GetBinWidth(b+1))
  	    if hist_mc[m].Integral()>0:
  	      hist_mc[m].Scale(hist.Integral()/hist_mc[m].Integral())
  	    hist_mc[m].SetLineColor(colors[m])
  	    hist_mc[m].SetStats(False)
  	    hist_mc[m].Draw("histesame")
  	    legend.AddEntry(hist_mc[m],mclegend[m],"l")

  	  hist.Draw("lesame")
  	  hist.Draw("axissame")

  	  legend.SetTextSize(0.04)
  	  legend.SetFillStyle(0)
  	  legend.Draw("same")

  	  canvas.cd(2)
  	  ratio=hist.Clone(hist.GetName()+"ratio")
  	  ratio.Divide(hist,hist)
  	  for b in range(hist.GetNbinsX()):
  	    if hist.GetBinContent(b+1)>0:
  	      ratio.SetBinError(b+1,hist.GetBinError(b+1)/hist.GetBinContent(b+1))
  	  ratio.SetTitle("")
  	  ratio.GetYaxis().SetTitle("Sim / Data")
  	  ratio.GetYaxis().SetTitleSize(0.11)
  	  ratio.GetYaxis().SetTitleOffset(0.5)
  	  ratio.SetMarkerSize(0.1)
  	  ratio.GetYaxis().SetLabelSize(0.13)
  	  ratio.GetYaxis().SetRangeUser(0,2)
  	  if var=="#chi":
  	    ratio.GetYaxis().SetRangeUser(0.5,1.5)
  	  ratio.GetXaxis().SetNdivisions(506)
  	  ratio.GetYaxis().SetNdivisions(503)
  	  ratio.GetXaxis().SetLabelColor(1)
  	  ratio.GetXaxis().SetTitle(label[variables.index(var)])
  	  ratio.GetXaxis().SetTitleSize(0.14)
  	  ratio.GetXaxis().SetTitleOffset(1.1)
  	  ratio.GetXaxis().SetLabelSize(0.12)
  	  ratio.Draw("histe")
          for m in range(1,len(mc.keys())+1):
  	    ratio_mc[m]=hist_mc[m].Clone(hist_mc[m].GetName()+"ratio")
  	    ratio_mc[m].Divide(hist_mc[m],hist)
  	    for b in range(hist.GetNbinsX()):
  	      if hist.GetBinContent(b+1)>0:
  	        ratio_mc[m].SetBinError(b+1,hist.GetBinError(b+1)/hist.GetBinContent(b+1))
  	    ratio_mc[m].Draw("histesame")
  	  ratio.Draw("axissame")

  	  canvas.cd(1)
  	  hist.GetYaxis().SetTitleOffset(1.2)

  	  canvas.SaveAs("chi_"+scenario+"_plots_"+var.replace("_","").replace("{","").replace("}","").replace("#","").replace("+","p").replace("-","m")+"_"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+postfix+".root")
  	  canvas.SaveAs("chi_"+scenario+"_plots_"+var.replace("_","").replace("{","").replace("}","").replace("#","").replace("+","p").replace("-","m")+"_"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+postfix+".pdf")
