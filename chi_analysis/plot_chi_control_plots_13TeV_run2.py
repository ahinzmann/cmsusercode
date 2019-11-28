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

   variables=["#chi","y_{boost}","p_{T1}","p_{T2}","y_{1}","y_{2}","METsumET","dPtsumPt","#Delta#phi"]
   label=["#chi","y_{boost}","p_{T1} [GeV]","p_{T2} [GeV]","y_{1}","y_{2}","missing E_{T} / #sum E_{T}","(p_{T1}-p_{T2})/(p_{T1}+p_{T2})","#Delta#phi"]
   variables=["#chi","y_{boost}","p_{T1}","p_{T2}","y_{1}","y_{2}","#Delta#phi"]
   label=["#chi","y_{boost}","p_{T1} [GeV]","p_{T2} [GeV]","y_{1}","y_{2}","#Delta#phi"]

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

   data=["data_2016",
         ]
   data2=["data_2017",
         ]
   data3=["data_2018",
         ]
   mc=[("2016_QCDmadgraph-HT200to300",1712000./56709875),
       ("2016_QCDmadgraph-HT300to500",347700./53096517),
       ("2016_QCDmadgraph-HT500to700",32100./52906552),
       ("2016_QCDmadgraph-HT700to1000",6831./36741540),
       ("2016_QCDmadgraph-HT1000to1500",1207./15210939),
       ("2016_QCDmadgraph-HT1500to2000",119.9/11839357),
       ("2016_QCDmadgraph-HT2000toInf",25.24/5947849),
       ]
   mc2=[("2017_QCDmadgraph-HT200to300",1545000./58990434),
       ("2017_QCDmadgraph-HT300to500",323300./58748739),
       ("2017_QCDmadgraph-HT500to700",30000./54366431),
       ("2017_QCDmadgraph-HT700to1000",6324./46924322),
       ("2017_QCDmadgraph-HT1000to1500",1090./16495598),
       ("2017_QCDmadgraph-HT1500to2000",101./11196479),
       ("2017_QCDmadgraph-HT2000toInf",20.43/5362513),
       ]
   mc3=[("2018_QCDmadgraph-HT200to300",1461000./54289442),
       ("2018_QCDmadgraph-HT300to500",311900./54512704),
       ("2018_QCDmadgraph-HT500to700",29070./53919811),
       ("2018_QCDmadgraph-HT700to1000",5962./48158738),
       ("2018_QCDmadgraph-HT1000to1500",1005./14945819),
       ("2018_QCDmadgraph-HT1500to2000",101.8/10707847),
       ("2018_QCDmadgraph-HT2000toInf",20.54/5329144),
       ]
   f_data=[]
   f_data2=[]
   f_data3=[]
   for name in data:
      f_data+=[TFile.Open(prefix+name+"_chi.root")]
   for name in data2:
      f_data2+=[TFile.Open(prefix+name+"_chi.root")]
   for name in data3:
      f_data3+=[TFile.Open(prefix+name+"_chi.root")]
   f_mc=[]
   f_mc2=[]
   f_mc3=[]
   for name,xsec in mc:
      f_mc+=[TFile.Open(prefix+name+"_chi.root")]
   for name,xsec in mc2:
      f_mc2+=[TFile.Open(prefix+name+"_chi.root")]
   for name,xsec in mc3:
      f_mc3+=[TFile.Open(prefix+name+"_chi.root")]

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

     hist=f_data[0].Get(prefix+data[0].replace("data_","")+var)
     for i in range(1,len(data)):
    	 hist.Add(f_data[i].Get(prefix+data[i].replace("data_","")+var))
     hist.SetLineColor(1)
     hist.SetMarkerStyle(24)
     hist.SetMarkerColor(1)
     hist.SetMarkerSize(0.2)
     hist.SetTitle("")
     hist.GetXaxis().SetLabelColor(0)
     hist.GetYaxis().SetTitle("N")
     hist.GetXaxis().SetRangeUser(1200,8400)
     hist.GetYaxis().SetRangeUser(0.5,hist.GetMaximum()*1.5)
     hist.GetXaxis().SetTitleOffset(1.1)
     hist.GetYaxis().SetTitleOffset(1.1)
     hist.GetXaxis().SetLabelSize(0.05)
     hist.GetYaxis().SetLabelSize(0.05)
     hist.GetXaxis().SetTitleSize(0.06)
     hist.GetYaxis().SetTitleSize(0.06)
     hist.SetStats(False)
     hist.Draw("pe")
     legend.AddEntry(hist,"Data (2016)","lpe")

     hist2=f_data2[0].Get(prefix+data2[0].replace("data_","")+var)
     for i in range(1,len(data2)):
    	 hist2.Add(f_data2[i].Get(prefix+data2[i].replace("data_","")+var))
     hist2.SetLineColor(2)
     hist2.SetMarkerStyle(25)
     hist2.SetMarkerColor(2)
     hist2.SetMarkerSize(0.2)
     hist2.SetStats(False)
     hist2.Draw("pesame")
     legend.AddEntry(hist2,"Data (2017)","lpe")
     
     hist3=f_data3[0].Get(prefix+data3[0].replace("data_","")+var)
     for i in range(1,len(data3)):
    	 hist3.Add(f_data3[i].Get(prefix+data3[i].replace("data_","")+var))
     hist3.SetLineColor(4)
     hist3.SetMarkerStyle(26)
     hist3.SetMarkerColor(4)
     hist3.SetMarkerSize(0.2)
     hist3.SetStats(False)
     hist3.Draw("pesame")
     legend.AddEntry(hist3,"Data (2018)","lpe")
     
     hist_mc=f_mc[0].Get(prefix+mc[0][0].replace("-HT200to300","")+var)
     for i in range(1,len(mc)):
         hist_mc.Add(f_mc[i].Get(prefix+mc[0][0].replace("-HT200to300","")+var),mc[i][1]/mc[0][1])
     hist_mc.Scale(hist.Integral(hist.FindBin(2400),hist.GetNbinsX())/hist_mc.Integral(hist_mc.FindBin(2400),hist_mc.GetNbinsX()))
     hist_mc.SetLineColor(1)
     hist_mc.SetStats(False)
     hist_mc.Draw("histsame")
     legend.AddEntry(hist_mc,"MG+Py QCD (2016)","l")

     hist_mc2=f_mc2[0].Get(prefix+mc2[0][0].replace("-HT200to300","")+var)
     for i in range(1,len(mc2)):
    	 hist_mc2.Add(f_mc2[i].Get(prefix+mc2[0][0].replace("-HT200to300","")+var),mc2[i][1]/mc2[0][1])
     hist_mc2.Scale(hist2.Integral(hist2.FindBin(2400),hist2.GetNbinsX())/hist_mc2.Integral(hist_mc2.FindBin(2400),hist_mc2.GetNbinsX()))
     hist_mc2.SetLineColor(2)
     hist_mc2.SetStats(False)
     hist_mc2.Draw("histsame")
     legend.AddEntry(hist_mc2,"MG+Py QCD (2017)","l")

     hist_mc3=f_mc3[0].Get(prefix+mc3[0][0].replace("-HT200to300","")+var)
     for i in range(1,len(mc3)):
    	 hist_mc3.Add(f_mc3[i].Get(prefix+mc3[0][0].replace("-HT200to300","")+var),mc3[i][1]/mc3[0][1])
     hist_mc3.Scale(hist3.Integral(hist3.FindBin(2400),hist3.GetNbinsX())/hist_mc3.Integral(hist_mc3.FindBin(2400),hist_mc3.GetNbinsX()))
     hist_mc3.SetLineColor(4)
     hist_mc3.SetStats(False)
     hist_mc3.Draw("histsame")
     legend.AddEntry(hist_mc3,"MG+Py QCD (2018)","l")

     hist.Draw("pesame")
     hist2.Draw("pesame")
     hist3.Draw("pesame")
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
     ratio.GetYaxis().SetTitleSize(0.18)
     ratio.GetYaxis().SetTitleOffset(0.3)
     ratio.SetMarkerSize(0.1)
     ratio.GetYaxis().SetLabelSize(0.2)
     ratio.GetYaxis().SetRangeUser(0,2)
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
     ratio_mc3=hist_mc3.Clone(hist_mc3.GetName()+"ratio")
     ratio_mc3.Divide(hist_mc3,hist3)
     for b in range(hist3.GetNbinsX()):
       if hist3.GetBinContent(b+1)>0:
    	 ratio_mc3.SetBinError(b+1,hist3.GetBinError(b+1)/hist3.GetBinContent(b+1))
     ratio_mc3.Draw("histsame")
     ratio.Draw("axissame")

     canvas.cd(3)
     ddratio=hist.Clone(hist.GetName()+"ddratio")
     ddratio.Divide(hist,hist)
     for b in range(hist.GetNbinsX()):
       if hist.GetBinContent(b+1)>0:
    	 ddratio.SetBinError(b+1,hist.GetBinError(b+1)/hist.GetBinContent(b+1))
     ddratio.SetTitle("")
     ddratio.GetYaxis().SetTitle("Data / Data")
     ddratio.GetYaxis().SetTitleSize(0.11)
     ddratio.GetYaxis().SetTitleOffset(0.5)
     ddratio.SetMarkerSize(0.1)
     ddratio.GetYaxis().SetLabelSize(0.13)
     ddratio.GetYaxis().SetRangeUser(0,2)
     ddratio.GetXaxis().SetNdivisions(506)
     ddratio.GetYaxis().SetNdivisions(503)
     ddratio.GetXaxis().SetLabelColor(1)
     ddratio.GetXaxis().SetTitle("dijet mass [GeV]")
     ddratio.GetXaxis().SetTitleSize(0.14)
     ddratio.GetXaxis().SetTitleOffset(1.1)
     ddratio.GetXaxis().SetLabelSize(0.12)
     ddratio.Draw("histe")
     ddratio_mc2=hist_mc2.Clone(hist_mc2.GetName()+"ddratio")
     ddratio_mc2.Divide(hist_mc2,hist_mc)
     for b in range(hist_mc.GetNbinsX()):
       if hist_mc.GetBinContent(b+1)>0:
    	 ddratio_mc2.SetBinError(b+1,hist_mc.GetBinError(b+1)/hist_mc.GetBinContent(b+1))
     ddratio_mc2.Draw("histsame")
     ddratio_mc3=hist_mc3.Clone(hist_mc3.GetName()+"ddratio")
     ddratio_mc3.Divide(hist_mc3,hist_mc)
     for b in range(hist_mc.GetNbinsX()):
       if hist_mc.GetBinContent(b+1)>0:
    	 ddratio_mc3.SetBinError(b+1,hist_mc.GetBinError(b+1)/hist_mc.GetBinContent(b+1))
     ddratio_mc3.Draw("histsame")
     ratio.Draw("axissame")

     canvas.cd(1)
     hist.GetYaxis().SetTitleOffset(1.2)
     
     canvas.SaveAs("chi_control_plots_mass"+postfix+".root")
     canvas.SaveAs("chi_control_plots_mass"+postfix+".pdf")

   for var in variables:
    log=(var=="p_{T1}" or var=="p_{T2}" or var=="METsumET" or var=="#Delta#phi")
    legends=[]
    for mass in range(len(massbins)):
        print var, str(massbins[mass][0])+"_"+str(massbins[mass][1])
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
    
        name=prefix+data[0].replace("data_","")+var+str(massbins[mass][0])+"_"+str(massbins[mass][1])
	if var=="#chi": name+="_rebin1"
        hist=f_data[0].Get(name)
        for i in range(1,len(data)):
            hist.Add(f_data[i].Get(name))
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
        legend.AddEntry(hist,"Data (2016)","lpe")

        name=prefix+data2[0].replace("data_","")+var+str(massbins[mass][0])+"_"+str(massbins[mass][1])
	if var=="#chi": name+="_rebin1"
        hist2=f_data2[0].Get(name)
        for i in range(1,len(data2)):
            hist2.Add(f_data2[i].Get(name))
	if var=="#chi":
            hist2=hist2.Rebin(len(chi_binnings[mass])-1,hist2.GetName()+"_rebin1",chi_binnings[mass])
 	    for b in range(hist2.GetNbinsX()):
	       hist2.SetBinContent(b+1,hist2.GetBinContent(b+1)/hist2.GetBinWidth(b+1))
	       hist2.SetBinError(b+1,hist2.GetBinError(b+1)/hist2.GetBinWidth(b+1))
      	hist2.SetLineColor(2)
        hist2.SetMarkerStyle(25)
        hist2.SetMarkerColor(2)
        hist2.SetMarkerSize(0.2)
	miny=0
	if hist2.Integral()>0:
            miny=log*0.1/hist2.Integral()
            hist2.Scale(1./hist2.Integral())
	hist2.SetTitle("")
        hist2.SetStats(False)
        hist2.Draw("pesame")
        legend.AddEntry(hist2,"Data (2017)","lpe")

        name=prefix+data3[0].replace("data_","")+var+str(massbins[mass][0])+"_"+str(massbins[mass][1])
	if var=="#chi": name+="_rebin1"
        hist3=f_data3[0].Get(name)
        for i in range(1,len(data3)):
            hist3.Add(f_data3[i].Get(name))
	if var=="#chi":
            hist3=hist3.Rebin(len(chi_binnings[mass])-1,hist3.GetName()+"_rebin1",chi_binnings[mass])
 	    for b in range(hist3.GetNbinsX()):
	       hist3.SetBinContent(b+1,hist3.GetBinContent(b+1)/hist3.GetBinWidth(b+1))
	       hist3.SetBinError(b+1,hist3.GetBinError(b+1)/hist3.GetBinWidth(b+1))
      	hist3.SetLineColor(4)
        hist3.SetMarkerStyle(25)
        hist3.SetMarkerColor(4)
        hist3.SetMarkerSize(0.2)
	miny=0
	if hist3.Integral()>0:
            miny=log*0.1/hist3.Integral()
	    hist3.Scale(1./hist3.Integral())
	hist3.SetTitle("")
        hist3.SetStats(False)
        hist3.Draw("pesame")
        legend.AddEntry(hist3,"Data (2018)","lpe")

	print "mass bin",mass,"data integral",hist.Integral()

        name=prefix+mc[0][0].replace("-HT200to300","")+var+str(massbins[mass][0])+"_"+str(massbins[mass][1])
	if var=="#chi": name+="_rebin1"
     	hist_mc=f_mc[0].Get(name)
     	for i in range(1,len(mc)):
     	    hist_mc.Add(f_mc[i].Get(name),mc[i][1]/mc[0][1])
	if var=="#chi":
            hist_mc=hist_mc.Rebin(len(chi_binnings[mass])-1,hist_mc.GetName()+"_rebin1",chi_binnings[mass])
 	    for b in range(hist_mc.GetNbinsX()):
	       hist_mc.SetBinContent(b+1,hist_mc.GetBinContent(b+1)/hist_mc.GetBinWidth(b+1))
	       hist_mc.SetBinError(b+1,hist_mc.GetBinError(b+1)/hist_mc.GetBinWidth(b+1))
	if hist_mc.Integral()>0:
            hist_mc.Scale(hist.Integral()/hist_mc.Integral())
      	hist_mc.SetLineColor(1)
        hist_mc.SetStats(False)
        hist_mc.Draw("histsame")
        legend.AddEntry(hist_mc,"MG+Py QCD (2016)","l")

        name=prefix+mc2[0][0].replace("-HT200to300","")+var+str(massbins[mass][0])+"_"+str(massbins[mass][1])
	if var=="#chi": name+="_rebin1"
     	hist_mc2=f_mc2[0].Get(name)
     	for i in range(1,len(mc2)):
     	    hist_mc2.Add(f_mc2[i].Get(name),mc2[i][1]/mc2[0][1])
	if var=="#chi":
            hist_mc2=hist_mc2.Rebin(len(chi_binnings[mass])-1,hist_mc2.GetName()+"_rebin1",chi_binnings[mass])
 	    for b in range(hist_mc2.GetNbinsX()):
	       hist_mc2.SetBinContent(b+1,hist_mc2.GetBinContent(b+1)/hist_mc2.GetBinWidth(b+1))
	       hist_mc2.SetBinError(b+1,hist_mc2.GetBinError(b+1)/hist_mc2.GetBinWidth(b+1))
	if hist_mc2.Integral()>0:
            hist_mc2.Scale(hist2.Integral()/hist_mc2.Integral())
      	hist_mc2.SetLineColor(2)
        hist_mc2.SetStats(False)
        hist_mc2.Draw("histsame")
        legend.AddEntry(hist_mc2,"MG+Py QCD (2017)","l")

        name=prefix+mc3[0][0].replace("-HT200to300","")+var+str(massbins[mass][0])+"_"+str(massbins[mass][1])
	if var=="#chi": name+="_rebin1"
     	hist_mc3=f_mc3[0].Get(name)
     	for i in range(1,len(mc3)):
     	    hist_mc3.Add(f_mc3[i].Get(name),mc3[i][1]/mc3[0][1])
	if var=="#chi":
            hist_mc3=hist_mc3.Rebin(len(chi_binnings[mass])-1,hist_mc3.GetName()+"_rebin1",chi_binnings[mass])
 	    for b in range(hist_mc3.GetNbinsX()):
	       hist_mc3.SetBinContent(b+1,hist_mc3.GetBinContent(b+1)/hist_mc3.GetBinWidth(b+1))
	       hist_mc3.SetBinError(b+1,hist_mc3.GetBinError(b+1)/hist_mc3.GetBinWidth(b+1))
	if hist_mc3.Integral()>0:
            hist_mc3.Scale(hist3.Integral()/hist_mc3.Integral())
      	hist_mc3.SetLineColor(4)
        hist_mc3.SetStats(False)
        hist_mc3.Draw("histsame")
        legend.AddEntry(hist_mc3,"MG+Py QCD (2018)","l")

        hist.Draw("lesame")
        hist2.Draw("lesame")
        hist3.Draw("lesame")
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
     	ratio.GetYaxis().SetTitleSize(0.18)
     	ratio.GetYaxis().SetTitleOffset(0.3)
     	ratio.SetMarkerSize(0.1)
     	ratio.GetYaxis().SetLabelSize(0.2)
     	ratio.GetYaxis().SetRangeUser(0,2)
	if var=="#chi":
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
     	ratio_mc3=hist_mc3.Clone(hist_mc3.GetName()+"ratio")
     	ratio_mc3.Divide(hist_mc3,hist3)
     	for b in range(hist3.GetNbinsX()):
     	  if hist3.GetBinContent(b+1)>0:
     	    ratio_mc3.SetBinError(b+1,hist3.GetBinError(b+1)/hist3.GetBinContent(b+1))
     	ratio_mc3.Draw("histsame")
     	ratio.Draw("axissame")

     	canvas.cd(3)
     	ddratio=hist.Clone(hist.GetName()+"ddratio")
     	ddratio.Divide(hist,hist)
     	for b in range(hist.GetNbinsX()):
     	  if hist.GetBinContent(b+1)>0:
     	    ddratio.SetBinError(b+1,hist.GetBinError(b+1)/hist.GetBinContent(b+1))
     	ddratio.SetTitle("")
     	ddratio.GetYaxis().SetTitle("Data / Data")
     	ddratio.GetYaxis().SetTitleSize(0.11)
     	ddratio.GetYaxis().SetTitleOffset(0.5)
     	ddratio.SetMarkerSize(0.1)
     	ddratio.GetYaxis().SetLabelSize(0.13)
     	ddratio.GetYaxis().SetRangeUser(0,2)
	if var=="#chi":
	  ddratio.GetYaxis().SetRangeUser(0.8,1.2)
     	ddratio.GetXaxis().SetNdivisions(506)
     	ddratio.GetYaxis().SetNdivisions(503)
     	ddratio.GetXaxis().SetLabelColor(1)
     	ddratio.GetXaxis().SetTitle(label[variables.index(var)])
     	ddratio.GetXaxis().SetTitleSize(0.14)
     	ddratio.GetXaxis().SetTitleOffset(1.1)
     	ddratio.GetXaxis().SetLabelSize(0.12)
     	ddratio.Draw("histe")
     	ddratio_mc2=hist_mc2.Clone(hist_mc2.GetName()+"ddratio")
     	ddratio_mc2.Divide(hist_mc2,hist_mc)
     	for b in range(hist_mc.GetNbinsX()):
     	  if hist_mc.GetBinContent(b+1)>0:
     	    ddratio_mc2.SetBinError(b+1,hist_mc.GetBinError(b+1)/hist_mc.GetBinContent(b+1))
     	ddratio_mc2.Draw("histsame")
     	ddratio_mc3=hist_mc3.Clone(hist_mc3.GetName()+"ddratio")
     	ddratio_mc3.Divide(hist_mc3,hist_mc)
     	for b in range(hist_mc.GetNbinsX()):
     	  if hist_mc.GetBinContent(b+1)>0:
     	    ddratio_mc3.SetBinError(b+1,hist_mc.GetBinError(b+1)/hist_mc.GetBinContent(b+1))
     	ddratio_mc3.Draw("histsame")
     	ratio.Draw("axissame")

        canvas.cd(1)
        hist.GetYaxis().SetTitleOffset(1.2)

        canvas.SaveAs("chi_control_plots_"+var.replace("_","").replace("{","").replace("}","").replace("#","").replace("+","p").replace("-","m")+"_"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+postfix+".root")
        canvas.SaveAs("chi_control_plots_"+var.replace("_","").replace("{","").replace("}","").replace("#","").replace("+","p").replace("-","m")+"_"+str(massbins[mass][0])+"_"+str(massbins[mass][1])+postfix+".pdf")
