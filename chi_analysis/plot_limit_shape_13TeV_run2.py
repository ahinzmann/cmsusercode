from ROOT import *
import ROOT
from math import log10,sqrt

#gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
#gROOT.Reset()
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

if __name__=="__main__":

 models=[3] #ADD
 #models+=[10] #QBH ADD6
 #models+=[11] #QBH RS1
 #models+=[60,61,62,63,64,65,66,67,68,69] #CI
 #models+=[70,71,72,73,74,75,76,77]
 #models+=[78,79,80,81,82,83,84,85,86,87]
 #models+=[30,31,32,33,34,35,36,37,38]
 #models+=[40,41,42,43,44,45,46,47,48]
 #models=[88,89]
 #models=[60,61]
 #models=[10,11]
 models+=[90,91,92,93,94] #alp
 models+=[95,96,97,98,99] #tripleG
 #models=[90,91]

 testStat="LHC"
 asym="a" #asymptotic CLS
 runs="2" # "2" or "3" or "23"

 limit_list={}

 for model in models:

    if model==1:
       signal="CIplusLL"    
    if model==2:
       signal="CIminusLL"    
    if model==3:
       signal="ADD"    
    if model==4:
       signal="cs_nn30nlo_0_"    
    if model==5:
       signal="cs_nn30nlo_0_"    
    if model==6:
       signal="cs_nn30nlo_0_"    
    if model==7:
       signal="cs_nn30nlo_0_"    
    if model==8:
       signal="cs_nn30nlo_0_"    
    if model==9:
       signal="cs_nn30nlo_0_"    
    if model==10:
       signal="QBH_"    
    if model==11:
       signal="QBH_"    

    if model>=18 and model<20:
       signal="cs_ct14nlo_"
    if model>=20 and model<30:
       signal="cs_nn30nlo_0_"
    if model>=30 and model<40:
       signal="cs_ct14nlo_"
    if model>=40 and model<50:
       signal="AntiCIplusLL"
    if model>=60 and model<70:
       signal="cs_ct14nlo_"
    if model>=70 and model<90:
       signal="cs_ct14nlo_"
    if model>=90 and model<95:
       signal="alp_QCD_fa"
    if model>=95 and model<100:
       signal="tripleG_QCD_CG"

    print signal,model

    f=file("limits"+testStat+asym+str(model)+"_"+signal+"_run"+runs+".txt")
    limits=eval(f.readline())
    #print limits

    canvas = TCanvas("","",0,0,300,300)
    #canvas.GetPad(0).SetLogy()
    mg=TMultiGraph()

    min_x=5000
    max_x=40000
    if "alp" in signal:
      min_x=1000
      max_x=5000
    if "tripleG" in signal:
      min_x=1000
      max_x=40000
    g0=TGraph(0)
    g0.SetPoint(0,min_x,0)
    g0.SetPoint(1,max_x,0)
    mg.Add(g0)
    
    g=TGraph(0)
    g_exp=TGraph(0)
    g_exp1m=TGraph(0)
    g_exp1p=TGraph(0)
    g_exp2m=TGraph(0)
    g_exp2p=TGraph(0)
    for mass,limit,error,exp,exp1m,exp1p,exp2m,exp2p in limits:
      if "tripleG" in signal: mass=1000./sqrt(mass)
      if limit==0: limit=1e-5
      if not testStat=="LHC":
       if exp>=1: exp=1e-5
       if exp1m>=1: exp1m=1e-5
       if exp1p>=1: exp1p=1e-5
       if exp2m>=1: exp2m=1e-5
       if exp2p>=1: exp2p=1e-5
      if limit>0:
        g.SetPoint(g.GetN(),mass,log10(limit))
      if exp>0:
        g_exp.SetPoint(g_exp.GetN(),mass,log10(exp))
      if exp1m>0:
        g_exp1m.SetPoint(g_exp1m.GetN(),mass,log10(exp1m))
      if exp1p>0:
        g_exp1p.SetPoint(g_exp1p.GetN(),mass,log10(exp1p))
      if exp2m>0:
        g_exp2m.SetPoint(g_exp2m.GetN(),mass,log10(exp2m))
      if exp2p>0:
        g_exp2p.SetPoint(g_exp2p.GetN(),mass,log10(exp2p))
    g.SetMarkerStyle(24)
    g.SetMarkerSize(0.5)
    g.SetLineColor(1)
    g.SetLineWidth(3)
    mg.Add(g)
    g_exp.SetMarkerStyle(24)
    g_exp.SetMarkerSize(0.5)
    g_exp.SetLineColor(2)
    g_exp.SetLineWidth(3)
    mg.Add(g_exp)
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
    mg.GetXaxis().SetTitle("scale [GeV]")
    if testStat=="LHC":
      mg.GetYaxis().SetTitle("log_{10}(signal strength)")
      miny=-1
      maxy=1
    else:
      mg.GetYaxis().SetTitle("log_{10}(CL_{S})")
      miny=-3
      maxy=0
    mg.GetYaxis().SetRangeUser(miny,maxy)
    
    if testStat=="LHC":
      cut=1
    else:
      cut=0.05
    l=TLine(min_x,log10(cut),max_x,log10(cut))
    l.SetLineColor(2)
    l.SetLineStyle(2)
    l.Draw("same")
    if testStat=="LHC":
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
    for i in range(max_x):
        mass=i*(max_x-limits[0][0])/max_x+limits[0][0]
        if mass<min_x or mass>max_x: continue
	if limit==0 and g.Eval(mass,0)>log10(cut):
	    limit=mass
	if exp==0 and g_exp.Eval(mass,0)>log10(cut):
	    exp=mass
	if exp1m==0 and g_exp1m.Eval(mass,0)>log10(cut):
	    exp1m=mass
	if exp1p==0 and g_exp1p.Eval(mass,0)>log10(cut):
	    exp1p=mass
	if exp2m==0 and g_exp2m.Eval(mass,0)>log10(cut):
	    exp2m=mass
	if exp2p==0 and g_exp2p.Eval(mass,0)>log10(cut):
	    exp2p=mass
    err=0
    if exp1m>0: err=exp-exp1m
    if exp1p>0 and exp1p-exp>exp-exp1m: err=exp1p-exp
     
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
    
    #if "tripleG" in signal:
    #  err=0
    #  if exp1m>0: err=100./exp-100./exp1m
    #  if exp1p>0 and 100./exp1p-100./exp>100./exp-100./exp1m: err=100./exp1p-100./exp
    #  print "limit: %.4f" % (100./max(1e-5,limit)), "& %.4f" % (100./max(1e-5,exp)), "$\pm$ %.4f" % (-err)
    #  print "limit: %.4f," % (100./max(1e-5,limit)), "%.4f," % (100./max(1e-5,exp)), "%.4f, %.4f, 0, 0" % ((100./max(1e-5,exp1m),(100./max(1e-5,exp1p))))
    #else:
    print "limit: %.1f" % (limit/1000.), "& %.1f" % (exp/1000.), "$\pm$ %.1f" % (err/1000.)
    print "limit: %.2f," % (limit/1000.), "%.2f," % (exp/1000.), "%.2f, %.2f, %.2f, %.2f" % ((exp1m)/1000.,(exp1p)/1000.,(exp2m)/1000.,(exp2p)/1000.)
    limit_list[model]=(limit/1000.,(exp/1000.),(exp1m)/1000.,(exp1p)/1000.,(exp2m)/1000.,(exp2p)/1000.)
    
    canvas.SaveAs('limits'+testStat+asym+str(model)+signal+'_run'+runs+'.pdf')
    #canvas.SaveAs('limits'+testStat+asym+str(model)+signal+'_run2.eps')
    
    
    
    
    
 for coupling in ["alp","tripleG"]:

    canvas = TCanvas("","",0,0,1800,1550)
    #canvas.GetPad(0).SetLogy()
    mg=TMultiGraph()

    points=[]
    if coupling=="alp":
      ymin=0
      ymax=4
      min_x=3.5
      max_x=10
      points=[(3.6,91),(4.2,92),(4.8,93),(5.4,94),(6.0,94),(7.0,94),(8.0,94),(9.0,94),(10.0,94),(12.0,94),(14.0,94),(16.0,94),(18.0,94),(20.0,94)] #(3.0,90),
    if coupling=="tripleG":
      ymin=0
      ymax=1.8
      min_x=3.5
      max_x=14
      points=[(3.6,96),(4.2,97),(4.8,98),(5.4,99),(6.0,99),(7.0,99),(8.0,99),(9.0,99),(10.0,99),(12.0,99),(14.0,99),(16.0,99),(18.0,99),(20.0,99)] #(3.0,95),
    g0=TGraph(0)
    g0.SetPoint(0,min_x,0)
    g0.SetPoint(1,max_x,0)
    mg.Add(g0)
    
    g=TGraph(0)
    g_exp=TGraph(0)
    g_band=TGraphAsymmErrors(0)
    g_band_2sigma=TGraphAsymmErrors(0)
    g_val=TGraph(0)
    g_val.SetPoint(0,3.6,0)
    g_val.SetPoint(1,3.6,10)
    
    for mass,model in points:
        obs=limit_list[model][0]
        exp=limit_list[model][1]
        exp1m=limit_list[model][2]
        exp1p=limit_list[model][3]
        exp2m=limit_list[model][4]
        exp2p=limit_list[model][5]
        if coupling=="alp":
          g.SetPoint(g.GetN(),mass,mass/obs)
          g_exp.SetPoint(g_exp.GetN(),mass,mass/exp)
          g_band.SetPoint(g_band.GetN(),mass,mass/exp)
          g_band_2sigma.SetPoint(g_band_2sigma.GetN(),mass,mass/exp)
          g_band.SetPointError(g_band.GetN()-1,0,0,mass/exp-min(mass/exp1p,mass/exp1m),max(mass/exp1p,mass/exp1m)-mass/exp)
          g_band_2sigma.SetPointError(g_band_2sigma.GetN()-1,0,0,mass/exp-min(mass/exp2p,mass/exp2m),max(mass/exp2p,mass/exp2m)-mass/exp)
        if coupling=="tripleG":
          g.SetPoint(g.GetN(),mass,pow(mass/obs,2))
          g_exp.SetPoint(g_exp.GetN(),mass,pow(mass/exp,2))
          g_band.SetPoint(g_band.GetN(),mass,pow(mass/exp,2))
          g_band_2sigma.SetPoint(g_band_2sigma.GetN(),mass,pow(mass/exp,2))
          g_band.SetPointError(g_band.GetN()-1,0,0,pow(mass/exp,2)-min(pow(mass/exp1p,2),pow(mass/exp1m,2)),max(pow(mass/exp1p,2),pow(mass/exp1m,2))-pow(mass/exp,2))
          g_band_2sigma.SetPointError(g_band_2sigma.GetN()-1,0,0,pow(mass/exp,2)-min(pow(mass/exp2p,2),pow(mass/exp2m,2)),max(pow(mass/exp2p,2),pow(mass/exp2m,2))-pow(mass/exp,2))
    g_band_2sigma.SetFillStyle(1001)
    g_band_2sigma.SetFillColor(kOrange)
    g_band_2sigma.SetLineColor(1)
    g_band_2sigma.SetLineStyle(1)
    g_band_2sigma.SetLineWidth(0)
    mg.Add(g_band_2sigma,"3")
    g_band.SetFillStyle(1001)
    g_band.SetFillColor(kGreen-3)
    g_band.SetLineColor(1)
    g_band.SetLineStyle(1)
    g_band.SetLineWidth(0)
    mg.Add(g_band,"3")
    g_exp.SetLineColor(1)
    g_exp.SetLineWidth(3)
    g_exp.SetLineStyle(kDashed)
    mg.Add(g_exp,"l")
    g.SetMarkerStyle(24)
    g.SetMarkerSize(0)
    g.SetLineColor(1)
    g.SetLineWidth(3)
    #mg.Add(g,"pl")
    g_val.SetLineColor(2)
    g_val.SetLineWidth(303)
    g_val.SetFillColor(2)
    g_val.SetFillStyle(3004)
    mg.Add(g_val,"l3")
    
    mg.Draw("apl")
    mg.SetTitle("")
    if coupling=="alp":
      mg.GetXaxis().SetTitle("f_{a} [TeV]")
      mg.GetYaxis().SetTitle("c_{g}")
    if coupling=="tripleG":
      mg.GetXaxis().SetTitle("#Lambda [TeV]")
      mg.GetYaxis().SetTitle("C_{G}")
    mg.GetYaxis().SetRangeUser(ymin,ymax)
    mg.GetYaxis().SetLabelSize(0.04)
    mg.GetYaxis().SetTitleSize(0.05)
    mg.GetYaxis().SetTitleOffset(0.95)
    mg.GetXaxis().SetLimits(min_x,max_x)
    #mg.GetXaxis().SetRangeUser(min_x_new,max_x_new)
    mg.GetYaxis().SetNdivisions(510)
    mg.GetXaxis().SetNdivisions(510)
    mg.GetXaxis().SetLabelSize(0.04)
    mg.GetXaxis().SetLabelOffset(0.015)
    mg.GetXaxis().SetTitleSize(0.05)
    mg.GetXaxis().SetTitleOffset(1.1)
   
    if coupling=="alp":
      lt=TLatex(4.5,3.6,"#splitline{#bf{ALP linear EFT}}{#bf{m_{a} = 1 MeV}}")
    if coupling=="tripleG":
      lt=TLatex(10,0.09,"#splitline{#bf{SMEFT}}{}")
    lt.SetTextSize(0.04)
    lt.SetTextFont(42)
    lt.Draw("same")

    #l=TLegend(0.13,0.5,0.43,0.72,"95% CL upper limits")
    if coupling=="alp":
      l=TLegend(0.5,0.25,0.8,0.54,"95% CL upper limits")
    if coupling=="tripleG":
      l=TLegend(0.2,0.61,0.5,0.90,"95% CL upper limits")
    l.SetFillColor(0)
    l.SetTextFont(42)
    l.SetFillStyle(0)
    l.SetTextSize(0.04)
    l.SetShadowColor(0)
    #l.AddEntry(g,"Observed","LP")
    l.AddEntry(g_exp,"Expected","LP")
    l.AddEntry(g_band,"Expected #pm 1 s.d.","F")
    l.AddEntry(g_band_2sigma,"Expected #pm 2 s.d.","F")
    l.AddEntry(g_val,"Validity limit of EFT","LF")
    l.Draw()

    # CMS
    leg2=TLatex(min_x,ymax+0.03,"#bf{CMS} #it{Preliminary}")
    #cmsPos=min_x+220/1000.
    #leg2=TLatex(cmsPos,ymax-0.17,"#bf{CMS}")
    leg2.SetTextFont(42)
    leg2.SetTextSize(0.05)
    # lumi
    lumiPos=max_x*(0.75 if coupling=="alp" else 0.7)
    leg3=TLatex(lumiPos,ymax+0.03,"138 fb^{-1} (13 TeV)")
    leg3.SetTextFont(42)
    leg3.SetTextSize(0.045)
    leg2.Draw("same")
    leg3.Draw("same")
    
    canvas.SaveAs('limits'+testStat+asym+coupling+'_coupling_run'+runs+'.pdf')
