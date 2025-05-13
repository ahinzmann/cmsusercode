from ROOT import *
import pdb, math, os
gROOT.SetBatch(True)

gROOT.ForceStyle()
gStyle.SetOptStat(0)
gStyle.SetOptFit(0)
gStyle.SetOptTitle(0)
gStyle.SetMarkerSize(0.5)
gStyle.SetHistLineWidth(1)
gStyle.SetMarkerStyle(20)
gStyle.SetLegendBorderSize(0)

from getVariations import getVariations
myVariations=getVariations()

chibins=myVariations.bins
styles=myVariations.styles

def calcHessianPDFUncert(myxsecs):
    myuncertup=0
    myuncertdown=0
    i=0
    xsc=[]
    for myxsec in myxsecs:
        #i+=1
        #if myxsec>0:
        #    myuncertup+=myxsec**2
        #if myxsec<0:
        #    myuncertdown+=myxsec**2
        #print i
        #if i==0: continue
        if i%2==0:
            xsc.append(myxsec)
        if i%2==1:
            xsc.append(myxsec)
            xsc.append(0)
            #print myxsecplus,myxsecminus
            #print xsc[0],xsc[1],xsc[2]
            #print xsc
            myuncertup+=((max(xsc))**2)
            myuncertdown+=((min(xsc))**2)
            xsc=[]
        i+=1
    myfinaluncertup=math.sqrt(myuncertup)
    myfinaluncertdown=(-1)*math.sqrt(myuncertdown)
    return myfinaluncertup,myfinaluncertdown

def calcMCPDFUncert(myxsecs):
    myuncertup=0
    myuncertdown=0
    for myxsec in myxsecs:
        if myxsec>0:
            myuncertup+=myxsec**2
        if myxsec<0:
            myuncertdown+=myxsec**2
    myfinaluncertup=math.sqrt(myuncertup/len(myxsecs))
    myfinaluncertdown=(-1)*math.sqrt(myuncertdown/len(myxsecs))
    return myfinaluncertup,myfinaluncertdown
    

def multiplyBinWidth(myhist):
    chibins=myVariations.bins
    mynewhist=TH1F(myhist.GetName()+"_norm",myhist.GetName()+"_norm",len(chibins)-1,chibins)
    for b in range(0,myhist.GetXaxis().GetNbins()):
        mynewhist.SetBinContent(b+1,myhist.GetBinContent(b+1)*myhist.GetBinWidth(b+1))
    return mynewhist

massbins=[[1200,1500],[1500,1900],[1900,2400],[2400,3000],[3000,3600],[3600,4200],[4200,4800],[4800,5400],[5400,6000],[6000,7000],[7000,13000]]
ciDir="./CIJET_fnl5662j_cs_001_ct14nlo_0_56_/"

order="NNLO"

if order=="NNLO":
  #PDF="ct14nnlo"
  PDF="nn31nnlo"
  mscale="m2"
  #mscale="pt12"
else:
  PDF="ct14nlo"
  mscale="pt12"

if PDF=="nn31nnlo":
  PDFmembers=100
else:
  PDFmembers=56

normalizeBack=False  

calcUncert=False        # calculate qcd+ci uncertainties
compareUncert=True      # compare qcd uncertainty with qcd+ci uncertainty

if calcUncert:
    
    scales=[[1.0,1.0],[0.5,0.5],[2.0,2.0],[0.5,1.0],[1.0,0.5],[1.0,2.0],[2.0,1.0]]
    
    qcdmufile=TFile.Open("2jet.NNLO.fnl5662j_mjj_chi_norm_v25_"+PDF+"_cppread_mu_"+mscale+".root")
    qcdmemfile=TFile.Open("2jet.NNLO.fnl5662j_mjj_chi_norm_v25_"+PDF+"_cppread_mem_"+mscale+".root")
    qcdstatfile=TFile.Open("2jet.NNLO.fnl5662j_mjj_chi_norm_v25_"+PDF+"_cppread_stat_"+mscale+".root")

    cimufiles=[]
    cimemfiles=[]
    for j in styles:
        for i in range(5,31):
            cimufile="CIJET_fnl5662j_cs_001_ct14nlo_0_56_"+str(int(i*1000))+"_"+j+"_mu.root"
            cimemfile="CIJET_fnl5662j_cs_001_ct14nlo_0_56_"+str(int(i*1000))+"_"+j+"_mem.root"
            cimufiles.append(cimufile)
            cimemfiles.append(cimemfile)
    
    scaleFactorsmu={}
    for f in cimufiles:
        print("read",f)
        scaleFactorsmu[f.replace(".root","")]=[]
        file=TFile(ciDir+f)
        newfile=TFile(f.replace("CIJET_","").replace("_mu","").replace("_0_56_","_").replace("ct14nlo",PDF+"_"+mscale).replace("_001_","_"),"RECREATE")
        for massbin in massbins:
            if massbins.index(massbin)<5:
                for scale in scales:
                    qcdhist=qcdmufile.Get("qcd_chi-"+str(massbin[0])+"-"+str(massbin[1])+"scale-"+str(scale[0])+"-"+str(scale[1]))
                    if scales.index(scale)==0:
                        scaleFactor=qcdhist.Integral()
                        scaleFactorsmu[f.replace(".root","")].append(scaleFactor)
                    qcdhistadd=qcdhist.Clone(f.replace(".root","")+"_"+qcdhist.GetName()+"_addmu")
                    qcdhistadd.Scale(1./qcdhistadd.Integral())
                    qcdhistadd.Write()
            else:
                for scale in scales:
                    n=f.replace("mu.root","")+"chi-"+str(massbin[0])+"-"+str(massbin[1])+"scale-"+str(scale[0])+"-"+str(scale[1])
                    cihist=file.Get(n.replace("-7000","-6600")) ## TODO NEW CONTACT CALCULATION
                    qcdhist=qcdmufile.Get("qcd_chi-"+str(massbin[0])+"-"+str(massbin[1])+"scale-"+str(scale[0])+"-"+str(scale[1]))
                    qcdhistadd=qcdhist.Clone(f.replace(".root","")+"_"+qcdhist.GetName()+"_addmu")
                    qcdhistadd.Add(cihist)
                    if scales.index(scale)==0:
                        scaleFactor=qcdhistadd.Integral()
                        scaleFactorsmu[f.replace(".root","")].append(scaleFactor)
                    #pdb.set_trace()
                    qcdhistadd.Scale(1./qcdhistadd.Integral())
                    qcdhistadd.Write()
        newfile.Close()
    
    scaleFactorsmem={}
    for f in cimemfiles:
        scaleFactorsmem[f.replace(".root","")]=[]
        file=TFile(ciDir+f)
        myfile=TFile(f.replace("CIJET_","").replace("_mem","").replace("_0_56_","_").replace("ct14nlo",PDF+"_"+mscale).replace("_001_","_"),"UPDATE")
        for massbin in massbins:
            if massbins.index(massbin)<5:
                for i in range(0,PDFmembers+1):
                    qcdhist=qcdmemfile.Get("qcd_chi-"+str(massbin[0])+"-"+str(massbin[1])+"PDF-"+str(i))
                    if i==0:
                        scaleFactor=qcdhist.Integral()
                        scaleFactorsmem[f.replace(".root","")].append(scaleFactor)
                    qcdhistadd=qcdhist.Clone(f.replace(".root","")+"_"+"qcd_chi-"+str(massbin[0])+"-"+str(massbin[1])+"PDF-"+str(i)+"_addmem")
                    qcdhistadd.Scale(1./qcdhistadd.Integral())
                    qcdhistadd.Write()
            else:
                for i in range(0,PDFmembers+1):
                    n=f.replace("mem.root","")+"chi-"+str(massbin[0])+"-"+str(massbin[1])+"PDF-"+(str(i) if PDFmembers==56 else "0") ## TODO NEW CONTACT CALCULATION WITH NNPDF
                    cihist=file.Get(n.replace("-7000","-6600")) ## TODO NEW CONTACT CALCULATION
                    qcdhist=qcdmemfile.Get("qcd_chi-"+str(massbin[0])+"-"+str(massbin[1])+"PDF-"+str(i))
                    qcdhistadd=qcdhist.Clone(f.replace(".root","")+"_"+"qcd_chi-"+str(massbin[0])+"-"+str(massbin[1])+"PDF-"+str(i)+"_addmem")
                    qcdhistadd.Add(cihist)
                    if i==0:
                        scaleFactor=qcdhistadd.Integral()
                        scaleFactorsmem[f.replace(".root","")].append(scaleFactor)
                    qcdhistadd.Scale(1./qcdhistadd.Integral())
                    qcdhistadd.Write()
                        
    for f in cimufiles:
        myfile=TFile(f.replace("CIJET_","").replace("_mu","").replace("_0_56_","_").replace("ct14nlo",PDF+"_"+mscale).replace("_001_","_"),"UPDATE")
        for massbin in massbins:
            histcentral=TH1F("chi-"+str(massbin[0])+"-"+str(massbin[1]),"chi-"+str(massbin[0])+"-"+str(massbin[1]),len(chibins)-1,chibins)
            histcentral.Sumw2()
            histscaleup=TH1F("chi-"+str(massbin[0])+"-"+str(massbin[1])+"ScaleUp","chi-"+str(massbin[0])+"-"+str(massbin[1])+"scaleUp",len(chibins)-1,chibins)
            histscaleup.Sumw2()
            histscaledown=TH1F("chi-"+str(massbin[0])+"-"+str(massbin[1])+"ScaleDown","chi-"+str(massbin[0])+"-"+str(massbin[1])+"scaleDown",len(chibins)-1,chibins)
            histscaledown.Sumw2()
            for b in range(1,len(chibins)):
                xsec=[]
                i=0
                for scale in scales:
                    myhist=myfile.Get(f.replace(".root","")+"_"+"qcd_chi-"+str(massbin[0])+"-"+str(massbin[1])+"scale-"+str(scale[0])+"-"+str(scale[1])+"_addmu")
                    xsec.append(myhist.GetBinContent(b))
                    i+=1
                    if i==1:
                        central=myhist.GetBinContent(b)
                        histcentral.SetBinContent(b,central)
                    if i==7:
                        scale_up=max(xsec)
                        scale_down=min(xsec)
                        histscaleup.SetBinContent(b,scale_up)
                        histscaledown.SetBinContent(b,scale_down)
    
            if normalizeBack:
                histcentral.Scale(scaleFactorsmu[f.replace(".root","")][massbins.index(massbin)])
                histscaleup.Scale(scaleFactorsmu[f.replace(".root","")][massbins.index(massbin)])
                histscaledown.Scale(scaleFactorsmu[f.replace(".root","")][massbins.index(massbin)])
            histcentral.Write()
            histscaleup.Write()
            histscaledown.Write()
    
    for f in cimemfiles:
        myfile=TFile(f.replace("CIJET_","").replace("_mem","").replace("_0_56_","_").replace("ct14nlo",PDF+"_"+mscale).replace("_001_","_"),"UPDATE")
        print("read",f)
        for massbin in massbins:
            #print massbin[0],massbin[1]
            histcentral=TH1F("chi-"+str(massbin[0])+"-"+str(massbin[1])+"backup","chi-"+str(massbin[0])+"-"+str(massbin[1]),len(chibins)-1,chibins)
            histcentral.Sumw2()
            histscaleup=TH1F("chi-"+str(massbin[0])+"-"+str(massbin[1])+"PDFUp","chi-"+str(massbin[0])+"-"+str(massbin[1])+"PDFUp",len(chibins)-1,chibins)
            histscaleup.Sumw2()
            histscaledown=TH1F("chi-"+str(massbin[0])+"-"+str(massbin[1])+"PDFDown","chi-"+str(massbin[0])+"-"+str(massbin[1])+"PDFDown",len(chibins)-1,chibins)
            histscaledown.Sumw2()
            for b in range(1,len(chibins)):
                #print "Bin number:",b
                xsec=[]
                i=0
                for i in range(0,PDFmembers+1):
                    myhist=myfile.Get(f.replace(".root","")+"_"+"qcd_chi-"+str(massbin[0])+"-"+str(massbin[1])+"PDF-"+str(i)+"_addmem")
                    i+=1
                    if i==1:
                        central=myhist.GetBinContent(b)
                        histcentral.SetBinContent(b,central)
                    else:
                        xsec.append(myhist.GetBinContent(b)-central)
                    if i==PDFmembers+1:
                        #print "central:",central
                        #print "xsec:",xsec
                        if PDFmembers==100: ## MC uncertainties
                          pdf_up,pdf_down=calcMCPDFUncert(xsec)
                        else: ## Hessian uncertainties
                          pdf_up,pdf_down=calcHessianPDFUncert(xsec)
                        #print pdf_up,pdf_down
                        histscaleup.SetBinContent(b,central+pdf_up)
                        histscaledown.SetBinContent(b,central+pdf_down)
    
            if normalizeBack:
                histcentral.Scale(scaleFactorsmem[f.replace(".root","")][massbins.index(massbin)])
                histscaledown.Scale(scaleFactorsmem[f.replace(".root","")][massbins.index(massbin)])
                histscaleup.Scale(scaleFactorsmem[f.replace(".root","")][massbins.index(massbin)])
            histcentral.Write()
            histscaleup.Write()
            histscaledown.Write()

    for f in cimemfiles:
        myfile=TFile(f.replace("CIJET_","").replace("_mem","").replace("_0_56_","_").replace("ct14nlo",PDF+"_"+mscale).replace("_001_","_"),"UPDATE")
        print("read",f)
        for massbin in massbins:
            #print massbin[0],massbin[1]
            histcentral=TH1F("chi-"+str(massbin[0])+"-"+str(massbin[1])+"backup","chi-"+str(massbin[0])+"-"+str(massbin[1]),len(chibins)-1,chibins)
            histcentral.Sumw2()
            histscaleup=TH1F("chi-"+str(massbin[0])+"-"+str(massbin[1])+"StatUp","chi-"+str(massbin[0])+"-"+str(massbin[1])+"StatUp",len(chibins)-1,chibins)
            histscaleup.Sumw2()
            histscaledown=TH1F("chi-"+str(massbin[0])+"-"+str(massbin[1])+"StatDown","chi-"+str(massbin[0])+"-"+str(massbin[1])+"StatDown",len(chibins)-1,chibins)
            histscaledown.Sumw2()
            for b in range(1,len(chibins)):
                #print "Bin number:",b
                xsec=[]
                myhistcentral=myfile.Get(f.replace(".root","")+"_"+"qcd_chi-"+str(massbin[0])+"-"+str(massbin[1])+"PDF-0"+"_addmem")
                central=myhistcentral.GetBinContent(b)
                histcentral.SetBinContent(b,central)
                myhist=qcdstatfile.Get("qcd_chi-"+str(massbin[0])+"-"+str(massbin[1])+"stat")
                histscaleup.SetBinContent(b,central*(1+myhist.GetBinContent(b)))
                histscaledown.SetBinContent(b,central*(1-myhist.GetBinContent(b)))
    
            if normalizeBack:
                histcentral.Scale(scaleFactorsstat[f.replace(".root","")][massbins.index(massbin)])
                histscaledown.Scale(scaleFactorsstat[f.replace(".root","")][massbins.index(massbin)])
                histscaleup.Scale(scaleFactorsstat[f.replace(".root","")][massbins.index(massbin)])
            histcentral.Write()
            histscaleup.Write()
            histscaledown.Write()

        myfile.Close()
    
    
if compareUncert:

    for style in styles:
        for i in range(5,31):
            ciqcdfilename="fnl5662j_cs_ct14nlo_".replace("ct14nlo",PDF+"_"+mscale)+str(i*1000)+"_"+style+".root"

            os.system("cp ../RunII/fnl5662i_v23_fix_CT14_ak4.root .")   # copy qcd file to current directory, this file will be used in the comparison
    
            hist1s=[]
            hist1scaleups=[]
            hist1scaledowns=[]
            hist1pdfups=[]
            hist1pdfdowns=[]
            hist2s=[]
            hist2scaleups=[]
            hist2scaledowns=[]
            hist2pdfups=[]
            hist2pdfdowns=[]
            hist2statups=[]
            hist2statdowns=[]
            hist3s=[]
            hist3scaleups=[]
            hist3scaledowns=[]
            hist3pdfups=[]
            hist3pdfdowns=[]
            hist3statups=[]
            hist3statdowns=[]
    
            file1=TFile("fnl5662i_v23_fix_CT14_ak4.root","UPDATE")
            for massbin in massbins[2:8]:
                histname="chi-"+str(massbin[0])+"-"+str(massbin[1])
                h1=file1.Get(histname)
                h1scaleup=file1.Get(histname+"scaleUp")
                h1scaledown=file1.Get(histname+"scaleDown")
                h1pdfup=file1.Get(histname+"PDFUp")
                h1pdfdown=file1.Get(histname+"PDFDown")
                hist1=multiplyBinWidth(h1)
                hist1scaleup=multiplyBinWidth(h1scaleup)
                hist1scaledown=multiplyBinWidth(h1scaledown)
                hist1pdfup=multiplyBinWidth(h1pdfup)
                hist1pdfdown=multiplyBinWidth(h1pdfdown)
                scaleFactor=hist1.Integral()
                hist1.Scale(1./scaleFactor)
                hist1scaleup.Scale(1./scaleFactor)
                hist1scaledown.Scale(1./scaleFactor)
                hist1pdfup.Scale(1./scaleFactor)
                hist1pdfdown.Scale(1./scaleFactor)
                hist1.Write()
                hist1scaleup.Write()
                hist1scaledown.Write()
                hist1pdfup.Write()
                hist1pdfdown.Write()
            file1.Close()
    
            file1new=TFile("fnl5662i_v23_fix_CT14_ak4.root")
            for massbin in massbins[2:8]:
                histname="chi-"+str(massbin[0])+"-"+str(massbin[1])
                hist1=file1new.Get(histname+"_norm")
                hist1scaleup=file1new.Get(histname+"scaleUp"+"_norm")
                hist1scaleup.Divide(hist1)
                hist1scaledown=file1new.Get(histname+"scaleDown"+"_norm")
                hist1scaledown.Divide(hist1)
                hist1pdfup=file1new.Get(histname+"PDFUp"+"_norm")
                hist1pdfup.Divide(hist1)
                hist1pdfdown=file1new.Get(histname+"PDFDown"+"_norm")
                hist1pdfdown.Divide(hist1)
                hist1s.append(hist1)
                hist1scaleups.append(hist1scaleup)
                hist1scaledowns.append(hist1scaledown)
                hist1pdfups.append(hist1pdfup)
                hist1pdfdowns.append(hist1pdfdown)
    
            file2=TFile.Open(ciqcdfilename.replace("m2","pt12"))
            for massbin in massbins:
                histname="chi-"+str(massbin[0])+"-"+str(massbin[1])
                hist2=file2.Get(histname)
                hist2scaleup=file2.Get(histname+"ScaleUp")
                hist2scaleup.Divide(hist2)
                hist2scaledown=file2.Get(histname+"ScaleDown")
                hist2scaledown.Divide(hist2)
                hist2pdfup=file2.Get(histname+"PDFUp")
                hist2pdfup.Divide(hist2)
                hist2pdfdown=file2.Get(histname+"PDFDown")
                hist2pdfdown.Divide(hist2)
                hist2statup=file2.Get(histname+"StatUp")
                hist2statup.Divide(hist2)
                hist2statdown=file2.Get(histname+"StatDown")
                hist2statdown.Divide(hist2)
                hist2s.append(hist2)
                hist2scaleups.append(hist2scaleup)
                hist2scaledowns.append(hist2scaledown)
                hist2pdfups.append(hist2pdfup)
                hist2pdfdowns.append(hist2pdfdown)
                hist2statups.append(hist2statup)
                hist2statdowns.append(hist2statdown)
        
            file3=TFile.Open(ciqcdfilename.replace("pt12","m2"))
            for massbin in massbins:
                histname="chi-"+str(massbin[0])+"-"+str(massbin[1])
                hist3=file3.Get(histname)
                hist3scaleup=file3.Get(histname+"ScaleUp")
                hist3scaleup.Divide(hist3)
                hist3scaledown=file3.Get(histname+"ScaleDown")
                hist3scaledown.Divide(hist3)
                hist3pdfup=file3.Get(histname+"PDFUp")
                hist3pdfup.Divide(hist3)
                hist3pdfdown=file3.Get(histname+"PDFDown")
                hist3pdfdown.Divide(hist3)
                hist3statup=file3.Get(histname+"StatUp")
                hist3statup.Divide(hist3)
                hist3statdown=file3.Get(histname+"StatDown")
                hist3statdown.Divide(hist3)
                hist3s.append(hist3)
                hist3scaleups.append(hist3scaleup)
                hist3scaledowns.append(hist3scaledown)
                hist3pdfups.append(hist3pdfup)
                hist3pdfdowns.append(hist3pdfdown)
                hist3statups.append(hist3statup)
                hist3statdowns.append(hist3statdown)
        
            canvas=TCanvas("MyCanvas","MyCanvas",0,0,600,300)
            saveName=ciqcdfilename.replace(".root","_xsec.pdf")
            saveName2=ciqcdfilename.replace(".root","_uncert.pdf")
            canvas.cd()
            canvas.Divide(4,2)
            
            legends=[]
            for j in range(3,11):
                
                legend=TLegend(0.2,0.6,0.8,0.9,(str(massbins[j][0])+"<m_{jj}<"+str(massbins[j][1])+" GeV"))
                legend.SetFillStyle(0)
                if j<8:
                  legend.AddEntry(hist1s[j-2],"2016 NLO mu=pTave","l")
                legend.AddEntry(hist2s[j],"Run2 NNLO mu=pTave","l")
                legend.AddEntry(hist3s[j],"Run2 NNLO mu=mjj","l")
                legends.append(legend)
        
            for j in range(3,11):
                canvas.cd(j-2)
                pad1=TPad("","",0, 0, 1, 1)
                pad1.Draw()
                pad1.cd()
                hist2s[j].SetLineColor(2)
                hist2s[j].SetLineStyle(3)
                hist2s[j].SetMaximum(hist2s[j].GetMaximum()*2)
                hist2s[j].GetXaxis().SetTitle("#chi_{dijet}")
                hist2s[j].Draw()
                hist3s[j].SetLineColor(6)
                hist3s[j].Draw("samehist")
                if j<8:
                  hist1s[j-2].SetLineColor(1)
                  hist1s[j-2].SetLineStyle(2)
                  hist1s[j-2].Draw("samehist")
                legends[j-3].Draw()
        
            canvas.SaveAs(saveName)
        
            canvas2=TCanvas("MyCanvas","MyCanvas",0,0,600,300)
            canvas2.cd()
            canvas2.Divide(4,2)
            
            legends2=[]
            for j in range(3,11):
                
                legend=TLegend(0.2,0.6,0.8,0.9,(str(massbins[j][0])+"<m_{jj}<"+str(massbins[j][1])+" GeV"))
                legend.SetFillStyle(0)
                if j<8:
                  legend.AddEntry(hist1pdfups[j-2],"2016 NLO mu=pTave PDF unc.","l")
                  legend.AddEntry(hist1scaleups[j-2],"2016 NLO mu=pTave scale unc.","l")
                legend.AddEntry(hist2pdfups[j],"Run2 NNLO mu=pTave PDF unc.","l")
                legend.AddEntry(hist2scaleups[j],"Run2 NNLO mu=pTave scale unc.","l")
                legend.AddEntry(hist2statups[j],"Run2 NNLO mu=pTave stat unc.","l")
                legend.AddEntry(hist3pdfups[j],"Run2 NNLO mu=mjj PDF unc.","l")
                legend.AddEntry(hist3scaleups[j],"Run2 NNLO mu=mjj scale unc.","l")
                legend.AddEntry(hist3statups[j],"Run2 NNLO mu=mjj stat unc.","l")
                legends2.append(legend)
        
            for j in range(3,11):
                canvas2.cd(j-2)
                pad2=TPad("","",0, 0, 1, 1)
                pad2.Draw()
                pad2.cd()
                pad2.SetLeftMargin(0.05)
                pad2.SetBottomMargin(0.08)
                pad2.SetTopMargin(0.03)
                pad2.SetRightMargin(0.05)
                gPad.SetTicky(2)

                hist2pdfups[j].SetLineColor(1)
                hist2pdfups[j].SetLineStyle(3)
                hist2pdfups[j].GetXaxis().SetTitle("#chi_{dijet}")
                hist2pdfups[j].SetMinimum(0.9)
                hist2pdfups[j].SetMaximum(1.2)
                hist2pdfups[j].Draw("same")
                hist2pdfdowns[j].SetLineColor(1)
                hist2pdfdowns[j].SetLineStyle(3)
                hist2pdfdowns[j].Draw("samehist")
                
                hist2scaleups[j].SetLineColor(2)
                hist2scaleups[j].SetLineStyle(3)
                hist2scaleups[j].GetXaxis().SetTitle("#chi_{dijet}")
                hist2scaleups[j].Draw("samehist")
                hist2scaledowns[j].SetLineColor(2)
                hist2scaledowns[j].SetLineStyle(3)
                hist2scaledowns[j].Draw("samehist")
                
                hist2statups[j].SetLineColor(8)
                hist2statups[j].SetLineStyle(3)
                hist2statups[j].GetXaxis().SetTitle("#chi_{dijet}")
                hist2statups[j].Draw("same")
                hist2statdowns[j].SetLineColor(8)
                hist2statdowns[j].SetLineStyle(3)
                hist2statdowns[j].Draw("samehist")
                
                if j<8:
                  hist1pdfups[j-2].SetLineColor(4)
                  hist1pdfups[j-2].SetLineStyle(2)
                  hist1pdfups[j-2].GetXaxis().SetTitle("#chi_{dijet}")
                  hist1pdfups[j-2].Draw("histsame")
                  hist1pdfdowns[j-2].SetLineColor(4)
                  hist1pdfdowns[j-2].SetLineStyle(2)
                  hist1pdfdowns[j-2].Draw("samehist")
                  
                  hist1scaleups[j-2].SetLineColor(6)
                  hist1scaleups[j-2].SetLineStyle(2)
                  hist1scaleups[j-2].GetXaxis().SetTitle("#chi_{dijet}")
                  hist1scaleups[j-2].Draw("samehist")
                  hist1scaledowns[j-2].SetLineColor(6)
                  hist1scaledowns[j-2].SetLineStyle(2)
                  hist1scaledowns[j-2].Draw("samehist")
                
                hist3pdfups[j].SetLineColor(7)
                hist3pdfups[j].SetLineStyle(1)
                hist3pdfups[j].GetXaxis().SetTitle("#chi_{dijet}")
                hist3pdfups[j].Draw("samehist")
                hist3pdfdowns[j].SetLineColor(7)
                hist3pdfdowns[j].SetLineStyle(1)
                hist3pdfdowns[j].Draw("samehist")
                
                hist3scaleups[j].SetLineColor(3)
                hist3scaleups[j].SetLineStyle(1)
                hist3scaleups[j].GetXaxis().SetTitle("#chi_{dijet}")
                hist3scaleups[j].Draw("samehist")
                hist3scaledowns[j].SetLineColor(3)
                hist3scaledowns[j].SetLineStyle(1)
                hist3scaledowns[j].Draw("samehist")
                
                hist3statups[j].SetLineColor(9)
                hist3statups[j].SetLineStyle(1)
                hist3statups[j].GetXaxis().SetTitle("#chi_{dijet}")
                hist3statups[j].Draw("samehist")
                hist3statdowns[j].SetLineColor(9)
                hist3statdowns[j].SetLineStyle(1)
                hist3statdowns[j].Draw("samehist")
                
                if j<8:
                  hist1s[j-2].SetLineColor(1)
                  hist1s[j-2].Draw("samehist")
                legends2[j-3].Draw()
        
            canvas2.SaveAs(saveName2)
