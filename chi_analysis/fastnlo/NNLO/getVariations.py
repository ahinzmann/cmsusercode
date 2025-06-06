from ROOT import *
from array import *
import subprocess,os

order="NNLO"
newci=True

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

if newci:
  ciscales=[10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30] # TeV
  massbins=[[5400,6000],[6000,7000],[7000,13000]]
else:
  ciscales=range(5,31) # TeV
  massbins=[[1200,1500],[1500,1900],[1900,2400],[2400,3000],[3000,3600],[3600,4200],[4200,4800],[4800,5400],[5400,6000],[6000,7000],[7000,13000]]

class getVariations:  # This class helps extract xsecs from NLOJET++ and CIJET++ program and store the xsecs into root histograms or python lists/dictionaries 
    
    def __init__(self):
        # Use version="i"for 2015 analysis and use version="j" for 2016 analysis
        # Default path for NLOJET++ cross section input (.log) files is ./
        # Defualt path for NLOJET++ cross section output files is ./
        # Default path for CIJET++ cross section input (.xsc) files is ./cixsecDir+style
        # Default path for CIJET++ cross section output files is ./cixsecDir
        self.version="j"
        self.BaseDir=os.getcwd()
        # other initiation
        self.bins=array('d',[1,2,3,4,5,6,7,8,9,10,12,14,16])
        self.styles=["LL-","LL+","RR-","RR+","VV-","VV+","AA-","AA+","V-A-","V-A+"]
        #self.styles=["AA-","AA+"]
        self.cixsecDir=("NEWCI" if newci else "CIJET_fnl5662"+self.version+"_cs_001_"+PDF+"_0_"+str(PDFmembers)+"_")

    def getqcdallmu(self):
      allmulist=[]
      for m in range(len(massbins)):
        #fileallmu=self.BaseDir+"/2jet.NNLO.fnl5662"+self.version+"_mjj"+str(m).replace("10","a")+"_chi_"+PDF+"_cppread_7_"+mscale+".log"
        fileallmu=self.BaseDir+"/2jet.NNLO.fnl5662"+self.version+"_mjj"+str(m).replace("10","a")+"_chi_norm_v25_"+PDF+"_allmu_"+mscale.replace("m2","m12")+".log"
        print("load",fileallmu)
        mu=["start"]
        muxsecs=["start"]
        with open(fileallmu) as f:
            for line in f:
                if "The scale factors xmur, xmuf chosen here are" in line:
                    if mu != ["start"] and muxsecs !=["start"]:
                        allmulist.append([mu,muxsecs])
                    mu=[float(line.split()[9].replace(',','')),float(line.split()[10])]
                    muxsecs=[]
                if "#" in line: continue
                if line.split()==[]: continue
                if line.split()[0]=="LHAPDF": continue
                if line.split()[0]=="CT14nlo": continue
                if line.split()[0]=="NNPDF31_nnlo_as_0118": continue
                if order=="NNLO":
                  x=float(line.split()[8])#/(massbins[m][1]-massbins[m][0])
                else:
                  x=float(line.split()[7])#/(massbins[m][1]-massbins[m][0])
                xsec=[massbins[m][0],massbins[m][1],float(line.split()[3]),float(line.split()[4]),x]
                muxsecs.append(xsec)
            allmulist.append([mu,muxsecs])
        print("found",len(allmulist),"elements")
      return allmulist

    def getqcdmup(self,points):
      muplist=[]
      for m in range(len(massbins)):
        filemup=self.BaseDir+"/2jet.NNLO.fnl5662"+self.version+"_mjj"+str(m).replace("10","a")+"_chi_norm_v25_"+PDF+"_"+points+"-norm_"+mscale.replace("m2","m12")+".log"
        print("load",filemup)
        muxsecscentral=[]
        muxsecsdown=[]
        muxsecsup=[]
        with open(filemup) as f:
            for line in f:
                if "#" in line: continue
                if line.split()==[]: continue
                if line.split()[0]=="LHAPDF": continue
                if line.split()[0]=="CT14nlo": continue
                if line.split()[0]=="NNPDF31_nnlo_as_0118": continue
                x=float(line.split()[1])/(massbins[m][1]-massbins[m][0])/(self.bins[int(line.split()[0])]-self.bins[int(line.split()[0])-1])
                xsec=[massbins[m][0],massbins[m][1],self.bins[int(line.split()[0])-1],self.bins[int(line.split()[0])],x]
                muxsecscentral.append(xsec)
                x=float(line.split()[1])*(1.+float(line.split()[2]))/(massbins[m][1]-massbins[m][0])/(self.bins[int(line.split()[0])]-self.bins[int(line.split()[0])-1])
                xsec=[massbins[m][0],massbins[m][1],self.bins[int(line.split()[0])-1],self.bins[int(line.split()[0])],x]
                muxsecsdown.append(xsec)
                x=float(line.split()[1])*(1.+float(line.split()[3]))/(massbins[m][1]-massbins[m][0])/(self.bins[int(line.split()[0])]-self.bins[int(line.split()[0])-1])
                xsec=[massbins[m][0],massbins[m][1],self.bins[int(line.split()[0])-1],self.bins[int(line.split()[0])],x]
                muxsecsup.append(xsec)
            muplist.append([[points,"central"],muxsecscentral])
            muplist.append([[points,"down"],muxsecsdown])
            muplist.append([[points,"up"],muxsecsup])
        print("found",len(muplist),"elements")
        print(muplist[0])
      return muplist

    def getqcdallmem(self):
      allmemlist=[]
      for m in range(len(massbins)):
        #fileallmem=self.BaseDir+"/2jet.NNLO.fnl5662"+self.version+"_mjj"+str(m).replace("10","a")+"_chi_"+PDF+"_cppread_0_"+mscale+".log"
        fileallmem=self.BaseDir+"/2jet.NNLO.fnl5662"+self.version+"_mjj"+str(m).replace("10","a")+"_chi_norm_v25_"+PDF+"_allmem_"+mscale.replace("m2","m12")+".log"
        print("load",fileallmem)
        mem=["start"]
        memxsecs=["start"]
        with open(fileallmem) as f:
            for line in f:
                if "The PDF member chosen here is" in line:
                    if mem != ["start"] and memxsecs !=["start"]:
                        allmemlist.append([mem,memxsecs])
                    mem=[int(line.split()[7])]
                    memxsecs=[]
                if "#" in line: continue
                if line.split()==[]: continue
                if line.split()[0]=="LHAPDF": continue
                if line.split()[0]=="CT14nlo": continue
                if order=="NNLO":
                  x=float(line.split()[8])#/(massbins[m][1]-massbins[m][0])
                else:
                  x=float(line.split()[7])#/(massbins[m][1]-massbins[m][0])
                xsec=[massbins[m][0],massbins[m][1],float(line.split()[3]),float(line.split()[4]),x]
                memxsecs.append(xsec)
            allmemlist.append([mem,memxsecs])
        print("found",len(allmemlist),"elements")
      return allmemlist

    def getqcdstat(self):
      statlist=[]
      for m in range(len(massbins)):
        filestat=self.BaseDir+"/2jet.NNLO.fnl5662"+self.version+"_mjj"+str(m).replace("10","a")+"_chi_norm_v25_"+PDF+"_stat_"+mscale.replace("m2","m12")+".log"
        print("load",filestat)
        stat=["start"]
        statxsecs=["start"]
        with open(filestat) as f:
            for line in f:
                if "Relative statistical uncertainties" in line:
                    if stat != ["start"] and statxsecs !=["start"]:
                        statlist.append([stat,statxsecs])
                    stat=["stat"]
                    statxsecs=[]
                if "#" in line: continue
                if line.split()==[]: continue
                if line.split()[0]=="LHAPDF": continue
                if line.split()[0]=="CT14nlo": continue
                if order=="NNLO":
                  x=float(line.split()[3])/(massbins[m][1]-massbins[m][0])/(self.bins[int(line.split()[0])]-self.bins[int(line.split()[0])-1])
                else:
                  print("WARNING: No log file for NLO stat uncertainty available, only for NNLO")
                  x=0
                xsec=[massbins[m][0],massbins[m][1],self.bins[int(line.split()[0])-1],self.bins[int(line.split()[0])],x]
                statxsecs.append(xsec)
            statlist.append([stat,statxsecs])
        print("found",len(statlist),"elements")
      return statlist

    def getciallmu(self):
##         ls = subprocess.Popen(['ls',self.BaseDir+'/'+self.cixsecDir],stdout=subprocess.PIPE)    
##         stdouts=[]
##         while True:
##             line = ls.stdout.readline()
##             stdouts.append(line)
##             if line == '' and ls.poll() != None:
##                 break        
##         files=[]
##         for stdout in stdouts:
##             if stdout=='': continue
##             if "_0_" not in stdout: continue
##             if ".root" in stdout: continue
##             file=stdout.replace("\n","")
##             files.append(file)
##     
##         cimudict={}
##         for file in files:
##             cimudict[file.replace('.xsc','')]=[]
##             with open(self.BaseDir+"/"+self.cixsecDir+"/"+file) as f:
##                 i=-9999
        cimudict={}
        for style in self.styles:
            for m in ciscales:
                name="CIJET_fnl5662"+self.version+"_cs_001_"+PDF+"_0_"+str(PDFmembers)+"_"+str(m*1000)+"_"+style
                cimudict[name]=[]
##                 for j in range(0,PDFmembers+1):
                if newci:
                  filename="CIJET_cs_fnl5662"+self.version+"_mjj8-a_"+PDF.replace("nlo","-nlo").replace("n-nlo","-nnlo-as-0118")+"_0_"+str(m*1000)+"_"+style+".xsc"
                else:
                  filename=name.replace("_0_"+str(PDFmembers)+"_","_0_")+".xsc"
                muxsecs=[]
                print(filename)
                with open(self.BaseDir+"/"+self.cixsecDir+("" if newci else style)+"/"+filename) as f:
                    i=-9999
                    for line in f:
                        i+=1
                        if "CT14nlo" in line: continue
                        if "Lambda" in line: continue
                        if "cj/4pi" in line: continue
                        if "Bin" in line:
                            i=0
                        if i==1:
                            mass_low=float(line.split()[2])
                            mass_high=float(line.split()[3])
                            chi_low=float(line.split()[0])
                            chi_high=float(line.split()[1])
                            #print(mass_low)
                        if i>2:
                            xsec=float(line.split()[3])
                            mur=float(line.split()[1])
                            muf=float(line.split()[0])
                            #print cimudict
                            if len(cimudict[name])<9:
                                mu=[mur,muf]
                                xsecs=[mass_low,mass_high,chi_low,chi_high,xsec]
                                muxsecs=[mu,xsecs]
                                cimudict[name].append(muxsecs)    
                            else:
                                xsecs=[mass_low,mass_high,chi_low,chi_high,xsec]
                                cimudict[name][i-3].append(xsecs)
        return cimudict

    def getciallmem(self):
        cimemdict={}
        for style in self.styles:
            for m in ciscales:
                name="CIJET_fnl5662"+self.version+"_cs_001_"+PDF+"_0_"+str(PDFmembers)+"_"+str(m*1000)+"_"+style
                cimemdict[name]=[]
                for j in range(0,PDFmembers+1):
                    if newci:
                      filename="CIJET_cs_fnl5662"+self.version+"_mjj8-a_"+PDF.replace("nlo","-nlo").replace("n-nlo","-nnlo-as-0118")+"_"+str(j)+"_"+str(m*1000)+"_"+style+".xsc"
                    else:
                      filename=name.replace('_0_'+str(PDFmembers)+'_','_'+str(j)+'_')+".xsc"
                    memxsecs=[]
                    print(filename)
                    with open(self.BaseDir+"/"+self.cixsecDir+("" if newci else style)+"/"+filename) as f:
                        #print f
                        mem=j
                        memxsecs.append(mem)
                        i=-9999
                        for line in f:
                            #print(line)
                            i+=1
                            if "CT14nlo" in line: continue
                            if "xsec" in line: continue
                            if "Lambda" in line: continue
                            if "cj/4pi" in line: continue
                            if "Bin" in line:
                                i=0
                            if i==1:
                                mass_low=float(line.split()[2])
                                mass_high=float(line.split()[3])
                                chi_low=float(line.split()[0])
                                chi_high=float(line.split()[1])
                                #print(mass_low)
                            if i>2:
                                xsec=float(line.split()[3])
                                mur=float(line.split()[1])
                                muf=float(line.split()[0])
                                if mur==1 and muf==1:
                                    xsecs=[mass_low,mass_high,chi_low,chi_high,xsec]
                                    memxsecs.append(xsecs)
                        cimemdict[name].append(memxsecs)
        return cimemdict
    
    def dictPrint(self,mydict):
        for key in mydict:
            print(key)
            print(mydict[key])

    def listPrint(self,mylist):
        for list in mylist:
            print(list)
            
    def listFill(self,mylist,uncerttype):
        myfile=TFile(self.BaseDir+"/2jet.NNLO.fnl5662"+self.version+"_mjj_chi_norm_v25_"+PDF+"_cppread_"+uncerttype+"_"+mscale+".root","RECREATE")
        for list in mylist:
            if uncerttype=="mu" or uncerttype=="mu6" or uncerttype=="mu30":
                mur=str(list[0][0])
                muf=str(list[0][1])
                i=0
                for point in list[1]:
                    #print point
                    if i==0:
                        hist=TH1F("qcd_chi-"+str(int(point[0]))+"-"+str(int(point[1]))+"scale-"+mur+"-"+muf,"chi-"+str(int(point[0]))+"-"+str(int(point[1]))+"scale-"+mur+"-"+muf,len(self.bins)-1,self.bins)
                        hist.Sumw2()
                    hist.Fill(point[2]+(point[3]-point[2])/2,point[4]*(point[3]-point[2])*(point[1]-point[0]))
                    i+=1
                    if i==12:
                        hist.Write()
                        i=0
            elif uncerttype=="mem":
                member=str(list[0][0])
                i=0
                for point in list[1]:
                    if i==0:
                        hist=TH1F("qcd_chi-"+str(int(point[0]))+"-"+str(int(point[1]))+"PDF-"+member,"chi-"+str(int(point[0]))+"-"+str(int(point[1]))+"PDF-"+member,len(self.bins)-1,self.bins)
                        hist.Sumw2()
                    hist.Fill(point[2]+(point[3]-point[2])/2,point[4]*(point[3]-point[2])*(point[1]-point[0]))
                    i+=1
                    if i==12:
                        hist.Write()
                        i=0
            elif uncerttype=="stat":
                i=0
                for point in list[1]:
                    if i==0:
                        hist=TH1F("qcd_chi-"+str(int(point[0]))+"-"+str(int(point[1]))+"stat","chi-"+str(int(point[0]))+"-"+str(int(point[1]))+"stat",len(self.bins)-1,self.bins)
                        hist.Sumw2()
                    hist.Fill(point[2]+(point[3]-point[2])/2,point[4]*(point[3]-point[2])*(point[1]-point[0]))
                    i+=1
                    if i==12:
                        hist.Write()
                        i=0
            else:
                print("Please specify style: mu or mem or stat.")
                    
    
    def dictFill(self,mydict,uncerttype):
        for key in mydict:
            myfile=TFile(self.BaseDir+"/"+self.cixsecDir+"/"+key+"_"+uncerttype+".root","RECREATE")
            histos=[]
            for list in mydict[key]:
                if uncerttype=="mu" or uncerttype=="mu6" or uncerttype=="mu30":
                    mur=str(list[0][0])
                    muf=str(list[0][1])
                    i=0
                    for j in range(1,len(list)):
                        if i==0:
                            hist=TH1F(key+"_"+"chi-"+str(int(list[j][0]))+"-"+str(int(list[j][1]))+"scale-"+mur+"-"+muf,"chi-"+str(int(list[j][0]))+"-"+str(int(list[j][1]))+"scale-"+mur+"-"+muf,len(self.bins)-1,self.bins)
                            hist.Sumw2()
                        hist.Fill(list[j][2]+(list[j][3]-list[j][2])/2,list[j][4])
                        i+=1
                        if i==12:
                            hist.Write()
                            i=0
                elif uncerttype=="mem":
                    member=str(list[0])
                    i=0
                    for j in range(1,len(list)):
                        if i==0:
                            h=key+"_"+"chi-"+str(int(list[j][0]))+"-"+str(int(list[j][1]))+"PDF-"+member
                            if h in histos:
                              print("skipping",h)
                              continue
                            histos+=[h]
                            hist=TH1F(h,h,len(self.bins)-1,self.bins)
                            hist.Sumw2()
                        hist.Fill(list[j][2]+(list[j][3]-list[j][2])/2,list[j][4])
                        i+=1
                        if i==12:
                            hist.Write()
                            i=0
                elif uncerttype=="stat":
                    i=0
                    for j in range(1,len(list)):
                        if i==0:
                            hist=TH1F(key+"_"+"chi-"+str(int(list[j][0]))+"-"+str(int(list[j][1]))+"stat","chi-"+str(int(list[j][0]))+"-"+str(int(list[j][1]))+"stat",len(self.bins)-1,self.bins)
                            hist.Sumw2()
                        hist.Fill(list[j][2]+(list[j][3]-list[j][2])/2,list[j][4])
                        i+=1
                        if i==12:
                            hist.Write()
                            i=0
                else:
                    print("Please specify style: mu or mem ot stat.")
                
            
if __name__ == "__main__":
    myVariations=getVariations()

    #qcdallmu=myVariations.getqcdallmu()
    #myVariations.listFill(qcdallmu,"mu")
    
    #qcdmu6p=myVariations.getqcdmup("6P")
    #myVariations.listFill(qcdmu6p,"mu6")
    
    #qcdmu30=myVariations.getqcdmup("30")
    #myVariations.listFill(qcdmu30,"mu30")
    
    #qcdallmem=myVariations.getqcdallmem()
    #myVariations.listFill(qcdallmem,"mem")

    #qcdstat=myVariations.getqcdstat()
    #myVariations.listFill(qcdstat,"stat")

    ciallmu=myVariations.getciallmu()
    myVariations.dictFill(ciallmu,"mu")
    
    ciallmem=myVariations.getciallmem()
    myVariations.dictFill(ciallmem,"mem")    
