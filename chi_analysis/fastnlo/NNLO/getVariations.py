from ROOT import *
from array import *
import subprocess,os

massbins=[[1200,1500],[1500,1900],[1900,2400],[2400,3000],[3000,3600],[3600,4200],[4200,4800],[4800,5400],[5400,6000],[6000,7000],[7000,13000]]

order="NNLO"

if order=="NNLO":
  PDF="ct14nnlo"
  mscale="m2"
else:
  PDF="ct14nlo"
  mscale="pt12"

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
        self.cixsecDir="CIJET_fnl5662"+self.version+"_cs_001_ct14nlo_0_56_"

    def getqcdallmu(self):
      allmulist=[]
      for m in range(len(massbins)):
        fileallmu=self.BaseDir+"/2jet.NNLO.fnl5662"+self.version+"_mjj"+str(m).replace("10","a")+"_chi_"+PDF+"_cppread_7_"+mscale+".log"
	print "load",fileallmu
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
		if order=="NNLO":
  		  x=float(line.split()[8])
  		else:
		  x=float(line.split()[7])
                xsec=[massbins[m][0],massbins[m][1],float(line.split()[3]),float(line.split()[4]),x]
                muxsecs.append(xsec)
            allmulist.append([mu,muxsecs])
        print "found",len(allmulist),"elements"
      return allmulist

    def getqcdallmem(self):
      allmemlist=[]
      for m in range(len(massbins)):
        fileallmem=self.BaseDir+"/2jet.NNLO.fnl5662"+self.version+"_mjj"+str(m).replace("10","a")+"_chi_"+PDF+"_cppread_0_"+mscale+".log"
	print "load",fileallmem
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
  		  x=float(line.split()[8])
  		else:
		  x=float(line.split()[7])
                xsec=[massbins[m][0],massbins[m][1],float(line.split()[3]),float(line.split()[4]),x]
                memxsecs.append(xsec)
            allmemlist.append([mem,memxsecs])
        print "found",len(allmemlist),"elements"
      return allmemlist

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
            for m in range(5,31):
                cimudict["CIJET_fnl5662"+self.version+"_cs_001_ct14nlo_0_56_"+str(m*1000)+"_"+style]=[]
##                 for j in range(0,57):
                filename="CIJET_fnl5662"+self.version+"_cs_001_ct14nlo_0_"+str(m*1000)+"_"+style+".xsc"
                muxsecs=[]
		print filename
                with open(self.BaseDir+"/"+self.cixsecDir+style+"/"+filename) as f:
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
                        if i>2:
                            xsec=float(line.split()[3])
                            mur=float(line.split()[1])
                            muf=float(line.split()[0])
                            #print cimudict
                            if len(cimudict[filename.replace('.xsc','').replace('_0_','_0_56_')])<9:
                                mu=[mur,muf]
                                xsecs=[mass_low,mass_high,chi_low,chi_high,xsec]
                                muxsecs=[mu,xsecs]
                                cimudict[filename.replace('.xsc','').replace('_0_','_0_56_')].append(muxsecs)    
                            else:
                                xsecs=[mass_low,mass_high,chi_low,chi_high,xsec]
                                cimudict[filename.replace('.xsc','').replace('_0_','_0_56_')][i-3].append(xsecs)
        return cimudict

    def getciallmem(self):
        cimemdict={}
        for style in self.styles:
            for m in range(5,31):
                cimemdict["CIJET_fnl5662"+self.version+"_cs_001_ct14nlo_0_56_"+str(m*1000)+"_"+style]=[]
                for j in range(0,57):
                    filename="CIJET_fnl5662"+self.version+"_cs_001_ct14nlo_"+str(j)+"_"+str(m*1000)+"_"+style+".xsc"
                    memxsecs=[]
		    print filename
                    with open(self.BaseDir+"/"+self.cixsecDir+style+"/"+filename) as f:
                        #print f
                        mem=j
                        memxsecs.append(mem)
                        i=-9999
                        for line in f:
                            #print line
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
                            if i>2:
                                xsec=float(line.split()[3])
                                mur=float(line.split()[1])
                                muf=float(line.split()[0])
                                if mur==1 and muf==1:
                                    xsecs=[mass_low,mass_high,chi_low,chi_high,xsec]
                                    memxsecs.append(xsecs)
                        cimemdict[filename.replace('.xsc','').replace('_'+str(j)+'_','_0_56_')].append(memxsecs)
        return cimemdict
    
    def dictPrint(self,mydict):
        for key in mydict:
            print key
            print mydict[key]

    def listPrint(self,mylist):
        for list in mylist:
            print list
            
    def listFill(self,mylist,uncerttype):
        myfile=TFile(self.BaseDir+"/2jet.NNLO.fnl5662"+self.version+"_mjj_chi_"+PDF+"_cppread_"+uncerttype+"_"+mscale+".root","RECREATE")
        for list in mylist:
            if uncerttype=="mu":
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
            else:
                print "Please specify tyle: mu or mem."
                    
    
    def dictFill(self,mydict,uncerttype):
        for key in mydict:
            myfile=TFile(self.BaseDir+"/"+self.cixsecDir+"/"+key+"_"+uncerttype+".root","RECREATE")
	    histos=[]
            for list in mydict[key]:
                if uncerttype=="mu":
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
			      print "skipping",h
			      continue
			    histos+=[h]
                            hist=TH1F(h,h,len(self.bins)-1,self.bins)
                            hist.Sumw2()
                        hist.Fill(list[j][2]+(list[j][3]-list[j][2])/2,list[j][4])
                        i+=1
                        if i==12:
                            hist.Write()
                            i=0
                else:
                    print "Please specify tyle: mu or mem."
                
            
if __name__ == "__main__":
    myVariations=getVariations()
    
    qcdallmu=myVariations.getqcdallmu()
    myVariations.listFill(qcdallmu,"mu")
    
    qcdallmem=myVariations.getqcdallmem()
    myVariations.listFill(qcdallmem,"mem")

    ciallmu=myVariations.getciallmu()
    myVariations.dictFill(ciallmu,"mu")
    
    ciallmem=myVariations.getciallmem()
    myVariations.dictFill(ciallmem,"mem")    