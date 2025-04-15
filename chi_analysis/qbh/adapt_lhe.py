names=[]
for signalMass in [4500,5000,5500,6000,6500,7000,7500,8000]:
  names+=["LHEfile_"+str(signalMass)+"_1"]
for signalMass in [7500,8000,8500,9000,9500,10000,11000,12000]:
  names+=["LHEfile_"+str(signalMass)+"_0"]

for name in names:
  infile=open(name+".lhe","r")
  outfile=open(name+"_unweighted.lhe","w")
  first_line_in_event=False
  for lin in infile.readlines():
    for l in lin.split("\\n"):
     #print(l)
     if first_line_in_event:
       splitl=l.split(" ")
       outl=" ".join(splitl[:2])+" 1.0 "+" ".join(splitl[3:])
       #outl=" ".join(splitl[:2])+" "+str(1000.*float(splitl[2]))+" "+" ".join(splitl[3:])
       first_line_in_event=False
     else:
       outl=l
     if "<event>" in l:
       first_line_in_event=True
     #print(outl)
     outfile.write(outl)
  infile.close()
  outfile.close()
