f=open("lhe_output.txt")

xsecs_rs1={}
xsecs_add={}
for l in f.readlines():
   if "Planck" in l and "=" in l:
      m=float(l.strip(" TeV\n").split("=")[-1].strip(" "))*1000
   if "dimensions = 10" in l:
      ADD=True
   if "dimensions = 5" in l:
      ADD=False
   if "cross section" in l and "pb" in l:
     if ADD:
      xsecs_add[m]=float(l.strip(" pb\n").split(" ")[-1].strip(" "))
     else:
      xsecs_rs1[m]=float(l.strip(" pb\n").split(" ")[-1].strip(" "))
print(xsecs_add, xsecs_rs1)
