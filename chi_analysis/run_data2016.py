counter=0
for i in range(11962):
   string="python plot_data_13TeV_desy_run2.py 0 "+str(i)
   if counter%5!=4:
       string+=" &"
   print string
   counter+=1
