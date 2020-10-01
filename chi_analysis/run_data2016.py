counter=0
for i in range(48):
   string="python plot_data_13TeV_desy_run2.py 1 "+str(i)+" > events2017_"+str(i)+".txt"
   if counter%5!=4:
       string+=" &"
   print string
   counter+=1
