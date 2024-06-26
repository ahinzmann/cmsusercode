from ROOT import *
from array import array
import os

def writeHistogram(name,points):
   chi_binning=array('d')
   for i in range(len(points)):
       chi_binning.append(points[i][0])
   chi_binning.append(points[len(points)-1][1])
   histogram=TH1F(name,name,len(chi_binning)-1,chi_binning)
   #histogram.Sumw2()
   for i in range(len(points)):
       histogram.Fill(points[i][0]+(points[i][1]-points[i][0])/2,points[i][2])
   histogram.Write()

if __name__=="__main__":
    paths=[("CIJET_MassBin","app1_CT10nlo_5000_LL-"),
           ("CIJET_MassBin","app1_CT10nlo_6000_LL-"),
           ("CIJET_MassBin","app1_CT10nlo_7000_LL-"),
           ("CIJET_MassBin","app1_CT10nlo_8000_LL-"),
           ("CIJET_MassBin","app1_CT10nlo_9000_LL-"),
           ("CIJET_MassBin","app1_CT10nlo_10000_LL-"),
           ("CIJET_MassBin","app1_CT10nlo_11000_LL-"),
           ("CIJET_MassBin","app1_CT10nlo_12000_LL-"),
           ("CIJET_MassBin","app1_CT10nlo_13000_LL-"),
           ("CIJET_MassBin","app1_CT10nlo_14000_LL-"),
           ("CIJET_MassBin","app1_CT10nlo_15000_LL-"),
           ("CIJET_MassBin","app1_CT10nlo_5000_LL+"),
           ("CIJET_MassBin","app1_CT10nlo_6000_LL+"),
           ("CIJET_MassBin","app1_CT10nlo_7000_LL+"),
           ("CIJET_MassBin","app1_CT10nlo_8000_LL+"),
           ("CIJET_MassBin","app1_CT10nlo_9000_LL+"),
           ("CIJET_MassBin","app1_CT10nlo_10000_LL+"),
           ("CIJET_MassBin","app1_CT10nlo_11000_LL+"),
           ("CIJET_MassBin","app1_CT10nlo_12000_LL+"),
           ("CIJET_MassBin","app1_CT10nlo_13000_LL+"),
           ("CIJET_MassBin","app1_CT10nlo_14000_LL+"),
           ("CIJET_MassBin","app1_CT10nlo_15000_LL+"),
           ("CIJET_MassBin","app1_CT10nlo_5000_RR-"),
           ("CIJET_MassBin","app1_CT10nlo_6000_RR-"),
           ("CIJET_MassBin","app1_CT10nlo_7000_RR-"),
           ("CIJET_MassBin","app1_CT10nlo_8000_RR-"),
           ("CIJET_MassBin","app1_CT10nlo_9000_RR-"),
           ("CIJET_MassBin","app1_CT10nlo_10000_RR-"),
           ("CIJET_MassBin","app1_CT10nlo_11000_RR-"),
           ("CIJET_MassBin","app1_CT10nlo_12000_RR-"),
           ("CIJET_MassBin","app1_CT10nlo_13000_RR-"),
           ("CIJET_MassBin","app1_CT10nlo_14000_RR-"),
           ("CIJET_MassBin","app1_CT10nlo_15000_RR-"),
           ("CIJET_MassBin","app1_CT10nlo_5000_RR+"),
           ("CIJET_MassBin","app1_CT10nlo_6000_RR+"),
           ("CIJET_MassBin","app1_CT10nlo_7000_RR+"),
           ("CIJET_MassBin","app1_CT10nlo_8000_RR+"),
           ("CIJET_MassBin","app1_CT10nlo_9000_RR+"),
           ("CIJET_MassBin","app1_CT10nlo_10000_RR+"),
           ("CIJET_MassBin","app1_CT10nlo_11000_RR+"),
           ("CIJET_MassBin","app1_CT10nlo_12000_RR+"),
           ("CIJET_MassBin","app1_CT10nlo_13000_RR+"),
           ("CIJET_MassBin","app1_CT10nlo_14000_RR+"),
           ("CIJET_MassBin","app1_CT10nlo_15000_RR+"),
           ("CIJET_MassBin","app1_CT10nlo_5000_AA-"),
           ("CIJET_MassBin","app1_CT10nlo_6000_AA-"),
           ("CIJET_MassBin","app1_CT10nlo_7000_AA-"),
           ("CIJET_MassBin","app1_CT10nlo_8000_AA-"),
           ("CIJET_MassBin","app1_CT10nlo_9000_AA-"),
           ("CIJET_MassBin","app1_CT10nlo_10000_AA-"),
           ("CIJET_MassBin","app1_CT10nlo_11000_AA-"),
           ("CIJET_MassBin","app1_CT10nlo_12000_AA-"),
           ("CIJET_MassBin","app1_CT10nlo_13000_AA-"),
           ("CIJET_MassBin","app1_CT10nlo_14000_AA-"),
           ("CIJET_MassBin","app1_CT10nlo_15000_AA-"),
           ("CIJET_MassBin","app1_CT10nlo_5000_AA+"),
           ("CIJET_MassBin","app1_CT10nlo_6000_AA+"),
           ("CIJET_MassBin","app1_CT10nlo_7000_AA+"),
           ("CIJET_MassBin","app1_CT10nlo_8000_AA+"),
           ("CIJET_MassBin","app1_CT10nlo_9000_AA+"),
           ("CIJET_MassBin","app1_CT10nlo_10000_AA+"),
           ("CIJET_MassBin","app1_CT10nlo_11000_AA+"),
           ("CIJET_MassBin","app1_CT10nlo_12000_AA+"),
           ("CIJET_MassBin","app1_CT10nlo_13000_AA+"),
           ("CIJET_MassBin","app1_CT10nlo_14000_AA+"),
           ("CIJET_MassBin","app1_CT10nlo_15000_AA+"),
           ("CIJET_MassBin","app1_CT10nlo_5000_VV-"),
           ("CIJET_MassBin","app1_CT10nlo_6000_VV-"),
           ("CIJET_MassBin","app1_CT10nlo_7000_VV-"),
           ("CIJET_MassBin","app1_CT10nlo_8000_VV-"),
           ("CIJET_MassBin","app1_CT10nlo_9000_VV-"),
           ("CIJET_MassBin","app1_CT10nlo_10000_VV-"),
           ("CIJET_MassBin","app1_CT10nlo_11000_VV-"),
           ("CIJET_MassBin","app1_CT10nlo_12000_VV-"),
           ("CIJET_MassBin","app1_CT10nlo_13000_VV-"),
           ("CIJET_MassBin","app1_CT10nlo_14000_VV-"),
           ("CIJET_MassBin","app1_CT10nlo_15000_VV-"),
           ("CIJET_MassBin","app1_CT10nlo_5000_VV+"),
           ("CIJET_MassBin","app1_CT10nlo_6000_VV+"),
           ("CIJET_MassBin","app1_CT10nlo_7000_VV+"),
           ("CIJET_MassBin","app1_CT10nlo_8000_VV+"),
           ("CIJET_MassBin","app1_CT10nlo_9000_VV+"),
           ("CIJET_MassBin","app1_CT10nlo_10000_VV+"),
           ("CIJET_MassBin","app1_CT10nlo_11000_VV+"),
           ("CIJET_MassBin","app1_CT10nlo_12000_VV+"),
           ("CIJET_MassBin","app1_CT10nlo_13000_VV+"),
           ("CIJET_MassBin","app1_CT10nlo_14000_VV+"),
           ("CIJET_MassBin","app1_CT10nlo_15000_VV+"),
           ("CIJET_MassBin","app1_CT10nlo_5000_V-A-"),
           ("CIJET_MassBin","app1_CT10nlo_6000_V-A-"),
           ("CIJET_MassBin","app1_CT10nlo_7000_V-A-"),
           ("CIJET_MassBin","app1_CT10nlo_8000_V-A-"),
           ("CIJET_MassBin","app1_CT10nlo_9000_V-A-"),
           ("CIJET_MassBin","app1_CT10nlo_10000_V-A-"),
           ("CIJET_MassBin","app1_CT10nlo_11000_V-A-"),
           ("CIJET_MassBin","app1_CT10nlo_12000_V-A-"),
           ("CIJET_MassBin","app1_CT10nlo_13000_V-A-"),
           ("CIJET_MassBin","app1_CT10nlo_14000_V-A-"),
           ("CIJET_MassBin","app1_CT10nlo_15000_V-A-"),
           ("CIJET_MassBin","app1_CT10nlo_5000_V-A+"),
           ("CIJET_MassBin","app1_CT10nlo_6000_V-A+"),
           ("CIJET_MassBin","app1_CT10nlo_7000_V-A+"),
           ("CIJET_MassBin","app1_CT10nlo_8000_V-A+"),
           ("CIJET_MassBin","app1_CT10nlo_9000_V-A+"),
           ("CIJET_MassBin","app1_CT10nlo_10000_V-A+"),
           ("CIJET_MassBin","app1_CT10nlo_11000_V-A+"),
           ("CIJET_MassBin","app1_CT10nlo_12000_V-A+"),
           ("CIJET_MassBin","app1_CT10nlo_13000_V-A+"),
           ("CIJET_MassBin","app1_CT10nlo_14000_V-A+"),
           ("CIJET_MassBin","app1_CT10nlo_15000_V-A+"),
           
           ("CIJET_MassBin","app1_cteq66_5000_LL-"),
           ("CIJET_MassBin","app1_cteq66_6000_LL-"),
           ("CIJET_MassBin","app1_cteq66_7000_LL-"),
           ("CIJET_MassBin","app1_cteq66_8000_LL-"),
           ("CIJET_MassBin","app1_cteq66_9000_LL-"),
           ("CIJET_MassBin","app1_cteq66_10000_LL-"),
           ("CIJET_MassBin","app1_cteq66_11000_LL-"),
           ("CIJET_MassBin","app1_cteq66_12000_LL-"),
           ("CIJET_MassBin","app1_cteq66_13000_LL-"),
           ("CIJET_MassBin","app1_cteq66_14000_LL-"),
           ("CIJET_MassBin","app1_cteq66_15000_LL-"),
           ("CIJET_MassBin","app1_cteq66_5000_LL+"),
           ("CIJET_MassBin","app1_cteq66_6000_LL+"),
           ("CIJET_MassBin","app1_cteq66_7000_LL+"),
           ("CIJET_MassBin","app1_cteq66_8000_LL+"),
           ("CIJET_MassBin","app1_cteq66_9000_LL+"),
           ("CIJET_MassBin","app1_cteq66_10000_LL+"),
           ("CIJET_MassBin","app1_cteq66_11000_LL+"),
           ("CIJET_MassBin","app1_cteq66_12000_LL+"),
           ("CIJET_MassBin","app1_cteq66_13000_LL+"),
           ("CIJET_MassBin","app1_cteq66_14000_LL+"),
           ("CIJET_MassBin","app1_cteq66_15000_LL+"),
           ("CIJET_MassBin","app1_cteq66_5000_RR-"),
           ("CIJET_MassBin","app1_cteq66_6000_RR-"),
           ("CIJET_MassBin","app1_cteq66_7000_RR-"),
           ("CIJET_MassBin","app1_cteq66_8000_RR-"),
           ("CIJET_MassBin","app1_cteq66_9000_RR-"),
           ("CIJET_MassBin","app1_cteq66_10000_RR-"),
           ("CIJET_MassBin","app1_cteq66_11000_RR-"),
           ("CIJET_MassBin","app1_cteq66_12000_RR-"),
           ("CIJET_MassBin","app1_cteq66_13000_RR-"),
           ("CIJET_MassBin","app1_cteq66_14000_RR-"),
           ("CIJET_MassBin","app1_cteq66_15000_RR-"),
           ("CIJET_MassBin","app1_cteq66_5000_RR+"),
           ("CIJET_MassBin","app1_cteq66_6000_RR+"),
           ("CIJET_MassBin","app1_cteq66_7000_RR+"),
           ("CIJET_MassBin","app1_cteq66_8000_RR+"),
           ("CIJET_MassBin","app1_cteq66_9000_RR+"),
           ("CIJET_MassBin","app1_cteq66_10000_RR+"),
           ("CIJET_MassBin","app1_cteq66_11000_RR+"),
           ("CIJET_MassBin","app1_cteq66_12000_RR+"),
           ("CIJET_MassBin","app1_cteq66_13000_RR+"),
           ("CIJET_MassBin","app1_cteq66_14000_RR+"),
           ("CIJET_MassBin","app1_cteq66_15000_RR+"),
           ("CIJET_MassBin","app1_cteq66_5000_AA-"),
           ("CIJET_MassBin","app1_cteq66_6000_AA-"),
           ("CIJET_MassBin","app1_cteq66_7000_AA-"),
           ("CIJET_MassBin","app1_cteq66_8000_AA-"),
           ("CIJET_MassBin","app1_cteq66_9000_AA-"),
           ("CIJET_MassBin","app1_cteq66_10000_AA-"),
           ("CIJET_MassBin","app1_cteq66_11000_AA-"),
           ("CIJET_MassBin","app1_cteq66_12000_AA-"),
           ("CIJET_MassBin","app1_cteq66_13000_AA-"),
           ("CIJET_MassBin","app1_cteq66_14000_AA-"),
           ("CIJET_MassBin","app1_cteq66_15000_AA-"),
           ("CIJET_MassBin","app1_cteq66_5000_AA+"),
           ("CIJET_MassBin","app1_cteq66_6000_AA+"),
           ("CIJET_MassBin","app1_cteq66_7000_AA+"),
           ("CIJET_MassBin","app1_cteq66_8000_AA+"),
           ("CIJET_MassBin","app1_cteq66_9000_AA+"),
           ("CIJET_MassBin","app1_cteq66_10000_AA+"),
           ("CIJET_MassBin","app1_cteq66_11000_AA+"),
           ("CIJET_MassBin","app1_cteq66_12000_AA+"),
           ("CIJET_MassBin","app1_cteq66_13000_AA+"),
           ("CIJET_MassBin","app1_cteq66_14000_AA+"),
           ("CIJET_MassBin","app1_cteq66_15000_AA+"),
           ("CIJET_MassBin","app1_cteq66_5000_VV-"),
           ("CIJET_MassBin","app1_cteq66_6000_VV-"),
           ("CIJET_MassBin","app1_cteq66_7000_VV-"),
           ("CIJET_MassBin","app1_cteq66_8000_VV-"),
           ("CIJET_MassBin","app1_cteq66_9000_VV-"),
           ("CIJET_MassBin","app1_cteq66_10000_VV-"),
           ("CIJET_MassBin","app1_cteq66_11000_VV-"),
           ("CIJET_MassBin","app1_cteq66_12000_VV-"),
           ("CIJET_MassBin","app1_cteq66_13000_VV-"),
           ("CIJET_MassBin","app1_cteq66_14000_VV-"),
           ("CIJET_MassBin","app1_cteq66_15000_VV-"),
           ("CIJET_MassBin","app1_cteq66_5000_VV+"),
           ("CIJET_MassBin","app1_cteq66_6000_VV+"),
           ("CIJET_MassBin","app1_cteq66_7000_VV+"),
           ("CIJET_MassBin","app1_cteq66_8000_VV+"),
           ("CIJET_MassBin","app1_cteq66_9000_VV+"),
           ("CIJET_MassBin","app1_cteq66_10000_VV+"),
           ("CIJET_MassBin","app1_cteq66_11000_VV+"),
           ("CIJET_MassBin","app1_cteq66_12000_VV+"),
           ("CIJET_MassBin","app1_cteq66_13000_VV+"),
           ("CIJET_MassBin","app1_cteq66_14000_VV+"),
           ("CIJET_MassBin","app1_cteq66_15000_VV+"),
           ("CIJET_MassBin","app1_cteq66_5000_V-A-"),
           ("CIJET_MassBin","app1_cteq66_6000_V-A-"),
           ("CIJET_MassBin","app1_cteq66_7000_V-A-"),
           ("CIJET_MassBin","app1_cteq66_8000_V-A-"),
           ("CIJET_MassBin","app1_cteq66_9000_V-A-"),
           ("CIJET_MassBin","app1_cteq66_10000_V-A-"),
           ("CIJET_MassBin","app1_cteq66_11000_V-A-"),
           ("CIJET_MassBin","app1_cteq66_12000_V-A-"),
           ("CIJET_MassBin","app1_cteq66_13000_V-A-"),
           ("CIJET_MassBin","app1_cteq66_14000_V-A-"),
           ("CIJET_MassBin","app1_cteq66_15000_V-A-"),
           ("CIJET_MassBin","app1_cteq66_5000_V-A+"),
           ("CIJET_MassBin","app1_cteq66_6000_V-A+"),
           ("CIJET_MassBin","app1_cteq66_7000_V-A+"),
           ("CIJET_MassBin","app1_cteq66_8000_V-A+"),
           ("CIJET_MassBin","app1_cteq66_9000_V-A+"),
           ("CIJET_MassBin","app1_cteq66_10000_V-A+"),
           ("CIJET_MassBin","app1_cteq66_11000_V-A+"),
           ("CIJET_MassBin","app1_cteq66_12000_V-A+"),
           ("CIJET_MassBin","app1_cteq66_13000_V-A+"),
           ("CIJET_MassBin","app1_cteq66_14000_V-A+"),
           ("CIJET_MassBin","app1_cteq66_15000_V-A+"),
           
	   ("CIJET_MassBin","app1_cteq6ll_5000_LL-"),
           ("CIJET_MassBin","app1_cteq6ll_6000_LL-"),
           ("CIJET_MassBin","app1_cteq6ll_7000_LL-"),
           ("CIJET_MassBin","app1_cteq6ll_8000_LL-"),
           ("CIJET_MassBin","app1_cteq6ll_9000_LL-"),
           ("CIJET_MassBin","app1_cteq6ll_10000_LL-"),
           ("CIJET_MassBin","app1_cteq6ll_11000_LL-"),
           ("CIJET_MassBin","app1_cteq6ll_12000_LL-"),
           ("CIJET_MassBin","app1_cteq6ll_13000_LL-"),
           ("CIJET_MassBin","app1_cteq6ll_14000_LL-"),
           ("CIJET_MassBin","app1_cteq6ll_15000_LL-"),
           ("CIJET_MassBin","app1_cteq6ll_5000_LL+"),
           ("CIJET_MassBin","app1_cteq6ll_6000_LL+"),
           ("CIJET_MassBin","app1_cteq6ll_7000_LL+"),
           ("CIJET_MassBin","app1_cteq6ll_8000_LL+"),
           ("CIJET_MassBin","app1_cteq6ll_9000_LL+"),
           ("CIJET_MassBin","app1_cteq6ll_10000_LL+"),
           ("CIJET_MassBin","app1_cteq6ll_11000_LL+"),
           ("CIJET_MassBin","app1_cteq6ll_12000_LL+"),
           ("CIJET_MassBin","app1_cteq6ll_13000_LL+"),
           ("CIJET_MassBin","app1_cteq6ll_14000_LL+"),
           ("CIJET_MassBin","app1_cteq6ll_15000_LL+"),
	   ("CIJET_MassBin","app1_cteq6ll_5000_RR-"),
           ("CIJET_MassBin","app1_cteq6ll_6000_RR-"),
           ("CIJET_MassBin","app1_cteq6ll_7000_RR-"),
           ("CIJET_MassBin","app1_cteq6ll_8000_RR-"),
           ("CIJET_MassBin","app1_cteq6ll_9000_RR-"),
           ("CIJET_MassBin","app1_cteq6ll_10000_RR-"),
           ("CIJET_MassBin","app1_cteq6ll_11000_RR-"),
           ("CIJET_MassBin","app1_cteq6ll_12000_RR-"),
           ("CIJET_MassBin","app1_cteq6ll_13000_RR-"),
           ("CIJET_MassBin","app1_cteq6ll_14000_RR-"),
           ("CIJET_MassBin","app1_cteq6ll_15000_RR-"),
           ("CIJET_MassBin","app1_cteq6ll_5000_RR+"),
           ("CIJET_MassBin","app1_cteq6ll_6000_RR+"),
           ("CIJET_MassBin","app1_cteq6ll_7000_RR+"),
           ("CIJET_MassBin","app1_cteq6ll_8000_RR+"),
           ("CIJET_MassBin","app1_cteq6ll_9000_RR+"),
           ("CIJET_MassBin","app1_cteq6ll_10000_RR+"),
           ("CIJET_MassBin","app1_cteq6ll_11000_RR+"),
           ("CIJET_MassBin","app1_cteq6ll_12000_RR+"),
           ("CIJET_MassBin","app1_cteq6ll_13000_RR+"),
           ("CIJET_MassBin","app1_cteq6ll_14000_RR+"),
           ("CIJET_MassBin","app1_cteq6ll_15000_RR+"),
           ("CIJET_MassBin","app1_cteq6ll_5000_AA-"),
           ("CIJET_MassBin","app1_cteq6ll_6000_AA-"),
           ("CIJET_MassBin","app1_cteq6ll_7000_AA-"),
           ("CIJET_MassBin","app1_cteq6ll_8000_AA-"),
           ("CIJET_MassBin","app1_cteq6ll_9000_AA-"),
           ("CIJET_MassBin","app1_cteq6ll_10000_AA-"),
           ("CIJET_MassBin","app1_cteq6ll_11000_AA-"),
           ("CIJET_MassBin","app1_cteq6ll_12000_AA-"),
           ("CIJET_MassBin","app1_cteq6ll_13000_AA-"),
           ("CIJET_MassBin","app1_cteq6ll_14000_AA-"),
           ("CIJET_MassBin","app1_cteq6ll_15000_AA-"),
           ("CIJET_MassBin","app1_cteq6ll_5000_AA+"),
           ("CIJET_MassBin","app1_cteq6ll_6000_AA+"),
           ("CIJET_MassBin","app1_cteq6ll_7000_AA+"),
           ("CIJET_MassBin","app1_cteq6ll_8000_AA+"),
           ("CIJET_MassBin","app1_cteq6ll_9000_AA+"),
           ("CIJET_MassBin","app1_cteq6ll_10000_AA+"),
           ("CIJET_MassBin","app1_cteq6ll_11000_AA+"),
           ("CIJET_MassBin","app1_cteq6ll_12000_AA+"),
           ("CIJET_MassBin","app1_cteq6ll_13000_AA+"),
           ("CIJET_MassBin","app1_cteq6ll_14000_AA+"),
           ("CIJET_MassBin","app1_cteq6ll_15000_AA+"),
           ("CIJET_MassBin","app1_cteq6ll_5000_VV-"),
           ("CIJET_MassBin","app1_cteq6ll_6000_VV-"),
           ("CIJET_MassBin","app1_cteq6ll_7000_VV-"),
           ("CIJET_MassBin","app1_cteq6ll_8000_VV-"),
           ("CIJET_MassBin","app1_cteq6ll_9000_VV-"),
           ("CIJET_MassBin","app1_cteq6ll_10000_VV-"),
           ("CIJET_MassBin","app1_cteq6ll_11000_VV-"),
           ("CIJET_MassBin","app1_cteq6ll_12000_VV-"),
           ("CIJET_MassBin","app1_cteq6ll_13000_VV-"),
           ("CIJET_MassBin","app1_cteq6ll_14000_VV-"),
           ("CIJET_MassBin","app1_cteq6ll_15000_VV-"),
           ("CIJET_MassBin","app1_cteq6ll_5000_VV+"),
           ("CIJET_MassBin","app1_cteq6ll_6000_VV+"),
           ("CIJET_MassBin","app1_cteq6ll_7000_VV+"),
           ("CIJET_MassBin","app1_cteq6ll_8000_VV+"),
           ("CIJET_MassBin","app1_cteq6ll_9000_VV+"),
           ("CIJET_MassBin","app1_cteq6ll_10000_VV+"),
           ("CIJET_MassBin","app1_cteq6ll_11000_VV+"),
           ("CIJET_MassBin","app1_cteq6ll_12000_VV+"),
           ("CIJET_MassBin","app1_cteq6ll_13000_VV+"),
           ("CIJET_MassBin","app1_cteq6ll_14000_VV+"),
           ("CIJET_MassBin","app1_cteq6ll_15000_VV+"),
           ("CIJET_MassBin","app1_cteq6ll_5000_V-A-"),
           ("CIJET_MassBin","app1_cteq6ll_6000_V-A-"),
           ("CIJET_MassBin","app1_cteq6ll_7000_V-A-"),
           ("CIJET_MassBin","app1_cteq6ll_8000_V-A-"),
           ("CIJET_MassBin","app1_cteq6ll_9000_V-A-"),
           ("CIJET_MassBin","app1_cteq6ll_10000_V-A-"),
           ("CIJET_MassBin","app1_cteq6ll_11000_V-A-"),
           ("CIJET_MassBin","app1_cteq6ll_12000_V-A-"),
           ("CIJET_MassBin","app1_cteq6ll_13000_V-A-"),
           ("CIJET_MassBin","app1_cteq6ll_14000_V-A-"),
           ("CIJET_MassBin","app1_cteq6ll_15000_V-A-"),
           ("CIJET_MassBin","app1_cteq6ll_5000_V-A+"),
           ("CIJET_MassBin","app1_cteq6ll_6000_V-A+"),
           ("CIJET_MassBin","app1_cteq6ll_7000_V-A+"),
           ("CIJET_MassBin","app1_cteq6ll_8000_V-A+"),
           ("CIJET_MassBin","app1_cteq6ll_9000_V-A+"),
           ("CIJET_MassBin","app1_cteq6ll_10000_V-A+"),
           ("CIJET_MassBin","app1_cteq6ll_11000_V-A+"),
           ("CIJET_MassBin","app1_cteq6ll_12000_V-A+"),
           ("CIJET_MassBin","app1_cteq6ll_13000_V-A+"),
           ("CIJET_MassBin","app1_cteq6ll_14000_V-A+"),
           ("CIJET_MassBin","app1_cteq6ll_15000_V-A+"),
           ]
    for path in paths:
     rootfile=TFile("_".join(path)+".root","RECREATE")
     for filename in os.listdir("."):
        skip=False
        if not ".xsc" in filename: skip=True
	for s in path:
	  if not s in filename: skip=True
        if skip: continue
        print filename
        lam=filename.split(".")[0].split("_")[-2]
        f=file(filename)
        pointsNLO=[]
        name=""
        histogram=None
	nextBin=False
	lowerChi=0
	upperChi=0
        for line in f.readlines():
            split=line.replace("D","e").strip(" ").replace("  "," ").replace("  "," ").replace("  "," ").split(" ")
            if nextBin:
		nextBin=False
                lowerChi=float(split[0])
                upperChi=float(split[1])
		name="chi-"+str(int(float(split[2])))+"-"+str(int(float(split[3])))
            elif "Bin information" in line:
                nextBin=True
            elif "Lambda" in line:
                pass
            elif "LHgrid" in line:
                pass
            elif "cj" in line:
                pass
            elif "muf" in line:
                pass
            elif len(split)==4 and float(split[0])==1 and float(split[1])==1:
	        if "cteq6ll":
                  pointsNLO+=[(lowerChi,upperChi,float(split[2])/(upperChi-lowerChi))]
                else:
		  pointsNLO+=[(lowerChi,upperChi,float(split[3])/(upperChi-lowerChi))]
        if pointsNLO!=[]:
            writeHistogram(name,pointsNLO)
        pointsNLO=[]
     rootfile.Close()
