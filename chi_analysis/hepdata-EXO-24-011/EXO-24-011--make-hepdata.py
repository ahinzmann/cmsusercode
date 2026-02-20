import hepdata_lib
from hepdata_lib import Submission
from hepdata_lib.helpers import *

def round_value_to_decimals(cont, decimals=3):
    """
    round all values in a dictionary to some decimals in one go
    default round to 3 digits after period
    possible use case: correlations where typical values are within -1,1
    : param cont : dictionary as returned e.g. by RootFileReader::read_hist_1d()
    : type  cont : dictionary
    : param decimals: how many decimals for the rounding
    : type  decimals: integer
    """

    decimals = int(decimals)

    for i, val in enumerate(cont):
        if isinstance(val, tuple):
            cont[i] = (round(val[0], decimals), round(val[1], decimals))
        else:
            cont[i] = round(val, decimals)
    return cont

submission = Submission()
submission.read_abstract("abstract.txt")
submission.add_link("Webpage with all figures and tables", "https://cms-results.web.cern.ch/cms-results/public-results/publications/EXO-24-011/")
#submission.add_link("arXiv", "not yet")
submission.add_record_id(0, "inspire")

input_dir="/data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/"
input_dir2="versions/run2ULNNLO_m2_NNPDF3/"

paper_plots=[]
paper_plots+=[("2.4 < M(JJ) < 3.0 TeV","Figure 2 upper left","datacard_shapelimit13TeV_combined_theory3_detector_run2")]
paper_plots+=[("3.0 < M(JJ) < 3.6 TeV","Figure 2 upper middle","datacard_shapelimit13TeV_combined_theory4_detector_run2")]
paper_plots+=[("3.6 < M(JJ) < 4.2 TeV","Figure 2 upper right","datacard_shapelimit13TeV_combined_theory5_detector_run2")]
paper_plots+=[("4.2 < M(JJ) < 4.8 TeV","Figure 2 middle left","datacard_shapelimit13TeV_combined_theory6_detector_run2")]
paper_plots+=[("4.8 < M(JJ) < 5.4 TeV","Figure 2 middle center","datacard_shapelimit13TeV_combined_theory7_detector_run2")]
paper_plots+=[("5.4 < M(JJ) < 6.0 TeV","Figure 2 middle right","datacard_shapelimit13TeV_combined_theory8_detector_run2")]
paper_plots+=[("6.0 < M(JJ) < 7.0 TeV","Figure 2 lower left","datacard_shapelimit13TeV_combined_theory9_detector_run2")]
paper_plots+=[("M(JJ) > 7.0 TeV","Figure 2 lower right","datacard_shapelimit13TeV_combined_theory10_detector_run2")]
paper_plots+=[("DM mediator","Figure 3","limitsDetLHCa_DMAxial_mdm1_v6b_run2")]
paper_plots+=[("ALP liniar EFT","Figure 4 left","limitsLHCaalp_coupling_run2")]
paper_plots+=[("SMEFT","Figure 4 right","limitsLHCatripleG_coupling_run2")]
paper_plots+=[("3.0 < M(JJ) < 3.6 TeV","Figure 5 upper left","datacard_shapelimit13TeV_combined_theory4_scales_pdfs_run2")]
paper_plots+=[("3.6 < M(JJ) < 4.2 TeV","Figure 5 upper right","datacard_shapelimit13TeV_combined_theory5_scales_pdfs_run2")]
paper_plots+=[("4.2 < M(JJ) < 4.8 TeV","Figure 5 middle left","datacard_shapelimit13TeV_combined_theory6_scales_pdfs_run2")]
paper_plots+=[("4.8 < M(JJ) < 5.4 TeV","Figure 5 middle center","datacard_shapelimit13TeV_combined_theory7_scales_pdfs_run2")]
paper_plots+=[("5.4 < M(JJ) < 6.0 TeV","Figure 5 middle right","datacard_shapelimit13TeV_combined_theory8_scales_pdfs_run2")]
paper_plots+=[("6.0 < M(JJ) < 7.0 TeV","Figure 5 lower left","datacard_shapelimit13TeV_combined_theory9_scales_pdfs_run2")]
paper_plots+=[("M(JJ) > 7.0 TeV","Figure 5 lower right","datacard_shapelimit13TeV_combined_theory10_scales_pdfs_run2")]
#print(paper_plots)

from hepdata_lib import Table
from hepdata_lib import RootFileReader
from hepdata_lib import Variable, Uncertainty
for massbin,figure,filename in paper_plots:
 
 if "Figure 2" in figure:
  table = Table("Fig. 2: Detector level CHI distribution ("+massbin+")")
  print(figure)
  table.description = "Normalized CHI = exp(|Y1-Y2|) distributions (1/N)*D(N)/DCHI in "+massbin+" with |Y1+Y2|/2< 1.1. The data distributions at detector-level are compared to NNLO predictions, corrected for the detector response. The error bars represent statistical and experimental systematic uncertainties combined in quadrature. Theoretical uncertainties are provided as well. Also given are the predictions for various BSM scenarios."
  table.location = "Data from "+figure
  table.keywords["reactions"] = ["P P --> JET JET X"]
  table.keywords["cmenergies"] = [13000]
  table.keywords["observables"] = ["DN/DCHI"]
  print(input_dir+input_dir2+filename)
  table.add_image(input_dir+input_dir2+filename+".pdf")

  reader = RootFileReader(input_dir+input_dir2+filename+".root")
  canvas = reader.retrieve_object("combined/main")
  print([p.GetName() for p in canvas.GetListOfPrimitives()])
  data = reader.read_graph("combined/main/DataSysStat")
  datasys = reader.read_graph("combined/main/DataSys")
  nnlomain = reader.read_hist_1d("combined/main/MainScale")
  nnlomainup = reader.read_hist_1d("combined/main/TheoryUp")
  nnlomaindown = reader.read_hist_1d("combined/main/TheoryDown")
  hsigs={}
  signames=[("hcib","$\Lambda^{-}_{LL}$ (CI) = 30 TeV"),
            ("hcic","$\Lambda^{+}_{VV}$ (CI) = 30 TeV"),
            ("hgrw","$\Lambda_{T}$ (GRW) = 13 TeV"),
            ("hqbh","$M_{QBH}$ ($n_{ED}$ = 6 ADD) = 9 TeV"),
            ("hdmb","$m_{Med}$ = 4 TeV ($g_q$ = 0.3)"),
            ("halp","$f_a/c_g$ (ALP) = 2.5 TeV"),
            ("htripleG","$\Lambda/\sqrt{C_G}$ (SMEFT) = 10 TeV")]
  sigcaptions=[]
  for signame,label in signames:
    if signame in canvas.GetListOfPrimitives():
     hsigs[signame]=reader.read_hist_1d("combined/main/"+signame)

  print(data)
  print(nnlomain)
  #print(nnlomainup)

  chi = Variable("CHI = exp(|Y1-Y2|)", is_independent=True, is_binned=True, units="")
  chi.values = [(round((([0.5]+data['x'])[i]+data['x'][i])/2.),round((data['x'][i]+(data['x']+[17])[i+1])/2.)) for i in range(len(data['x']))]

  obs = Variable("Measured detector level", is_independent=False, is_binned=False, units="")
  obs.values = round_value_to_decimals(data['y'],5)
  obs.add_qualifier("LUMINOSITY", 138.0, "fb$^{-1}$")
  total = Uncertainty("total", is_symmetric=False)
  total.values = round_value_to_decimals(data['dy'],5)
  obs.add_uncertainty(total)
  
  thnnlomain = Variable("QCD NNLO + EW NLO", is_independent=False, is_binned=False, units="")
  thnnlomain.values = round_value_to_decimals(nnlomain['y'],5)
  theory = Uncertainty("theory", is_symmetric=False)
  theory.values = round_value_to_decimals([(nnlomaindown['y'][i]-nnlomain['y'][i],nnlomainup['y'][i]-nnlomain['y'][i]) for i in range(len(nnlomain['y']))],4)
  thnnlomain.add_uncertainty(theory)
  table.add_variable(chi)
  table.add_variable(obs)
  table.add_variable(thnnlomain)
  for signame,label in signames:
    if signame in canvas.GetListOfPrimitives():
     sig = Variable(label, is_independent=False, is_binned=False, units="")
     sig.values = round_value_to_decimals(hsigs[signame]['y'],5)
     table.add_variable(sig)
  table.add_additional_resource("ROOT file", input_dir+input_dir2+filename+".root", copy_file=True)  # optional
  submission.add_table(table)

 if "Figure 3" in figure:
  table = Table("Fig. 3: Exclusion limit ("+massbin+")")
  print(figure)
  table.description = "The 95% CL upper limits on the universal quark coupling $g_{q}$ as a function of a vector or axial-vector mediator mass, with $g_{DM}=1.0$ and $m_{DM}=1$ GeV. The observed and expected limits and their variation at the 1 and 2 standard deviation levels are given."
  table.location = "Data from "+figure
  table.keywords["reactions"] = ["P P --> JET JET X"]
  table.keywords["cmenergies"] = [13000]
  print(input_dir+filename)
  table.add_image(input_dir+filename+".pdf")

  reader = RootFileReader(input_dir+filename+".root")
  #canvas = reader.retrieve_object("main")
  #print([p.GetName() for p in canvas.GetListOfPrimitives()])
  observed = reader.read_graph("main/Observed")
  expected = reader.read_graph("main/Expected")
  expected1down = reader.read_hist_1d("main/low")
  expected1up = reader.read_hist_1d("main/high")
  expected2down = reader.read_hist_1d("main/low2")
  expected2up = reader.read_hist_1d("main/high2")
  print(observed)
  print(expected)
  print(expected1up)
  print(expected1down)
  #print(nnlomainup)

  m = Variable("$m_{Med}$", is_independent=True, is_binned=False, units="TeV")
  m.values = observed["x"]

  obs = Variable("Observed", is_independent=False, is_binned=False, units="")
  obs.values = round_value_to_decimals(observed['y'],3)
  obs.add_qualifier("LUMINOSITY", 138.0, "fb$^{-1}$")
 
  exp = Variable("Expected 95% CL upper limit on $g_{q}$", is_independent=False, is_binned=False, units="")
  exp.values = round_value_to_decimals(expected['y'],3)
  sd1 = Uncertainty("1 s.d.", is_symmetric=False)
  sd1.values = round_value_to_decimals([(expected1down['y'][i]-expected['y'][i],expected1up['y'][i]-expected['y'][i]) for i in range(len(expected['y']))],3)
  sd2 = Uncertainty("2 s.d.", is_symmetric=False)
  sd2.values = round_value_to_decimals([(expected2down['y'][i]-expected['y'][i],expected2up['y'][i]-expected['y'][i]) for i in range(len(expected['y']))])
  exp.add_uncertainty(sd1)
  exp.add_uncertainty(sd2)
  
  table.add_variable(m)
  table.add_variable(obs)
  table.add_variable(exp)
  table.add_additional_resource("ROOT file", input_dir+filename+".root", copy_file=True)  # optional
  submission.add_table(table)

 if "Figure 4" in figure:
  if "ALP" in massbin:
    table = Table("Fig. 4 left: Exclusion limit ("+massbin+")")
  else:
    table = Table("Fig. 4 right: Exclusion limit ("+massbin+")")
  print(figure)
  if "ALP" in massbin:
    table.description = "The 95% CL upper limits on the ALP gluon coupling, $c_{g}$, as a function of the characteristic energy scale, $f_{a}$, assuming $m_{a}=1$ MeV. Only dijet angular distributions with $m_{jj}<f_{a}$ are used to obtain the limits on $c_{g}$. The observed limits, expected limits and their variation at the 1 and 2 standard deviation levels are given."
  elif "SMEFT" in massbin:
    table.description = "The 95% CL upper limits on the anomalous triple gluon coupling, $C_{G}$, in SMEFT as a function of the BSM physics energy scale, $\Lambda$. Only dijet angular distributions with $m_{jj}<\Lambda$ are used to obtain the limits on $C_{G}$. The observed limits, expected limits and their variation at the 1 and 2 standard deviation levels are given."
  table.location = "Data from "+figure
  table.keywords["reactions"] = ["P P --> JET JET X"]
  table.keywords["cmenergies"] = [13000]
  print(input_dir+filename)
  table.add_image(input_dir+filename+".pdf")

  reader = RootFileReader(input_dir+filename+".root")
  #canvas = reader.retrieve_object("main")
  #print([p.GetName() for p in canvas.GetListOfPrimitives()])
  observed = reader.read_graph("main/Observed")
  expected = reader.read_graph("main/Expected")
  expected1 = reader.read_graph("main/Expected1sd")
  expected2 = reader.read_graph("main/Expected2sd")
  print(observed)
  print(expected)
  print(expected1)
  print(expected2)
  #print(nnlomainup)

  if "ALP" in massbin:
    m = Variable("$f_{a}$", is_independent=True, is_binned=False, units="TeV")
  else:
    m = Variable("$\Lambda$", is_independent=True, is_binned=False, units="TeV")
  m.values = observed["x"]

  obs = Variable("Observed", is_independent=False, is_binned=False, units="")
  obs.values = round_value_to_decimals(observed['y'],(2 if "ALP" in massbin else 3))
  obs.add_qualifier("LUMINOSITY", 138.0, "fb$^{-1}$")
 
  if "ALP" in massbin:
    exp = Variable("Expected 95% CL upper limit on $c_{g}$", is_independent=False, is_binned=False, units="")
  else:
    exp = Variable("Expected 95% CL upper limit on $C_{G}$", is_independent=False, is_binned=False, units="")
  exp.values = round_value_to_decimals(expected['y'],(2 if "ALP" in massbin else 3))
  sd1 = Uncertainty("1 s.d.", is_symmetric=False)
  sd1.values = round_value_to_decimals(expected1['dy'],(2 if "ALP" in massbin else 3))
  sd2 = Uncertainty("2 s.d.", is_symmetric=False)
  sd2.values = round_value_to_decimals(expected2['dy'],(2 if "ALP" in massbin else 3))
  exp.add_uncertainty(sd1)
  exp.add_uncertainty(sd2)
  
  table.add_variable(m)
  table.add_variable(obs)
  table.add_variable(exp)
  table.add_additional_resource("ROOT file", input_dir+filename+".root", copy_file=True)  # optional
  submission.add_table(table)

 if "Figure 5" in figure:
  table = Table("Fig. 5: Particle level CHI distribution ("+massbin+")")
  print(figure)
  table.description = "Normalized CHI = exp(|Y1-Y2|) distributions (1/SIG)*D(SIG)/DCHI in "+massbin+" with |Y1+Y2|/2< 1.1 for 138 fb$^{-1}$ of integrated luminosity at $\sqrt{s} = 13$ TeV. The data distributions corrected for detector effects are compared to NNLO predictions with central scale choices $\mu_{F}$ = $\mu_{R}$ = $m_{jj}$ and $\mu_{F}$ = $\mu_{R}$ = $<p_{T}>$. Prediction with the alternative CT14 PDF set, and a comparison with the NLO prediction used in the previous publication are also given. Theoretical uncertainties with the central scale $\mu_{F}$ = $\mu_{R}$ = $m_{jj}$ and using NNPDF3.1 are provided. We also provide prediction without EW corrections applied."
  table.location = "Data from "+figure
  table.keywords["reactions"] = ["P P --> JET JET X"]
  table.keywords["cmenergies"] = [13000]
  table.keywords["observables"] = ["DSIG/DCHI"]
  print(input_dir+input_dir2+filename)
  table.add_image(input_dir+input_dir2+filename+".pdf")

  reader = RootFileReader(input_dir+input_dir2+filename+".root")
  #canvas = reader.retrieve_object("combined/main")
  #print([p.GetName() for p in canvas.GetListOfPrimitives()])
  data = reader.read_graph("combined/main/DataSysStat")
  datasys = reader.read_graph("combined/main/DataSys")
  nnlomain = reader.read_hist_1d("combined/main/MainScale")
  nnlomainup = reader.read_hist_1d("combined/main/TheoryUp")
  nnlomaindown = reader.read_hist_1d("combined/main/TheoryDown")
  nnloalt = reader.read_hist_1d("combined/main/AltScale")
  nloold = reader.read_hist_1d("combined/main/OldNLO")

  reader = RootFileReader(input_dir+input_dir2+filename.replace("pdfs","pdfs_noewk")+".root")
  noewk = reader.read_hist_1d("combined/main/NoEWK")
  print(data)
  print(nnlomain)
  #print(nnlomainup)

  chi = Variable("CHI = exp(|Y1-Y2|)", is_independent=True, is_binned=True, units="")
  chi.values = [(round((([0.5]+data['x'])[i]+data['x'][i])/2.),round((data['x'][i]+(data['x']+[17])[i+1])/2.)) for i in range(len(data['x']))]

  obs = Variable("Measured", is_independent=False, is_binned=False, units="")
  obs.values = round_value_to_decimals(data['y'],5)
  obs.add_qualifier("LUMINOSITY", 138.0, "fb$^{-1}$")
  total = Uncertainty("total", is_symmetric=False)
  total.values = round_value_to_decimals(data['dy'],5)
  sys = Uncertainty("sys", is_symmetric=False)
  sys.values = round_value_to_decimals(datasys['dy'],5)
  obs.add_uncertainty(sys)
  obs.add_uncertainty(total)
  
  thnnlomain = Variable("QCD NNLO + EW NLO ($\mu=m_{jj}$)", is_independent=False, is_binned=False, units="")
  thnnlomain.values = round_value_to_decimals(nnlomain['y'],5)
  theory = Uncertainty("theory", is_symmetric=False)
  theory.values = round_value_to_decimals([(nnlomaindown['y'][i]-nnlomain['y'][i],nnlomainup['y'][i]-nnlomain['y'][i]) for i in range(len(nnlomain['y']))],5)
  thnnlomain.add_uncertainty(theory)
  thnnloalt = Variable("QCD NNLO + EW NLO ($\mu=<p_{T}>$)", is_independent=False, is_binned=False, units="")
  thnnloalt.values = round_value_to_decimals(nnloalt['y'],5)
  thnloold = Variable("QCD NLO + EW NLO ($\mu=<p_{T}>$)", is_independent=False, is_binned=False, units="")
  thnloold.values = round_value_to_decimals(nloold['y'],5)
  thnoewk = Variable("QCD NNLO ($\mu=m_{jj}$)", is_independent=False, is_binned=False, units="")
  thnoewk.values = noewk['y']
  
  table.add_variable(chi)
  table.add_variable(obs)
  table.add_variable(thnnlomain)
  table.add_variable(thnnloalt)
  table.add_variable(thnloold)
  table.add_variable(thnoewk)
  table.add_additional_resource("ROOT file", input_dir+input_dir2+filename+".root", copy_file=True)  # optional
  submission.add_table(table)

if True:
   figure="Add. Fig. 1"
   table = Table("Correlation matrix")
   table.description = "Correlation matrix of the maximum likelihood estimators of the signal strength modifiers. The matrix is obtained after the fit to the data. The bin numbers correspond to 11 times the index (starting at 0) of the $m_{jj}$ (3.0,3.6,4.2,4.8,5.4,6.0,7.0) bin plus the index (starting at 0) of the CHI bin (1,2,3,4,5,6,7,8,9,10,12,14)."
   table.location = ""#"Data from "+figure
   table.keywords["reactions"] = ["P P --> JET JET X"]
   table.keywords["cmenergies"] = [13000]

   corr_file="/data/dust/user/hinzmann/dijetangular/CMSSW_8_1_0/src/cmsusercode/chi_analysis/versions/run2ULNNLO_m2_NNPDF3/datacard_shapelimit13TeV_unfold_withUncertainties_correlationMatrix_run2.root"
   reader = RootFileReader(corr_file)
   data = reader.read_hist_2d("c1/correlation")
   from hepdata_lib import Variable
   x = Variable("Bin number 1 (11 m$_{jj}$-bin + $\chi$-bin)", is_independent=True, is_binned=False)
   x2 = Variable("Bin number 2 (11 m$_{jj}$-bin + $\chi$-bin)", is_independent=True, is_binned=False)
   y = Variable("Correlation coefficient", is_independent=False, is_binned=False)
   x.values = data["x"]
   x2.values = data["y"]
   y.values = round_value_to_decimals(data["z"],3)
   #print(x.values,x2.values,y.values)
   table.add_variable(x)
   table.add_variable(x2)
   table.add_variable(y)
   submission.add_table(table)

outdir = "output"
submission.create_files(outdir, remove_old=True)
