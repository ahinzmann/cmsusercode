Here is the list of all scripts to run the dijet angular analysis.
The order of the execution matters in some cases.

-------- MC production in CMSSW_7_1_X:

submit_ci_samples_13TeV.py # produce GEN samples for LO QCD and LO QCD+CI
#submit_add_samples_13TeV.py # produce GEN samples for LO QCD+ADD CUETP8M1 tune
submit_herwig_samples_13TeV.py # produce GEN samples for QCD

#calculate_crosssections_uzh.py # extract LO QCD, LO QCD+CI and LO QCD+ADD from production Pythia log files, Herwig cross sections have to be taken from LHC.log file by hand.
xsdj_table.py # extract LO DM cross sections from LHC headers

copy-samples.py # copy GEN samples from PSI to CERN
copy-tree.sh # copy data and full simulation QCD from PSI to CERN

--------- MC production in CMSSW_10_6_19 -------:

submit_add_samples_CP5.py # produce GEN samples for LO QCD and QCD+ADD CP5 tune
submit_add_samples_CUETP8M1.py # produce GEN samples for LO QCD and QCD+ADD CUETP8M1 tune
submit_ci_samples_CP5.py # produce GEN samples for LO QCD+CI CP5 tune
calculate_xsecs_desy.py > xsecs_13p6TeV.txt # extract LO QCD, LO QCD+CI and LO QCD+ADD from production Pythia log files

source /cvmfs/grid.desy.de/etc/profile.d/grid-ui-env.sh
source /cvmfs/cms.cern.ch/cmsset_default.sh
# setup https://twiki.cern.ch/twiki/bin/view/CMS/QuickGuideMadGraph5aMCatNLO#Quick_tutorial_on_how_to_produce
cp tripleG*.dat /nfs/dust/cms/user/hinzmann/dijetangular/madgraphProduction/genproductions/bin/MadGraph5_aMCatNLO
nohup ./gridpack_generation.sh tripleG_QCD_HT2000toInf ./ # create gridpacks for tripleG
nohup ./gridpack_generation.sh tripleG_QCD_HT4000toInf ./ # create gridpacks for tripleG
cp alp*.dat /nfs/dust/cms/user/hinzmann/dijetangular/madgraphProduction/genproductions/bin/MadGraph5_aMCatNLO
nohup ./gridpack_generation.sh alp_QCD_HT2000toInf ./ # create gridpacks for alp
nohup ./gridpack_generation.sh alp_QCD_HT4000toInf ./ # create gridpacks for alp
submit_madgraph_samples_13TeV.py # produce samples for alp and tripleG

-------- Data analysis in CMSSW_10_6_X:

plotSignal_13TeV_desy_run2_loop.py # call plotSignal_13TeV_desy_run2.py # produce dijet angular histograms from GEN samples for QCD, CI and ADD samples
plotSignal_13TeV_desy_dm_run2.py # call plotSignal_13TeV_desy_run2.py # produce dijet angular histograms from GEN samples for DM samples
plotSignal_13TeV_desy_alp_run2.py # call plotSignal_13TeV_desy_run2.py # produce dijet angular histograms from GEN samples for alp samples
plotSignal_13TeV_desy_tripleG_run2.py # call plotSignal_13TeV_desy_run2.py # produce dijet angular histograms from GEN samples for tripleG samples

run_data.py run_dataSingleMuon.sh merge_data.sh # call plot_data_13TeV_desy_run2.py # produce dijet angular histograms for data and fullsim QCD samples

plotSignal_jes_13TeV_run2_loop.py # calls plotSignal_jes_13TeV_run2.py # produce JES uncertainty histograms from CI and QCD GEN samples
plot_chi_jes_plots_13TeV_combine_run2.py # plot JES uncertainty histograms
plotSignal_jer_13TeV_run2_loop.py # calls plotSignal_jer_13TeV_run2.py # produce JER uncertainty histograms from CI and QCD GEN samples
plot_chi_jer_plots_13TeV_combine_run2.py # plot JER uncertainty histograms
plotSignal_pu_13TeV_run2.sh # calls plotSignal_pu_13TeV_run2.py # produce PU uncertainty histograms from QCD Fullsim samples
plot_chi_pu_plots_13TeV_combine_run2.py # plot PU uncertainty histograms
plot_chi_prefire_plots_13TeV_combine_run2.py # plot prefire uncertainty histograms

plot_chi_control_plots_13TeV_run2.py # plot data-MC comparisons
plot_chi_model_plots_13TeV_run2.py # plot smearing vs fullsim comparisons
plot_chi_trigger_plots_13TeV_run2.py # plot trigger efficiencies

plot_nonpert_13TeV.py # plot non perturbative correction histograms
plot_dm_pdf_plots_13TeV_run2.py  # compute PDF uncertainties for DM

add_systematics_13TeV_run2_loop.py # call add_systematics_13TeV_run2.py # add systematic shift histograms, NLOQCD and data histograms in the datacards for each CI, ADD, QBH, alp, tripleG and DM signal
plot_chi_uncertainties_13TeV_run2.py # plot summary of all systematic uncertainties

plot_chi_combined_data_13TeV_dm_run2.py # DM comparison plots
plot_chi_combined_theory_13TeV_run2.py # final result plot

-------- Limit calculation in CMSSW_8_1_X:

calculate_limit_shape_13TeV_run2.py # calculate CLS for each CI, ADD, DM and QBH, alp, tripleG model
calculate_limit_shape_13TeV_dm_run2.py # calculate CLS for each DM model
plot_limit_shape_13TeV_run2.py # compute CLS limit for each CI, ADD and QBH, alp, tripleG model
plot_limit_shape_13TeV_dm_run2.py # compute CLS limit for each DM model
plot_limit_summary.py #  summary of limits on various models (still 8 TeV version)

summarize_pvalues_13TeV_run2.py # make table of p-values of the data in each mass bin
calculate_goodness_13TeV_run2.py # calculate goodness of fit measure of the data in each mass bin
summarize_goodness_13TeV_run2.py # make table of goodness of fit measure
calculate_bestfit.py # make table of bestfit JER (still 8 TeV version)

python unfold_13TeV_run2.py # run unfolding with combine
python plot_unfolded_distribution_13TeV_run2.py # plot unfolded data distribution
python plot_correlation_matrix_13TeV_run2.py #  plot correlation matrix
