******************************

Simulation experiments for m2a

******************************

General idea: starting from a set of genes, a tree and a taxon list
- analyse data with m2a model, using either codeml or bayescode, single-gene (m2a) and multi-gene (mm2a) version,
  and under various priors
- based on these analyses, re-simulate data
  (using either m2a or mm2a,
  possibly modulating prevalence of pos-sel and effect sizes for single-gene re-simulations,
  and either linking or shrinking branch lengths and nuc rates across genes for multi-gene re-simulations)
- re-analyse simulated data
- comparing estimates under the three approaches and the alternative priors with true values
  tabulate fdr, etc.

This multi-step analysis is conducted through the following sequence of 5 python scripts
(+ a few interventions inbetween, for running analyses on cluster):

# step1_empruns.py <data_folder> <gene_file> <tree_file> <exp_folder>

starting from a list of single-gene alignments contained in a data folder, a tree and a target name for the experiment:
- creates an experiment folder, with subexperiment folder (called 'empirical'), containing reformatted data
- in 'empirical', prepare data (formatted for codeml, m2a, mm2a) and scripts for running analyses 
- in the case of mm2a: prepares script for the posterior predictive analysis (to run after mcmc)

-> submit jobs for all analyses under codeml, m2a, m2a;
   all scripts (.sh) are in the subfolders named codeml, m2a and mm2a in <exp_folder>/empirical/
-> once runs are completed, submit jobs for ppred analyses for mm2a
   (scripts in mm2a subfolder in <exp_folder>/empirical)


# step2_emppost.py <exp_folder>

post-analysis of all runs on empirical data. Produces summary files:
- .sortedparams  : main results tabulated across methods for all genes sorted by decreasing DlnL (codeml)
- .hyper         : median and 95CI for hyperparameters
- .genefdr       : posterior estimate of FDR


# step3_ppred.py <exp_folder>

re-simulation based on parameter values inferred at step 2, using a modified post pred approach

- single gene resimulations (using m2a and based on codeml and m2a output)
  modulating number of genes under pos sel (cutoff applied on list of genes sorted by decreasing codeml dlnl)
  -> 30%, 10%, 3% of all genes
  modulating effect sizes: shrinking posw or dposom (or both) by a factor 0.3 or 0.1

- multi-gene re-simulations: simple post pred under mm2a (either shared or shrunken runs)

each simulation settings creates its own subfolder in <exp_folder>, each containing simulated data and tree:
simu03                      : 3% of genes under positive selection, effect sizes (posw dposom) based on m2a estimate
simu10                      : 10% of genes under positive selection (effect sizes unmodified)
simu30                      : 30% of genes under positive selection (effect sizes unmodified)
simu30_shrink_dposom01      : 30% of genes under positive selection (dposom x 0.1)
simu30_shrink_dposom03      : 30% of genes under positive selection (dposom x 0.3)
simu30_shrink_posw01        : 30% of genes under positive selection (posw x 0.1)
simu30_shrink_posw03        : 30% of genes under positive selection (posw x 0.3)
simushared                  : post pred under mm2a, with all genes sharing same branch lengths and nucleotide rates
simushrunken                : post pred under mm2a, with genes having different but similar branch lengths and nuc rates

# step4_simruns.py <exp_folder>

prepare scripts for running all analyses of re-simulated data (under all simulation settings)

-> again submit jobs for all analyses under codeml, m2a, m2a;
   all scripts (.sh) are in the subfolders named codeml, m2a and mm2a in <exp_folder>/<simu_name>/
   where <simu_name> is the name of a given re-simulation experiment
-> once runs are completed, submit jobs for ppred analyses for mm2a
   (scripts in mm2a subfolder in <exp_folder>/<simu_name>)

# step5_simpost.py <exp_folder>

post analysis: comparing estimates under the three approaches with true values
most important output:
- .sortedparams
- .hyper
- .genefdr

