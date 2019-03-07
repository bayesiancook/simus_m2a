
# based on externally given data and tree,
# creates a new folder for the experiment: exp_folder
# a subfolder called original_data, for the sub-experiment (analyses on original data)
# within it, a list of genes and a tree

new_experiment data_folder gene_file taxon_file tree_file exp_folder

# given an experiment/subexperiment path, with data already formatted in it
# create codeml data
# create codeml scripts
# create bayes single and multi gene data folders
# calls prepare_codeml prepare_m2adata prepare_m2arun
# or do it at once

prepare_subexperiment experiment/subexperiment

run_m2a experiment/subexperiment options basename (batchmode?)
run_codeml
run_mm2a


starting point:

a data folder with a set of single-gene alignments
a taxon list
a tree

a base folder for the project, with all scripts in it

# make_experiment.py data_path gene_list taxon_list tree target_folder
for a given gene list and a target folder:
    make target folder
    make a data folder in target folder
    in this folder, and for each gene in genelist,
    make minimal single-gene alignments, pruning out taxa not in taxon list
    make subtree based on taxon list

this gives a typical experiment folder with a data folder in it
    in experiment folder
        all.list: the list of genes
        all.tree: a complete tree (spanning the union of all taxa across gene list)
    in experiment/data folder:
        the corresponding single gene alignments, codeml-readable

(1) given and experiment folder and data folder in it, conduct following analyses:

# prepare_codeml.py experiment_name
codeml analyses
    make a folder for codeml analyses
    in this folder,
    make minimal single-gene alignments and single-gene trees
    make codeml ctl file and batch files
    run codeml analyses

single-gene bayescode analyses
    make a folder for bayescode single-gene analyses
    in this folder,
    make codeml minimal single-gene alignments and single-gene trees, pruning out the taxa not in taxon list
    under alternative priors (in same folder), make batches for single-gene m2a analyses
    run those analyses

multi-gene bayescode analyses
    make a folder for bayescode multi-gene analyses
    in this folder,
    make maximal single-gene alignments (all with complete set of taxa)
    add complete tree in folder
    under alternative priors, make batches for multi-gene m2a analyses
    run those analyses

post analyses
    compare each bayescode single-gene and multi-gene analysis with codeml
    (global comparison?)
    with an option if true parameter values and true siteom are available
    list genes in decreasing delta_lnl under codeml: this will be a reference list for simulations

simulations

    parameters:
        the list of genes with decreasing delta_lnL
        a threshold for delta_lnL
        a basename

    simulations based on codeml parameter values
    output: a list of minimal single-gene alignments + complete tree + list + truesiteom and trueparam
    -> in basename_codeml

    post pred simulations based on bayescode single-gene analyses (under uninformative prior)
    output: a list of minimal single-gene alignments + complete tree + list + truesiteom and trueparam
    -> in basename_bayes_singlegene_options

    post pred simulations based on multi-gene analyses (bl and nucrates jointly ind / shrunken / shared)
    output: a list of maximal single-gene alignments, which should then be further reduced
    -> in basename_bayse_multigene_options

    each simulation this conducted produces a typical data folder, which can then be used to run the typical analysis
    post analyses: can be compared with true values

