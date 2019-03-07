#! /usr/bin/python3.5

from prepare_subexperiment import *
from runcodeml import * 
from runm2a import * 
from runmm2a import * 
from simuparams import * 

path2codeml_occigen = "/home/nlartillot/paml4.8/src/"
path2codeml_pbil = "/beegfs/home/lartillot/othersofts/paml4.8/src/"
path2bayescode_occigen = "/scratch/cnt0027/bbe0449/nlartillot/bayescode/data/"

exp_folder = sys.argv[1] 

(pi, purw_mean, purw_var, purw_invconc, posw_mean, posw_var, posw_invconc, purom_mean, purom_var, purom_invconc, dposom_mean, dposom_var, dposom_invshape) = get_empirical_moments(exp_folder)

infm2a_options= "-x 1 200 -pi {0} -purom {1} {2} -dposom {3} {4} -purw {5} {6} -posw {7} {8}".format(pi, purom_mean, purom_invconc, dposom_mean, dposom_invshape, purw_mean, purw_invconc, posw_mean, posw_invconc)

infmm2a_options= "-x 1 200 -nucrates ind -bl ind +G -pi {0} 0 -purom {1} {2} -dposom {3} {4} -purw {5} {6} -posw {7} {8}".format(pi, purom_mean, purom_invconc, dposom_mean, dposom_invshape, purw_mean, purw_invconc, posw_mean, posw_invconc)

print(infm2a_options)
print(infmm2a_options)

runm2a(exp_folder, infm2a_options, "infm2a", machine="occigen", queue="none", core = 24, njobs_per_batch = 24, time = 2, mem=16, path2batch = exp_folder + "/singlegene/", path2run = path2bayescode_occigen)

runmm2a(exp_folder, infmm2a_options, "infmm2a", machine="occigen", queue="none", core = 24, time = 2, mem=16, path2batch = exp_folder + "/multigene/", path2run = path2bayescode_occigen)

