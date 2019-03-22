#! /usr/bin/python3.5

from runcodeml import * 

path2codeml_occigen = "/home/nlartillot/paml4.8/src/"
path2codeml_pbil = "/beegfs/home/lartillot/othersofts/paml4.8/src/"

exp_folder = sys.argv[1] 

# codeml

reruncodeml(exp_folder, machine="pbil", queue="none", core = 1, njobs_per_batch = 1, time = 24, mem=2, path2batch = exp_folder + "/codeml/", path2run = path2codeml_pbil)
# reruncodeml(exp_folder, machine="occigen", queue="none", core = 24, njobs_per_batch = 24, time = 24, mem=16, path2batch = exp_folder + "/codeml/", path2run = path2codeml_occigen)


