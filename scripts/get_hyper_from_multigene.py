#! /usr/bin/python3.5

import sys
import os
import parsemm2a
from makebatch import *

def get_hyper_from_multigene(chain_name, burnin = 100):
    
    res = parsemm2a.parse_list(chain_name, burnin = burnin, with_sites = False)
    est = res[9]
    hyper = dict()

    with open(chain_name + ".hyper", 'w') as outfile:
        for (i,name) in enumerate(parsemm2a.hypernamelist):
            (mean, min, max) = est[i]
            hyper[name] = mean
            outfile.write("{0}\t{1}\n".format(name, mean))

    print("hyper-param estimates in {0}.hyper\n".format(chain_name))

    return hyper


if __name__ == "__main__":

    import sys
    if len(sys.argv) == 1:
        print("get_hyper_from_multigene chain_name burnin")
        sys.exit()

    chain_name = sys.argv[1]
    burnin = int(sys.argv[2])

    hyper = get_hyper_from_multigene(chain_name, burnin = burnin)

