import sys
import numpy

# to be further elaborated 

# take a bayescode m2a analysis (single gene) 
# get alignment size from param file
# then open .sitepp file, computes posterior probabilities for positive selection
# select sites above cutoff, output number of selected sites and estimated fdr

chain_name = sys.argv[1]
burnin = int(sys.argv[2])
cutoff = float(sys.argv[3])

with open(chain_name + ".param",'r') as param_file:
    param_file.readline()
    ali_name = param_file.readline().replace("\n","").split("\t")[0]
    with open(ali_name,'r') as ali_file:
        (ntaxa,nnuc) = ali_file.readline().replace("\n","").split()

nsite = int(nnuc) // 3

print(ali_name,nnuc)

with open(chain_name + ".sitepp",'r') as sitepp_file:
    full_list = [[float(word) for word in line.replace("\n","").split()] for i,line in enumerate(sitepp_file) if i>=burnin ]
    transposed_list = [list(i) for i in zip(*full_list)]

if (len(transposed_list) != nsite):
    print("error: non matching number of sites")
    raise

postpp = [numpy.mean(sample) for sample in transposed_list]
selected = {i:pp for i,pp in enumerate(postpp) if pp > cutoff}
tdr = numpy.mean([selected[i] for i in selected])

print()
print("#selected : ",len(selected))
print("fdr       : ",1-tdr)
print()

print("site\tpp\tfdr")
tot = 0
totpp = 0
for (pp,i) in sorted([(selected[i],i) for i in selected],reverse=True):
    tot = tot + 1
    totpp = totpp + pp
    fdr = totpp / tot
    print("{0}\t{1:2f}\t{2:2f}".format(i,pp,fdr))

