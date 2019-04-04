import sys
import numpy
import re

class CodemlParseError(Exception):
    pass

def parse(ali_name, force=True) :

    match_nsite = r"^ns\s+=\s+(\d+)\s+ls\s+=\s+(\d+)\s*$"
    match_model1 = r"^Model\s1:\sNearlyNeutral\s\(2\scategories\)$"
    match_model2 = r"Model\s2:\sPositiveSelection\s\(3\scategories\)$"
    match_lnL = r"^lnL\(ntime:\s*\d+\s*np:\s*\d+\s*\):\s*(\-\d+\.\d+)\s+\+0\.\d+$"
    match_length = r"^tree\slength\s=\s+(\d+\.\d+)$"
    match_kappa = r"^kappa\s\(ts/tv\)\s=\s+(\d+\.\d+)$"
    match_prop = r"^p:\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)$"
    match_omega = r"^w:\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)$"
    match_neb = r"^Naive Empirical Bayes"
    match_beb = r"^Bayes Empirical Bayes"
    match_site = r"^\s+(\d+)\s.\s+(\d+\.\d+)\*{0,2}\s+..\s+(\d+\.\d+)"
    match_grid = r"^The grid"

    lnL2 = 0
    lnL1 = 0
    p2 = 0
    om2 = 0
    beb_selected = dict()
    neb_selected = dict()
    beb_sitepp = []
    neb_sitepp = []
    kappa = 0
    length = 0
    p0 = 0
    p1 = 0
    p2 = 0
    om0 = 0
    om1 = 0
    om2 = 0
    tree = ""

    try:
        with open(ali_name,'r') as outfile:
            
            for line in outfile :
                line = line.rstrip('\n')
                m = re.match(match_nsite, line)
                if m :
                    nsite = int(m.group(2))
                    neb_sitepp = [0 for i in range(nsite)]
                    beb_sitepp = [0 for i in range(nsite)]
                    break
            else:
                print("did not find nsite", ali_name)
                raise CodemlParseError

            for line in outfile :
                line = line.rstrip('\n')
                m = re.match(match_model1, line)
                if m :
                    break
            else:
                print("did not find model 1", ali_name)
                raise CodemlParseError

            for line in outfile :
                line = line.rstrip('\n')
                m = re.match(match_lnL, line)
                if m :
                    lnL1 = float(m.group(1))
                    lnL2 = lnL1
                    break
            else :
                print("did not find lnL1", ali_name)
                raise CodemlParseError

            for line in outfile :
                line = line.rstrip('\n')
                m = re.match(match_model2, line)
                if m :
                    break
            else:
                print("did not find model 2", ali_name)
                raise CodemlParseError

            for line in outfile :
                line = line.rstrip('\n')
                m = re.match(match_lnL, line)
                if m :
                    lnL2 = float(m.group(1))
                    break
            else :
                print("did not find lnL2", ali_name)
                raise CodemlParseError

            for line in outfile :
                line = line.rstrip('\n')
                m = re.match(match_length, line)
                if m :
                    length = float(m.group(1))
                    for i in range(4) :
                        line = outfile.readline()
                    tree = line.rstrip('\n')
                    break
            else :
                print("did not find length", ali_name)
                raise CodemlParseError

            for line in outfile :
                line = line.rstrip('\n')
                m = re.match(match_kappa, line)
                if m :
                    kappa = float(m.group(1))
                    break
            else :
                print("did not find kappa", ali_name)
                raise CodemlParseError

            for line in outfile :
                line = line.rstrip('\n')
                m = re.match(match_prop, line)
                if m :
                    p0 = float(m.group(1))
                    p1 = float(m.group(2))
                    p2 = float(m.group(3))
                    break
            else :
                print("did not find proportions", ali_name)
                raise CodemlParseError

            for line in outfile :
                line = line.rstrip('\n')
                m = re.match(match_omega, line)
                if m :
                    om0 = float(m.group(1))
                    om1 = float(m.group(2))
                    om2 = float(m.group(3))
                    break
            else :
                print("did not find omega", ali_name)
                raise CodemlParseError

            for line in outfile :
                line = line.rstrip('\n')
                m = re.match(match_neb, line)
                if m :
                    break

            else :
                print("did not find Naive Empirical Bayes", ali_name)
                raise CodemlParseError

            neb_selected = dict()

            for line in outfile :
                line = line.rstrip('\n')
                m2 = re.match(match_beb, line)
                if m2 :
                    break

                m = re.match(match_site, line)
                if m :
                    site = int(m.group(1))
                    pp = float(m.group(2))
                    om = float(m.group(3))
                    neb_selected[site-1] = pp

            else :
                print("did not find Bayes Empirical Bayes", ali_name)
                raise CodemlParseError

            neb_sitepp = [0 for i in range(nsite)]
            for (i,pp) in neb_selected.items():
                neb_sitepp[i] = pp

            beb_selected = dict()

            for line in outfile:
                line = line.rstrip('\n')
                m2 = re.match(match_grid, line)
                if m2 :
                    break

                m = re.match(match_site, line)
                if m :
                    site = int(m.group(1))
                    pp = float(m.group(2))
                    om = float(m.group(3))
                    beb_selected[site-1] = pp

            else :
                print("did not find grid (end of parsing)", ali_name)
                raise CodemlParseError

            beb_sitepp = [0 for i in range(nsite)]
            for (i,pp) in beb_selected.items():
                beb_sitepp[i] = pp

    except  CodemlParseError:
        print("default values for codeml results")

    except  FileNotFoundError:
        print("default values for codeml results")

    return [lnL2 - lnL1, p2, om2, om2, om2, beb_selected, beb_sitepp, neb_selected, neb_sitepp, kappa, length, p0, p1, p2, om0, om1, om2, tree]


def parse_list(codeml_dir, gene_list, write_output = False):

    score = dict()
    posw = dict()
    posom = dict()
    minposom = dict()
    maxposom = dict()
    selectedsites = dict()
    sitepp = dict()

    for gene in gene_list:
        res = parse(codeml_dir + gene + ".codeml")
        [score[gene], posw[gene], posom[gene], minposom[gene], maxposom[gene], selectedsites[gene], sitepp[gene]] = res[0:7]

    if write_output:
        with open(chain_name + ".postanalysis", 'w') as outfile:
            outfile.write("gene\tpp\tposw\tom\tmin\tmax")
            if with_sites:
                outfile.write("\tsites")
            outfile.write("\n")
            for gene in gene_list:
                outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(gene, score[gene], posw[gene], posom[gene], minposom[gene], maxposom[gene]))
                if with_sites:
                    outfile.write("\t{0}".format(sitepp[gene]))
                outfile.write("\n")

    return [score, posw, posom, minposom, maxposom, selectedsites, sitepp]

if __name__ == "__main__":

    import sys
    ali_name = sys.argv[1]
    # [lnL1, lnL2, kappa, length, p0, p1, p2, om0, om1, om2, tree, selected, sitepp] = parse(ali_name)
    # print(lnL1,lnL2, lnL2-lnL1)
    # print(kappa,length)
    # print(p0,p1,p2)
    # print(om0,om1,om2)
    # print(tree)
    # print(selected)


