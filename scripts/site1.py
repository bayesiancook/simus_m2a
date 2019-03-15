
    # for each method and each cutoff
    # record number of selected sites, mean pp, and fdr at the site level, and this, across all genes 
    # should do it only for those genes that are above a given threshold?
    nsel = dict()
    etdr = dict()
    if fromsimu:
        tdr = dict()

    for cutoff in cutoff_list:
        cnsel = dict()
        cetdr = dict()
        if fromsimu:
            ctdr = dict()
        for name in namelist:
            gcnsel = dict()
            gcetdr = dict()
            if fromsimu:
                gctdr = dict()
            for gene in genelist:
                gcnsel[gene] = sum([(pp>cutoff) for (i,pp) in selectedsites[name][gene].items()])
                gcetdr[gene] = sum([pp for (i,pp) in selectedsites[name][gene].items() if pp > cutoff])
                if fromsimu:
                    gctdr[gene] = sum([ ((pp > cutoff) and (truesiteom[gene][i] > 1.0)) for (i,pp) in selectedsites[name][gene].items()])
            cnsel[name] = gcnsel
            cetdr[name] = gcetdr
            if fromsimu:
                ctdr[name] = gctdr
        nsel[cutoff] = cnsel
        etdr[cutoff] = cetdr
        if fromsimu:
            tdr[cutoff] = ctdr

