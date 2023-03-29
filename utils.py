def partition_bins(dt, keys, bins, keyName):
    n=len(keys)
    idx = sorted(list(zip(range(0, n), keys)), key=lambda k : k[1])
    remain = n
    binning = [0]*n
    bin=0
    while remain > 0 and bin < bins:
        per_bin = remain//(bins-bin)
        for i in range(0,per_bin):
            binning[idx[n-remain+i][0]] = bin
        i = per_bin
        while n-remain+i<n and idx[n-remain+i][1]==idx[n-remain+per_bin-1][1]:
            binning[idx[n-remain+i][0]] = bin
            i += 1
        remain -= i
        bin += 1
    return dt.assign(**{ keyName : binning })

    
