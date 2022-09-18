import pandas as pd, sys, glob, numpy as np
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

import concurrent.futures, itertools
max_workers = None #nProcesses - leave as default on HPC [ie None]; set as no. processors on your machine -1 for laptops/Desktops 
chunksize = 1 #break up of parallelisation chunks; large values may be faster for larger datasets

data = sys.argv[1] # 
fldr = sys.argv[2]

files = glob.glob1(fldr, '*best') #get snpEff output files
dat = pd.read_csv(data, header=0) #get sample name/phenotype data
dat = pd.concat([dat.acc, dat.trait], axis=1) #core of df (which gets modified) must comprise only two columns

#main processing function
def fnc(msu, dat):
    try: #load files sequentially (i.e. not parallelised but not time consuming)
        snps = pd.read_csv(f"{fldr}/{msu}", header = None, sep='\t')
        posn_samps = pd.read_csv(f"{fldr}/{msu.replace('.ann.best','.posn_samps.csv')}", header=0)
        vcf = pd.read_csv(f"{fldr}/{msu.replace('.ann.best','.vcf')}", header=0, sep='\t')
    except: return

    #id SNPs with multi ALT alleles
    pos_cnts = vcf.groupby(['POS', 'ALT']).size().reset_index().rename(columns={0:'cnt'})
    pos_cnts2 = pos_cnts.groupby(['POS']).size().reset_index().rename(columns={0:'cnt'})
    multi = pos_cnts2[pos_cnts2.cnt >1].POS.unique()

    for mult in multi:
        alts = vcf[vcf.POS == mult].sort_values('ALT').ALT.values
        cnt = 1
        for alt in alts:
            try:
                snps[1][(snps[1] == mult) & (snps[4] == alt)] = f"{mult}_{cnt}"
            except:
                pass
            cnt += 1
            
    snps[1] = snps[1].astype(str)
    posn_samps.marker = posn_samps.marker.astype(str)

    #build haplotypes at separate columns in df
    mks = snps[1].sort_values().tolist() #NB no headers with snpEff outputs
    for mk in mks:
        dat[mk] = np.nan
        tmp = posn_samps[posn_samps.marker == mk]
        samps = tmp.samples.tolist()[0][2:-2].split("', '") #bit ugly but works from outputs
        dat[mk][dat.acc.isin(samps)] = snps[4][snps[1] == mk].values[0] #alt alleles
        dat[mk][~dat.acc.isin(samps)] = snps[3][snps[1] == mk].values[0] #ref alleles

    dat['hap'] = dat[mks].agg(''.join, axis=1) #aggregate haplotype columns into single column
    haps = dat.hap.value_counts().keys().tolist() #list haplotypes (is ordered by freq)
    dat['val'] = None
    for h in haps: dat.val[dat.hap == h] = dat.trait[dat.hap == h].mean() #column of means by haplotype

    #id haplos that have mean trait value above or below thresh
    sign = dat.copy() #[dat.groupby('hap')['hap'].transform('size') > 0] #[dat.val <= q1] #lower percentiles
    dat = dat.iloc[:, :2] #strip dat back to original two columns for next input file

    #tidy up df
    del sign['trait'], sign['hap'], sign['val']
    for mk in mks: #rm invariant SNPs in sign loci
        try:
            if sign[mk].unique().size == 1: del sign[mk]
        except: continue

    #write file if any interesting haplotypes left
    if not sign.empty: sign.to_csv(f"{fldr}/{msu.replace('.ann.best', '.haps.csv')}", index = False)

def main():
    #run parallelised function
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        executor.map(fnc, files, itertools.repeat(dat), chunksize=chunksize)

if __name__ == '__main__': main()

