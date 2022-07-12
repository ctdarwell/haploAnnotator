import pandas as pd, sys, glob, os
import warnings
warnings.filterwarnings("ignore")

vcfs = sys.argv[1] #e.g. LOC*vcf or rgd_*vcf
fldr = sys.argv[2]

#out_fldr = 'tmp' #!!!!!!!!!!!!!!!!!!!!!!!!

panel = sys.argv[3] #eg IRIS W00
files = glob.glob1(fldr, vcfs) #PRINT LEN FILES!!!!!!!!!!!!!!!!!!!!??????????????

#files = ['gene_Bph6.vcf'] #!!!!!!!!!!!!!!!!!!!!!!

for vcf in files:
    cols = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','tbr']

    samps, dat = [], []
    for line in open(f"{fldr}/{vcf}"):
        if line.startswith(panel): samp = line.strip()
        else:
            try:
                samps.append(samp)
            except:
                continue
            dat.append(line.strip().split('\t'))

    if len(dat) == 0: #THIS MAYBE INEFFECTIVE
        try:
            del cols, line, samp
            os.remove(f"{fldr}/{vcf}")
        except: continue
        continue

    df = pd.DataFrame(dat)
    df = pd.concat([df, pd.DataFrame(samps)], axis=1)
    cols.append('samp')
    df.columns = cols
    del df['FORMAT']
    del df['tbr']

    df.INFO = '.'
    df.QUAL = '.'
    df.FILTER = '.'
    df.POS = df.POS.astype(int)
    df = df.sort_values('POS')

    #id SNPs with multi ALT alleles
    pos_cnts = df.groupby(['POS','REF','ALT']).size().reset_index().rename(columns={0:'cnt'})
    pos_cnts2 = pos_cnts.groupby(['POS']).size().reset_index().rename(columns={0:'cnt'})
    multi = pos_cnts2[pos_cnts2.cnt >1].POS.unique()
    df2 = df.copy()

    for mult in multi:
        alts = pos_cnts[pos_cnts.POS == mult].sort_values('ALT')
        cnt = 1
        for i in alts.index:
            df2.POS[(df2.POS == mult) & (df2.REF == alts.REF.loc[i]) & (df2.ALT == alts.ALT.loc[i])] = f"{mult}_{cnt}"
            cnt += 1
            


    #write csv of posn samples
    posn_samps = []
    for pos in df2.POS.unique(): posn_samps.append([pos, df2.samp[df2.POS == pos].tolist()])
    posn_samps = pd.DataFrame(posn_samps)
    posn_samps.columns = ['marker', 'samples']
    posn_samps.to_csv(f"{fldr}/{vcf.replace('.vcf','.posn_samps')}.csv", index = False)

    del df['samp']
    df = df.drop_duplicates() #MAY NEED TO THINK ABOUT THIS!!!

    df.to_csv(f"{fldr}/{vcf}", index = False, sep='\t')
    del cols, line, samp

