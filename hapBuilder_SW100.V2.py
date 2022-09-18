import matplotlib.pyplot as plt
import pandas as pd, sys, glob, numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
import xlsxwriter
import concurrent.futures, itertools
max_workers = None #nProcesses - leave as default on HPC [ie None]; set as no. processors on your machine -1 for laptops/Desktops 
chunksize = 1 #break up of parallelisation chunks; large values may be faster for larger datasets


#write colour coded XL sheet matching PDF
def writeXL(tab, colours, gene, out_dir): 
    colours = pd.DataFrame(colours)
    writer = pd.ExcelWriter(f"{out_dir}/{gene}.xlsx", engine='xlsxwriter')

    tab.to_excel(writer, sheet_name='Sheet1', index=False)

    workbook = writer.book
    worksheet = writer.sheets['Sheet1']

    #map color df to hap df - done col by col because different vals (e.g. AGCT) may have diff colours in diff cols
    #a col range maybe A2:A10 (here generated as var 'txt3') BUT if more than 26 columns then eg AC2:AC10 (thus complicated)
    for i in range(tab.shape[1]):
        dik = {}
        for j in range(tab.shape[0]):
            if tab.iloc[j, i] not in dik: dik[f'"{tab.iloc[j, i]}"']  = colours.iloc[j, i]

        if i < 26: txt3 = f"{chr(i+65)}2:{chr(i+65)}{tab.shape[0] + 1}" #cells to colour
        else:
            first = chr((i) // 26 + 64)
            m = ((i // 26) - 1) * 26
            txt3 = f"{first}{chr(i - m - 25 + 64)}2:{first}{chr(i - m - 25 + 64)}{tab.shape[0] + 1}" #cells to colour

        for val, color in dik.items():
            fmt = workbook.add_format({'bg_color': color})
            worksheet.conditional_format(txt3, {'type': 'cell',
                                               'criteria': '=',
                                               'value': val,
                                               'format': fmt})

    writer.save()

#START
fldr = sys.argv[1]
out_dir = sys.argv[2]
files = glob.glob1(fldr, '*haps.csv')
print(files)

def mkFigs(file):
    df = pd.read_csv(f"{fldr}/{file}", header=0)

    #check for nans - different runs may have different snpEff positions called depending on inputted samples - nans should result ONLY fron REF alleles
    for col in df.columns[1:]:
        if np.isin(True, pd.isna(df[col])):
            x = df[col].unique() #id all variants
            x = x[~pd.isnull(x)] #rm nans
            for qw in x: #rm ALT alleles
                if ',<NON_REF>' in qw: x = np.setdiff1d(x, qw)
            if x.shape[0] == 1: df[col].fillna(x[0], inplace=True) #normally should be one REF allele left
            if x.shape[0] > 1: #occasionally >1 REF - convert all REFs and nans to a single REF
                for qw in x:
                    df[col][df[col] == qw] = x[0]
                    df[col].fillna(x[0], inplace=True)
                print(f"WHOA: {file} {col}\nThis is a sanity checker!!!\nThere is more than one remaining allele in this column where we should expect only a REF allele - maybe check your *ann.best files\nA single REF was called throughout this column\n")
            if x.shape[0] == 0: #only ALTs are remaining - therefore call an ambiguous REF elsewhere in the column
                df[col].fillna('N', inplace=True)

    df['hap'] = df[df.columns[1:].tolist()].agg(' '.join, axis=1) #compress haps
    df.hap = df.hap.str.replace(',<NON_REF>', '[ALT]', regex=True)
    df['hap_no'] = 'x'
    c = 1
    for hap in df.hap.unique():
        df['hap_no'][df.hap == hap] = f"Hap{str(c).zfill(3)}" 
        c += 1

    df = df.sort_values('acc')
    df.to_csv(f"{out_dir}/{file.replace('.haps.csv','.haplotypes.csv')}", index=False) #write files showing hap_no. against each sample
    df2 = df.iloc[:, 1:]
    df3 = df2.groupby(['hap']).size().reset_index().rename(columns={0:'cnt'})
    df3 = df3.astype(object)

    #count no.s of each hap
    df2['counts'] = 'x'
    for hap in df3.hap.unique():
        try: df2.counts[df2.hap == hap] = df3.cnt[df3.hap == hap].squeeze()
        except: df2.counts[df2.hap == hap] = df3.cnt[df3.hap == hap].iloc[0]

    df4 = df2.drop_duplicates()
    
    cols = []
    for col in df4.columns[:-1]:
        if df4[df4[col].str.contains('NON_REF')].shape[0] > 0:
            cols.append(col)

    snpsDF = df4[cols + ['hap_no','counts']]
    snpsDF = snpsDF.reset_index(drop=True)
    snpsDF = snpsDF.sort_values('counts', ascending=False)


    #NOW make colours array to code REF V ALT variants in table
    colours = []
    booleanDictionary = {True: 'orange', False: 'b'}
    indelDictionary = {False: 'silver', True: 'yellow'}
    for col in snpsDF.columns[:-2]: #cycle thru rows (i.e. haps)        
        tmp = snpsDF[col]
        colors = list(tmp.str.contains('NON_REF').map(booleanDictionary)) #alt alleles get orange

        #now deal with indels and multiSNP ALTs
        tmp2 = tmp.str.replace(',<NON_REF>', '') #allows assessment of alleles as Indels (w/ len >1) or otherwise
        tmp3 = tmp2.str.replace(',', '_') #after NON_REF txt rmvd convert any remaining comma separating two or more ALT snps so can be id'd
        qw = tmp2.str.len() > 1 #id in tmp2 if element lengths > 1 - if so an indel
        if (qw.unique().size == 2 or tmp2.str.len().max() > 1) and False in tmp3.str.contains('_').unique(): #qw unq booleans will be one if no indels
            colors = qw.map(indelDictionary) #indels get yellow integer repesenting length or silver 0
            snpsDF[col] = np.core.defchararray.add((np.array(tmp2.str.len()) - 1).astype(str), 'bp')

        colours.append(colors)
        snpsDF[col] = snpsDF[col].str.replace(',<NON_REF>','')

    colours.extend([snpsDF.shape[0] * ['pink']] * 2)
    colours = np.array(colours).T
    
    #tidy up
    chrom = str(int(file[6:8]))
    new_cols = []
    for col in snpsDF.columns[:-2]: new_cols.append(f"{chrom}:{col}")
    new_cols.extend(['HapNo.', 'Counts'])
    snpsDF.columns = new_cols
    
    #make colour coded table & write as PDF
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.axis('tight')
    ax.axis('off')
    the_table = ax.table(cellText=snpsDF.values, colColours=['lime'] * snpsDF.shape[1],
                         cellColours=colours, colLabels=snpsDF.columns, loc='center',
                         cellLoc='center')
    [cell.set_text_props(fontproperties=FontProperties(weight='bold')) for key, cell in the_table.get_celld().items()] #make bold font 
    pp = PdfPages(f"{out_dir}/{file.replace('.haps.csv', '')}.altVars.pdf")
    pp.savefig(fig, bbox_inches='tight')
    pp.close()
    plt.close()
    
    #write an Excel version
    writeXL(snpsDF, colours, file.replace('.haps.csv', ''), out_dir)


def main(): #run parallelised function
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        executor.map(mkFigs, files, chunksize=chunksize)

if __name__ == '__main__': main()


