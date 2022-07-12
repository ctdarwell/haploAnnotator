#!/bin/bash
#SBATCH --cpus-per-task=24

#The next line should be altered to specifically identify that versions of BCFtools you are using
#module load BCFtools

source BPH_Annot_vars #input parameters file
python installer.py xlsxwriter pandas scipy matplotlib

#ADD ''' HERE IF NEEDED TO PART RUN WORFLOW:


### *** START OF WORKFLOW ***
mkdir $fldr 
mkdir $figs_out

rm $fldr/* #should be empty fldr as we write using '>>'


bcftask(){
while read p; #samples from samps
do echo $(echo $p | rev | cut -d'/' -f 1 | rev | cut -d'.' -f 1) >> $3/$(echo $1 | cut -d',' -f $7).vcf
bcftools view --min-ac=1  $p $(echo $1 | cut -d',' -f $4):$(echo $1 | cut -d',' -f $5)-$(echo $1 | cut -d',' -f $6 | tr -d $'\n') | sed -n '/AC=/,$p' >> $3/$(echo $1 | cut -d',' -f $7).vcf
done < $2
}

export -f bcftask
xargs -a $genes -rd'\n' -P20 -I'{}' bash -c 'bcftask "$1" "$2" "$3" "$4" "$5" "$6" "$7"' bash {} "$samps" "$fldr" "$chrom" "$first" "$end" "$loc_col" #syntax note: {} represents $genes

echo 1. DONE: bcftools-IRRI
wait

#rm repeat entries (i.e. from diff accs), format for snpEff & write allele X sample file: *.posn_samps.csv; sys.argv[2] denotes identifying prefix of sample name (to distinguish lines from rest of data)
python vcf4snpEffLoop_SW.V3.py $gene_id*vcf $fldr $db_id #file search term, folder, sample prefix
echo 2. DONE: vcf4snpEffLoop.py-1
wait



#run snpEff - id potentially influential SNPs
snpefftask(){
java -Xmx4g -jar /fs/home/cliveterence.dar/snpEff/snpEff.jar $2 $1 > $(echo $1 | cut -d'.' -f 1).ann.vcf
#select only loci with significant predicted impact
grep 'MODERATE' $(echo $1 | cut -d'.' -f 1).ann.vcf >> $(echo $1 | cut -d'.' -f 1).ann.best
grep 'HIGH' $(echo $1 | cut -d'.' -f 1).ann.vcf >> $(echo $1 | cut -d'.' -f 1).ann.best
}

export -f snpefftask
find $fldr/ -name $gene_id*vcf | xargs -P20 -I'{}' -d'\n' bash -c 'snpefftask "$1" "$2"' bash {} $db #searches for all LOC*vcf files

#rm $fldr/*.ann.vcf
#find $fldr/*best -size 0 -delete
echo 3. DONE: snpEff
wait



#evaluate which IRRI haplotypes (based on snpEff influential SNPs), have extremes phenotypes 
python hapXphenoPredictor_BPH.V2.py $data $fldr
echo 4. DONE: hapXphenoPredictorSNPs.py - checked haps!
wait


python hapBuilder_BPH.py $fldr $figs_out
wait

echo DONE $fldr $(date +%F_%T) >> fin.txt


