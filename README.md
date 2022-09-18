"# haploAnnotator" 

A variant of riceExplorer

This pipeline annotates all functional variants in a stipulated gene region giving graphical and tabular output to illustrate where SNPs and Indels are found relative to a reference genome. It also evaluates mean phenotypic values against each haplotype in a crop (or other species) WGS panel. 


It differs from riceExplorer as it uses a single crop panel with phenotypic data against each accession. See the CSV files for formatting. Run as `./BPH_Annot.bash` after modifying the `BPH_Annot_vars` file.

There is also the option of running the pipeline without pheontypic data to simply out a graphical illustration of functional variants within genes without associated trait means (please run `reduci.bash`).


You must have BCFTOOLS and SNPEFF software installed on your system (see riceExplorer docs)

The  `BPH_Annot_vars` file, requires the following variable set:

"fldr" - a folder to be created for processed output files<br/>
"figs_out" - a folder to be created for graphical outputs<br/>
"genes" - a file of gene regions to be evaluated e.g. `bph_genes.csv`<br/>
"samps" - a file listing gvcf paths e.g., `kpp_AUC_paths.csv`<br/>
"data" - trait data with columns 'acc' and 'trait' e.g., `kpp_AUC.csv` # NB if you use `reduci.bash` you need a dummy 'trait' column<br/>
"gene_id" - suffix of annotated genes e.g. 'gene'<br/>
"db_id" - prefices of database gvcf files (e.g. 'W_OS_' - ALL files must start with same suffix) <br/>
"db" - Name of snpEff database required e.g. Oryza_sativa <br/>
"chrom" - Chromosomes column number in "genes" file<br/>
"loc_col" - Annotated gene names column number in "genes"  file<br/>
"first" - Column number of first base positions in "genes" file<br/>
"end" - Column number of end base positions in "genes" file<br/>
"nChroms" - No. of chromosomes for the organism<br/>




