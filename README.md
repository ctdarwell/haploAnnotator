"# haploAnnotator" 

A variant of riceExplorer

This pipeline annotates all functional variants in a stipulated gene region giving graphical and tabular output to illustrate where SNPs and Indels are found relative to a reference genome. It also evaluates mean phenotypic values against each haplotype in a crop (or other species) WGS panel. There is also the option of running the pipeline without pheontypic data to simply out a graphical illustration of functional variants within genes without associated trait means (please download and unzip the `reduci.zip` utility).

It differs from riceExplorer as it uses a single crop panel with phenotypic data against each accession. See the CSV files for formatting. Run as `./BPH_Annot.bash` after modifying the `BPH_Annot_vars` file.

You must have BCFTOOLS and SNPEFF software installed on your system (see riceExplorer docs)

The Hap

