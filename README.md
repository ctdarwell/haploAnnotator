"# haploAnnotator" 

A variant of riceExplorer

This pipeline annotates all functional variants in a stipulated gene region giving graphical and tabular output to illustrate where SNPs and Indels are found relative to a reference genome. It also evaluates mean phenotypic values against each haplotype in a crop (or other species) WGS panel. 

The difference is that it uses a single crop panel with phenotypic data against each accession. See the CSV files for formatting. Run as `./BPH_Annot.bash` after modifying the `BPH_Annot_vars` file.
