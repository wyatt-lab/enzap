# Copy number

Copy number analyses included two distinct analyses. 
1) Genome-wide copy number profile were generated, ploidy models were fit algorithmically and then manually inspected and adjusted if needed to estimated ctDNA%. This step also informs chromosome-arm and whole-gene level copy number states. 
2) Focal intragenic copy number assessment. This analysis zooms into key genes with dense probe tiling (_TP53_, _RB1_, _PTEN_) and identifies focal, putatively deleterious/pathogenic copy number changes relative to the base ploidy identified in the genome-wide analysis above.  

## 1) Genome-wide copy number profiles 
Coverage at probes across the genome were calculated with the sam count tool in Seqkit v 0.8.1 (https://github.com/annalam/seqkit). The probes_hg38.bed file is a tsv with 4 columns: chromosome, probe start genomic location, probe end genomic location, feature annotation (gene name or "SNP grid"). An example can be found below.   
```
sam count sample.bam /path/to/reference/probes_hg38.bed > /project/coverage/sample.tsv
```
```
chr1	26696003	26696123	ARID1A
chr1	26696115	26696235	ARID1A
chr1	26696227	26696347	ARID1A
```

Then, using copynum2 (included in the julia folder), create a copy number profile background file with the ctDNA-negative samples in the cohort. These are inputted as a list of sample names, one sample per row, in cna_controls.txt.  
```
copynum2 background /path/to/reference/probes_hg38.bed $(cat /path/to/reference/cna_controls.txt) > path/to/reference/copynum_panel_background.jls
```

Generate Generate copy number plots for visualization of whole genome 





