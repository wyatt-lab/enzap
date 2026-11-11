# Copy number

Copy number analyses included two distinct analyses. 
1) Genome-wide copy number profile were generated, ploidy models were fit algorithmically and then manually inspected and adjusted if needed to estimated ctDNA%. This step also informs chromosome-arm and whole-gene level copy number states. 
2) Focal intragenic copy number assessment. This analysis zooms into key genes with dense probe tiling (_TP53_, _RB1_, _PTEN_) and identifies focal, putatively deleterious/pathogenic copy number changes relative to the base ploidy identified in the genome-wide analysis above.  

## 1) Genome-wide copy number profiles 
Coverage at probes across the genome were calculated with the sam count tool in Seqkit v 0.8.1 (https://github.com/annalam/seqkit). 
```
sam count sample.bam /path/to/reference/probes_hg38.bed > /project/coverage/sample.tsv
```
The probes_hg38.bed file is a tsv with 4 columns: chromosome, probe start genomic location, probe end genomic location, feature annotation (gene name or "SNP grid"). An example:
```
chr1	26696003	26696123	ARID1A
chr1	26696115	26696235	ARID1A
chr1	26696227	26696347	ARID1A
```

Then, using copynum2 (included in the julia folder), create a copy number profile background file with the ctDNA-negative samples in the cohort. These are inputted as a list of sample names, one sample per row, in cna_controls.txt. Repeat this with a subset of WBC samples that have no germline copy number alterations. 
```
copynum2 background /path/to/reference/probes_hg38.bed $(cat /path/to/reference/cna_controls_cfdna.txt) > path/to/reference/copynum_panel_background_cfdna.jls
copynum2 background /path/to/reference/probes_hg38.bed $(cat /path/to/reference/cna_controls_wbc.txt) > path/to/reference/copynum_panel_background_wbc.jls
```

Generate files containing logratio values across the probes, one file per sample, including for WBC samples.   
Generate copy number plots for visualization of whole genome by switching the --output-format to "PDF" instead of "IGV". These will show genome-wide copy number profile as well as violin plots for coverage logratio over targeted genes.   
   
The cn_features bed file contains the location of the probed region for each target gene.
```
copynum2 call --report-dir=/project/logratio_tracks --output-format=IGV --features=/path/to/reference/cn_features.bed --snps_dir=/path/to/mutato_analyze/mut/files --segment=0.2 /path/to/reference/copynum_panel_background_cfdna.jls $(ls -1 | grep "cfDNA" ) &

copynum2 call --report-dir=/project/logratio_tracks --output-format=IGV --features=/path/to/reference/cn_features.bed --snps_dir=/path/to/mutato_analyze/mut/files --segment=0.2 /path/to/reference/path/to/reference/copynum_panel_background_wbc.jls $(ls -1 | grep "WBC" ) &
```
An example line in the cn_features.bed file (do not include headers in file):
```
chr1	26696003	26780786	ARID1A
```

Then, in the julia REPL, run subclones_windows for copy number model fitting from the /project/logratio_tracks folder. 
```
include("/path/to/tools/subclones_windows")
truncal_model("sample_name", interactive=true, adjust_depth=0.567, vcf_path="-")
```
An interactive window will load, with genome-wide probs and lines indicating ploidy levels. The pre-calculated ctDNA% and diploid line can be adjusted until the model best fits the observed data.

2) Focal copy number analyses
This was performed by plotting logratio and B allele frequency per probe across the gene, overlaying annotation of exon locations and ploidy lines derived from step 1. Focal copy number changes are visually identified from these plots. The logratio tracks from patient-matched WBC are plotted below the cfDNA track, to ensure that we do not misinterpret germline copy number differences as somatic.









