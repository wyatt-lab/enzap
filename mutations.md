## Somatic mutations
Somatic mutations were called in parallel with Vardict and Mutato Analyze. Only mutations detected by both variant callers were then further curated, filtered, and annotated.

### Vardict
Vardict is run in paired sample mode, where FILENAME is the cfDNA bam file, and FILENAME_WBC is the patient-matched white blood cell germline control. BASENAME_FILE amd BASENAME_WBC refer to the cfDNA and WBC sample names without the ".bam" extension.  
target_bed_file refers to a bed file of probe locations in our targeted panel, including the genome-wide SNP grid.  
  
We required a minimum of 4 reads and 0.2% VAF for initial Vardict calling. These thresholds are made more stringent in downstream filtering steps (see step 5 below).  
Note that testsomatic.R, var2vcf_paired.pl, and  are standard helper scripts in the vardict-java workflow.  

```
vardict-java \
-b "${FILENAME}|${FILENAME_WBC}" \
-f 0.002 -r 4 --nosv -k 0 -N "${BASENAME_FILE}" \
-c 1 -S 2 -E 3 ${target_bed_file}  \
-G /path/to/reference/hg38.fa \
| testsomatic.R \
| var2vcf_paired.pl -A -N "${BASENAME_FILE}|${BASENAME_WBC}" -f 0.002 > "/path/to/output/folder/${BASENAME_FILE}.primary.vcf"
```

Mutations that fall within the targeted regions (ie within targeted genes, not the genome-wide SNP grid) are then 
1) Sorted  
2) Decomposed with vt decompose_blocksub (Hall lab, Atks/vt on GitHub)
3) Annotated using table_annovar.pl from ANNOVAR, incorporating annotations from refGene, COSMIC, GNOMAD, bdSNP, and ClinVar
4) Converted from VCF into TSV with one row per variant using GATK VariantsToTable
5) Further filtered and processed with the code in scripts/vardict_somatic_filter.py. For our analyses, a minimum of 5 reads and 0.25% VAF was required, along with additional filters.

### Mutato Analyze
Mutato analyze was run as below. The sample_sheet.tsv contains patient and sample IDs, along with annotation of sample type: "Cancer" for cfDNA samples and "WBC" for white blood cell controls. "background"  
```
mutato analyze --alt-reads=5 --alt-frac=0.0025 --min-mapq=45 \
               --min-end-dist=15 --bg-ratio=20 --germline-ratio=3 \
               --hotspots="${ref_matti}/mutation_hotspots.tsv" \
               --output-dir="../mutations/mutato_analyze/" \
               /path/to/reference/hg38.fa "/path/to/reference/sample_sheet.tsv"
```
