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
Mutato analyze was run as below. The sample_sheet.tsv contains patient and sample IDs, along with annotation of sample type: "Cancer" for cfDNA samples and "WBC" for white blood cell controls. The additional presence of the string "background" (case-insensitive) marks samples to be used as control samples to calculate cohort-wide background error for each mutation called. In the ENZA-p analyses, the set of controls were 21 ENZA-p samples where no evidence of ctDNA was observed (no called somatic mutations, copy number changes both focal or chromosomal arm level, or structural variants).  
```
mutato analyze --alt-reads=5 --alt-frac=0.0025 --min-mapq=45 \
               --min-end-dist=15 --bg-ratio=20 --germline-ratio=3 \
               --hotspots="${ref_matti}/mutation_hotspots.tsv" \
               --output-dir="../mutations/mutato_analyze/" \
               /path/to/reference/hg38.fa "/path/to/reference/sample_sheet.tsv"
```
The sample sheet should be formatted like the following example:
```
PATIENT	SAMPLE	TYPE
Patient-001	Patient-001_cfDNA	Cancer
Patient-001	Patient-001_WBC	WBC
Patient-002	EnzaP-002_cfDNA  Cancer background
Patient-002	Patient-002_WBC	WBC
```

The mutations can then be annotated with mutato_analyze_tools annotate, inputting GNOMAD, COSMIC, and ClinVar database files. The tool is included here in the julia folder.  
Then, somatic and germline mutations are called separately.
'''
mutato_analyze_tools somatic --protein_altering=false --ox-filter=false --min-germline-reads=20 > "mutato_analuze_somatic_include_noncoding.tsv"
mutato_analyze_tools germline defects --min-germline-reads=50 > "mutato_analyze_germline.tsv"
'''
The output is then further filtered and annotated as described in the publication, using custom python scripts.   

The two sets of mutations from Vardict and R are intersected to identify mutations called by both, and those are manually curated (facilitated with IGV screenshots) to remove artefact and putative variants derived from clonal hematopoiesis.

