## Somatic mutations
Somatic mutations were called in parallel with Vardict and Mutato Analyze. Only mutations detected by both variant callers were then further curated, filtered, and annotated.

### Vardict
Vardict is run in paired sample mode, where FILENAME is the cfDNA bam file, and FILENAME_WBC is the patient-matched white blood cell germline control.  
We required a minimum of 4 reads and 0.2% VAF for Vardict calling.  
BASENAME_FILE amd BASENAME_WBC refer to the cfDNA and WBC sample names without the ".bam" extension. 
target_bed_file refers to a bed file of probe locations in our targeted panel, including the genome-wide SNP grid. 

Note that testsomatic.R and var2vcf_paired.pl are standard helper scripts in the vardict-java workflow.  

```
vardict-java \
-b "${FILENAME}|${FILENAME_WBC}" \
-f 0.002 -r 4 --nosv -k 0 -N "${BASENAME_FILE}" \
-c 1 -S 2 -E 3 ${target_bed_file}  \
-G /path/to/reference/hg38.fa \
| testsomatic.R \
| var2vcf_paired.pl -A -N "${BASENAME_FILE}|${BASENAME_WBC}" -f 0.002 > "/path/to/output/folder/${BASENAME_FILE}.primary.vcf"
```
