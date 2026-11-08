## Somatic mutations
Somatic mutations were called in parallel with Vardict and Mutato Analyze. Only mutations detected by both variant callers were then further curated, filtered, and annotated.

### Vardict
```
        vardict-java \
        -b "${FILENAME}|${FILENAME_WBC}" \
        -f $min_vaf \
        -r $min_reads \
        --nosv \
        -k 0 \
        -N "${BASENAME_FILE}" \
        -th "${SLURM_CPUS_PER_TASK}" \
        -c 1 -S 2 -E 3 ${bed}  \
        -G $reference 2>"$workingdir/logs/vardict/${BASENAME_FILE}.err" \
        | testsomatic.R 2>>"$workingdir/logs/vardict/${BASENAME_FILE}.err" \
        | var2vcf_paired.pl -A -N "${BASENAME_FILE}|${BASENAME_WBC}" -f $min_vaf > "$output/${BASENAME_FILE}.primary.vcf" 2>>"$workingdir/logs/vardict/${BASENAME_FILE}.err" &&

        bedtools intersect \
        -a "$output/${BASENAME_FILE}.primary.vcf" \
        -b /groups/wyattgrp/users/jbacon/reference/bed_files/SNP_spikein_snps.bed \
        -header \
        -wa > "$output/mut_snp_grid/${BASENAME_FILE}.snp.vcf" 2>>"$workingdir/logs/vardict/${BASENAME_FILE}.err" &&

        bedtools intersect \
        -a "$output/${BASENAME_FILE}.primary.vcf" \
        -b /groups/wyattgrp/users/jbacon/reference/bed_files/SNP_spikein_snps.bed \
        -header \
        -v > "$output/mut_somatic/${BASENAME_FILE}.somatic.vcf" 2>>"$workingdir/logs/vardict/${BASENAME_FILE}.err"

```
