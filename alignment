Note that in the code below, multiple steps to calculate quality metrics are interspersed through alignment steps, but have been omitted here for brevity.

## Alignment 
Step 1: Paired end sequencing fastq files were trimmed (adapters, poly-G, low quality 3' ends) with fastp. Create a temp directory to hold intermediate files. The FILENAME variable holds the name of your sample file name, which will stay constant throughout processing.
```
fastp \
-i $fastqdir/${FILENAME}_1.fq.gz \ -I $fastqdir/${FILENAME}_2.fq.gz \
-o $temp_dir/trimmed_fastq/${FILENAME}.trimmed_1.fq.gz \ -O $temp_dir/trimmed_fastq/${FILENAME}.trimmed_2.fq.gz \
--dont_eval_duplication --cut_tail --trim_poly_g --length_required 32 \
--html $workingdir/seq_metrics/html/${FILENAME}.fastq_metrics.html \
--json $workingdir/seq_metrics/json/${FILENAME}.fastq_metrics.json \
--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
--adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
--thread 8
```

Step 2: uBAMs were constructed with fqbio and sorted with sambamba.
```
fgbio --compression 0 --async-io FastqToBam \
  --input $temp_dir/trimmed_fastq/${FILENAME}.trimmed_1.fq.gz $temp_dir/trimmed_fastq/${FILENAME}.trimmed_2.fq.gz \
  --output /dev/stdout --read-structure 3M2S+T 3M2S+T --read-group-id A --umi-tag RX \
  --sample ${FILENAME} --library lib1 --platform illumina \
| sambamba sort -t 4 -m 16G --sort-picard -o $temp_dir/uBAM/${FILENAME}.unmapped.IDT.bam /dev/stdin
```

Step 3: Convert the uBAM back to fastq with picard SamToFastq to align with BWA-MEM2, then picard MergeBamAlignment to combine aligned reads with read information from the uBAM. Sort with sambamba. This is the first alignment, pre-consensus and pre-local realignment steps. 
```
picard SamToFastq \
--INPUT $temp_dir/uBAM/${FILENAME}.unmapped.IDT.bam \
--VALIDATION_STRINGENCY SILENT --COMPRESSION_LEVEL 0 \
--TMP_DIR $temp_dir/SamToFastq/${FILENAME}/ \
--INTERLEAVE true --FASTQ /dev/stdout 2>"$workingdir/logs/${FILENAME}.err" \
| bwa-mem2 mem -t 14 -K 120000000 -p -Y /path/to/reference/genome/hg38 - 2>>"$workingdir/logs/${FILENAME}.err" \
| picard -Xmx16G MergeBamAlignment \
--UNMAPPED $temp_dir/uBAM/${FILENAME}.unmapped.IDT.bam \
--ALIGNED /dev/stdin \
--OUTPUT $temp_dir/mapped_pre_consensus/${FILENAME}.mapped_pre_consensus_unsorted.bam \
--REFERENCE_SEQUENCE /path/to/reference/hg38.fa
--TMP_DIR $temp_dir/MergeBamAlignment/${FILENAME}/ \
--SORT_ORDER unsorted --VALIDATION_STRINGENCY SILENT --CLIP_OVERLAPPING_READS false \
--CLIP_ADAPTERS false --ALIGNER_PROPER_PAIR_FLAGS true --EXPECTED_ORIENTATIONS FR \
--MAX_INSERTIONS_OR_DELETIONS -1 2>>"$workingdir/logs/${FILENAME}.err"

sambamba sort \
--tmpdir $temp_dir/sambamba/${FILENAME}/ --out $temp_dir/mapped_pre_consensus/${FILENAME}.mapped_pre_consensus.bam \
$temp_dir/mapped_pre_consensus/${FILENAME}.mapped_pre_consensus_unsorted.bam 2>>"$workingdir/logs/idt/1-align/${FILENAME}.err"
```

Step 4: Perform base quality score recalibration with GATK and index the recalibrated bam. The known-sites files are a part of the Broad Institute's resource bundle for GATK and can be downloaded here: https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle. 

```
gatk BaseRecalibratorSpark \
-I $temp_dir/mapped_pre_consensus/${FILENAME}.mapped_pre_consensus.bam \
-R /path/to/reference/hg38.fa \
--known-sites /path/to/resources/folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf \
--known-sites /path/to/resources/folder/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
--known-sites /path/to/resources/folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz \
--tmp-dir $temp_dir/BaseRecalibratorSpark/${FILENAME}/ \
--spark-master local[16] \
-O $temp_dir/recal_table/${FILENAME}.recaldata.table

gatk ApplyBQSR  \
--input $temp_dir/mapped_pre_consensus/${FILENAME}.mapped_pre_consensus.bam \
--bqsr-recal-file $temp_dir/recal_table/${FILENAME}.recaldata.table \
--reference /path/to/reference/hg38.fa \
--tmp-dir $temp_dir/ApplyBQSR/${FILENAME}/ \
--output $temp_dir/recal_table/${FILENAME}.recal.bam

sambamba index $temp_dir/recal_table/${FILENAME}.recal.bam
```

Step 5: Perform local realignment with ABRA2 to improve accuracy of indel calling and indel VAF. Then sort, fix read mate information that ABRA2 may have rendered no longer correct, and group reads by position and UMI.
```
abra2 \
--in $temp_dir/recal_table/${FILENAME}.recal.bam \
--out $temp_dir/unmapped_consensus/${FILENAME}.abra2.bam \
--ref /path/to/reference/hg38.fa \
--threads 16 --no-edge-ci --nosort \
--targets "/path/to/reference/exons_of_targeted_panel.bed" \
--tmpdir $temp_dir/ABRA2/${FILENAME}/

sambamba sort \
--sort-picard --compression-level 0 \
--tmpdir $temp_dir/sambamba/${FILENAME}/ --out /dev/stdout \
$temp_dir/unmapped_consensus/${FILENAME}.abra2.bam \
| picard FixMateInformation \
--INPUT /dev/stdin --OUTPUT /dev/stdout --COMPRESSION_LEVEL 0 --VALIDATION_STRINGENCY SILENT \
--TMP_DIR $temp_dir/FixMateInformation/${FILENAME}/ \
| fgbio GroupReadsByUmi \
--input /dev/stdin --output $temp_dir/unmapped_consensus/${FILENAME}.grouped.bam \
--strategy Paired --edits 1 --min-map-q 20 \
--family-size-histogram $temp_dir/seq_metrics/group_reads/${FILENAME}.consensus_metrics
```

Step 6: Collapse reads into duplex consensus families with fgbio, then sort. 
```
fgbio --compression 0 CallDuplexConsensusReads \
--input $temp_dir/unmapped_consensus/${FILENAME}.grouped.bam --output /dev/stdout \
--min-reads 1 1 0  --min-input-base-quality 20 \
--read-name-prefix duplex --consensus-call-overlapping-bases true \
| sambamba sort \
--tmpdir $temp_dir/sambamba/${FILENAME}/ --sort-picard \
-o $temp_dir/unmapped_consensus/${FILENAME}.consensus.unmapped.bam /dev/stdin 
```
Step 7: Realign the consensus reads in the same way as in Step 3, and perform local realignment on the consensus mapped bams as in Step 5, with the addition of the following to filter consensus reads and remove some read tags to decrease storage size. Sort with sambamba afterwards. This is the final bam for downstream analyses. 
```
fgbio -Xmx16G --compression 0 FilterConsensusReads \
--input /dev/stdin --output /dev/stdout \
--ref /path/to/reference/hg38.fa  --min-reads 1 1 0 \
--require-single-strand-agreement true --reverse-per-base-tags true \
--max-read-error-rate 0.025 --min-base-quality 20 --max-no-calls 0.15 --max-base-error-rate 0.1 \
| fgbio -Xmx16G RemoveSamTags \
--input /dev/stdin --output $temp_dir/mapped_consensus/${FILENAME}.filtered.bam \
--tags-to-remove ad bd ae be ac bc aq bq cd ce
```
