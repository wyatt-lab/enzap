# Enza-P circulating tumour DNA analyses

This repository contains the custom code used to analyse copy number profiles, somatic variants, and structural variants in the paper "Circulating tumour DNA fraction and androgen receptor genotype predict benefit for addition of 177Lu-PSMA-617 to enzalutamide in advanced prostate cancer in the ENZA-p trial". 

Access to patient sequencing data in the EnzaP trial is restricted but may be shared upon request (see Data Availability section of paper).

## Pre-requisites
Our analysis software and scripts are written using the Julia and Rust programming languages. To use our code, you will first need to:

1. Install Julia version 1.8.5: [https://julialang.org/downloads/](https://julialang.org/downloads/oldreleases/)  
2. Install the Rust compiler and toolchain: https://www.rust-lang.org/tools/install  
3. Add the custom Julia modules in subfolder `julia_modules` included here to a folder on your system
4. Add your folder to the LOAD_PATH, and install required Julia packages from inside the Julia terminal:

```
push!(LOAD_PATH, "/path/to/your/folder")
using Pkg;
Pkg.add(split("Printf Random StatsBase DelimitedFiles Statistics Distributions CSV Plot2 PyPlot PDFIO FastaIO HypothesisTests BlackBoxOptim PyCall KernelDensity Interpolations Combinatorics Suppressor"))
```
5. Install Mutato: https://www.github.com/annalam/mutato
6. Install Breakfast: https://www.github.com/annalam/breakfast
7. Install ANNOVAR: https://annovar.openbioinformatics.org/en/latest/user-guide/download/
8. Install VarDictJava, GRIDSS, and GRIPSS in your conda environment

Note that in the code below, mkdir commands to create referenced folders are omitted for brevity.

## Alignment 
Paired end sequencing fastq files were trimmed (adapters, poly-G, low quality 3' ends) with fastp. Create a temp directory to hold intermediate files. The FILENAME variable holds the name of your sample file name, which will stay constant throughout processing.
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

uBAMs were constructed with fqbio and sorted with sambamba.
```
fgbio --compression 0 --async-io FastqToBam \
  --input $temp_dir/trimmed_fastq/${FILENAME}.trimmed_1.fq.gz $temp_dir/trimmed_fastq/${FILENAME}.trimmed_2.fq.gz \
  --output /dev/stdout --read-structure 3M2S+T 3M2S+T --read-group-id A --umi-tag RX \
  --sample ${FILENAME} --library lib1 --platform illumina \
| sambamba sort -t 4 -m 16G --sort-picard -o $temp_dir/uBAM/${FILENAME}.unmapped.IDT.bam /dev/stdin
```

Convert the uBAM back to fastq with picard SamToFastq to align with BWA-MEM2, then picard MergeBamAlignment to combine aligned reads with read information from the uBAM. Sort with sambamba. This is the first alignment, pre-consensus and pre-local realignment steps. 
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

Perform base quality score recalibration with GATK and index the recalibrated bam. The known-sites files are a part of the Broad Institute's resource bundle for GATK and can be downloaded here: https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle. 

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

Perform local realignment with ABRA2 to improve accuracy of indel calling and indel VAF. Then sort, fix read mate information that ABRA2 may have rendered no longer correct, and group reads by position and UMI.
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

Collapse reads into duplex consensus families with fgbio, then sort. 
```
fgbio --compression 0 CallDuplexConsensusReads \
--input $temp_dir/unmapped_consensus/${FILENAME}.grouped.bam --output /dev/stdout \
--min-reads 1 1 0  --min-input-base-quality 20 \
--read-name-prefix duplex --consensus-call-overlapping-bases true \
| sambamba sort \
--tmpdir $temp_dir/sambamba/${FILENAME}/ --sort-picard \
-o $temp_dir/unmapped_consensus/${FILENAME}.consensus.unmapped.bam /dev/stdin 

fgbio CollectDuplexSeqMetrics \
--input $temp_dir/unmapped_consensus/${FILENAME}.grouped.bam \
--output $temp_dir/seq_metrics/consensus/${FILENAME}/${FILENAME}
```























