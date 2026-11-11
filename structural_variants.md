# Structural variants 

Structural variants (SVs) were detected with two different callers, GRIDSS and Breakfast (https://github.com/annalam/breakfast). All SVs called by either tool were manually inspected to remove artefacts and assess pathogenicity as described in the publication.  

## Breakfast

Detect variants, filter for minimum read depth, and annotate. A description and example of the sample_sheet.tsv can be found in the mutations.md file. 
```
breakfast detect --max-frag-len=1000 --anchor-len=30 sample.bam /path/to/reference/hg38.fa > /project/sv/detect/sample.sv

cd /project/sv/detect/
breakfast filter \
  --min-reads=4 \
  --output-dir="../annotated" \
  /path/to/reference/sample_sheet.tsv

cd /project/sv/annotated/
breakfast annotate sample.sv /groups/wyattgrp/users/jbacon/reference/matti_hg38/ensembl84_genes.bed > sample.sv
```

## GRIDSS 
GRIDSS and GRIPSS were run with parameters as described in Methods.

AR variants were additionally scrutinized to identify the ligand-binding domain truncating structural variants. 
