# Blacklist regions

*Summary*

To remove signal-artefact regions from the sequencing data for downstream analysis, we generated blacklist files for the Atlantic salmon and rainbow trout genomes. First, uint8 mappability files were created for 75, 100, and 150 k-mers across each chromosome in the genome using Umap then all k-mer mappability results were combined. Finally, blacklist regions were defined using Blacklist taking in all available ChIP-seq input BAM files (32 and 29 respective input controls for Atlantic salmon and rainbow trout) and combined uint8 files. Blacklists were generated separately for classic ChIP (BodyMap) and μChIPmentation (DevMap) input controls, with the resulting blacklisted regions merged.
