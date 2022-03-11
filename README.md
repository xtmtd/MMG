# MMG —— a tool for evolutionary and ecological studies using Illumina sequencers for mixed-species samples

## Introduction

Mitochondrial metagenomics (MMG) pipeline for the rapid mitogenome assembly, integrating a fast, accurate read mapper for filtering non-mitochondrial reads, a seed-and-extend assembler for assembling species-specific mitogenomes while detecting ‘noisy’ species/sequences potentially obstructing target assembly.

MMG assembly procedure for each dataset was completed in a few hours on desktop PCs while maintaining high accuracy and completeness except for some very closely related taxa. Excluding ‘noisy’ reads including chimera of non-targeted species could improve the target assembly, par-ticularly for those closely-related species. Read abundance was correctly assessed for all tested species except for very closely related ones (COI divergence 1.5%). Short barcodes as the reference can have almost identical detection power but require at least an order of magnitude greater sequencing depth than mitogenomes. Sequencing amount of 1 Gbp per bulk sample is usually sufficient to detect species richness and abundance against mitogenome reference.

![image](https://user-images.githubusercontent.com/45136134/157005857-7e00689b-0d7a-4009-993b-9162a634420a.png)


## Requirements

Some bioinformatic tools are neccessary for above scripts. Most of them are recommended to be added into the environmental paths. Softwares, versions and source ate listed below.

NextGenMap v0.5.5 (https://github.com/Cibiv/NextGenMap/)

SAMtools v1.9 (https://github.com/samtools/samtools)

NOVOPlasty v4.3.1 (https://github.com/ndierckx/NOVOPlasty)

Seqkit v2.0.0 (https://github.com/shenwei356/seqkit)

Vsearch v2.14.2(https://github.com/torognes/vsearch)


## User manual

Type 'bash MMG.sh'

Chcek these four notes after the assembly carefully although more details have been described in this script: 

  1. Check the file '1-assembly/seeds_vs_mitogenome.txt', which aligned seed sequecnes to the assemblies. The identity (the third column) and the seed cover ratio (the fourth column) may be useful for the final determination.

  2. For the species failed to assemble (listed in 1-assembly/round2/list.fail), carefully scan their NOVOPlasty assembly progress in the folder 1-assembly/roundX/ and encourage to assemble them separately by adjusting parameters, such as k-mer values. Some 'failed' species may be due to their shorter length (< 10,000 bp, see the file 'list.incomplete'), and you can recover them to the 'correct' ones depending on your aims.

  3. Check the file '2-chimera/possible_chimera.list', which listed the possible chimera sequences. Carefully curate them! Aligning results '2-chimera/XXX/blast.out' and new assembly results in the same folder may provide some hints.

  4. All assembled sequences are deposited in 1-assembly/novoplasty.fasta. Don't forget to delete chimera or replace the incorrect sequences!


## Contact

Please send emails to Dr. Feng Zhang (xtmtd.zf@gmail.com).
