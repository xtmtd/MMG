# MMG.sh —— a script for evolutionary and ecological studies using Illumina sequencers for mixed-species samples

## Introduction

Mitochondrial metagenomics (MMG) pipeline for the rapid mitogenome assembly, integrating a fast, accurate read mapper for filtering non-mitochondrial reads, a seed-and-extend assembler for assembling species-specific mitogenomes while detecting ‘noisy’ species/sequences potentially obstructing target assembly.

MMG assembly procedure for each dataset was completed in a few hours on desktop PCs while maintaining high accuracy and completeness except for some very closely related taxa. Excluding ‘noisy’ reads including chimera of non-targeted species could improve the target assembly, particularly for those closely related-species.

![image](https://user-images.githubusercontent.com/45136134/164358722-9d3a569a-223a-4efe-aca7-8b1a494ff201.png)


## Requirements

Some bioinformatic tools are neccessary for above scripts. Most of them are recommended to be added into the environmental paths. Softwares, versions and source ate listed below.

NextGenMap v0.5.5 (https://github.com/Cibiv/NextGenMap/)

SAMtools v1.9 (https://github.com/samtools/samtools)

NOVOPlasty v4.3.1 (https://github.com/ndierckx/NOVOPlasty)

Seqkit v2.0.0 (https://github.com/shenwei356/seqkit)

Vsearch v2.14.2 (https://github.com/torognes/vsearch)


## User manual

 ● Merge all the forward or reverse reads：

    cat *.1.fq | pigz -p $THREADS -c > 1.fq.gz

    cat *.2.fq | pigz -p $THREADS -c > 2.fq.gz

 ● Assemble the mitochondrial genomes using MMG.sh:
    
  Type 'bash MMG.sh'
  
  √ Some tips:
   
   1. The script will check the installation directory of NextGenMap, NOVOPlasty, Samtools, Seqkit and Vsearch automatically. Please input the installation directory (absolute path, e.g. /usr/local/bin) as prompted if is not found.
   
   2. Input the name (with its path) of the reference mitogenome FASTA file used for filtering non-mitochondrial reads. They are usually downloaded from NCBI based on target taxa.
   
   3. Input the name (with its path) of forward Fastq file (e.g. illumina.R1.fq.gz).
   
   4. Input the name (with its path) of reverse Fastq file (e.g. illumina.R2.fq.gz).
   
   5. Input the name (with its path) of seed file, e.g. COI barcode sequences of each species.
   
   6. Input the read length of sequencing data (e.g. 150, 250).
   
   7. Input the number of threads/cores (e.g. 8).
   

  √ Chcek these four notes after the assembly carefully although more details have been described in this script: 

   1. Check the file '1-assembly/seeds_vs_mitogenome.txt', which aligned seed sequecnes to the assemblies. The identity (the third column) and the seed cover ratio (the fourth column) may be useful for the final determination.

   2. For the species failed to assemble (listed in 1-assembly/round2/list.fail), carefully scan their NOVOPlasty assembly progress in the folder 1-assembly/roundX/ and encourage to assemble them separately by adjusting parameters, such as k-mer values. Some 'failed' species may be due to their shorter length (< 10,000 bp, see the file 'list.incomplete'), and you can recover them to the 'correct' ones depending on your aims.

   3. Check the file '2-chimera/possible_chimera.list', which listed the possible chimera sequences. Carefully curate them! Aligning results '2-chimera/XXX/blast.out' and new assembly results in the same folder may provide some hints.

   4. All assembled sequences are deposited in 1-assembly/novoplasty.fasta. Don't forget to delete chimera or replace the incorrect sequences!

## Citation

Du, S.; Dong, J.; Godeiro, N.N.; Wu, J.; Zhang, F. Advancing Mitochondrial Metagenomics: A New Assembly Strategy and Validating the Power of Seed-Based Approach. Diversity 2022, 14, 317. https://doi.org/10.3390/d14050317

## Contact

Please send emails to Dr. Feng Zhang (xtmtd.zf@gmail.com).
