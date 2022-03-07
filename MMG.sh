#/usr/bin/bash
#by ZF 20210818

#filter non-mitochondrial reads from the raw clean data
#tools NextGenMap, samtools, NOVOPlasty, Seqkit, Vsearch may be used
#type 'bash MMG.sh'


#check packages
echo "Checking packages ......"

#check NextGenMap
if [ $(which ngm) ]
    then
      echo "NextGenMap ...... OK"
      EXE_NGM=$(which ngm)
      DIR_NGM=${EXE_NGM%/*}
    else
      until [ -x $DIR_NGM/ngm ]
        do
          read -p "NextGenMap is not found. Please input its installation directory (absolute path, e.g. /usr/local/bin):      " DIR_NGM_TEMP
          DIR_NGM=$(realpath $(echo $DIR_NGM_TEMP | sed "s/'//g"))
        done
      echo "NextGenMap ...... OK"
fi


#check NOVOPlasty
echo
if [ $(which NOVOPlasty4.3.1.pl) ]
    then
      echo "NOVOPlasty ...... OK"
      EXE_NOVOPLASTY=$(which NOVOPlasty4.3.1.pl)
      DIR_NOVOPLASTY=${EXE_NOVOPLASTY%/*}
    else
      until [ -x $DIR_NOVOPLASTY/NOVOPlasty4.3.1.pl ]
        do
          read -p "NOVOPlasty v4.3 is not found. Please input its installation directory (absolute path, e.g. /home/zf/install/NOVOPlasty-4.3.1):      " DIR_NOVOPLASTY_TEMP
          DIR_NOVOPLASTY=$(realpath $(echo $DIR_NOVOPLASTY_TEMP | sed "s/'//g"))
        done
      echo "NOVOPlasty ...... OK"
fi


#check Samtools
echo
if [ $(which samtools) ]
    then
      echo "SAMTOOLS ...... OK"
      EXE_SAMTOOLS=$(which samtools)
      DIR_SAMTOOLS=${EXE_SAMTOOLS%/*}
    else
      until [ -x $DIR_SAMTOOLS/samtools ]
        do
          read -p "SAMTOOLS is not found. Please input its installation directory (absolute path, e.g. /usr/local/bin):      " SAMTOOLS_TEMP
          DIR_SAMTOOLS=$(realpath $(echo $DIR_SAMTOOLS_TEMP | sed "s/'//g"))
        done
      echo "SAMTOOLS ...... OK"
fi


#check Seqkit
echo
if [ $(which seqkit) ]
    then
      echo "SEQKIT ...... OK"
      EXE_SEQKIT=$(which seqkit)
      DIR_SEQKIT=${EXE_SEQKIT%/*}
    else
      until [ -x $DIR_SEQKIT/seqkit ]
        do
          read -p "SEQKIT is not found. Please input its installation directory (absolute path, e.g. /usr/bin):      " DIR_SEQKIT_TEMP
          DIR_SEQKIT=$(realpath $(echo $DIR_SEQKIT_TEMP | sed "s/'//g"))
        done
      echo "SEQKIT ...... OK"
fi


#check Vsearch
echo
if [ $(which vsearch) ]
    then
      echo "Vsearch ...... OK"
      EXE_VSEARCH=$(which vsearch)
      DIR_VSEARCH=${EXE_VSEARCH%/*}
    else
      until [ -x $DIR_VSEARCH/vsearch ]
        do
          read -p "VSEARCH is not found. Please input its installation directory (absolute path, e.g. /usr/local/bin):      " DIR_VSEARCH_TEMP
          DIR_VSEARCH=$(realpath $(echo $DIR_VSEARCH_TEMP | sed "s/'//g"))
        done
      echo "VSEARCH ...... OK"
fi


#input the name (with its path) of reference mitogenome file
echo
read -p "Please input the name (with its path) of reference mitogenome file (FASTA format) from closely-repated species:      " REF_MITO_TEMP
REF_MITO0=$(realpath $(echo $REF_MITO_TEMP | sed "s/'//g"))
cat $REF_MITO0 | col -b > REFERENCE.fasta
REF_MITO=$(realpath REFERENCE.fasta)
echo ref_mitogenome="$REF_MITO0" >> parameters.txt


#input the name (with its path) of forward Fastq file
echo
read -p "Please input the name (with its path) of forward Fastq file (.fq, .fastq, or gzipped):      " FORWARD_TEMP
FORWARD=$(realpath $(echo $FORWARD_TEMP | sed "s/'//g"))
echo forward_reads="$FORWARD" >> parameters.txt


#input the name (with its path) of reverse Fastq file
echo
read -p "Please input the name (with its path) of reverse Fastq file (.fq, .fastq, or gzipped):      " REVERSE_TEMP
REVERSE=$(realpath $(echo $REVERSE_TEMP | sed "s/'//g"))
echo reverse_reads="$REVERSE" >> parameters.txt


#input the name (with its path) of seed file
echo
read -p "Please input the name (with its path) of seed file (FASTA format) containing all species. The frequently used seed sequences are COIs, e.g. standard 658 bp-COI. Actually, shorter COI sequences often reduce assembly errors, for example, shorten COI to the 300 bp of 3' end. Sequence name in seed FASTA file should be simplfied, such as species name '>Sinella_curviseta':      " SEED_TEMP
SEED0=$(realpath $(echo $SEED_TEMP | sed "s/'//g"))
cat $SEED0 | col -b > SEED_RENAME.fasta
SEED=$(realpath SEED_RENAME.fasta)
echo seeds="$SEED0" >> parameters.txt


#Check the read length
echo
read -p "Please input the read length of sequencing data (e.g. 150, 250):      " READ_LENGTH
  until [ $READ_LENGTH -gt 0 ]
    do
      read -p "Please input the read length of sequencing data (e.g. 150, 250):      " READ_LENGTH
    done
echo read_length="$READ_LENGTH" >> parameters.txt


#Check the threads can be used
echo
read -p "Please input the number of threads/cores (e.g. 8):      " THREADS
  until [ $THREADS -gt 0 ]
    do
      read -p "Please input the correct integer for the number of threads/cores (e.g. 8):      " THREADS
    done
echo threads="$THREADS" >> parameters.txt


DIR_CURR=$(pwd)
rm -rf 0-mitoreads

mkdir 0-mitoreads
cd 0-mitoreads

#map raw sequencing reads to the mitogenome database (e.g. filter_RefSeq.fasta) and transform the resulting file into BAM format
echo
echo "Mapping raw reads to the mitogenome reference ......"
$DIR_NGM/ngm -r $REF_MITO -1 $FORWARD -2 $REVERSE -t $THREADS -i 0.3 -b -o mapping.bam 1>>ngm.log 2>&1


#extract the reads with at least one (forward, reverse or both directions) of a pair of reads which can be aligned to the mitogenome database
$DIR_SAMTOOLS/samtools view -b -f 1 -F 12 -@ $THREADS mapping.bam > map_map.bam
$DIR_SAMTOOLS/samtools view -b -f 4 -F 264 -@ $THREADS mapping.bam > unmap_map.bam
$DIR_SAMTOOLS/samtools view -b -f 8 -F 260 -@ $THREADS mapping.bam > map_unmap.bam


#merge these BAM files
$DIR_SAMTOOLS/samtools merge map.mito.bam map_map.bam unmap_map.bam map_unmap.bam


#transform the BAM format of the candidate mitochondrial reads into interleaved FASTQ format
$DIR_SAMTOOLS/samtools collate -uO map.mito.bam | $DIR_SAMTOOLS/samtools fastq - -1 1.fq -2 2.fq -@ $THREADS
#$DIR_SAMTOOLS/samtools bam2fq map.mito.bam > mito.fq


cd ..
rm *ngm


cd $DIR_CURR
echo "Assembling mitogenomes ......"
rm -rf config.txt 1-assembly species.list SEEDS.fasta

#generate the species list
cat $SEED | grep ">" | sed "s/>//g" > species.list

#assembly of round 1
mkdir -p 1-assembly/round1
cd 1-assembly/round1
SPECIES1=$(head -n 1 ../../species.list)
mkdir $SPECIES1
#modify configure for NOVOPlasty
cp $DIR_NOVOPLASTY/config.txt .
cat $SEED | seqkit grep -r -p "$SPECIES1" -w 0 > $SPECIES1/seed.fasta
sed -i "3s/Test/$SPECIES1/" config.txt
sed -i "10c Seed Input            = $SPECIES1/seed.fasta" config.txt
sed -i "12c Reference sequence    = " config.txt
sed -i "14c Chloroplast sequence  = " config.txt
sed -i "18c Read Length           = "$READ_LENGTH"" config.txt
sed -i "23c Forward reads         = "$DIR_CURR"/0-mitoreads/1.fq" config.txt
sed -i "24c Reverse reads         = "$DIR_CURR"/0-mitoreads/2.fq" config.txt
sed -i "s/\/path\/to\/reads\/reads_[0-9].fastq//g" config.txt
sed -i "25c Store Hash            = yes" config.txt
sed -i "37c Output path           = $SPECIES1/" config.txt

#separately assemble the first species to store HASH files, which can be skipped in subsequent NOVOPasty assemblings
perl $DIR_NOVOPLASTY/NOVOPlasty*pl -c config.txt
echo
sed -i "22c Combined reads        = " config.txt
sed -i "23c Forward reads         = $SPECIES1/HASH2B_$SPECIES1.txt" config.txt
sed -i "24c Reverse reads         = $SPECIES1/HASH2C_$SPECIES1.txt" config.txt
sed -i "25c Store Hash            = $SPECIES1/HASH_$SPECIES1.txt" config.txt

cat ../../species.list | sed '1d' > temp.list
for SPECIES in $(cat temp.list)
  do
    echo
    mkdir $SPECIES
    cat $SEED | seqkit grep -r -p "$SPECIES" -w 0 > $SPECIES/seed.fasta
    sed -i "3c Project name          = "$SPECIES"" config.txt
    sed -i "10c Seed Input            = "$SPECIES"/seed.fasta" config.txt
    sed -i "37c Output path           = $SPECIES/" config.txt
    perl $DIR_NOVOPLASTY/NOVOPlasty*pl -c config.txt
  done
rm temp.list


#deposite assemblies longer than 10000 bp and highly matching seed sequences (identity >=0.95, query cover 100%)

#define a function used for assembly quality assessment
ASSESS_ASSEMBLY_fun() {
echo -e Species"\t""Length""\t""Average organelle coverage\tCircularized" > mitogenome.statistics
for SPECIES in $(cat $1)
  do
    cd $SPECIES
    if [ -s C*fasta ]; then
      cat C*fasta | seqkit replace -p .+ -r "novoplasty_{nr}" --nr-width 3 > temp.fa
      SEED_LENGTH=$(cat seed.fasta | ../seqkit stat -T | sed '1d' | cut -f5)
      ../vsearch --usearch_global seed.fasta --db temp.fa --blast6out - --maxaccepts 0 --maxrejects 0 --id 0.2 > vsearch.out
      LINE=$(cat vsearch.out | wc -l)
      for seq in $(seq $LINE)
        do
          seq_name=$(cut -f2 vsearch.out | sed -n "$seq"p)
          seq_length=$(cut -f10 vsearch.out | sed -n "$seq"p)
          seq_similarity=$(cut -f3 vsearch.out | sed -n "$seq"p)
          seed_cover=$(cut -f4 vsearch.out | sed -n "$seq"p)
          num1=$(echo "scale=1;($seq_similarity-95)"|bc)
          num2=`echo "$num1 >= 0" |bc`
          test "$num2" = 1 -a "$seed_cover" -ge "$SEED_LENGTH" && echo -e "$seq_name""\t""$seq_length" >> candidates     ##这一步仅把相似度大与95%的物种放到candidates中
        done
      if [ -s candidates ]; then
        cat candidates | sort -k2nr | head -n 1 | cut -f1 > seq.list
        cat temp.fa | ../seqkit grep -f seq.list -w 0 | sed "s/*//g" | sed "1c >"$SPECIES"" > novoplasty_"$SPECIES".fasta
        LENGTH=$(cat novoplasty*.fasta | ../seqkit stat -T | sed '1d' | cut -f5)
        COVERAGE=$(cat log_*.txt | grep "Average organelle coverage" | awk '{print $5}')
        test -s Cir*fasta && CIR=C || CIR=""
        test "$LENGTH" -ge 10000 && echo -e "$SPECIES""\t""$LENGTH""\t""$COVERAGE""\t""$CIR" >> ../mitogenome.statistics
        test "$LENGTH" -ge 10000 && cat novoplasty*fasta >> ../novoplasty.fasta
        test "$LENGTH" -lt 10000 && echo "$SPECIES" >> ../list.incomplete
        rm candidates seq.list temp.fa
      else
        echo "$SPECIES" >> ../list.incorrect
      fi
    else
      echo "$SPECIES" >> ../list.incorrect
      rm temp.fa
    fi
    cd ..
  done
cat list.incomplete list.incorrect > list.fail
}
export -f ASSESS_ASSEMBLY_fun

ln -s $DIR_VSEARCH/vsearch .
ln -s $DIR_SEQKIT/seqkit .
cat ../../species.list | ASSESS_ASSEMBLY_fun
rm vsearch seqkit


#check those assembled duplicated sequences and re-assemble failed species
if [ -s list.fail ]; then
    cd .. && mkdir round2 && cd round2
    #check noisy sequences possibly affecting incorrect assemblies
    cat $SEED | $DIR_SEQKIT/seqkit grep -v -f ../round1/list.incorrect -w 0 > temp.seeds.fa
    for species in $(cat ../round1/list.incorrect); do cat ../round1/$species/C*fasta >> incorrect.fa; done
    $DIR_VSEARCH/vsearch --usearch_global incorrect.fa --db temp.seeds.fa --blast6out - --maxaccepts 0 --maxrejects 0 --id 0.85 --threads $THREADS > blast.out
    cat blast.out | cut -f2 | sort -u > list.noisy
    
    #filter noisy reads and re-assemble failed species
    for species in $(cat list.noisy); do cat ../round1/$species/C*fasta >> noisy.fasta; done
    $DIR_NGM/ngm -r noisy.fasta -1 $DIR_CURR/0-mitoreads/1.fq -2 $DIR_CURR/0-mitoreads/2.fq -t $THREADS -i 0.99 --skip-mate-check -b -o mapping.bam 1>>ngm.log 2>&1
    $DIR_SAMTOOLS/samtools view -b -f 12 -F 256 -@ $THREADS mapping.bam > unmap_unmap.bam
    $DIR_SAMTOOLS/samtools collate -uO unmap_unmap.bam | $DIR_SAMTOOLS/samtools fastq - -1 1.fq -2 2.fq -@ $THREADS
    #$DIR_SAMTOOLS/samtools bam2fq unmap_unmap.bam > mito.fq

    
    SPECIES1=$(head -n 1 ../round1/list.fail)
    mkdir $SPECIES1
    cp ../round1/config.txt .
    cat $SEED | $DIR_SEQKIT/seqkit grep -r -p "$SPECIES1" -w 0 > $SPECIES1/seed.fasta
    sed -i "3c Project name          = "$SPECIES1"" config.txt
    sed -i "6c K-mer                 = 23" config.txt
    sed -i "10c Seed Input            = $SPECIES1/seed.fasta" config.txt
    sed -i "11c Extend seed directly  = yes" config.txt
    sed -i "23c Forward reads         = 1.fq" config.txt
    sed -i "24c Reverse reads         = 2.fq" config.txt
    sed -i "25c Store Hash            = yes" config.txt
    sed -i "37c Output path           = $SPECIES1/" config.txt

    perl $DIR_NOVOPLASTY/NOVOPlasty*pl -c config.txt
    echo

    sed -i "22c Combined reads        = " config.txt
    sed -i "23c Forward reads         = $SPECIES1/HASH2B_$SPECIES1.txt" config.txt
    sed -i "24c Reverse reads         = $SPECIES1/HASH2C_$SPECIES1.txt" config.txt
    sed -i "25c Store Hash            = $SPECIES1/HASH_$SPECIES1.txt" config.txt

    cat ../round1/list.fail | sed '1d' > temp.list
    if [ -s temp.list ]; then
      for SPECIES in $(cat temp.list)
        do
          echo
          mkdir $SPECIES
          cat $SEED | seqkit grep -r -p "$SPECIES" -w 0 > $SPECIES/seed.fasta
          sed -i "3c Project name          = "$SPECIES"" config.txt
          sed -i "10c Seed Input            = "$SPECIES"/seed.fasta" config.txt
          sed -i "37c Output path           = $SPECIES/" config.txt
          perl $DIR_NOVOPLASTY/NOVOPlasty*pl -c config.txt
        done
     else
      echo
    fi


    #keep assemblies longer than 10000 bp, seed identity >=0.95 and seed cover 100%
    ln -s $DIR_VSEARCH/vsearch .
    ln -s $DIR_SEQKIT/seqkit .
    cat ../round1/list.fail | ASSESS_ASSEMBLY_fun
    rm -rf vsearch seqkit temp.seeds.fa incorrect.fa *ngm temp.list
    cd ..
  
else
  echo
fi
      
cat round*/novoplasty.fasta > novoplasty.fasta
$DIR_VSEARCH/vsearch --usearch_global $SEED --db novoplasty.fasta --blast6out - --maxaccepts 0 --maxrejects 0 --id 0.9 --threads $THREADS > seeds_vs_mitogenome.txt
cp round1/mitogenome.statistics .
cat round2/mitogenome.statistics | sed '1d' >> mitogenome.statistics


cd $DIR_CURR

#detect chimera
mkdir 2-chimera
cd 2-chimera

cat $DIR_CURR/1-assembly/novoplasty.fasta | grep "^>" | sed "s/>//g" > novoplasty.list
for SPECIES in $(cat novoplasty.list)
  do
    mkdir $SPECIES
    cd $SPECIES
    cat $DIR_CURR/1-assembly/novoplasty.fasta | $DIR_SEQKIT/seqkit grep -r -p "$SPECIES" -w 0 | tr -d '\n' | fold -w 200 | nl -n rz -s " " | sed 's/^/>fragment_/g' | sed 's/ /\n/g' > $SPECIES.fragments.fasta
    cat $DIR_CURR/1-assembly/novoplasty.fasta | $DIR_SEQKIT/seqkit grep -v -r -p "$SPECIES" -w 0 > target.fasta
    $DIR_VSEARCH/vsearch --usearch_global $SPECIES.fragments.fasta --db target.fasta --blast6out - --maxaccepts 0 --maxrejects 0 --id 0.97 --threads $THREADS --blast6out blast.out 1>>vsearch.log 2>&1
    #the number of highly similar fragments more than one indicates the sequences may be a chimera
    num=$(tail -n 1 vsearch.log | awk '{print $5}')
    test "$num" -ge 2 && echo "$SPECIES" >> ../possible_chimera.list
    cd ..
  done

#extracted chimera sequences
if [ -s possible_chimera.list ]; then
  for SPECIES in $(cat possible_chimera.list)
    do
      cd $SPECIES
      cat blast.out | cut -f1 | sort -u > list.fragments
      cat "$SPECIES".fragments.fasta | $DIR_SEQKIT/seqkit grep -f list.fragments -w 0 | $DIR_SEQKIT/seqkit seq -s | sed "1i >"$SPECIES"" | $DIR_SEQKIT/seqkit seq -w 0 >> ../chemera.fasta
      cd ..
    done

  #filter noisy reads and re-assemble failed species
  $DIR_NGM/ngm -r chemera.fasta -1 $DIR_CURR/1-assembly/round2/1.fq -2 $DIR_CURR/1-assembly/round2/2.fq -t $THREADS -i 0.99 --skip-mate-check -b -o mapping.bam 1>>ngm.log 2>&1  $DIR_SAMTOOLS/samtools view -b -f 12 -F 256 -@ $THREADS mapping.bam > unmap_unmap.bam
  $DIR_SAMTOOLS/samtools collate -uO unmap_unmap.bam | $DIR_SAMTOOLS/samtools fastq - -1 1.fq -2 2.fq -@ $THREADS
  #$DIR_SAMTOOLS/samtools bam2fq unmap_unmap.bam > mito.fq

  #assemble
    SPECIES1=$(head -n 1 possible_chimera.list)
    cp ../1-assembly/round1/config.txt .
    cat $SEED | $DIR_SEQKIT/seqkit grep -r -p "$SPECIES1" -w 0 > $SPECIES1/seed.fasta
    sed -i "3c Project name          = "$SPECIES1"" config.txt
    sed -i "6c K-mer                 = 23" config.txt
    sed -i "10c Seed Input            = $SPECIES1/seed.fasta" config.txt
    sed -i "11c Extend seed directly  = yes" config.txt
    sed -i "23c Forward reads         = 1.fq" config.txt
    sed -i "24c Reverse reads         = 2.fq" config.txt
    sed -i "25c Store Hash            = yes" config.txt
    sed -i "37c Output path           = $SPECIES1/" config.txt

    perl $DIR_NOVOPLASTY/NOVOPlasty*pl -c config.txt
    echo

    sed -i "22c Combined reads        = " config.txt
    sed -i "23c Forward reads         = $SPECIES1/HASH2B_$SPECIES1.txt" config.txt
    sed -i "24c Reverse reads         = $SPECIES1/HASH2C_$SPECIES1.txt" config.txt
    sed -i "25c Store Hash            = $SPECIES1/HASH_$SPECIES1.txt" config.txt

    cat possible_chimera.list | sed '1d' > temp.list
    if [ -s temp.list ]; then
      for SPECIES in $(cat temp.list)
        do
          echo
          cat $SEED | seqkit grep -r -p "$SPECIES" -w 0 > $SPECIES/seed.fasta
          sed -i "3c Project name          = "$SPECIES"" config.txt
          sed -i "10c Seed Input            = "$SPECIES"/seed.fasta" config.txt
          sed -i "37c Output path           = $SPECIES/" config.txt
          perl $DIR_NOVOPLASTY/NOVOPlasty*pl -c config.txt
        done
     else
      echo
    fi


    #keep assemblies longer than 10000 bp, seed identity >=0.95 and seed cover 100%
    ln -s $DIR_VSEARCH/vsearch .
    ln -s $DIR_SEQKIT/seqkit .
    cat possible_chimera.list | ASSESS_ASSEMBLY_fun
    rm vsearch seqkit *ngm temp.list
    cd ..

else
  echo
  echo "Congratulations. No chimeras were detected."
fi


echo
echo "Check the file '1-assembly/seeds_vs_mitogenome.txt', which aligned seed sequecnes to the assemblies. The identity (the third column) and the seed cover ratio (the fourth column) may be useful for the final determination."
echo
echo "For the species failed to assemble (listed in 1-assembly/round2/list.fail), carefully scan their NOVOPlasty assembly progress in the folder 1-assembly/roundX/ and encourage to assemble them separately by adjusting parameters, such as k-mer values. Some 'failed' species may be due to their shorter length (< 10,000 bp, see the file 'list.incomplete'), and you can recover them to the 'correct' ones depending on your aims."
echo
echo "Check the file '2-chimera/possible_chimera.list', which listed the possible chimera sequences. Carefully curate them! Aligning results '2-chimera/XXX/blast.out' and new assembly results in the same folder may provide some hints."
echo
echo "All assembled sequences are deposited in 1-assembly/novoplasty.fasta. Don't forget to delete chimera or replace the incorrect sequences!!!"

