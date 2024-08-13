# Paired End Read Hakai Adapted Holi Pipeline - SC samples with just refseq mitochondrial database
##### holi-pipeline_SC-MT.md
##### Authors: Libby Natola, Evan Morien
##### Last modified: 25 July 2024

I took out all the explanations as they were distracting for this rerun. refer to holi-pipeline.md for step descriptions. I already ran the read preparation steps in holi-pipeline_SC-NT.md. In lieu of duplicating the work I have linked to the cleaned files for this analysis. Refer to that file for details on read processing. 

### Start a screen

```
screen -S holi-pipeline-SC-MT
```


### Set up conda environment

```
conda activate holi-pipeline
```


### Set up your workspace

```
run='SC-holi-MT'
mkdir /data2/ancient_eDNA/SC_runs/"$run" && cd /data2/ancient_eDNA/SC_runs/"$run"
mkdir raw_fastq
```

Gather your input fastq files

```
ln -s /data2/ancient_eDNA/SC_runs/raw_fastq/*fastq.gz raw_fastq/
bpath=$(pwd)
raw_files_R1="$bpath/raw_fastq/*R1*fastq.gz"
threads=32
mem=50GB
```


### Map the reads

Here I cut out the read processing steps as those are long and we've already done them in holi-pipeline_SC-NT.md. Instead I've symlinked the adap2_kmer2$bname files and run the mapping step to the refseq mitochondrial database only. 


```bash
for infileR1 in $raw_files_R1; do
  infileR2=$(echo $infileR1 | sed 's/_R1_/_R2_/') # get filepath for R2 of same sample
  bnameR1=$(basename "$infileR1" .fastq.gz) # set the R1 basename eg CLO-35-EAP_S107_L001_R1_001
  bnameR2=$(echo $bnameR1 | sed 's/_R1_/_R2_/') # set the R2 basename eg CLO-35-EAP_S107_L001_R2_001
  bname=$(echo $bnameR1 | sed 's/_R1_001//') # set the sample basename eg CLO-35-EAP_S107_L001
  mkdir "$bpath/$bname"_holi && cd "$bpath/$bname"_holi # create and move to the working directory, eg CLO-35-EAP_S107_L001_holi

  # link the cleaned adap2kmer file from nt to here
  ln -s /data2/ancient_eDNA/SC_runs/SC-holi-NT/$bname*/adap2_kmer2_$bname*.pp.rmdup.fastq .

  # map to mitochondrial references only
  for DB in /mnt/Genomics/Working/databases/holi_bowtie2_indexedDBs/refseq/mtDNA/mitochondrion.?
  do
  echo Mapping $bname against $DB
  bowtie2 --threads $threads -k 1000 -x $DB -U adap2_kmer2_$bname*.pp.rmdup.fastq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
  done
    
  cd $bpath
done
```

### Sort and merge bam files

```bash

temp='/data2/scratch-ln'
mkdir $temp
export TMPDIR=$temp

for infileR1 in $raw_files_R1; do
  infileR2=$(echo $infileR1 | sed 's/_R1_/_R2_/') # get filepath for R2 of same sample
  bnameR1=$(basename "$infileR1" .fastq.gz) # set the R1 basename eg CLO-35-EAP_S107_L001_R1_001
  bnameR2=$(echo $bnameR1 | sed 's/_R1_/_R2_/') # set the R2 basename eg CLO-35-EAP_S107_L001_R2_001
  bname=$(echo $bnameR1 | sed 's/_R1_001//') # set the sample basename eg CLO-35-EAP_S107_L001
  cd "$bpath/$bname"_holi # move to the working directory, eg CLO-35-EAP_S107_L001_holi

  for file in *.bam; do
  	echo sorting $file 
    sambamba sort -n -o $file.sorted.bam --tmpdir=$temp -t $threads -m $mem $file
  done

  cd $bpath
done

rm -r $temp
```

Next, run the merging loop. 

```bash
for infileR1 in $raw_files_R1; do
  infileR2=$(echo $infileR1 | sed 's/_R1_/_R2_/') # get filepath for R2 of same sample
  bnameR1=$(basename "$infileR1" .fastq.gz) # set the R1 basename eg CLO-35-EAP_S107_L001_R1_001
  bnameR2=$(echo $bnameR1 | sed 's/_R1_/_R2_/') # set the R2 basename eg CLO-35-EAP_S107_L001_R2_001
  bname=$(echo $bnameR1 | sed 's/_R1_001//') # set the sample basename eg CLO-35-EAP_S107_L001
  cd "$bpath/$bname"_holi # move to the working directory, eg CLO-35-EAP_S107_L001_holi

  echo Merging $bname bam files
  samtools merge --verbosity 5 -f -c -n -@ $threads -O SAM -o $bname.merged.sam.gz "$bpath/$bname"_holi/*sorted.bam 

  cd $bpath
done
```

### metaDMG analysis

metaDMG is a software that estimates, quantifies, and visualizes postmortem damage of aDNA reads AND classifies reads to taxonomy. First get set up by activating the new conda environment (we created this up at the top of the document) and setting some filepaths to access third-party libraries and database files. 

```bash
conda activate metaDMG
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib"
nam=/data/taxonomyDBs/NCBI_taxonomy/2023-03-02/names.dmp
nod=/data/taxonomyDBs/NCBI_taxonomy/2023-03-02/nodes.dmp
acc=/mnt/Genomics/Working/databases/holi_bowtie2_indexedDBs/arctic_animal_genomes/combined_taxid_accssionNO_20230302.gz
```

Run metaDMG.

```bash
mkdir metaDMG/ && cd metaDMG
metaDMG config ../*/*merged.sam.gz --names $nam --nodes $nod --acc2tax $acc --metaDMG-cpp /data2/programs/metaDMG-cpp/metaDMG-cpp --config-file config.metaDMG.$run.yaml --custom-database --bayesian --parallel-samples 1 --cores-per-sample $threads
metaDMG compute config.metaDMG.$run.yaml
metaDMG convert --output ancient_eDNA.$run.csv config.metaDMG.$run.yaml
metaDMG plot config.metaDMG.$run.yaml --output $run.plots.pdf
```

Your output files are now ready for R analysis! 


