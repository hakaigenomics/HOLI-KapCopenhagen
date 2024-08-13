# Paired End Read Hakai Adapted Holi Pipeline - SC samples with just NT database
##### holi-pipeline_SC-NT.md
##### Authors: Libby Natola, Evan Morien
##### Last modified: 19 July 2024

Authorship note: this file was made by Libby Natola using and modifying the file notes.processing.KWY_cave_sample.sh by Evan Morien, which used and modified the KapCopenhagen version of the Holi ancient DNA pipeline, as detailed here: https://github.com/miwipe/KapCopenhagen

The basic procedure detailed in this file is to prepare and quality check fastq sequences, map reads to references, sort and combine mapped reads, and identify DNA damage and classify taxonomy using LCA. We begin with paired end ancient DNA shotgun sequencing reads and end with metaDMG damage data per taxon. 

I've written the script in a way that it will require less adaptation for unique runs, <span style="color:red">***pay attention to red, bolded, italicized text***</span>, which indicates lines you will have to reformat for your own particular run.

### Start a screen

This code will take a long time to run, and if you're running it in terminal without a screen open you will probably end up disconnecting from fucus and killing the run. I recommend starting a screen and naming it so you will recognize it tomorrow, or in a couple weeks, when you go to check on your job or run the next step. 

```
screen -S holi-pipeline-SC-NT
```

At any point you can detach this screen using `ctl-a` (press these keys at the same time) then `d` (let go of control and a, then press this key alone), and the code will keep running in the background just as you left it. The `-S` argument allows you to name your screens, so if you want to reconnect to your screen after you've detached, run `screen -r`, find the name for your screen (it will have a string of numbers followed by '.holi-pipeline' if you have used the same name for your screen as I have here), and reattach with `screen -r holi-pipeline`. 

### Set up conda environment

**Run this the first time you run this script only**! Once you've set up your conda environment on fucus it should be good to go for the future. This step activates the conda environment as well so you can skip the second code block in this section. You will be prompted to okay some processes, just type 'y'. MetaDMG has some conflicting dependencies with the environment we just set up, so create a second environment for that program, which we will activate later when we get to that step

Set up the conda environments like so:

```
conda create --name holi-pipeline
conda activate holi-pipeline
conda install bioconda::fastq-tools
conda install bioconda::sga
conda activate holi-pipeline
```

```
cat <<EOF > holi-metadmg-environment.yaml
name: metaDMG
channels:
  - conda-forge
dependencies:
  - python==3.9
  - pip
  - pip:
    - metaDMG[all]
EOF

mamba env create --file holi-metadmg-environment.yaml
```


**If you've already set up the conda environments** once before, you still need to activate the environment you made every time you log out of/into fucus. Run this one line each time you begin the script to change to this environment with these packages/versions

```
conda activate holi-pipeline
```


### Set up your workspace

<span style="color:red">***Set a name for your run/analysis***</span>. Here we are setting it as a bash variable so we can add the run name into directories and filenames without having to go in and edit all the file names manually every time we run the loop script below. If you're ever unsure you can test the variable you set with the echo command, eg `echo $run` should print the name you set for the run variable to the screen. Try to limit file names to these characters which aren't special characters as those can have unintended behaviours. These should all be fair game: [a-zA-Z0-9,._+:@%/-]. You want this name to be descriptive of the full run but concise, for example, when running the a51 cave samples, we could use `run='a51_holi'`. 

```
run='SC-holi-NT'
```

<span style="color:red">***Specify the pathway where you want the files for your analysis saved.***</span> Following our example above, we made the directory with $run title inside `/data2/ancient_eDNA/$run`, then moved into it. Make sure to change the filepath in front of `/$run` in both of these lines. In addition, make a directory for the raw_fastq files.

```
mkdir /data2/ancient_eDNA/SC_runs/"$run" && cd /data2/ancient_eDNA/SC_runs/"$run"
mkdir raw_fastq
```

<span style="color:red">***Gather your input fastq files.***</span> These input files can be large, so copying them to a more convenient place on fucus takes up lots of space. I only recommend copying your fastq files to this workspace if they were not sequenced by Hakai, and if other people have access to them and might accidentally delete or overwrite your raw files. If your raw data are on the NAS (most cases at Hakai), you can use a symlink so it's easier to navigate to the files without taking up the diskspace to have those in two places. To do so, replace `/path/to/fastqs/` with your absolute path to the raw fastq files on the NAS. 

```
ln -s /data2/ancient_eDNA/SC_runs/raw_fastq/*fastq.gz raw_fastq/
```


Set the variable for the pathways to the raw R1 files. List just the R1s (this is arbitrary, it could just as easily be the R2s). We need to run the loop once for every SAMPLE (every sample has two files) because within the loop we need to concatenate the lists from both forward and reverse reads, so we need to have already processed both the forward and reverse fastq files for that sample in the same iteration of the loop. The R2s will be called in the loop as well by modifying the filepaths we set here for the R1s.


```
bpath=$(pwd)
raw_files_R1="$bpath/raw_fastq/*R1*fastq.gz"
```

Set your max computational resources, note that fucus has 48 threads and 1 TB of RAM. Always leave a minimum of 10GB of RAM and 1 CPU (i.e. 2 threads) free, so the operating system has resources to use while your job is running. Failure to do so will result in considerable system slowdown, and possible fatal error (loss of computation time for your job) and system crash (loss of computation time for everyone else). 

```
threads=32
mem=50GB
```


### Clean the data

This step runs a loop to produce fastq files without bacterial reads, polyA tails, adapter sequence, low-complexity sequences, short sequences, or duplicates. It cycles through each of the raw R1 files, sets variables for the input file pathways, sets basenames without filepaths for the F and R reads, and sets a basename that is applicable to the sample (no R1/R2 text), and makes new output directories for that sample. It is important that your files follow the normal Hakai naming conventions wherein forward read fastq basenames end in "R1_001" and reverse read fastq basenames end in "R2_001" or the search and replace code here won't work correctly. The loop then starts by converting F and R fastqs to fasta, which is the file format required for the next step, blast-n. 

The loop then uses blast-n to align reads against bacterial databases, it uses a variable to call the number of threads to specify for blastn, so make sure you've set that variable above. It runs with the "megablast" task. It limits output to the first 5 matches with a maximum e-value of 0.00001, a percent identity of at least 98%, and a query coverage of at least 90%. We also specify the bacterial reference database to use, the output formatting, and the output filepath. This step is run first on the R1 file and then on R2. 

For the next step in the loop we compile a list of the reads that the blast step identified as bacterial, then concatenate the bacterial reads from the F and R files together into a master list of reads to exclude from our final fastqs. Next, new versions of the raw files are made without the reads identified as bacteria. The next step is to remove unneccesary intermediate files that will clutter up your workspace. 

Then the loop moves on with the read trimming. It uses the program fastq-grep to remove poly-A tails (and polyTs for reverse complements) and adapter sequences. With SGA, we remove sequences with a DUST complexity score above 1 and reads shorter than 30 bp. SGA is also used to index and filter the reads for duplicates. Finally, the read length distribution is calculated. 

By copying and pasting this full loop into your terminal, you're running the bash loop interactively (not from a script file). This means you don't have to specify the shebang (#!) at the beginning of the code or make anything executable. The lines that begin with 'echo' are designed to print to the screen so you can tell what stages of the loop are running, which can be helpful if you need to do any troubleshooting.


```bash
for infileR1 in $raw_files_R1; do
  infileR2=$(echo $infileR1 | sed 's/_R1_/_R2_/') # get filepath for R2 of same sample
  bnameR1=$(basename "$infileR1" .fastq.gz) # set the R1 basename eg CLO-35-EAP_S107_L001_R1_001
  bnameR2=$(echo $bnameR1 | sed 's/_R1_/_R2_/') # set the R2 basename eg CLO-35-EAP_S107_L001_R2_001
  bname=$(echo $bnameR1 | sed 's/_R1_001//') # set the sample basename eg CLO-35-EAP_S107_L001

  # make a folder for the sample's output files to go into, cd into it, also make a filt/ folder for the filtered samples to go into
  mkdir "$bpath/$bname"_holi && cd "$bpath/$bname"_holi # create and move to the working directory, eg CLO-35-EAP_S107_L001_holi
  mkdir filt

  echo Preprocessing. Removing bacterial reads
  # convert to fasta
  seqtk seq -a $infileR1 > filt/$bnameR1.fasta 
  seqtk seq -a $infileR2 > filt/$bnameR2.fasta
  
  echo Preprocessing. blastn: sample against refseq bacterial genomes
  # blast to find reads that map to bacteria, R1 and R2
  blastn -task megablast -num_threads $threads -evalue 1e-5 -max_target_seqs 5 -perc_identity 98 -qcov_hsp_perc 90 -db /mnt/Genomics/Working/databases/holi_bowtie2_indexedDBs/refseq/bacteria_blastDB/ref_prok_rep_genomes -outfmt '6 qseqid stitle sacc pident qcovs evalue bitscore' -query filt/$bnameR1.fasta -out filt/$bnameR1.blast_out
  blastn -task megablast -num_threads $threads -evalue 1e-5 -max_target_seqs 5 -perc_identity 98 -qcov_hsp_perc 90 -db /mnt/Genomics/Working/databases/holi_bowtie2_indexedDBs/refseq/bacteria_blastDB/ref_prok_rep_genomes -outfmt '6 qseqid stitle sacc pident qcovs evalue bitscore' -query filt/$bnameR2.fasta -out filt/$bnameR2.blast_out
  
  echo Preprocessing. Compiling list of reads to remove.
  # make list of the names of those reads
  awk -F"\t" '{print $1}' filt/$bnameR1.blast_out | sort | uniq > filt/$bnameR1.bacterial_reads
  awk -F"\t" '{print $1}' filt/$bnameR2.blast_out | sort | uniq > filt/$bnameR2.bacterial_reads
  
  # concatenate the R1 and R2 lists, AdapterRemoval won't work if the read is missing from only one of the fastqs
  cat filt/$bnameR1.bacterial_reads filt/$bnameR2.bacterial_reads | sort | uniq > filt/$bname.bacterial_reads

  echo "Preprocessing. removing reads with >=98% identity with refseq bacterial genomes"
  #filter fastq based on identified bacterial reads above
  /data2/programs/bbmap/filterbyname.sh in=$infileR1 out=filt/$bnameR1.fastq names=filt/$bname.bacterial_reads overwrite=true
  /data2/programs/bbmap/filterbyname.sh in=$infileR2 out=filt/$bnameR2.fastq names=filt/$bname.bacterial_reads overwrite=true
  
  # clean up unnecessary files
  # rm filt/*.fasta # optionally run this if you don't want to hold on to any of the bacterial reads data 
  gzip filt/*.blast_out &
  gzip filt/*.bacterial_reads &

  # Read Trimming
  echo Step 1. Removing poly A tails
  fastq-grep -v "AAAAA$" filt/$bnameR1.fastq > kmer_$bnameR1.fastq
  fastq-grep -v "AAAAA$" filt/$bnameR2.fastq > kmer_$bnameR2.fastq

  echo Step 2. Removing reverse complemented A tails
  fastq-grep -v "^TTTTT" kmer_$bnameR1.fastq > kmer2_$bnameR1.fastq
  fastq-grep -v "^TTTTT" kmer_$bnameR2.fastq > kmer2_$bnameR2.fastq

  echo Step 3. Removing rememnants adapter sequence 1 = AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
  fastq-grep -v "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" kmer2_$bnameR1.fastq > adap1_kmer2_$bnameR1.fastq
  fastq-grep -v "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" kmer2_$bnameR2.fastq > adap1_kmer2_$bnameR2.fastq
  
  echo Step 4. Removing remnants adapter sequence 2 = ATCTCGTATGCCGTCTTCTGCTTG
  fastq-grep -v "ATCTCGTATGCCGTCTTCTGCTTG" adap1_kmer2_$bnameR1.fastq > adap2_kmer2_$bnameR1.fastq
  fastq-grep -v "ATCTCGTATGCCGTCTTCTGCTTG" adap1_kmer2_$bnameR2.fastq > adap2_kmer2_$bnameR2.fastq

  echo Step 6. sga preprocessing: remove low complexity and short reads
  sga preprocess --dust-threshold=1 -m 30 adap2_kmer2_$bnameR1.fastq -o adap2_kmer2_$bnameR1.pp.fastq 
  sga preprocess --dust-threshold=1 -m 30 adap2_kmer2_$bnameR2.fastq -o adap2_kmer2_$bnameR2.pp.fastq 

  echo Step 7. sga index reads
  sga index --algorithm=ropebwt --threads=30 adap2_kmer2_$bnameR1.pp.fastq 
  sga index --algorithm=ropebwt --threads=30 adap2_kmer2_$bnameR2.pp.fastq

  echo Step 8. sga filter to remove duplicates
  sga filter --threads=30  --no-kmer-check adap2_kmer2_$bnameR1.pp.fastq -o adap2_kmer2_$bnameR1.pp.rmdup.fastq 
  sga filter --threads=30  --no-kmer-check adap2_kmer2_$bnameR2.pp.fastq -o adap2_kmer2_$bnameR2.pp.rmdup.fastq 

  echo Step 9. Calculating read length distribution and outputting file
  cat adap2_kmer2_$bnameR1.pp.rmdup.fastq | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > adap2_kmer2_$bnameR1.pp.rmdup.fastq.read_length.txt
  cat adap2_kmer2_$bnameR2.pp.rmdup.fastq | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > adap2_kmer2_$bnameR2.pp.rmdup.fastq.read_length.txt

  # clean up temporary files
  rm kmer*
  rm adap1*
  rm *wt
  rm *sai
  rm *rmdup.discard.fa
  rm adap2_kmer2_*001.fastq
  rm adap2_kmer2_*001.pp.fastq
  rm filt/*.fasta

  # return to the base path so the directory for the next sample goes to the correct path
  cd $bpath

done
```

### Map the reads

Now that we have cleaned the bacteria, polyA-tails, adapters, low quality reads, and duplicates from our data, we can start mapping the remaining reads to databases in order to determine possible sources for the sequences. We have several databases to choose from, and if you're interested in a particular taxon it will reduce your running time by quite a bit to map your reads only to the databases containing your taxa of choice (eg, if you're only interested in mammals select the vert_mammal database) or organelle of your choice (eg, you can use a mitochondrial DNA database if you have reads enriched for mitochondrial reads).  

Again, this is performed through a large loop that runs several processes for each sample, looping through all the samples present in the raw data directory. Within that loop, we prepare the environment by setting the variables for the file currently being handled and moving into the sample's directory. Then, because each database is composed of multiple different fasta files, we loop through each database file, print the sample name and fasta file for debugging help, and map both the forward and reverse cleaned reads files to the current fasta with the software bowtie2. We set the threads to our specifications (remember to set the $threads variable as detailed above), report up to 1000 valid, distinct alignments for each read (in descending order by alignment score), and discard reads with no alignments. The output sam file is then from standard input and saved as a bamfile. When the loop has completed mapping reads to all available fasta files for a specified database the loop for the following database is begun. This continues until reads are mapped to all databases, then the outermost loop returns to the base directory and continues to the next sample. 

<span style="color:red">***To change the databases used***</span> remove blocks calling any unwanted databases below (delete from "for DB in..." to the next "done"), or substitute the pathways to the fasta files of your choice at the beginning of the mapping loop ("for DB in <file-path>"). Be sure to use only bowtie2 indexed files. 

```bash
for infileR1 in $raw_files_R1; do
  infileR2=$(echo $infileR1 | sed 's/_R1_/_R2_/') # get filepath for R2 of same sample
  bnameR1=$(basename "$infileR1" .fastq.gz) # set the R1 basename eg CLO-35-EAP_S107_L001_R1_001
  bnameR2=$(echo $bnameR1 | sed 's/_R1_/_R2_/') # set the R2 basename eg CLO-35-EAP_S107_L001_R2_001
  bname=$(echo $bnameR1 | sed 's/_R1_001//') # set the sample basename eg CLO-35-EAP_S107_L001
  cd "$bpath/$bname"_holi # move to the working directory, eg CLO-35-EAP_S107_L001_holi

  # map to NT references only
  for DB in /mnt/Genomics/Working/databases/holi_bowtie2_indexedDBs/NT/nt.?
  do
  echo Mapping $bname against $DB
  bowtie2 --threads $threads -k 1000 -x $DB -U adap2_kmer2_$bname*.pp.rmdup.fastq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
  done
  
  for DB in /mnt/Genomics/Working/databases/holi_bowtie2_indexedDBs/NT/nt.??
  do
  echo Mapping $bname against $DB
  bowtie2 --threads $threads -k 1000 -x $DB -U adap2_kmer2_$bname*.pp.rmdup.fastq --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
  done
  
  cd $bpath
done
```

### Sort and merge bam files

The next step is to sort all bam files with the read alignments we have just created and merge them into one single sam file. Sorting the files allows downstream analyses to access data for specific regions of the genome without having to read the entire file which increases efficiency. Merging the bams into one sam allows us to compare all mappings for each read in each sample across multiple databases.


The program sambamba requires a lot of temporary file space, and can sometimes exceed the fucus /tmp space. <span style="color:red">***Make a variable with the pathway for your own temp directory***</span> in /data2 with your name or initials on it. For example, I've used  `temp='/data2/scratch-ln'`. Then create the directory and export it so any processes in this environment will use it as the temp directory. 

```
temp='/data2/scratch-<your-initials-here>'
mkdir $temp
export TMPDIR=$temp
```

Run the sorting loop. It will proceed through each sample, each time running another loop that sorts each bam file by read name using the specified threads and memory settings. Be sure to set the `--tmpdir` parameter as well as the `-m` parameter taking into account the size of the bam files you are sorting. `-m 100MB` for block size seems adequate for bam files of 500MB or less in size. 10MB was only adequate up to ~70MB bam files. 

```bash
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

Next, run the merging loop. It will read in the read name-sorted bam files and output a gzipped sam file with detailed output to the terminal screen. 

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

metaDMG works by first creating a config file containing filepaths to the input and database files, which it will use in subsequent steps when you're running the software using the compute command. Here we specify the NCBI names.dmp and nodes.dmp files as well as the file mapping NCBI accession numbers to taxon ids. We add the pathway for the metaDMG-cpp files, name the config file with our run name (set earlier in this document), and specify that we are using a custom database, we are running the full Bayesian model, one sample at a time, with the number of threads set in our variable. Then this config file is used to run the analysis. Finally, the results are converted to a combined csv for all samples, and generates many plots in one single pdf file. 

```bash
mkdir metaDMG/ && cd metaDMG
metaDMG config ../*/*merged.sam.gz --names $nam --nodes $nod --acc2tax $acc --metaDMG-cpp /data2/programs/metaDMG-cpp/metaDMG-cpp --config-file config.metaDMG.$run.yaml --custom-database --bayesian --parallel-samples 1 --cores-per-sample $threads
metaDMG compute config.metaDMG.$run.yaml
metaDMG convert --output ancient_eDNA.$run.csv config.metaDMG.$run.yaml
metaDMG plot config.metaDMG.$run.yaml --output $run.plots.pdf
```

Your output files are now ready for R analysis! 


