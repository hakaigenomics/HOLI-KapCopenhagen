# Bacterial read removal for metagenomic pipelines
##### Author: Libby Natola
##### Last modified: 7 June, 2024

This file inclues instructions and code to remove bacteria from shotgun reads prior to a bioinformatic pipeline designed for metagenomics. We have done this for two main reasons. The first is that bacterial contamination is often present in non-bacterial reference sequences, which causes false species assingments. For example, there is lots of Achromobacter xylosoxidans contamination in the turkey genome on RefSeq. If you have A. xylosoxidans in your sample it would map to reads in the turkey genome and your taxonomic binner will tell you you have turkey DNA. If you're not interested in the bacteria at all you can run this step to remove that risk. Second, removing the bacterial taxa we aren't interested in will pare down the file sizes and run your subsequent analyses more efficiently. 

This step is taken from the holi pipeline script prepared by Evan Morien, which he adapted from KapCopenhagen scripts. Note that this step is included in the holi pipeline so you don't need to run this prior to that pipeline. I've written the script in a way that it will require less adaptation for unique runs, <span style="color:red">***pay attention to red, bolded, italicized text***</span>, which indicates lines you will have to reformat for your own particular run.

### Start a screen

This code can take a long time to run, and if you're running it in terminal without a screen open you will probably end up disconnecting from fucus and killing the run. I recommend starting a screen and naming it so you will recognize it tomorrow when you go to check on your run. 

```
screen -S remove-bacteria
```

At any point you can detach this screen using `ctl-a` (press these keys at the same time) then `d` (let go of control and a, then press this key alone), and the code will keep running in the background just as you left it. The `-S` argument allows you to name your screens, but the screen program will still append a pid to it, so if you want to reconnect to your screen after you've detached, run `screen -r`, find the full name for your screen (it will have a string of numbers followed by '.remove-bacteria'), and reattach with `screen -r remove-bacteria`. 

### Set up conda environment

**Run this the first time you run this script only**! Once you've set up your conda environment on fucus it should be good to go for the future. This step activates the conda environment as well so you can skip the second code block in this section. You will be prompted to okay some processes, just type 'y'. 

Set up a conda environment like so:

```
conda create --name remove-bacteria
conda activate remove-bacteria
conda install bioconda::fastq-tools
```

**If you've already set up the conda environment** once before, you still need to activate the environment you made. Run this one line each time you begin the script to change to this environment with these packages/versions

```
conda activate remove-bacteria
```

### Set up your workspace

<span style="color:red">***Set a name for your run/analysis***</span>. Here we are setting it as a bash variable so we can add the run name into directories and filenames without having to go in and edit all the file names manually every time we run the loop script below. If you're ever unsure you can test the variable you set with the echo command, eg `echo $run` should print the name you set for the run variable to the screen. Try to limit file names to these characters which aren't special characters as those can have unintended behaviours. These should all be fair game: [a-zA-Z0-9,._+:@%/-]. You want this name to be descriptive of the full run but concise, for example, when running the a51 cave samples, we could use `run='a51'`. 

```
run='test_analysis'
```

<span style="color:red">***Specify the pathway where you want the files for your analysis saved.***</span> Following our example above, we made the directory with $run title inside `/data2/ancient_eDNA/$run`, then moved into it. Make sure to change the filepath in front of `/$run` in both of these lines. In addition, make a directory for the raw_fastq files.

```
mkdir /<desired/filepath>/$run && cd /<desired/filepath>/$run
mkdir raw_fastq
```

<span style="color:red">***Gather your input fastq files.***</span> These input files can be large, so copying them to a more convenient place on fucus takes up lots of space. I only recommend copying your fastq files to this workspace if they were not sequenced by Hakai, and if other people have access to them and might accidentally delete or overwrite your raw files. If your raw data are on the NAS (most cases at Hakai), you can use a symlink so it's easier to navigate to the files without taking up the diskspace to have those in two places. To do so, replace `/path/to/fastqs/` with your absolute path to the raw fastq files on the NAS. 

```
ln -s /path/to/fastqs/*.fastq.gz raw_fastq/
```

From here on out, if you followed the steps above the following steps will be universal so you shouldn't have to edit the filenames or filepaths and you should be able to copy and paste the code. 

Set the variable for the pathways to the raw R1 files. List just the R1s (this is arbitrary, it could just as easily be the R2s). We need to run the loop once for every SAMPLE (every sample has two files) because within the loop we need to concatenate the lists from both forward and reverse reads, so we need to have already processed both the forward and reverse fastq files for that sample in the same iteration of the loop. The R2s will be called in the loop as well by modifying the filepaths we set here for the R1s.


```
raw_files_R1='raw_fastq/*R1*.fastq.gz' 
```

Make an output directory for the resulting fastq files with the bacteria removed.

```
mkdir bacteria-removed_fastqs
```

### Bacterial removal loop

This step runs a loop to produce fastq files with bacterial reads removed. It cycles through each of the raw R1 files, sets variables for the input file pathways, sets basenames without filepaths for the F and R reads, and sets a basename that is applicable to the sample (no R1/R2 text). It is important that your files follow the normal Hakai naming conventions wherein forward read fastq basenames end in "R1_001" and reverse read fastq basenames end in "R2_001" or the search and replace code here won't work correctly. The loop then starts by converting F and R fastqs to fasta, which is the file format required for the next step, blast-n. 

The loop then uses blast-n to align reads against bacterial databases, it uses 38 threads which is most of what is available on fucus, so make sure no one else is running lots of threads at the moment. It runs with the "megablast" task. It limits output to the first 5 matches with a maximum e-value of 0.00001, a percent identity of at least 98%, and a query coverage of at least 90%. We also specify the reference database to use, the output formatting, and the output filepath. This step is run first on the R1 file and then on R2. 

For the next step in the loop we compile a list of the reads that the blast step identified as bacterial, then concatenate the bacterial reads from the F and R files together into a master list of reads to exclude from our final fastqs. Next, new versions of the raw files are made without the reads identified as bacteria. The last step is to remove unneccesary intermediate files that will clutter up your workspace. 

By copying and pasting this full loop into your terminal, you're running the bash loop interactively (not from a script file). This means you don't have to specify the shebang (#!) at the beginning of the code or make anything executable. The lines that begin with 'echo' are designed to print to the screen so you can tell what stages of the loop are running, which can be helpful if you need to do any troubleshooting.

```bash
for infileR1 in $raw_files_R1
do
  infileR2=$(echo $infileR1 | sed 's/_R1_/_R2_/') # get filepath for R2 of same sample
  bnameR1=$(basename "$infileR1" .fastq.gz) # set the R1 basename eg CLO-35-EAP_S107_L001_R1_001
  bnameR2=$(echo $bnameR1 | sed 's/_R1_/_R2_/') # set the R2 basename eg CLO-35-EAP_S107_L001_RR_001
  bname=$(echo $bnameR1 | sed 's/_R1_001//') # set the sample basename eg CLO-35-EAP_S107_L001

  echo Preprocessing. Removing bacterial reads
  # convert to fasta
  seqtk seq -a $infileR1 > bacteria-removed_fastqs/$bnameR1.fasta 
  seqtk seq -a $infileR2 > bacteria-removed_fastqs/$bnameR2.fasta
  
  echo Preprocessing. blastn: sample against refseq bacterial genomes
  # blast to find reads that map to bacteria, R1 and R2
  blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 5 -perc_identity 98 -qcov_hsp_perc 90 -db /mnt/Genomics/Working/databases/holi_bowtie2_indexedDBs/refseq/bacteria_blastDB/ref_prok_rep_genomes -outfmt '6 qseqid stitle sacc pident qcovs evalue bitscore' -query bacteria-removed_fastqs/$bnameR1.fasta -out bacteria-removed_fastqs/$bnameR1.blast_out
  blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 5 -perc_identity 98 -qcov_hsp_perc 90 -db /mnt/Genomics/Working/databases/holi_bowtie2_indexedDBs/refseq/bacteria_blastDB/ref_prok_rep_genomes -outfmt '6 qseqid stitle sacc pident qcovs evalue bitscore' -query bacteria-removed_fastqs/$bnameR2.fasta -out bacteria-removed_fastqs/$bnameR2.blast_out
  
  echo Preprocessing. Compiling list of reads to remove.
  # make list of the names of those reads
  awk -F"\t" '{print $1}' bacteria-removed_fastqs/$bnameR1.blast_out | sort | uniq > bacteria-removed_fastqs/$bnameR1.bacterial_reads
  awk -F"\t" '{print $1}' bacteria-removed_fastqs/$bnameR2.blast_out | sort | uniq > bacteria-removed_fastqs/$bnameR2.bacterial_reads
  
  # concatenate the R1 and R2 lists, AdapterRemoval won't work if the read is missing from only one of the fastqs
  cat bacteria-removed_fastqs/$bnameR1.bacterial_reads bacteria-removed_fastqs/$bnameR2.bacterial_reads | sort | uniq > bacteria-removed_fastqs/$bname.bacterial_reads

  echo "Preprocessing. removing reads with >=98% identity with refseq bacterial genomes"
  #filter fastq based on identified bacterial reads above
  /data2/programs/bbmap/filterbyname.sh in=$infileR1 out=bacteria-removed_fastqs/$bnameR1.fastq names=bacteria-removed_fastqs/$bname.bacterial_reads overwrite=true
  /data2/programs/bbmap/filterbyname.sh in=$infileR2 out=bacteria-removed_fastqs/$bnameR2.fastq names=bacteria-removed_fastqs/$bname.bacterial_reads overwrite=true
  
  # clean up unnecessary files
  # rm bacteria-removed_fastqs/*.fasta # optionally run this if you don't want to hold on to any of the bacterial reads data 
  gzip bacteria-removed_fastqs/*.blast_out &
  gzip bacteria-removed_fastqs/*.bacterial_reads &
done
```

Now your files should have removed any bacteria that aligned well to what was available in our copy of the RefSeq bacterial database. There could still be bacteria that aren't included in that database, or perhaps ancient bacterial reads with lots of damage that now have < 98 % identity to references in the database, but this should give you much cleaner taxonomies. 
