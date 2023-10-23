# Basic Introduction to Bioinformatics for Undergraduate Teaching

__Main author:__  Trevor T. Bringloe  
__Affiliation:__  University of New Brunswick (UNB)   
__Group:__        Biology Department   
__Location:__     New Brunswick, Canada  
__Affiliated publication:__  
__Contact:__      e-mail: tbringloe@gmail.com | tel: (506)-259-2288


- [Objective](#objective)
- [Project Summary](#project-summary)
- [Disclaimer](#Disclaimer)
- [Tutorial 1 Linux based commands](#tutorial-1-linux-based-commands)
- [Tutorial 2 High Performance Computing](#tutorial-2-high-performance-computing)
- [Lab 1 Assembling LEGO K-mers](#lab-1-assembling-lego-k-mers)
- [Lab 2 Assembling and annotating organellar genomes](#lab-2-assembling-and-annotating-organellar-genomes)
- [Tutorial 3 Distilling Norwegian algal turf read datasets](#tutorial-3-distilling-norwegian-algal-turf-read-datasets)
- [Lab 3 Assignment](#lab-3-assignment)
- [Final Word](#final-word)
- [Acknowledgements](#acknowledgements)
- [References](#references)


## Objective
The repository provides a basic introduction to bioinformatics intended for upper level undergraduate students. The content introduces concepts related to read quality control, assembly of short reads, gene annotations, read mapping, and similarity searches of DNA barcodes against NCBI ([BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_PROGRAMS=megaBlast&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on&LINK_LOC=blasthome) searches) and the Barcode of Life Data Systems ([BOLD](http://v4.boldsystems.org/)). The content required to complete the project and lab assignments detailed here are not provided in their entirety due to data storage limitations on github, or otherwise because materials are physical (i.e. LEGO assembly lab). Note also, steps to complete are not verbatim, but require some level of user input to modify directory pathways, file names, ect. The content can be modified according to the level of immmersion expected of students; as drafted here, computational steps are provided but not performed by students. Given the opportunity (i.e. a full semester course in bioinformatics) students should be expected to (eventually) complete steps themselves to ensure engagement on a deeper level. The ultimate objective of a bioinformatics course should be to impart on students the capacity to develop a workflow de novo to achieve a specified task(s).

## Project Summary
Macroalgae (seaweeds) represent a conglomerate of species, particularly turf forming species that grow entangled with one another close to rocky coastal substrate. Moreover, the seaweeds themselves are home to a host of bacterial and eukaryotic taxa, many of which are unknown to science. Seaweed holobiomes are not simply a grouping of species, they represent an ecosystem of relationships, some obligate, some passive, some mutualistic, and some parasitic. The complexity of seaweed holobiomes (and holobiomes in general) in terms of the number of species and their relation to one another simply outscales human comprehension. Nonetheless, we can use big-data approaches to gleam some insight into these complex biological worlds.

Bioinformatics offers one such glimpse into the structure of holobiomes. Between August 7-12, 2022, I was fortunate enough to join colleagues on a research cruise to collect seaweeds in Tjongspollen Fjord, south of Bergen, Norway, as part of a taxonomic survey funded by the Norwegian Taxonomic Inititiative. Prior to the cruise, I also dove several sites in the area of Bergen with colleagues from France. We were able to collect and document algal turfs, including photo evidence of the turf samples as seen from a dissecting microscope. These samples were subsequently preserve in silica, had DNA extracted at the University of New Brunswick, Canada, and were sequenced at Genome Quebec using the NovaSeq6000 platform, targeting 100 million 150 bp paired-reads reads per sample in 25 samples (though closer to 3 billion reads were generated in total). Our objective was to infer the presence of inconspicous species using molecular data that were otherwise escaping labour intensive sorting under the microscope. Because we were interested in species level inferences, we used a whole genome approach rather than traditional metabarcoding (see introduction to bioinformatics lecture).

Using the tutorials provided here, students will understand how to distil these whole genome sequencing datasets into information that would allow inferences regarding species present in algal turf samples. Students should consult the files provided here, including sample metadata, and follow links to relevant sites for more information. The computationally intensive steps have been completed a priori, but labs nonetheless guide students through the concepts underpinning the bioinformatics methods used here, which were carried out in a high performance computing environment. A project report detailing species found within a particular algal turf sample is expected, including introduction, methods, results, and discussion.

Students can access specimen photos via a google drive link provided in class.

## Disclaimer
The material presented here is intended as a basic introduction to bioinformatics. The topics covered are by no means comprehensive, and potentially more efficient coding options are certainly available. The intention is to provide students with exposure to bioinformatics. The information is accurate to the best of my knowledge, and targetted toward biological students with a desire to incorporate computational methods into their learning. The material is not intended to foster formal learning in coding, scripting, understandning computional theory, ect. It is simply meant to get students introduced to bioinformatics and gain some momentum applying principals of read quality assessment, assembly, mapping, and extracting biological information from concensus sequences. The material rather serves as an initial framework for expanding into a fully fledged course.

## Tutorial 1 Linux based commands

Most students are used to operating systems that are highly visual and interactive (i.e. the user clicks icons to issue commands). Linux based computing, however, allows the user to execute desired actions using the command line interface. This provides the user with much greater flexibility for manipulating files in an efficient manner, while still enabling more complex functions to be executed using various built-in and external scripting languages. This tutorial provides basic commands used in linux. Exercises are provided to test commands and reinforce comprehension, which evidently require students to have access to a linux based environment (e.g. a virtual box or software such as Ubuntu). This tutorial is basically meant to supplement learning and allow students to better understand the commands provided in labs. It is by no means comprehensive; the intention is to put these tools in the hands of forthcoming biologists. Most people learn how to drive without ever knowing exactly how their vehicle works. Students should seek out a computer science course or linux-based workshop to further learning in this area.

A couple additional things to note. You cannot visualize or interact with linux file systems the same way we click and see directory heirarchies in windows or apple OS. Rather, the user must create and familiarize themselves with directory heirarchies, and use pathways (e.g. /home/trevor/data) to move about and interact with their environment and the files therein. Some third party programs can help, however, with connecting and navigating this environment, specifically if the environment is remote and you must login to it (e.g. [Filezilla](https://filezilla-project.org/). A good SSH client for remote work is [PuTTy](https://www.putty.org/). 
A hashtag indicates to not execute proceeding text as commands. They are therefore a useful tool for annotating scripts.
You can use the TAB key to finish a filename or pathway. If multiple options are present, they will be displayed. This can save on typing.
The up arrow key can be used to recall previously used commands, in order of most recent.
Many commands will display output at the command line interface. End any command with '> file.name' to output command to a the specified file (i.e. something more informative than file.name).
You can set environment variables, which sub in text in specified locations (see examples below). This can cut down on text, and facilitate more advanced commands such as looping through sample IDs
Watch out for double and single quotations, they must be closed, otherwise your command/script will break. Single quotes tell the sytem to interpret text as written.
You can pipe the output of one command into a new command using |. This helps avoid bulky commands/scripts and writing intermediate files.
You can paste copied text into a command window by right clicking

*list of common linux-based commands*

```
pwd # list present working directory, useful to get oriented in your environment
pwd
/home/tbringloe/Norway_turfs

ls # list files in present directory, or in specified directory. Can be used in combination with file names and wildcards.
ls *_R1.fastq.gz
TTB000600_R1.fastq.gz
TTB000601_R1.fastq.gz

cd # change directory, use in combination with downstream file pathway, or use full pathway to exit and enter new pathways
cd /home/tbringloe/Norway_turfs/
cd raw_reads # changes to file raw_reads folder in present working directory

.. # move up one directory

nano # this will open and display contents of any specified file (providing it is text). This is useful for editing files. Be sure to save edits when exiting with CTRL+X.

cat # list, combine or write file contents, useful for combining files or listing contents to be used in another command

cp # copy specified file
cp TTB000600_R1.fastq.gz raw_reads # copies fastq file to raw_read folder in present working directory

mkdir # make a new specified directory
mkdir trimmed_reads # creates new directory trimmed_reads in present working directory

mv # move specified file to new location or rename file
mv raw_reads/TTB000600_R1.fastq.gz trimmed_reads # moves file TTB000600_R1.fastq.gz in raw_reads directory to trimmed_reads folder
mv trimmed_reads/TTB000600_R1.fastq.gz trimmed_reads/TTB000600_trimmed_R1.fastq.gz # renames file moved to trimmed reads above

rm trimmed_reads/TTB000600_trimmed_R1.fastq.gz # removes file TTB000600_trimmed_R1.fastq.gz from trimmed_reads directory, can also use rm -R to delete files recursively
rmdir trimmed_reads # removes directory trimmed_reads, along with its contents

grep # prints lines matching criteria in a specified file. Can also use -v flag to pring inverse (any lines not matching criteria)
grep TTB sample.list
TTB000600
TTB000601
TTB000602

df -H # displays disk usage in human readable format

head -n 10 # displays first ten lines in a specified file
tail -n 10 # displays last ten lines in a specified file

top # displays running processes on the system, along with resource usage in terms of cores and RAM
htop # similar to top, but with different display
ps -x --forest # will display a tree of processes and subprocesses

echo # will display string to standard output, useful to output text in scripts or try commands without running them, particularly if variables are set

bash script.sh # execute lines of code in script.sh

# Some more advanced commands
nohup # when working on a remote server, use nohup to run command so that it continues to run if you exit the session or connection drops
screen # an alternative to nohup, which allows users to run tasks in multiple windows. Will also prevent command from hanging if connection drops of session is exited.

CTRL+Z # will suspend the currently running process
bg # run suspended process in the background
fg # bring background process into the foreground (i.e. commmand line)
```

*Some example command strings that are useful in bioinformatics*

```
# Create a list of sample IDs
cd /home/tbringloe/raw_reads
ls *_R1.fastq.gz > sample.IDs
# since we only want the prefix for downstream commmands, remove _R1.fastq.gz using sed
sed -i 's/_R1.fastq.gz//g' sample.IDs

# Loop command through set of sample IDs stored in text file sample.list. Cat reads sample.IDs file, then for each line execute the command
cat sample.IDs | while read line
do
fastqc "$line"_R1.fastq.gz "$line"_R2.fastq.gz -o .
done

# Can also use xargs command to perform command in parallel using multiple system threads.
cat sample.list | xargs -I {} -n 1 -P $threads sh -c "fastqc {}_R1.fastq.gz {}_R2.fastq.gz -o ."

# Set variables so you can quickly swap out relevant information througout lines of code
marker=coxI
echo ls /home/tbringloe/assembled_reads/"$marker"
ls /home/tbringloe/assembled_reads/coxI

# create forks in your workflow based on specified criteria. You can create elaborate bash scripts this way
trim_raw_reads=1 # 1 to trim raw reads, 0 to skip step
echo running workflow
if [$trim_raw_reads -eq 1]
then
  echo run trimmmomatic
  echo <insert trimmomatic code here>
else
  echo User has specified to skip raw read trimming
fi
exit 0
```

## Tutorial 2 High Performance Computing
Bioinformatics involves working with large datasets containing many GBs or even TBs of sequence information. Because there is so much information to work with, we cannot use normal laptops or desktops. We need computers with a lot of storage (i.e. disk) space for storing and writing files, and a lot of RAM (Random Access Memory) to work with lots of data at a given moment. Furthermore, resources can be distributed across nodes and cores in a system, allowing users to run many commands at the same time, or divide tasks into many smaller tasks, thus allowing users to complete commands in a shorter time frame. These needs are met by High Performance Computers (i.e. High Performance Computing; HPC). Canadian Universities share access to Compute Canada servers, but many institutions host private servers. One of the challenges on a shared system is tasks must be executed in an efficient manner; the system cannot be overloaded with tasks (otherwise this would quickly crash the system) and computational resources must be used in the most efficient manner possible to ensure the most amount of computation gets completed in the shortest amount of linear time. To facilitate this, shared servers use slurm scripts, which users use to submit their tasks to perform. Tasks are then queued and run when resources become available.

Here is an example slurm script one can submit to Compute Canada servers. The slurm script may look slightly different depending on how the shared system has been orchestrated.
```
#!/bin/bash

#SBATCH --account=def-xxxxxxx

#SBATCH --time=0-03:59:59

#SBATCH --job-name=fastqc_raw_check

#SBATCH -n 4

#SBATCH --mem=10G

#SBATCH -N 1

#SBATCH --mail-type=END

#SBATCH --mail-user=tbringloe@gmail.com

#module to load
module load StdEnv/2020
module load fastqc/0.11.9

#job commands
threads=20
cat sample.list | xargs -I {} -n 1 -P $threads sh -c "fastqc {}_R1.fastq.gz {}_R2.fastq.gz -o ."
```
In this example, the user has specified the account under which to run tasks (computation resources are carefully allocated and monitored across accounts to ensure resources are utilized in an equitable manner). The maximum amount of time to run the task(s) (i.e. wall time) has been specified at just under 4 hours. The job name has been specified as 'fastqc_raw_check', which will appear when monitoring the job status. -n specifies to use 4 threads or tasks when running the command, and --mem specifies the amount of extra memory or RAM to allocate to the job (usually there is a maximum amount of memory provided per node, so users will need to request additional memory if needed). -N indicates to run on a single node (most actions can be performed on a single node, except in more advanced tasks that require running seperate tasks on multiple nodes). --mail-user specifies an email where notifications can be sent when the job enters a new status (e.g. --mail-type=BEGIN, END, FAIL). Shared environments typically have staff who can install software of interest; as such, programs needed to run commands are loaded as modules, and become available for use once loaded. Oftentimes, system dependencies must also be loaded to run particular programs. Job commands are the commands you want to run. When working in a HPC environment, some specific and additional commands are available for use.

**NEVER RUN TASKS DIRECTLY ON A SHARED HIGH PERFORMANCE COMPUTER. DEPENDING ON THE SIZE OF THE TAKS, YOU RUN THE RISK OF CRASHING THE ENTIRE SYSTEM. THIS WILL RESULT IN POTENTIAL PENALTIES FROM STAFF, AND A LOT OF GRUMPY USERS WHOSE JOBS YOU KILLED IN ONE FELL SWOOP. ALWAYS USE SALLOC IF AN INTERACTIVE SESSION IS REQUIRED**

*list of common HPC commands*

```
# These commands work on compute canada servers, but may differ on different HPC systems

sbatch # submit slurm script to the system queue
sbatch yourscript.slurm

squeue -u tbringlo # check status of jobs submitted by specified user

scancel <jobID> # cancel a job that is queued or running

module spider <program.of.interest> # lists modules matching specified text. Useful to determine availability, including versions and dependencies

module load <program.of.interest> # load specified module

salloc --time 60 --mem 20 # launch an interactive session with a time limit of 60 minutes and with 20 GB of RAM. This is particularly useful for troubleshooting commands without repeatedly submitting jobs that fail, or smaller tasks that require more user input

exit # leave interactive or global session

diskusage_report # shows current storage and storage limits, both in terms of disk space and number of files. Projects are allocated a set amount of disk space.

sshare -l -A <user_account>_cpu -u tbringlo # check account and user usage, including levelFS, which proxies queue priority; <1 and priority is low because the account has used more resources than allocated, >1 and priority is higher because the account has used less resources than allocated
```

When working in an HPC environment, users must be mindful of the resources they are requesting and using. Large jobs take longer to process in the queue. The user is therefore incentivized to request only the required amount of resources to complete a task; estimating the appropriate amount of resources needed is gained with experience. When an account utilizes a lot of resources, users under that account receive less priority in the queue. Users are therefore also incentivized to be efficient with their resources, or risk long wait times to complete tasks. Use the above command sshare to monitor resource usage and queue priority.

When considering usage, users must consider two factors, the number of tasks to run and the job wall time. For instance, a single tasks running on a single thread for 24 hours=24 hours of computation time. But 32 tasks on 32 threads (or one task divided among 32 threads) running for 24 hours=768 (32 x 24) computation hours. Regardless of whether your tasks uses the available threads requested, the user account is "charged" for the time resources are tied up performing a given job. It is therefore critical to ensure the resources requested are being used efficiently.


## Lab 1 Assembling LEGO k-mers
Complete the following lab after the lecture introducing concepts related to bioinformatics, and the tutorials above. Note, if completing the computional steps in lab, the following commands require some user input to run, i.e. defining particular variables and establishing sample lists.

*Exercise 1: Evaluating read quality*

Before working with read files, it is critical to evaluate quality and remove error prone data. You have been provided [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) reports generated for several read datasets, before and after read trimming. They are under the GithHub folder [fastqc](https://github.com/tbringloe/Bioinformatics-CrashCourse/tree/main/fastqc). The following command was performed using [TRIMMOMATIC](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf):
```
# variables to declare for trimmomatic

trailing=20
headcrop=15
avgqual=25 
minlen=75
adapter=ill_adap_uni.fa

# Run trimmomatic, looping command over sample IDs contained in sample.list (text file)
cat sample.list | while read line
do
$java -jar $trimmomatic PE -threads 32 raw_reads/"$line"_R1.fastq.gz raw_reads/"$line"_R2.fastq.gz trimmed_reads/"$line"_R1_trimmed.fastq.gz trimmed_reads/"$line"_R1_tup.fastq.gz trimmed_reads/"$line"_R2_trimmed.fastq.gz trimmed_reads/"$line"_R2_tup.fastq.gz ILLUMINACLIP:raw_reads/$adapter:2:30:10 TRAILING:$trailing HEADCROP:$headcrop AVGQUAL:$avgqual MINLEN:$minlen
done
```

```
# Run fastqc, using xargs command to run fastqc over all samples simultaneously using multiple threads (i.e. in parallel to reduce linear computation time)
cat sample.list | xargs -I {} -n 1 -P $threads sh -c "fastqc {}_R1.fastq.gz {}_R2.fastq.gz -o ."
# Combine reports into a single html using multiQC
$multiqc . -o .
```

Investigate the FASTQC reports (individually or using multiQC) and answer the following questions.
1.	Briefly define the read trimming variables specified above by consulting the [TRIMMOMATIC](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) manual.
trailing=
headcrop=
avgqual= 
minlen=
adapter=

3.	Fill the following table based on the fastqc reports:

| File name | # of read before trimming | # of reads after trimming | Average read quality before trimming | Average read quality after trimming | GC content before trimmming | GC content after trimming |
| --- | --- | --- | --- | --- | --- | --- |
| TTB000601 |
| TTB000606 |
| TTB000611 |

3.	List notable differences between the before and after trimmming reports pertaining to the quality figures.

4.	Visit this [bad quality report](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html). What are potential issues with the read data and what trimming could be done to improve the quality of the specified read dataset?

*Exercise 2: Build an assembly graph* 

You have been provided long LEGO pieces with sequences written on them. These represent k-mers. The pieces have been attached in the true sequence using k-mers of 4 base pairs. In your group, find recurring (i.e. same instances of) k-mers. Where this occurs, break the LEGO sequence and attach string to connect the recurring k-mers.

1. Draw a diagram of your assembly graph (an example can be provided in lab, or see [BANDAGE](https://rrwick.github.io/Bandage/) for conceptual examples, though ours will be much simpler than what is shown here). You can also optionally take a picture and provide it with your report. Provide a figure title, indicating the k-mer size used.

Next, repeat the exercise, but flip all the LEGO pieces to the side where the k-mer sequences are 8 base pairs. Start overlapping and connecting the pieces at your k-mer size-1. Continue until all pieces are connected.

2. As before, draw a diagram of your assembly graph. Provide a figure title, indicating the k-mer size used.

3. Here we have assumed unlimited coverage for the assemblies. What was the effect of k-mer size on the assembly graphs? How did repetitive k-mers impact the graph?

4. Had coverage been limited, how might the assembly at the larger k-mer size been impacted?

*Exercise 3: Blast assembled sequences*

You have now an assembled sequence. Compare the sequence to NCBI's database using [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_PROGRAMS=megaBlast&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on&LINK_LOC=blasthome) and to [BOLD](http://v4.boldsystems.org/) records under the "identification" tab to assign taxonomic information to the sequences.

1.	What are the top hits for each sequence, and to what percentage do the sequences match existing records? Using the information in BOLD, can you describe the global distribution of the identified species?

## Lab 2 Assembling and annotating organellar genomes
Now that our reads are trimmed, and we are confident our data is of high quality, we can begin assembling reads. One option is to assemble the entire dataset. This would produce millions of contigs. It would also be computationally demanding (potentially impossible depending on the complexity of the dataset and the resources required to assemble). Rather than assemble the entire dataset, we will assemble the organellar genomes. This can be achieved using NOVOPlasty, a clever program that leverages a seed sequence and coverage information to infer circular organellar genomes. The following tutorial will produce mitochondrial and chloroplast sequences in two samples. The remainder of the lab will be dedicated to understanding annotations.

Here are the original sample metadata:

|Sample ID | Species | Sampling Location | Lat | Long | Date Sampled | Collector | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| TTB000538 | Coccotylus brodiei | Hakonsund | 60.17599 | 5.11074 | Aug-5-2022 | T.T. Bringloe| NA |
| TTB000539 | Meredithia microphylla | Hakonsund | 60.17599 | 5.110174 | Aug-5-2022 | T.T. Bringloe| NA |

[NOVOPlasty](https://github.com/ndierckx/NOVOPlasty) is a perl script that elongates a specified seed file. The seed file represents a sequence, either in the target organism or a closely related organism. For instance, DNA barcodes coxI and rbcL are good seed candidates for assembling algal organellar genomes (because that reference data is abundant). The script will then map reads to the seed sequence, determine high coverage k-mers, and start elongating the assembly graph from these k-mers while incorporating coverage information to hopefully avoid any breaks in the assembly. This works because organellar genomes occur as many copies within cells, meaning they are disproportionately represented in read datasets. Other assemblers such as [SPAdes](https://github.com/ablab/spades) also incorporates coverage information into the assembly graph. Once the program elongates enough such that it begins overlapping the other side of the same sequence, it confirms the sequence is circular and outputs the organellar genome as a fasta file.

```
# clone github files into the working directory
git clone https://github.com/ndierckx/NOVOPlasty.git
cd NOVOPlasty

# Update the config file, ideally saving to a seperate file; annotations have been provided for relevant fields, but the full config file has descriptors at the end:
Project:
-----------------------
Project name          = SAMPLEID_TAXA_mito # give your project a meaningful name
Type                  = mito # mito of chloro, depending on whether you are fetching a mitochondrial or chloroplast genome
Genome Range          = 20000-50000 # approximate expected range in lengths (base pairs)
K-mer                 = 57 # k-mer at which to assemble
Max memory            =
Extended log          = 0
Save assembled reads  = no # optional flag to save mapped reads
Seed Input            = TAXA_coxI.fasta # the fasta file with the seed sequence
Extend seed directly  = no # you can incorporate the seed sequence directly if you are confident it represents a true sequence from the target organellar genome, but it is safer to simply leave as no. NOVOPlasty will then incorporate any differences in the seed sequence.
Reference sequence    =
Variance detection    =
Chloroplast sequence  =

Dataset 1:
-----------------------
Read Length           = 135 # length of the reads, note, we performed a head crop of 15 bp, making read lengths 150-15=135
Insert size           = 300 # approximate insert size of the libraries, i.e. the approximate distance between paired reads
Platform              = illumina # platform used for sequencing
Single/Paired         = PE # PE=paired-end reads
Combined reads        = 
Forward reads         = SAMPLEID_R1_trimmed.fastq # forward reads dataset, must be uncompressed
Reverse reads         = SAMPLEID_R2_trimmed.fastq # reverse reads dataset, must be uncompressed
Store Hash            =

# Be sure to add the seed files to the NOVOPlasty directory and uncompress the read files, also in the NOVOPlasty directory
# Example command to uncompress .gz file : gzip -d TTB000632_R1_trimmed.fastq.gz
# Run NOVOPlasty, either via the command line interface if on a private computer, or via a slurm script as detailed above when working with HPC
perl NOVOPlasty4.3.3.pl -c TTB000539_Meredithia_mito_config.txt
```
It is possible to run NOVOPlasty over many samples simultaneously, we simply need more advanced code to swap in relevant sample details and create new output directories. See [this](https://github.com/tbringloe/WGS-NOVAC) tutorial and bash script for a potential solution. Note, NOVOPlasty won't necessarily produce a confirmed circular sequence. Depending on the complexity of the dataset, coverage, seed, k-mer size, ect. the output may occur as a fragmented or partial genome. Users will have to troubleshoot.

Now that we have a circular sequence, we can annotate the genomes. There are several approaches to annotating genomes. Because we are looking at relatively simple organellar genomes, we can use [MFannot](https://megasun.bch.umontreal.ca/RNAweasel/), an annotation database maintained at the Univeristy of Montreal. Visit the [MFannot](https://megasun.bch.umontreal.ca/apps/mfannot/) page and submit the mitochondrial and plastid sequences, selecting (1) standard and (11) bacterial, Archael, and plant plast genetic codes, respectively. Once the submission completes running, download the new annotation files.

The annotations are not usable in their current format. In order to view annotations in Geneious, we need to conver the MFannot to gff3 file format. Run a simple perl script to convert the file to one we can use.

```
# The third field is the MFannot output, the fourth is the sequence name in geneious, the fifth is the file name to write output to
perl MFannotSQL2GGF3.pl mfannot_75ed117a7fa2.fasta.new.sqn TTB000605_Coccotylus_truncatus_plastid > TTB000605_Coccotylus_plastid.gff3
```

*Exercise 1 observing annotations in geneious*

Make sure the organellar genomes are now in Geneious ([files here](https://github.com/tbringloe/Bioinformatics-CrashCourse/tree/main/Organellar_genomes)), and named accordingly. Drag the gff3 file into geneious and the annotations should automatically apply to the appropriate genome. Have a look at the annotations and answer the following questions.

1. How many coding sequences (labelled in yellow) were annotated for each genome? RNA sequences (labelled red)? Transfer RNA (leballed in pink)? Are the numbers the same or different between the two species?

2. Find and [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_PROGRAMS=megaBlast&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on&LINK_LOC=blasthome) the coxI (mitochondrial gene) and rbcL (chloroplast gene) sequences for both species. Were the original species identifications correct? Elaborate.

*Exercise 2 align genomes*

Now download the [MAUVE](https://darlinglab.org/mauve/mauve.html) alignment plugin for geneious. This alignment algorithm detects homologous "blocks" between two genomes. This is handy because organellar genomes are not always [co-linear](https://academic.oup.com/gbe/article/13/7/evab124/6290714); genes and intergenic regions can be present/absent between genomes, and some regions may appear in different orders. One cannot expect an alignment between genomes then, without the option of "chopping" the alignment up. Select both mitochondrial genomes and perform the MAUVE alignment option. Do the same for the chloroplast sequences.

3. Are the mitochondrial and chloroplast genomes co-linear, or do genomic regions appear out of order across genomes?

*Bonus exercise Group II introns*

4. Have a closer look at the coxI sequence of Meredithia. Investigate the group II intron by [blasting](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_PROGRAMS=megaBlast&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on&LINK_LOC=blasthome) the sequence. Try also a [blastx](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Tdranslations&PROGRAM=blastx&PAGE_TYPE=BlastSearch&BLAST_SPEC=) search (searching the translated sequence, i.e. amino acid level). What is the closest match? Blast the flanking coxI sequence and determine its closest match. Can you explain the origin of the group II intron and provide hypotheses as to why it might appear specifically in the coxI sequence?

5. What purpose might it serve to have fully resolved organellar genomes as opposed to DNA barcoding genes?

## Tutorial 3 Distilling Norwegian algal turf read datasets

For this lab, students will review the various steps to distill large raw read files to curated BLAST reports for DNA barcode regions, which can be used to interpret species present within the algal turf samples. The tutorial will cover steps taken to process the following three samples: TTB000601, TTB000606, TTB000611. Specimen photos can be found here: https://drive.google.com/drive/folders/1EX7cczCrXq7ODOZHp4ddxP6jAzai07zQ?usp=sharing

|Sample ID | Species | Sampling Location | Lat | Long | Date Sampled | Collector | BioSample Accession | Mitogenome resolved? | Chloroplast resolved? | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TTB000535 | Odonthalia dentata | Kleppesjoen | 60.1847 | 5.14938 | Aug-5-2022 | T.T. Bringloe, L. LeGall, J. Utge | SRR25569791 | Yes | Yes | NA |
| TTB000538 | Phyllophora crispa | Hakonsund | 60.17599 | 5.110174 | Aug-5-2022 | T.T. Bringloe, L. LeGall, J. Utge | SRR25569790 | Yes | Yes | NA |
| TTB000539 | Mededithia microphylla | Hakonsund | 60.17599 | 5.110174 | Aug-5-2022 | T.T. Bringloe, L. LeGall, J. Utge | SRR25569779 | Yes | Yes | NA |
| TTB000600 | Cladostephus spongiosus | Tjongspollen, Station 2 | 59.67424 | 5.233623 | Aug-9-2022 | T.T. Bringloe | SRR25569773 | No | Yes | NA |
| TTB000601 | Gelidium spinosum | Tjongspollen, Station 2 | 59.67424 | 5.233623 | Aug-9-2022 | T.T. Bringlo e| SRR25569772 | No | Yes | NA |
| TTB000605 | Coccotylus brodiei | Tjongspollen, Station 2 | 59.67424 | 5.233623 | Aug-9-2022 | T.T. Bringloe | SRR25569771 | Yes | Yes | NA |
| TTB000606 | Ascophyllum nodosum | Tjongspollen, Station 2 | 59.67424 | 5.233623 | Aug-9-2022 | T.T. Bringloe | SRR25569770 | Yes | Yes | NA |
| TTB000609 | Gelidium spinosum | Tjongspollen, Station 2 | 59.67424 | 5.233623 | Aug-9-2022 | T.T. Bringloe | SRR25569769 | Yes | No | NA |
| TTB000610 | Chorda filum | Tjongspollen, Station 2 | 59.67424 | 5.233623 | Aug-9-2022 | T.T. Bringloe | SRR25569768 | No | No | NA |
| TTB000611 | Laminaria hyperborea | Tjongspollen, Station 4 | 59.69405 | 5.246691 | Aug-10-2022 | L. LeGall, J. Utge, Franz Goecke | SRR25569767 | Yes | Yes | Stipe scrapes |
| TTB000612 | Laminaria hyperborea | Tjongspollen, Station 4 | 59.69405 | 5.246691 | Aug-10-2022 | L. LeGall, J. Utge, Franz Goecke | SRR25569789 | Yes | Yes | Mid-blade |
| TTB000614 | Palmaria palmata | Tjongspollen, Station 4 | 59.69405 | 5.246691 | Aug-10-2022 | L. LeGall, J. Utge, Franz Goecke | SRR25569788 | Yes | No | NA |
| TTB000615 | Codium vermilara | Tjongspollen, Station 3 | 59.66574 | 5.225333 | Aug-10-2022 | L. LeGall, J. Utge, Franz Goecke | SRR25569787 | No | No | NA |
| TTB000617 | Chondrus crispus | Tjongspollen, Station 3 | 59.66574 | 5.225333 | Aug-10-2022 | T.T. Bringloe | SRR25569786 | Yes | Yes | NA |
| TTB000618 | Cladostephus spongiosus | Tjongspollen, Station 3 | 59.66574 | 5.225333 | Aug-10-2022 | T.T. Bringloe | SRR25569785 | No | No |
| TTB000620 | Phycodrys sp.1NB | Tjongspollen, Station 3 | 59.66574 | 5.225333 | Aug-10-2022 | L. LeGall, J. Utge, Franz Goecke | SRR25569784 | No | Yes | NA |
| TTB000621 | Phyllophora crispa | Tjongspollen, Station 3 | 59.66574 | 5.225333 | Aug-10-2022 | L. LeGall, J. Utge, Franz Goecke | SRR25569783 | Yes | Yes | NA |
| TTB000629 | Alaria esculenta | Dredge 2 | 59.36036 | 5.09016 | Aug-11-2022 | T.T. Bringloe | SRR25569782 | Yes | Yes | NA |
| TTB000631 | Membranoptera alata | Dredge 2 | 59.36036 | 5.09016 | Aug-11-2022 | T.T. Bringloe | SRR25569781 | No | Yes | NA |
| TTB000632 | Phycodrys rubens | Dredge 2 | 59.36036 | 5.09016 | Aug-11-2022 | T.T. Bringloe | SRR25569780 | No | Yes | NA |
| TTB000634 | Ptilota gunneri | Dredge 2 | 59.36036 | 5.09016 | Aug-11-2022 | T.T. Bringloe | SRR25569778 | No | Yes | NA |
| TTB000636 | Saccorhiza polysiches| Dredge 2 | 59.36036 | 5.09016 | Aug-11-2022 | T.T. Bringloe | SRR25569777 | Yes | Yes | NA |
| TTB000637 | Phycodrys rubens | Dredge 2 | 59.36036 | 5.09016 | Aug-11-2022 | T.T. Bringloe | SRR25569776 | No | Yes | NA |
| TTB000639 | Rhodomela confervoides | Dredge 2 | 59.36036 | 5.09016 | Aug-11-2022 | T.T. Bringloe | SRR25569775 | Yes | Yes | NA |
| TTB000641 | Coccotylus brodiei | Hakonsund | 60.17599 | 5.110174 | Aug-5-2022 | T.T. Bringloe, L. LeGall, J. Utge | SRR25569774 | Yes | Yes | NA |

Our goal is to identify species in the algal turfs. For this, we need taxonomically informative sequence data. This is present in the read datasets, but must be extracted. Specifically, we can use well-established DNA barcode sequences to identify species (see lecture material). One strategy is to assemble the entire read datasets and "fish" for sequences corresponding to DNA barcodes. As noted above, however, this is computationally intensive. Another approach is to "fish" for DNA barcode sequences at the read level, then assemble reads corresponding to DNA barcodes to retrieve taxonomically informative sequences. The first step in this process is to create a reference database of DNA barcode sequences. In order to reduce computational time during read mapping, using [mmseqs2](https://github.com/soedinglab/MMseqs2) (Steinegger & Söding, 2017) we will cluster sequences at 90% similarity, and keep a single representative sequence from each cluster (reducing our database from millions of sequences to ~200,000). We will use this to map reads and identify candidate DNA barcode reads. The following code is meant to illustrate the steps taken to arrive at a reference database. Students are not expected to understand these steps, they are meant to supplement learning. Students will engage with data downstream, post mapping.

```
# We can retrieve DNA barcode sequences from BOLD. This is a better option, as the data is curated towards taxonomically informative regions, whereas NCBI would have full length DNA barcode genes that would create interpretation problems downstream when we start blasting sequences.
# taxon.list is a text file of Phyla I have chosen to include in the database. It represents all the phyla on BOLD, plus some lower level designations. This was done strategically to aid in curating barcode alignments downstream.
marker=COI-5P # Also ran declaring rbcL and tufA
cat taxon.list | while read line
do
wget http://v3.boldsystems.org/index.php/API_Public/sequence?taxon="$line"
# download will represent all sequence data for given phyla, so we need to use grep to retrieve lines with relevant marker
grep "$marker" -A 1 sequence\?taxon\="$line" > sequence\?taxon\="$line"_"$marker".fasta
# another reverse grep to remove dashed introduced with previous command
grep -vx -e "--" sequence\?taxon\="$line"_"$marker".fasta > sequence\?taxon\="$line"_"$marker"_final.fasta
done

# Use mmseqs2 to cluster at 90% similarity and retain a representative sequence for the database
module load mmseqs2/13-45111
cat taxon.list | while read line
do
mmseqs createdb sequence\?taxon\="$line"_"$marker"_final.fasta sequence\?taxon\="$line"_"$marker"_final.db
#Store temporary files in temporary directory
mkdir tmp_dir
mmseqs cluster sequence\?taxon\="$line"_"$marker"_final.db "$line"_db_clu tmp_dir --threads 4 --cov-mode 0 -c 0.90 --min-seq-id 0.90
mmseqs result2repseq sequence\?taxon\="$line"_"$marker"_final.db sequence\?taxon\="$line"_"$marker"_final.db_clu sequence\?taxon\="$line"_"$marker"_final.db_clu_rep
mmseqs result2flat sequence\?taxon\="$line"_"$marker"_final.db sequence\?taxon\="$line"_"$marker"_final.db sequence\?taxon\="$line"_"$marker"_final.db_clu_rep sequence\?taxon\="$line"_"$marker"_final.rep.fasta --use-fasta-header
done

# Sequences were aligned using clustal-omega alignment
module load clustal-omega/1.2.4
cat taxon.list | while read line
do
clustalo --threads=32 -i sequence\?taxon\="$line"_"$marker"_final.rep.fasta -o sequence\?taxon\="$line"_"$marker"_final.rep.aligned.fasta
done
```

The [clustal-omega](https://www.ebi.ac.uk/Tools/msa/clustalo/) alignments were curated in geneious to ensure the database consisted of only DNA barcode regions (some sequences hang off the end of these regions; we don't want reads mapping outside). Sequences with a lot of missing data were also deleted, and alignments were manually edited where possible.

Now that we have reasonably comprehensive files of available DNA barcode sequences (we started with 1,291,373 sequences downloaded from BOLD, down to 181,554 representative sequences), we can create our database and map reads using [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) (Langmead & Salzberg 2012). We will use local alignments (i.e. read sequences can be clipped) matching to within 25% similarity to our reference database.

```
# concatenate all the barcodes across phyla to create marker specific databases
marker=coxI
threads=32
cat *_"$marker"_final.rep.fasta > "$marker"_barcodes_final.fasta
# index the database so bowtie2 can retrieve sequences for mapping. Here we build in directory IDX (mkdir if needed)
bowtie2-build "$marker"_barcodes_final.fasta IDX/BOLD_"$marker"

# Map reads from sample datasets. Prefix is stored in file sample.list. --al-gz saves mapped reads into the mapped_reads directory and into a new fastq.gz file
cat sample.list | while read line
do
bowtie2 --local --score-min G,1,1 --al-conc-gz mapped_reads/"$marker"/"$line"_"$marker".fastq.gz --al-gz mapped_reads/"$marker"/"$line"_"$marker".fastq.gz -p $threads -x IDX/BOLD_"$marker" -1 trimmed_reads/"$line"_R1_trimmed.fastq.gz -2 trimmed_reads/"$line"_R2_trimmed.fastq.gz -S "$line"_"$marker".sam
# reads are also stored in an alignment file. We won't need this
rm "$line"_"$marker".sam
done

# repeat analysis for other markers of interest
```

Next, the saved mapped reads (which correspond to putative DNA barcode reads) are assembled using [SPAdes](https://github.com/ablab/spades) short read assembler (Bankevich et al. 2012). SPAdes can assemble at different k-mer sizes specified with the -k flag, and produce an assembly concensus across k-mers.

```
module load spades/3.15.4
marker=coxI
cat sample.list | while read line
do
spades.py -1 mapped_reads/"$marker"/"$line"_"$marker".R1.fastq.gz -2 mapped_reads/"$marker"/"$line"_"$marker".R2.fastq.gz -k 21,33,55,77,109 -t $threads -m 50 -o assembled_reads/"$marker"/"$line"
# move and rename assemblies
mv assembled_reads/"$marker"/"$line"/contigs.fasta barcode_results/"$marker"/"$line"_"$marker".fasta
done
```

Now we have assembled barcode sequences from our read datasets, among many other types of sequences (by catch). Now we need to assign taxonomic information and curate the dataset into something manageable and easy to interpret. We start this process by creating a local blast database of our barcode sequences, and blasting the assembled barcodes to this database. Sequences with positive hits, and the beginning and end position of alignments of positive hits, are then extracted from the initial assemblies. This way we can weed out assembled by-catch (reads that mapped to our database but don't represent true barcodes). The curated sequences are then [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_PROGRAMS=megaBlast&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on&LINK_LOC=blasthome) compared to the full NCBI nt database.

```
# NBCI build-db function does not like duplicate names in the input fasta, so this one liner cleared that up
awk '/^>/{f=!d[$1];d[$1]=1}f' in.fa > out.fa

# Make a custom nucleotide blast database
module load blast+/2.12.0
marker=coxI
makeblastdb -in "$marker"_barcodes_final.fasta -dbtype nucl -parse_seqids -out blastdb/BOLD_"$marker"

# Blast assembly output to custom database
cat sample.list | while read line
do
blastn -query barcode_results/"$marker"/"$line"_"$marker".fasta -db BLAST/BOLD_"$marker" -num_threads 32 -max_target_seqs 10 -outfmt "7 qseqid staxids bitscore std stitle evalue" -out blast_results/"$marker"/"$line"_"$marker"_blastout.txt
done

# This series of commmands will curate the assemblies so we are left with putative barcode sequences, and only barcode regions
module load seqtk/1.3
module load bedtools/2.30.0
mkdir sequence_extraction
cat sample.list | while read line
do
# Use grep to extract sequences with hits to the barcode database. We can reverse grep the hashtag symbol, since this represents sequences without a hit
grep -v -E "#" blast_results/"$marker"/"$line"_"$marker"_blastout.txt > sequence_extraction/"$line".hits.txt
# We want to exclude positive hits with an alignment length less than 250 (more artifacts at smaller alignment legnths). Start by sorting lines by sequence names followed by alignment length
sort -k 1,1 -k 7,7rn sequence_extraction/"$line".hits.txt > sequence_extraction/"$line".hits.sorted.txt
# Since the blast output contains 10 lines of positive hits, we remove redundancy by extracting the first instance of unique sequence names (i.e. extracts hit with longest alignment)
sort -k 1,1 -u sequence_extraction/"$line".hits.sorted.txt > sequence_extraction/"$line".hits.sorted.unique.txt
# Now we can remove sequences with alignment lengths less than 250
awk '$7 >= 250' sequence_extraction/"$line".hits.sorted.unique.txt > sequence_extraction/"$line".hits.sorted.unique.250.txt
# Extract sequence IDs for seqtk
awk '{print $1}' sequence_extraction/"$line".hits.sorted.unique.250.txt > sequence_extraction/"$line".hits.sorted.unique.250.IDs.txt
# Run seqtk to extract relevant sequences from assembly fasta file
seqtk subseq barcode_results/"$marker"/"$line"_"$marker".fasta sequence_extraction/"$line".hits.sorted.unique.250.IDs.txt > sequence_extraction/"$line".hits.sorted.unique.250.IDs.fasta
# Extract start and stop positions for relevant sequences
awk 'BEGIN { OFS="\t" }{print $1,$10,$11}' sequence_extraction/"$line".hits.sorted.unique.250.txt > sequence_extraction/"$line".hits.sorted.unique.250.regions.bed
# Extract barcode regions from fasta subset
bedtools getfasta -fi sequence_extraction/"$line".hits.sorted.unique.250.IDs.fasta -bed sequence_extraction/"$line".hits.sorted.unique.250.regions.bed -fo sequence_extraction/"$marker"/"$line".hits.sorted.unique.250.IDs."$marker".barcodes.fasta
# Clean up environment by removing intermediate files
rm sequence_extraction/"$line".*
done
```
Now we have curated files of DNA barcode regions for our samples. We can then BLAST these sequences against the full NCBI non-redundant nucleotide (nt) database. Users might consider BLASTX against the protein (nr) database to ID highly divergent sequences. Note, the nt database is typically downloaded and regularly updated on shared servers. So we can specify the directory path to the database (blocked out here)

```
marker=coxI
cat file.list | while read line
do
blastn -query sequence_extraction/"$marker"/"$line".hits.sorted.unique.250.IDs."$marker".barcodes.fasta -db <directory.pathway>/blast_dbs/2022_03_23/nt -num_threads 32 -max_target_seqs 10 -outfmt "7 qseqid staxids bitscore std stitle evalue" -out blast_results/"$marker"/"$line"_"$marker"_blastout_regions_nr_24viii23.txt
done
```

Now we have curated DNA barcode sequences from our samples, and BLAST output, which can be used to infer species present in the algal turfs :)

## Lab 3 Assignment

Because bioinformatic workflows can be complicated and non-intuitive, it can be helpful to produce a flow diagram of the steps taken to achieve results. [Mermaid live](https://mermaid.live/) is a solid resouce for generating such diagrams. The language is somewhat unique, but a basic flow diagram could look like this:

```
# Mermaid lines of code to generate a generic workflow diagram
graph TD
A((Input data 1 <br> file format)) --> B(Step 1 <br> programs used)
B --> C(Step 2 <br> programs used)
C --> D(Step 3 <br> programs used)
D --> DJ
DE((Input data 2 <br> file format)) --> DF(Step 4 <br> programs used)
DF --> DJ(Step 5 + Programs used)
DJ --> F(Step 6 <br> Programs used)
G((Input data 3 <br> file format)) --> H(Step 7 <br> Programs used)
H --> J
F --> J(Step 8 <br> Program used)
J --> JI(Final dataset)
style A fill:#f9f,stroke:#333,stroke-width:4px
style DE fill:#f9f,stroke:#333,stroke-width:4px
style G fill:#f9f,stroke:#333,stroke-width:4px
style JI fill:#bbf,stroke:#f66,stroke-width:2px,stroke-dasharray: 5 5
```
[Online FlowChart & Diagrams Editor - Mermaid Live Editor.pdf](https://github.com/tbringloe/Bioinformatics-CrashCourse/files/12441754/Online.FlowChart.Diagrams.Editor.-.Mermaid.Live.Editor.pdf)
![Online FlowChart   Diagrams Editor - Mermaid Live Editor](https://github.com/tbringloe/Bioinformatics-CrashCourse/assets/79611349/85222c03-1b62-4563-89f4-75e494c06e11)

For this assignment, plug in the above code into Mermaid live edit, and create a workflow diagram of the steps taken to arrive at the curated DNA barcode datasets of Norwegian algal turf samples. You do not necessarily need to follow the above code; feel free to research other coding options to make your diagram striking and unique.

## Final word

This is the end of the crash-course in bioinformatics. I hope the material was instructive. Feel free to email me at tbringloe@gmail.com if you have questions of concerns.

## Acknowledgements
All the inspirational students, postdocs, mentors, DFO research scientists, and forum junkies across the globe who contributed to my own learning journey in bioinformatics.

## References
**Other potentially useful tutorials:**

Variant calling for phylogenomic and population genomics: [https://github.com/tbringloe/WGS-NOVAC](https://github.com/tbringloe/WGS-NOVAC/tree/DFO_workflow_2023)

Reference genome assembly and annotation: [https://github.com/tbringloe/Monodontid_assemblies_2023](https://github.com/tbringloe/Monodontid_assemblies_2023)

**Citations for programs:**

Bankevich, A., et al. (2012). SPAdes: A new genome assembly algorith and its applications to single-cell sequencing. Journal of Computational Biology. 19: 455-477.

Bolger, A. M., Lohse, M., & Usadel, B. (2015). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics. 30, 2114-20.

Danacek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham. A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. Gigascience, 10, giab008. https://github.com/samtools/bcftools.

Darling, A. C. E., Mau, B., Blattner, F. R., & Perna, N. T. (2004). Mauve: multiple alignment of conserved genomic sequence with rearrangements. Genome Research, 14, 1394-1403.

Dierckxsens, N., Mardulyn, P., & Smits, G. (2017). NOVOPlasty: de novo assembly of organelle genomes from whole genome data. Nucleic Acids Research, 45, e18. https://github.com/ndierckx/NOVOPlasty

Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9, 357-9. http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup. (2009). The sequence alignment/map format and SAMtools. Bioinformatics 25: 2078-9. http://www.htslib.org/

Steinegger, M, Söding, J. (2017) MMseqs2 enables sensitive protein sequence searching for the analysis of massive datasets. Nature Biotechnology, 35, 1026-1028.
