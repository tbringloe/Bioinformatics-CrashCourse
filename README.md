# Basic Introduction to Bioinformatics for Undergraduate Teaching

__Main author:__  Trevor T. Bringloe  
__Affiliation:__  University of New Brunswick (UNB)   
__Group:__        Biology Department   
__Location:__     New Brunswick, Canada  
__Affiliated publication:__  
__Contact:__      e-mail: tbringloe@gmail.com | tel: (506)-259-2288


- [Objective](#objective)
- [Project Summary](#project-summary)
- [Tutorial 1 Linux based commands](#tutorial-1-linux-based-commands)
- [Tutorial 2 High Performance Computing](#tutorial-2-high-performance-computing)
- [Lab 1 Assembling LEGO K-mers](#lab-1-assembling-lego-k-mers)
- [Lab 2 Assembling organellar genomes](#lab-2-assembling-organellar-genomes)
- [Lab 3 Distilling Norwegian algal turf read datasets](#lab-3-distilling-norwegian-algal-turf-read-datasets)
- [Acknowledgements](#acknowledgements)
- [References](#references)


## Objective
The repository provides a basic introduction to bioinformatics intended for upper level undergraduate students. The content covers read quality control, assembly of short reads, read mapping, and similarity searches of DNA barcodes against NCBI ([BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_PROGRAMS=megaBlast&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on&LINK_LOC=blasthome) searches) and the Barcode of Life Data Systems ([BOLD](http://v4.boldsystems.org/)). The content required to complete the project and lab assignments detailed here are not provided in their entirety due to data storage limitations on github, or otherwise because materials are physical (i.e. lego assembly lab).

## Project Summary
Macroalgae (seaweeds) represent a conglomerate of species, particularly turf forming species that grow entangled with one another close to rocky coastal substrate. Moreover, the seaweeds themselves are home to a host of bacterial and eukaryotic taxa, many of which are unknown to science. Seaweed holobiomes are not simply a grouping of species, they represent an ecosystem of relationships, some obligate, some passive, some mutualistic, and some parasitic. The complexity of seaweed holobiomes (and holobiomes in general) in term of the number of species and their relation to one another simmply outscales human comprehension. Nonetheless, we can use big-data approaches to gleam some insight into these complex biological worlds.

Bioinformatics offers one such glimpse into the structure of holobiomes. August 7-12, 2022. I was fortunate enough to join colleagues on a research cruise to collect seaweeds in Tjongspollen Fjord, south of Bergen, Norway as part of a taxonomic survey funded by the Norwegian Taxonomic Inititiative. Prior to the cruise, I also dove several sites in the area of Bergen with colleagues from France. We were able to collect and document algal turfs, including photo evidence of the turf samples as seen from a dissecting microscope. These samples were subsequently preserve in silica, had DNA extracted at the University of New Brunswick, Canada, and were sequenced at Genome Quebec using the NovaSeq6000 platform, targeting 100 million 150 bp paired-reads reads in 25 samples (though closer to 3 billion reads were generated in total). Our objective was to infer the presence of inconspicous species using molecular data that were otherwise escaping labour intensive sorting under the microscope. Because we were interested in species level inferences, we used a whole genome approach rather than traditional metabarcoding (see introduction to bioinformtics lecture).

Using the tutorials provided here, students are expected to take these whole genome sequencing datasets, and distil them into usable information that would allow inferences regarding species present in algal turf samples. Students should consult the files provided here, including sample metadata, and follow links to relevant sites for more information. Some of the computationally intensive steps have been completed a priori, but labs nonetheless guide students through the concepts underpinning bioinformatics otherwise carried out in a high performance computing environment. A project report detailing species found within a particular algal turf sample is expected, including introduction, methods, results, and discussion, along with supplemental code and one figure generated in R.

## Tutorial 1 Linux based commands

Most students are used to operating systems that are highly visual and interactive (i.e. the user clicks icons to issue commands). Linux based computing, however, allows the user to execute desired actions using the command line interface. This provides the user with much greater flexibility for manipulating files in an efficient manner, while still enabling more complex functions to be executed using various built-in and external scripting languages. This tutorial provides basic commands used in linux. Students are not required to do anything for this tutorial, though ideally, students would have access to a linux based environment (e.g. a virtual box or software such as Ubuntu) and test commands themselves to foster comprehension. This tutorial is basically meant to supplement learning and allow students to better understand the commands provided in labs. It is by no means comprehensive; the intention is to put these tools in the hands of up and coming biologists. Most people learn how to drive without ever knowing exactly how their vehicle works.

A couple further things to note. You cannot visualize or interact with linux file systems the same way we click and see file heirarchies in windows or apple OS. Rather, the user must create and familiarize themselves with file heirarchies, and use file pathways (e.g. /home/trevor/data) to move about and interact with their environment and the files therein. Some third party programs can help, however, with connecting and navigating this environment, specifically if the environment is remote and you must login to it (e.g. [Filezilla](https://filezilla-project.org/). A good SSH client for remote work is [PuTTy](https://www.putty.org/). 
A hashtag indicates to not execute proceeding text as commands. They are therefore a useful tool for annotating scripts.
You can use the TAB key to finish a filename or pathway. If multiple options are present, they will be displayed. This can save typing out everything.
The up arrow key can be used to recall previously used commands, in order of more recent
Many commands will display output at the command line interface. End any command with '> file.name' to output command to a the specified file (i.e. something more informative than file.name.
You can set environment variables, which sub in text in specified locations. This can cut down on text, and facilitate more advanced commands such as looping through sample IDs
Watch out for double and single quotations, they must be closed, otherwise your command/script will break. Single quotes tell the sytem to interpret text as written.
You can pipe the output of one command into a new command using |. This helps avoid bulky writing and writing intermediate files.

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
mkdir trimmed_reads # creates new directoru trimmed_reads in present working directory

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
  echo insert trimmomatic code here
else
  echo User has specified to skip raw read trimming
fi
exit 0
```

## Tutorial 2 High Performance Computing
Bioinformatics usually involves working with large datasets containing many GBs or even TBs of sequence information. Because there is so much information to work with, we cannot use normal laptops or desktops. We need computers with a lot of storage (i.e. disk) space for storing and writing files, and a lot of RAM (Random Access Memory) to work with lots of data at a given moment. Furthermore, resources can be distributed across nodes and cores in a system, allowing users to run many commands at the same time, or divide tasks into many smaller tasks, thus allowing users to complete commands in a shorter time frame. These needs are met by High Performance Computers. Canadian Universities share access to Compute Canada servers, but many institutions host private servers. One of the challenges on a shared system is tasks must be executed in an efficient manner; the system cannot be overloaded with tasks (otherwise this would quickly crash the system) and computational resources must be used in the most efficient manner possible to ensure the most amount of computation gets completed in the shortest amount of linear time. To facilitate this, shared servers use slurm scripts, which users use to submit their tasks to perform. Tasks are then queued and run when resources become available.

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
In this example, the user has specified the account under which to run tasks (computation resources are carefully allocated and monitored across accounts to ensure resources are utilized in an equitable manner). The maximum amount of time to run the task(s) (i.e. wall time) has been specified at just under 4 hours. The job name has been specified as 'fastqc_raw_check', which will appear when monitoring the job status. -n specifies to use 4 threads or tasks when running the command, and --mem specifies the amount of extra memory or RAM to allocate to the job (usually there is a maximum amount of memory provided per node, so users will need to request additional memory if needed). -N indicates to run on a single node (most actions can be performed on a single node, except in more advanced tasks that require running seperate tasks on multiple nodes). --mail-user specified an email where notifications can be sent when the job enters a new status (e.g. --mail-type=BEGIN, END, FAIL). Shared environments typically have staff who can install software of interest; as such, programs needed to run commands are loaded as modules, and become available for use once loaded. Oftentimes, system dependencies must also be loaded to run particular programs. Job commands are the commands you want to run.



## Lab 1 Assembling LEGO k-mers
Following the lecture introducing concepts related to bioinformatics, and a tutorial on linux based High Performance Computing, complete the following lab. Note, the following commands require some user input to run, i.e. defining particular variables and establishing sample lists.

*Exercise 1: Evaluating read quality*

Before working with read files, it is critical to evaluate quality and remove error prone data. You have been provided [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) reports generated for several read datasets, before and after read trimming. The following command was performed using [TRIMMOMATIC](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf):
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

Investigate the FASTQC reports and answer the following questions.
1.	Briefly define the variables specified stated above by consulting the [TRIMMOMATIC](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) manual.

2.	Fill the following table:

| File name | # of read before trimming | # of reads after trimming | Average read quality before trimming | Average read quality after trimming | GC content before trimmming | GC content after trimming |
| --- | --- | --- | --- | --- | --- | --- |
| TTB000601 |
| TTB000606 |
| TTB000611 |

3.	List notable differences between the before and after trimmming reports pertaining to the quality figures.

4.	Visit this [bad quality report](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html). What are potential issues with the read data and what trimming could be done to improve the quality of the specified read dataset?

*Exercise 2: Build an assembly graph* 

You have been provided long LEGO pieces with sequences labelled on them. These represent k-mers. Start by turning all pieces to the side labelled with the red circle. This represents your read data as k-mers at length XX. Start overlapping and connecting the pieces at your k-mer size-1 (i.e. XX, your pieces should overlap almost entirely, extending the graph one node at a time). Note when multiple pathways occur, and whether m there are multiple unconnected pathways. 

1. Draw a diagram of your assembly graph (an example can be provided in lab, or see [BANDAGE](https://rrwick.github.io/Bandage/) for conceptual examples, though ours will be much simpler than what is shown here). Provide a figure title, indicating the k-mer size used.

Next, repeat the exercise, but flip all the LEGO pieces to the site labelled with a blue square. This represents your read data as k-mers at length XX. As before, start overlapping and connecting the pieces at your k-mer size-1.

2. As before, draw a diagram of your assembly graph. Provide a figure title, indicating the k-mer size used.

3. Here we have assumed unlimited coverage for the assemblies. What was the effect of k-mer size on the assembly graphs?

4. Had coverage been limited, how might the assembly at the larger k-mer size been impacted?

*Exercise 3: Blast assembled sequences*

You have now assembled several sequences. Compare the sequences to NCBI's database using [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_PROGRAMS=megaBlast&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on&LINK_LOC=blasthome) and to [BOLD](http://v4.boldsystems.org/) records under the "identification" tab to assign taxonomic information to the sequences.

1.	What are the top hits for each sequence, and to what percentage to the sequences match existing records?

## Lab 2 Assembling organellar genomes

## Lab 3 Distilling Norwegian algal turf read datasets

For this lab, students will go through the various steps to distill large raw read files to curated BLAST reports for DNA barcode regions, which can be used to interpret species present within the algal turf samples. The lab will process the following three samples:

## Acknowledgements
All the inspirational students, postdocs, mentors, DFO research scientists, and forum junkies across the globe who contributed to my own learning journey in bioinformatics, of which this workflow is a direct result.

## References
**Other potentially useful tutorials:**

Variant calling for phylogenomic and population genomics: [https://github.com/tbringloe/WGS-NOVAC](https://github.com/tbringloe/WGS-NOVAC/tree/DFO_workflow_2023)

Reference genome assembly and annotation: https://github.com/tbringloe/Monodontid_assemblies_2023

**Citations for programs:**

Bolger, A. M., Lohse, M., & Usadel, B. (2015). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics. 30, 2114-20.

Danacek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham. A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. Gigascience, 10, giab008. https://github.com/samtools/bcftools.

Dierckxsens, N., Mardulyn, P., & Smits, G. (2017). NOVOPlasty: de novo assembly of organelle genomes from whole genome data. Nucleic Acids Research, 45, e18. https://github.com/ndierckx/NOVOPlasty

Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9, 357-9. http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup. (2009). The sequence alignment/map format and SAMtools. Bioinformatics 25: 2078-9. http://www.htslib.org/
