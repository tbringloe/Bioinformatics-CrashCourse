# Basic Introduction to Bioinformatics for Undergraduate Teaching

__Main author:__  Trevor T. Bringloe  
__Affiliation:__  Fisheries and Oceans Canada (DFO)   
__Group:__        NA   
__Location:__     New Brunswick, Canada  
__Affiliated publication:__  
__Contact:__      e-mail: tbringloe@gmail.com | tel: (506)-259-2288


- [Objective](#objective)
- [Project Summary](#project-summary)
- [Lab 1 Assembling LEGO K-mers](#lab-1-assembling-lego-k-mers)
- [Acknowledgements](#acknowledgements)
- [References](#references)


## Objective
The repository provides a basic introduction to bioinformatics intended for upper level undergraduate students. The content covers read quality control, assembly of short reads, read mapping, and similarity searches of DNA barcodes against NCBI ([BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_PROGRAMS=megaBlast&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on&LINK_LOC=blasthome) searches) and the Barcode of Life Data Systems ([BOLD](http://v4.boldsystems.org/)). The content required to complete the project and lab assignments detailed here are not provided in their entirety due to data storage limitations on github, or otherwise because materials are physical (i.e. lego assembly lab).

## Project Summary
Macroalgae (seaweeds) represent a conglomerate of species, particularly turf forming species that grow entangled with one another close to rocky coastal substrate. Moreover, the seaweeds themselves are home to a host of bacterial and eukaryotic taxa, many of which are unknown to science. Seaweed holobiomes are not simply a grouping of species, they represent an ecosystem of relationships, some obligate, some passive, some mutualistic, and some parasitic. The complexity of seaweed holobiomes (and holobiomes in general) in term of the number of species and their relation to one another simmply outscales human comprehension. Nonetheless, we can use big-data approaches to gleam some insight into these complex biological worlds.

Bioinformatics offers one such glimpse into the structure of holobiomes. August 7-12, 2022. I was fortunate enough to join colleagues on a research cruise to collect seaweeds in Tjongspollen Fjord, south of Bergen, Norway as part of a taxonomic survey funded by the Norwegian Taxonomic Inititiative. Prior to the cruise, I also dove several sites in the area of Bergen with colleagues from France. We were able to collect and document algal turfs, including photo evidence of the turf samples as seen from a dissecting microscope. These samples were subsequently preserve in silica, had DNA extracted at the University of New Brunswick, Canada, and were sequenced at Genome Quebec using the NovaSeq6000 platform, targeting 100 million 150 bp paired-reads reads in 25 samples (though closer to 3 billion reads were generated in total). Our objective was to infer the presence of inconspicous species using molecular data that were otherwise escaping labour intensive sorting under the microscope. Because we were interested in species level inferences, we used a whole genome approach rather than traditional metabarcoding (see introduction to bioinformtics lecture).

Using the tutorials provided here, students are expected to take these whole genome sequencing datasets, and distil them into usable information that would allow inferences regarding species present in algal turf samples. Students should consult the files provided here, including sample metadata, and follow links to relevant sites for more information. Some of the computationally intensive steps have been completed a priori, but labs nonetheless guide students through the concepts underpinning bioinformatics otherwise carried out in a high performance computing environment. A project report detailing species found within a particular algal turf sample is expected, including introduction, methods, results, and discussion, along with supplemental code and one figure generated in R.

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
```

Investigate the FASTQC reports and answer the following questions.
1.	Briefly define the variables specified stated above by consulting the [TRIMMOMATIC](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) manual.

2.	Fill the following table:

| File name | # of read before trimming | # of reads after trimming | Average read quality before trimming | Average read quality after trimming | GC content before trimmming | GC content after trimming |
| --- | --- | --- | --- | --- | --- | --- |
| TTB000600 |
| TTB000606 |
| TTB000611 |

3.	List notable differences between the before and after trimmming reports pertaining to the quality figures.

4.	Visit this [bad quality report](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html). What are potential issues with the read data and what trimming could be done to improve the quality of the specified read dataset?

Exercise 2: Build an assembly graph 

and answer the following questions.

Exercise 3: Blast assembled sequences
1.	What are the top BLAST hits for each sequence assembled with the lego k-mers?


## Acknowledgements
All the inspirational students, postdocs, mentors, DFO research scientists, and forum junkies across the globe who contributed to my own learning journey in bioinformatics, of which this workflow is a direct result.

## References
**Other potentially useful tutorials:**

Variant calling for phylogenomic and population genomics: [https://github.com/tbringloe/WGS-NOVAC](https://github.com/tbringloe/WGS-NOVAC/tree/DFO_workflow_2023)

Reference genome assembly and annotation: https://github.com/tbringloe/Monodontid_assemblies_2023

**Citations for programs:**

Danacek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham. A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. Gigascience, 10, giab008. https://github.com/samtools/bcftools.

Dierckxsens, N., Mardulyn, P., & Smits, G. (2017). NOVOPlasty: de novo assembly of organelle genomes from whole genome data. Nucleic Acids Research, 45, e18. https://github.com/ndierckx/NOVOPlasty

Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9, 357-9. http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup. (2009). The sequence alignment/map format and SAMtools. Bioinformatics 25: 2078-9. http://www.htslib.org/
