# Basic Introduction to Bioinformatics for Undergraduate Teaching

__Main author:__  Trevor T. Bringloe  
__Affiliation:__  Fisheries and Oceans Canada (DFO)   
__Group:__        NA   
__Location:__     New Brunswick, Canada  
__Affiliated publication:__  
__Contact:__      e-mail: tbringloe@gmail.com | tel: (506)-259-2288


- [Objective](#objective)
- [Project Summary](#project-summary)

- [Acknowledgements](#acknowledgements)
- [References](#references)


## Objective
The repository provides a basic introduction to bioinformatics intended for upper level undergraduate students. The content covers read quality control, assembly of short reads, read mapping, and similarity searches of DNA barcodes against NCBI (BLAST searches) and the Barcode of Life Data Systems (BOLD). The content required to complete the project and lab assignments detailed here are not provided in their entirety due to data storage limitations on github. However, you may request missing content via email.

## Project Summary
Macroalgae (seaweeds) represent a conglomerate of species, particularly turf forming species that grow entangled with one another close to rocky coastal substrate. Moreover, the seaweeds themselves are home to a host of bacterial and eukaryotic taxa, many of which are unknown to science. Seaweed holobiomes consist not just of an ecosystem of species, but also relationships, some obligate, some passive, some mutualistic, some parasitic. The complexity of the number of species living in seaweed holobiomes and their relation to one another simmply outscales human comprehension. Nonetheless, we can use big-data approaches to at least gleam some insight into these complex biological worlds.

Bioinformatics offers one such glimpse into the structure of holobiomes. August 7-12, 2022. I was fortunate enough to join colleagues on a research cruise to collect seaweeds in Tjongspollen Fjord, south of Bergen, Norway as part of a taxonomic survey funded by the Norwegian Taxonomic Inititiative. Prior to the cruise, I also dove several sites in the area of Bergen with colleagues from France. We were able to collect and document algal turfs, including photo evidence of the turf samples as seen from a dissecting microscope. These samples were subsequently preserve in silica, had DNA extracted at the University of New Brunswick, Canada, and were sequenced at Genome Quebec using the NovaSeq6000 platform, targeting 100 million 150 bp paired-reads reads in 25 samples (though closer to 3 billion reads were generated in total). 

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
