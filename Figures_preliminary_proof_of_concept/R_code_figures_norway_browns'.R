library(ggplot)
library(viridis)
library(tidyverse)

#read in data on OTUS
OTUS <- read.csv("OTU_plot22xii23.txt", sep="\t")
#Order data
OTUS$Classification <- factor(OTUS$Classification, levels=c("Ectocarpales", "Phaeophyceae", "Haptophyceae", "Pelagophycidae", "Chrysophyceae", "Bacillariophyceae", "Eukaryota", "Chlorophyta", "Coscinodiscophyceae", "Prasinophyceae", "Stylonematophyceae", "Compsopogonophyceae", "Bangiophyceae", "Florideophyceae"))
OTUS$SpecimenID <- factor(OTUS$SpecimenID, levels=c("TTB000618", "TTB000600", "TTB000610", "TTB000612", "TTB000611", "TTB000629", "TTB000606", "TTB000636", "TTB000615", "TTB000609", "TTB000614", "TTB000601", "TTB000539", "TTB000617", "TTB000605", "TTB000641", "TTB000621", "TTB000538", "TTB000631", "TTB000620", "TTB000637", "TTB000632", "TTB000634", "TTB000639", "TTB000535"))

#Plot data
gg.OTUS <- ggplot(OTUS, aes(y=OTUname,x=SpecimenID, fill=Classification)) + geom_tile(group="Classification") + theme_bw() + facet_grid(Classification ~ ., space="free", scale="free") + scale_fill_manual(values=c("#fde725", "#fde725", "#c8e020", "#90d743", "#5ec962", "#35b779", "#20a486", "#21918c", "#287c8e", "#31688e", "#3b528b", "#443983", "#481f70", "#440154"))
gg.OTUS

#Code to produce stacked bar charts
#Import and reoder as above
gg.rbcL.classes <- ggplot(counts, aes(y=SpecimenID, x=SequencePercentage, fill=Classification)) + geom_bar(position="stack", stat="identity") + theme_bw() + scale_fill_manual(values=c("#fde725", "#c8e020", "#90d743", "#5ec962", "#35b779", "#20a486", "#21918c", "#287c8e", "#31688e", "#3b528b", "#443983", "#481f70", "#440154"))
gg.rbcL.classes

#Save to pdf
ggsave(filename="OTUS_22xii23.pdf", plot=gg.OTUS, height=24, width=12)