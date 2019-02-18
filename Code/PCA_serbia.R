library(stats)
library(ade4)
library(ape)
library(adegenet)
library(phangorn)
library(hierfstat)
library(vcfR)
library(pegas)
library(poppr)
library(radiator)
library(ggplot2)

#Convert from vcf to my_dnabin
vcf_file <- ("~/Documents/Resultados/Apodemus_Poland/Demultiplex/concatenated/apodemus/Best_parameter_Paris/test_m2_M4_n5/divergence/populations/serbia/filtered.recode.vcf")
vcf <- read.vcfR(vcf_file, verbose = FALSE)
vcf
record <- 130

#Upload population map
pop_apodemus<- read.csv ("~/Documents/Resultados/Apodemus_Poland/Demultiplex/concatenated/apodemus/Best_parameter_Paris/test_m2_M4_n5/divergence/populations/serbia/pop_map", sep = "\t", header = FALSE)
dim(as.data.frame(pop_apodemus))
pop_map_ordered <- pop_apodemus[order(pop_apodemus$V1),]
pop_apodemus
pop_map_ordered
pop_apodemus[order(pop_apodemus$V3),]

pop_map_ordered

my_genind<-vcfR2genind(vcf)
indNames(my_genind)

names <- as.vector(indNames(my_genind))
names_order<-data.frame(samples=names, row.names = NULL)
names_order
pop_apodemus_ordered<- pop_apodemus[ order(match(pop_apodemus$V1, names_order$samples)),]
as.data.frame(pop_apodemus_ordered)
stopifnot(all(indNames(my_genind) == pop_apodemus_ordered$V1))
indNames(my_genind)
row.names(pop_apodemus_ordered)=pop_apodemus_ordered$V1


levels(pop_apodemus_ordered$V3)<-c("A. flavicollis", "A.sylvaticus", "Undetermined")

my_genind@pop<- as.factor(pop_apodemus_ordered$V3)
my_genind
my_genind@pop


#my_genind@strata
#Barplot number of samples per population
barplot(table(pop(my_genind)), col = funky(5), las=3, ylab = "Sample size")

#pca
x.apo<- tab(my_genind, freq=TRUE, NA.method="mean")
pca.apo <- dudi.pca(x.apo, center = TRUE, scale = FALSE)
s.label(pca.apo$li)
eig.perc <- 100*pca.apo$eig/sum(pca.apo$eig)
head(eig.perc)
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank())
p<-ggplot(pca.apo$li, aes(x=Axis1,y=Axis2))+
  theme+xlab("PCA1=9.81")+
  theme+ylab("PCA2=4.47")+
  geom_point((aes(color=my_genind$pop)), size=2.5)+
  theme(legend.title = element_blank())+
  theme(legend.text=element_text(size=12))
p