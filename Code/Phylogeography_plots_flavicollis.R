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

#Convert from vcf to my_dnabin
vcf_file <- ("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/batch_301.vcf_final_filtered.recode.vcf")
vcf <- read.vcfR(vcf_file, verbose = FALSE)
vcf
record <- 130

#Upload population map
pop_apodemus<- read.csv ("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/pop_map_flav_filtered_pop.csv", sep = "\t", header = FALSE)
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

my_genind@pop<- as.factor(pop_apodemus_ordered$V3)
my_genind
my_genind@pop


#my_genind@strata
#Barplot number of samples per population
barplot(table(pop(my_genind)), col = funky(5), las=3, ylab = "Sample size")
#map with samples for analys
library(ggmap)
register_google(key = "AIzaSyBTBr0ek5JXmoQ0VvWahWhaSwMkKcvHdJI")  
map <- get_map(location = 'Europe', zoom = 3)
ggmap(map)
map_colors<- c( "red", "orange", "purple", "blue", "green", "dark green", "grey", "brown", "black")
map_coordinates<- read.csv("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/pop_coordinates.csv", sep= "\t", header= FALSE)

head(map_coordinates)
map_coordinates
names(map_coordinates)[5]<-paste("Number of Samples")
mapPoints <- ggmap(map) +
  geom_point(aes(x = V4, y = V3), data = map_coordinates , size= 3)

mapPoints
mapPointsLegend <- mapPoints +
  #   scale_size_area(breaks = sqrt(c(1, 2, 3, 4, 5, 10)), labels = c(1, 2, 3, 4, 5, 10), name = "Number of samples")
  geom_point(data= map_coordinates, aes(x = V4, y = V3, colour=as.factor(map_coordinates$`Number of Samples`)))+
  labs(colour ="Number of Samples")+
  scale_colour_manual(values=map_colors)
mapPointsLegend
par(mfrow=c(1,1))

#pca
x.apo<- tab(my_genind, freq=TRUE, NA.method="mean")
pca.apo <- dudi.pca(x.apo, center = TRUE, scale = FALSE)
s.label(pca.apo$li)

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank())


p<-ggplot(pca.apo$li, aes(x=Axis1,y=Axis2, label=indNames(my_genind)))+
  theme+xlab("PCA1=4.39")+
  theme+ylab("PCA2=3.30")+
  geom_text(label=indNames(my_genind))+
  geom_text((aes(color=my_genind$pop)))+
  theme(legend.position="none")
p
p
eig.perc <- 100*pca.apo$eig/sum(pca.apo$eig)
head(eig.perc)



barplot(pca.apo$eig)+
  title(main="Eigenvalues")
loadingplot(pca.apo$c1^2)
pca.apo


#dapc
grp <- find.clusters(my_genind, max.n.clust=20)
#180, 3
grp
dapc1 <- dapc(my_genind, grp$grp)
#190,3

print<- as.data.frame(grp$grp)
dapc1
scatter(dapc1)
scatter(dapc1, posi.da="bottomright", bg="white", pch=17:22)
scatter(dapc1, scree.da=FALSE, bg="white", pch=20,  cell=0, cstar=0, col=myCol, solid=.4,
        cex=1,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:3))
myCol <- c("orange","red", "green")
scatter(dapc1, posi.da="bottomright",  bg="white",
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="bottomleft")
scatter(dapc1, scree.da=FALSE, bg="white", pch=20,  cell=0, cstar=0, col=myCol,
        cex=3 ,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:3), posi.leg = "none")
dapc1

#######compoplot
dapc.results <- as.data.frame(dapc1$posterior)
dapc.results$pop <- pop(my_genind)
dapc.results$indNames <- rownames(dapc.results)

library(reshape2)
dapc.results <- melt(dapc.results)
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = myCol) 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
p

####map
setwd("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/")
#prepare geno file
geno_sylv<-vcf2geno("batch_301.vcf_final_filtered.recode.vcf", output.file = "flav_05_filtered")

geno_sylv<- read.geno("flav_05_filtered.geno")
genotype_sylv<-geno_sylv

#CHange 9 by NA
genotype_sylv[genotype_sylv ==9]<- NA
genotype_sylv
#run analysis
bj.at = snmf("flav_05_filtered.geno", K = 1:18, ploidy = 2, entropy = T,
             CPU = 8, project = "new")
plot(bj.at, col = "blue4", cex = 1.4, pch = 19)

qmatrix2 = Q(bj.at, K=2)
qmatrix3 = Q(bj.at, K=3)
qmatrix4 = Q(bj.at, K=4)
qmatrix5 = Q(bj.at, K=5)
qmatrix6 = Q(bj.at, K=6)
qmatrix7 = Q(bj.at, K=7)
qmatrix8 = Q(bj.at, K=8)
qmatrix9= Q(bj.at, K=9)
qmatrix10= Q(bj.at, K=10)
qmatrix11= Q(bj.at, K=11)
qmatrix12= Q(bj.at, K=12)
#coordenadas
pop_map_sylv<- read.csv("~/Documents/Resultados/Phylogeography_apodemus/sylvaticus/Pop_map/pop_map_europe_sylv_filtered_coordinates_populations.csv", header=FALSE, sep = "\t")
pop_map_sylv

coordinates_sylv<- pop_map_sylv[, 4:5]
coordinates_sylv
str(coordinates_sylv)
coord_sylv<- coordinates_sylv
str(coord_sylv)
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")
asc.raster="~/Downloads/etopo1(1).asc"
grid=createGridFromAsciiRaster(asc.raster)
constraints=getConstraintsFromAsciiRaster(asc.raster, cell_value_min=0)

lColorGradients = list( 
  c("gray95",brewer.pal(9,"Reds")),
  c("gray95",brewer.pal(9,"Greens")),
  c("gray95",brewer.pal(9,"Blues")),
  c("gray95",brewer.pal(9,"YlOrBr")),
  c("gray95",brewer.pal(9,"RdPu")),
  c("gray95",brewer.pal(9,"Purples")),
  c("gray95",brewer.pal(9,"Oranges")),
  c("gray95",brewer.pal(9,"Greys"))
)
# Use display.brewer.all() to display color gradients provided by RColorBrewer
display.brewer.all()
k2_matrix<-as.matrix(K2_ordered[-1])
map2<-maps(matrix =k2 , coord_sylv, grid, constraints, method = "max",
           main = "Ancestry coefficients", xlab = "Longitude", ylab = "Latitude", cex = .5, colorGradientsList=lColorGradients)
k3_matrix<-as.matrix(K3_ordered[-1])
map3<-maps(matrix = qmatrix3, coord_sylv, grid, constraints, method = "max",
           main = "Ancestry coefficients", xlab = "Longitude", ylab = "Latitude", cex = .5, colorGradientsList=lColorGradients)
k4_matrix<-as.matrix(K4_ordered[-1])
map4<-maps(matrix = qmatrix4, coord_sylv, grid, constraints, method = "max",
           main = "Ancestry coefficients", xlab = "Longitude", ylab = "Latitude", cex = .5, colorGradientsList=lColorGradients)
k5_matrix<-as.matrix(K5_ordered[-1])
map5<-maps(matrix = qmatrix5, coord_sylv, grid, constraints, method = "max",
           main = "Ancestry coefficients", xlab = "Longitude", ylab = "Latitude", cex = .5, colorGradientsList=lColorGradients)
k6_matrix<-as.matrix(k6_ordered[-1])
map6<-maps(matrix = qmatrix6, coord_sylv, grid, constraints, method = "max",
           main = "Ancestry coefficients", xlab = "Longitude", ylab = "Latitude", cex = .5, colorGradientsList=lColorGradients)
k7_matrix<-as.matrix(k7_ordered[-1])
map7<-maps(matrix = qmatrix7, coord_sylv, grid, constraints, method = "max",
           main = "Ancestry coefficients", xlab = "Longitude", ylab = "Latitude", cex = .5, colorGradientsList=lColorGradients)
map8<-maps(matrix = qmatrix8, coord_sylv, grid, constraints, method = "max",
           main = "Ancestry coefficients", xlab = "Longitude", ylab = "Latitude", cex = .5, colorGradientsList=lColorGradients)


#####structure
Admixcol<- c("green","orange","red","blue","purple", "yellow", "brown", "black")
tble=read.csv("~/Desktop/structure_sylv2.csv", header=TRUE, sep="\t")
barplot(t(as.matrix(tble[3:4])), col=Admixcol,
        ylab="Ancestry", border=NA, names.arg = tble$X.1,las =2, cex.names = 0.8)

tble=read.csv("~/Desktop/structure_sylv2.csv", header=TRUE, sep="\t")
barplot(t(as.matrix(tble[6:8])), col=Admixcol,
        ylab="Ancestry", border=NA, names.arg = tble$X,las =2, cex.names = 0.8)
tble=read.csv("~/Desktop/structure_sylv2.csv", header=TRUE, sep="\t")
barplot(t(as.matrix(tble[10:13])), col=Admixcol,
        ylab="Ancestry", border=NA, names.arg = tble$X,las =2, cex.names = 0.8)
tble=read.csv("~/Desktop/structure_sylv2.csv", header=TRUE, sep="\t")
barplot(t(as.matrix(tble[15:19])), col=Admixcol,
        ylab="Ancestry", border=NA, names.arg = tble$X,las =2, cex.names = 0.8)
tble=read.csv("~/Desktop/structure_sylv2.csv", header=TRUE, sep="\t")
barplot(t(as.matrix(tble[21:26])), col=Admixcol,
        ylab="Ancestry", border=NA, names.arg = tble$X,las =2, cex.names = 0.8)


#####admixture
cv_admixture<- read.csv("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/Admixture/cv.csv", sep = "\t", header = FALSE)
cv_admixture
cv_admixture$V1<- order(levels(cv_admixture$V1))
cv_admixture$V1
str(cv_admixture)
ggplot()+
  geom_point(data=cv_admixture, aes(x=V1, y=V2), color="blue")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank())+
  xlab("K")+ ylab("CV error")+
  theme(axis.title.x = element_text(color="black", size=12),
        axis.title.y = element_text(color="black", size=12),
        text=element_text((size=14)))


######k2
k2<-read.table("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/Admixture/1/plink.2.Q")
pop_map<-read.csv("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/pop_map_flav_filtered_pop.csv", sep="\t", header = FALSE)
K2<-cbind(pop_map$V1, k2)
K2
order_list<-read.csv("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/pop_map_fla_filtered_admix.csv", sep = "\t", header = FALSE)
K2_ordered<-K2[ order(match(K2$`pop_map$V1`, order_list$V1)),]
K2_ordered
K2_ordered<- K2_ordered[c(1,2,3)]
bark2<-barplot(t(as.matrix(K2_ordered[,2:3])), col=Admixcol,
               ylab="Ancestry",names.arg = order_list$V1, border=NA,las =2, cex.names = 0.8)

#####k3
k3<-read.table("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/Admixture/1/plink.3.Q")
pop_map<-read.csv("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/pop_map_flav_filtered_pop.csv", sep="\t", header = FALSE)
K3<-cbind(pop_map$V1, k3)
order_list<-read.csv("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/pop_map_fla_filtered_admix.csv", sep = "\t", header = FALSE)
K3_ordered<-K3[ order(match(K3$`pop_map$V1`, order_list$V1)),]
K3_ordered
K3_ordered<- K3_ordered[c(1,4,2,3)]
bark3<-barplot(t(as.matrix(K3_ordered[,2:4])), col=Admixcol,
               ylab="Ancestry",names.arg = order_list$V1, border=NA,las =2, cex.names = 0.8)

#####k4
k4<-read.table("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/Admixture/1/plink.4.Q")
pop_map<-read.csv("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/pop_map_flav_filtered_pop.csv", sep="\t", header = FALSE)
pop_map
K4<-cbind(pop_map$V1, k4)
K4
order_list<-read.csv("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/pop_map_fla_filtered_admix.csv", sep = "\t", header = FALSE)
order_list
K4_ordered<-K4[ order(match(K4$`pop_map$V1`, order_list$V1)),]
K4_ordered
K4_ordered<- K4_ordered[c(1,5,4,3,2)]
bark4<-barplot(t(as.matrix(K4_ordered[,2:5])), col=Admixcol,
               ylab="Ancestry",names.arg = order_list$V1, border=NA,las =2, cex.names = 0.8)


####k5
k5<-read.table("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/Admixture/1/plink.5.Q")
pop_map<-read.csv("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/pop_map_flav_filtered_pop.csv", sep="\t", header = FALSE)
K5<-cbind(pop_map$V1, k5)
order_list<-read.csv("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/pop_map_fla_filtered_admix.csv", sep = "\t", header = FALSE)
K5_ordered<-K5[ order(match(K5$`pop_map$V1`, order_list$V1)),]
K5_ordered
K5_ordered<- K5_ordered[c(1,3,5,6,2,4)]
bark5<-barplot(t(as.matrix(K5_ordered[,2:6])), col=Admixcol,
               ylab="Ancestry",names.arg = order_list$V1, border=NA,las =2, cex.names = 0.8)

#####k6
k6<-read.table("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/Admixture/1/plink.6.Q")
pop_map<-read.csv("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/pop_map_flav_filtered_pop.csv", sep="\t", header = FALSE)
K6<-cbind(pop_map$V1, k6)
order_list<-read.csv("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/pop_map_fla_filtered_admix.csv", sep = "\t", header = FALSE)
K6_ordered<-K6[ order(match(K6$V1, order_list$V1)),]
K6_ordered
K6_ordered<- K3_ordered[c(1,3,2)]
bark6<-barplot(t(as.matrix(K6_ordered[,2:7])), col=Admixcol,
               ylab="Ancestry",names.arg = order_list$V1, border=NA,las =2, cex.names = 0.8)


####maps polimorphic loci

d <- read.csv("~/Documents/Resultados/Phylogeography_apodemus/flavicollis/p5_without_FR15/filtered/genetic_diversity.csv", header = TRUE, sep="\t")
str(d)
library(sp)
dsp <- SpatialPoints(d[,3:4], proj4string=CRS("+proj=longlat +datum=NAD83"))
dsp <- SpatialPointsDataFrame(dsp, d)
dsp



require(maptools)

getinfo.shape("~/Downloads/TM_WORLD_BORDERS_SIMPL-0.3/TM_WORLD_BORDERS-0.3.shp")
world.map <- readShapeSpatial("~/Downloads/TM_WORLD_BORDERS_SIMPL-0.3/TM_WORLD_BORDERS-0.3.shp")
world.map@data = world.map@data[world.map@data$AREA > 30000,]
plot(world.map)


# define groups for mapping
cuts <- c(0,10,20,30,40,50,60,70,80,90,100)
# set up a palette of interpolated colors
blues <- colorRampPalette(c('yellow', 'orange', 'blue', 'dark blue'))
pols <- list("sp.polygons", world.map, fill = "lightgray")
a<-spplot(dsp, "Polymorphic.Loci", cuts=cuts, col.regions=blues(5), sp.layout=pols, pch=20, cex=2)
a
####maps heterozygosity


# define groups for mapping
cuts <- c(0,0.05,0.10,0.150,0.20,0.25,0.30)
# set up a palette of interpolated colors
blues <- colorRampPalette(c('yellow', 'orange', 'purple','blue', 'dark blue'))
pols <- list("sp.polygons",world.map, fill = "lightgray")
b<-spplot(dsp, "Exp.Het", cuts=cuts, col.regions=blues(5), sp.layout=pols, pch=20, cex=2)
b

####maps Pi

# define groups for mapping
cuts <- c(0,0.05,0.10,0.150,0.20,0.25,0.30,0.35)
# set up a palette of interpolated colors
blues <- colorRampPalette(c('yellow', 'orange', 'purple','blue', 'dark blue'))
pols <- list("sp.polygons", world.map, fill = "lightgray")
c<-spplot(dsp, "Pi", cuts=cuts, col.regions=blues(5), sp.layout=pols, pch=20, cex=2)
c
