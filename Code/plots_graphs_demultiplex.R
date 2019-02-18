#Code to build piecharts for different demultiplex options
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(cowplot)
          
#STACKED BARPLOT
                     
summary<-read.csv("~/Desktop/DEMULTIPLEX_SUMMARY.csv", sep = ",", header = TRUE)
summary 
levels(summary$read)<-c("Ambiguous Barcodes","Ambiguous RAD-cutsite", "Reads containing adapter sequence","Retained Reads")

p4 <- ggplot(data=summary, aes(y = Percentage, x = type, fill = read, label=Percentage))+
  geom_bar(stat="identity")+
  theme(axis.text = element_text(size = 15)) + 
  theme(axis.title = element_text(size = 15)) +
  theme(text = element_text(size = 15)) + 
  geom_text(aes(label = paste0(summary$Percentage/1,"%")), 
                 position = position_stack(vjust = 0.5), size=7)+
  xlab("")+
  theme(legend.title=element_blank())+
  ylab("Percentage of reads")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p4

ggsave("demultiplex.pdf", p4,
       scale = 1, width = 180, height = NA, units = c("mm"),
       dpi = "retina", limitsize = TRUE, device = "pdf", path = NULL)
#barplots
perfect<-read.csv("~/Desktop/Demultiplex_perfect_barcode_rad.csv", sep = ",", header = TRUE)
b1<-read.csv("~/Desktop/Demultiplex_0_mismatches.csv", sep = ",", header = TRUE)
b2<-read.csv("~/Desktop/Demultiplex_1_mismatches.csv", sep = ",", header = TRUE)
b3<-read.csv("~/Desktop/Demultiplex_2_mismatches.csv", sep = ",", header = TRUE)
barplot(perfect)
p1 <- ggplot(data=perfect[1:96,], aes(y = Retained, x = Barcode))+
  geom_bar(stat="identity") + theme(axis.text.x = element_text(size=6, angle=90),
          axis.text.y = element_text(size=6)) + coord_flip()+
  theme(axis.title.x=element_blank())
p2 <- ggplot(data=b1[1:96,], aes(y = Retained, x = Barcode))+
  geom_bar(stat="identity") + theme(axis.text.x = element_text(size=6, angle=90),
                                    axis.text.y = element_text(size=6)) + coord_flip()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(axis.title.x=element_blank())
p3 <- ggplot(data=b2[1:96,], aes(y = Retained, x = Barcode))+
  geom_bar(stat="identity") + theme(axis.text.x = element_text(size=6, angle=90),
                                    axis.text.y = element_text(size=6)) + coord_flip()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(axis.title.x=element_blank())
p4 <- ggplot(data=b3[1:96,], aes(y = Retained, x = Barcode))+
  geom_bar(stat="identity") + theme(axis.text.x = element_text(size=6, angle=90),
                                    axis.text.y = element_text(size=6)) + coord_flip()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(axis.title.x=element_blank())
combo<- plot_grid(p1+theme(legend.position = "none"),
                  p2+theme(legend.position = "none"),
                  p3+theme(legend.position = "none"),
                  p4+theme(legend.position = "none"),
                  align = 'v',
                  labels = c("A", "B", "C", "D"),
                  hjust = -0.1,
                  nrow = 1)
combo

