#Code used to build the graphs used in the Discussion:
#Barplot retained reads
quaddrad_test<- read.csv("~/Desktop/quaddrad_test.csv", sep = ",", header = TRUE)
quaddrad_test$Samples<- as.factor(quaddrad_test$Samples)
p4 <- ggplot(data=quaddrad_test, aes(y = Percentage, x = Samples, fill = Type, label=Percentage))+
  geom_bar(stat="identity")+
  theme_bw()+
  geom_text(aes(label = paste0(quaddrad_test$Percentage,"%"),y=Percentage),size = 4, 
position = position_stack(vjust = 0.5))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p4<-p4 + theme(legend.position="bottom")
p4


#Barplot Assembled loci, polymorphic loci and SNPs per sample
quaddrad_assemble<-read.csv("~/Desktop/Suplementary_materials_thesis/quaddrad_test_samples_assembled.csv", sep = ",", header = TRUE)
ggplot(data=quaddrad_assemble, aes(x=order, y=value, fill=type)) +
  geom_col(stat="identity", color="black", position=position_dodge())+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
theme(legend.position="bottom")
        
# Barplot coverage samples
quaddrad_cov<-read.csv("~/Desktop/Suplementary_materials_thesis/coverage_quaddrad_Test.csv", sep = ",", header = FALSE)
quaddrad_cov$V1<-factor(quaddrad_cov$V1, levels=quaddrad_cov$V1)
quaddrad_cov
ggplot(data=quaddrad_cov, aes(x=V1, y=V2, fill=V3)) +
  geom_col(stat="identity", color="black", position=position_dodge())+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())+
  xlab("Sample")+ylab("Coverage")+
  theme(legend.title = element_blank())+
  theme(legend.position="bottom")

