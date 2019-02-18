#Fst/(1-fst) vs Physical distances flavicollis
setwd("~/Documents/thesis/")
distances<-read.csv("genetic_diversity_flav_R.csv", sep = ",", header=FALSE)
levels(distances$V6)<-c("Normal-Normal", "Normal-Historical", "Historical-Historical")

Normal<-distances[distances$V6=='Normal-Normal',]
Skin<-distances[distances$V6=='Historical-Historical',]
SkinNormal <-distances[distances$V6=='Normal-Historical',]

my.formula <- y ~ x
p<-ggplot(data=distances, aes(x=V3, y=V5)) + 
  geom_point(aes(color=(distances$V6)))+
  geom_smooth(method='lm')+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)+
  xlab("Distance")+
  ylab("Fst/(1-Fst)")+
  theme(legend.title=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  theme(legend.text=element_text(size=14))+
  guides(color=guide_legend(override.aes=list(fill=NA)))
p
p1<-ggplot(data=distances, aes(x=V3, y=V5)) + 
  geom_point(aes(color=(distances$V6)))+
  geom_smooth(method='lm')+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)+
  xlab("Distance")+
  ylab("Fst/(1-Fst)")+
  theme(legend.title=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  theme(legend.text=element_text(size=14))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme(legend.position="none")
p1

a<-ggplot(data=Normal, aes(x=V3, y=V5)) + 
  geom_point(color="#F8766D")+
  xlab("Distance")+
  ylab("Fst/(1-Fst)")+
  geom_smooth(method='lm')+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)+
  theme(legend.title=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  theme(legend.text=element_text(size=14))+
  guides(color=guide_legend(override.aes=list(fill=NA)))
a

b<-ggplot(data=Skin, aes(x=V3, y=V5)) + 
  geom_point(color="#619CFF")+
  geom_smooth(method='lm')+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)+
  xlab("Distance")+
  ylab("Fst/(1-Fst)")+
  theme(legend.title=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  theme(legend.text=element_text(size=14))+
  guides(color=guide_legend(override.aes=list(fill=NA)))
b

c<-ggplot(data=SkinNormal, aes(x=V3, y=V5)) + 
  geom_point( color= "#00BA38")+
  geom_smooth(method='lm')+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)+
  xlab("Distance")+
  ylab("Fst/(1-Fst)")+
  theme(legend.title=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  theme(legend.text=element_text(size=14))+
  guides(color=guide_legend(override.aes=list(fill=NA)))
c


q<-plot_grid(p1, a, b, c, labels = c("A", "B", "C", "D"), align = "vh",hjust = -1, nrow = 2)
q
legend<-get_legend(p)


r <- plot_grid( q,legend, rel_widths = c(2, .5))
r


#sylvaticus
setwd("~/Documents/thesis/")
distances<-read.csv("genetic_diversity_sylv_R3.csv", sep = ",", header=FALSE)

my.formula <- y ~ x
sylv<-ggplot(data=distances, aes(x=V9, y=V10)) + 
  geom_point(color="#F8766D")+
  geom_smooth(method='lm')+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)+
  xlab("Distance")+
  ylab("Fst/(1-Fst)")+
  theme(legend.title=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  theme(legend.text=element_text(size=14))+
  guides(color=guide_legend(override.aes=list(fill=NA)))
sylv
