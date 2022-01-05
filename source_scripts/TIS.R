TIS <- c("CCL5", "CD27", "CD274", "CD276", 
         "CD8A", "CMKLR1", "CXCL9", "HLA-DQA1", 
         "HLA-DRB1","HLA-E", "IDO1", "LAG3", 
         "NKG7", "PDCD1LG2","PSMB10", "STAT1", "TIGIT")
library(GSVA)
library(ggpubr)
GSVA <- gsva(expr=as.matrix(CGGA_exp1), 
             gset.idx.list=list(TIS), 
             method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE
             )

GSVA <- as.data.frame(t(GSVA))

CGGA_exp[1:2,1:3]
labels_tis <- merge(CGGA_exp[,c("PDCD1","SLC11A1")],GSVA, by.x="row.names", by.y="row.names")

labels_tis$group=ifelse(labels_tis$mRNA>median(labels_tis$mRNA),"High","Low")
labels_tis$group <- factor(labels_tis$group, levels=c("High", "Low"))

plotTIS <- ggplot(labels_tis, aes(x=group, y=V1)) +
  geom_boxplot(color="black", outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.25), size=0.7, 
              aes(x=group, y=V1, color =group)) +
  labs( y = "TIS Score", x= "") +
  theme_classic() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values=c( "pink", "lightblue")) +
  geom_hline(yintercept=0.075, linetype="dashed", 
             color = "red", size=0.5) +
  theme(plot.title = element_text(hjust = 0.5, size=8, face="bold"), panel.spacing = unit(1, "lines"),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10), 
        legend.text=element_text(size=10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "white", fill=NA, size=0.3),
        strip.background = element_rect(colour = "white", fill = "white")) +
  stat_compare_means(method = "kruskal.test", label="p.format", size=3)+
  theme(legend.position = "top")
plotTIS