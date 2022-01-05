CGGA=readRDS(file = "CGGA.Rds")
CGGA_exp=CGGA$expr
CGGA_phe=CGGA$pData

CGGA_exp=CGGA_exp[,9:ncol(CGGA_exp)]

CGGA_exp1=as.data.frame(t(CGGA_exp))
CGGA_phe$group=ifelse(CGGA_exp$SLC11A1>median(CGGA_exp$SLC11A1),"High","Low")
library(IOBR)
library(dplyr)
library(IOBR)
library(EPIC)
library(estimate) 
library(MCPcounter)
library(tidyverse)
library(tidyHeatmap)
library(maftools)
library(ggpubr)
library(ggplot2)
library(survival)
library(data.table)
data<-deconvo_tme(eset = CGGA_exp1,
                  method = "mcpcounter")

CGGA_sub=CGGA_exp[,c("PDCD1","SLC11A1")]
CGGA_sub$Sample=rownames(CGGA_sub)
data=merge(CGGA_sub,data,by.x = "Sample",by.y = "ID")

colnames(data)
TumorPurity=ggscatter(data, 
          x= "SLC11A1", 
          y = "TumorPurity_estimate",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray")
          )+stat_cor(method = "pearson", label.x = 0,label.y = 1.05)  # Add correlation coefficient

library(patchwork)
ESTIMATES+Stromal+Immune+TumorPurity+plot_layout(ncol = 4)


library(survminer)
library(survival)
data.back=data
data$cell=ifelse(data$CD8_T_cells_MCPcounter>
                    median(data$CD8_T_cells_MCPcounter),
                  "High","Low")
data$SLC11A1=CGGA_exp$SLC11A1
data$SLC11A1=ifelse(data$SLC11A1>median(data$SLC11A1),"High","Low")
data$survival=CGGA_phe$survival
data$status=CGGA_phe$status

fit<- survfit(Surv(survival, status) ~ cell+SLC11A1,
              data = data)

# Basic survival curves
ggsurvplot(fit, 
           data = data,
           pval = T,
           legend = c(0.7,0.9),
           palette = "jco")

data_immune=read.csv("data_XX.csv",header = T)

data_immune=merge(data_immune,data,by.x = "X",by.y = "ID")

table(data_immune$Response,data_immune$group2)

data_immune$group=ifelse(data_immune$CD4_T>median(data_immune$CD4_T),
                  "High","Low")

library(plyr)
data_immune2=data_immune[,c("group2","Response")]
data_immune2$v=1
data_immune2=ddply(
                   data_immune2,
                   "group2",
                   transform,
                   percent=v/sum(v)*100
                   )


library(ggsci)
data_immune2$group2=factor(data_immune2$group2,levels = c("Low","High"))
ggplot(data_immune2,
       aes(group2,
           percent,
           fill=Response))+
  geom_bar(stat = "identity")+
  theme_minimal()+
  scale_fill_nejm()+
  xlab("SLC11A1")+
  ylab("Percent(%)")+
  coord_flip()
ggsave("barplot.pdf",height = 1.5,width = 10)

fit<- survfit(Surv(survival, status) ~ group+group2,
              data = data_immune)
ggsurvplot(fit, data = data_immune)
