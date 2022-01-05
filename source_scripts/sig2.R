library(tidyverse)
library(data.table)

sig<-fread("sig_TAM.csv",header = T,data.table = F)
sig[sig==""]<-NA
sig.list<-list()
for (i in 6:7) {
  name<-colnames(sig)[i]
  sig.list[[name]]<-na.omit(as.character(sig[1:99,i]))
}

data_expr<-CGGA_exp1
library(GSVA)
es<-gsva(as.matrix(data_expr),
         sig.list,method="ssgsea",
         abs.ranking=T,
         ssgsea.norm = TRUE)
es<-as.data.frame(es)

x=as.numeric(data_expr["SLC11A1",])
data_all<-data.frame(gs=NULL,
                     cor=NULL,
                     p=NULL)

y<-as.numeric(es[1,]-es[2,])
corres<-cor.test(x,y,method = "pearson")
data_tmp<-data.frame(cor=corres$estimate,
                     p=corres$p.value,
                     dataset="CGGA")
data_all<-cbind(data_all,data_tmp)
