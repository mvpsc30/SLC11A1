library(data.table)
sig<-fread("sig_TAM.csv",header = T,data.table = F)
sig[sig==""]<-NA
sig.list<-list()
for (i in 6:7) {
  name<-colnames(sig)[i]
  sig.list[[name]]<-na.omit(as.character(sig[1:99,i]))
}

library(GSVA)
es<-gsva(as.matrix(CGGA_exp1),
         sig.list,method="ssgsea",
         abs.ranking=T,
         ssgsea.norm = TRUE)
es<-as.data.frame(es)

x=as.numeric(CGGA_exp1["SLC11A1",])
data_all<-data.frame(gs=NULL,
                     cor=NULL,
                     p=NULL)
for (i in 1:7) {
  y<-as.numeric(es[i,])
  corres<-cor.test(x,y,method = "pearson")
  data_tmp<-data.frame(gs=rownames(es)[i],
                       cor=corres$estimate,
                       p=corres$p.value)
  data_all<-rbind(data_all,data_tmp)
}

y<-as.numeric(es[1,]-es[2,])
corres<-cor.test(x,y,method = "pearson")
data_tmp<-data.frame(cor=corres$estimate,
                     p=corres$p.value)




