library(survminer) 
library(survival)
mRNA_exprSet=read.csv("ITGB2_Vital_clinical.csv",
                      header = T)
mRNA_exprSet$group=ifelse(mRNA_exprSet$mRNA>median(mRNA_exprSet$mRNA),"High","Low")

fit=survfit(Surv(survival,status)~group,
              data = mRNA_exprSet)

p1=ggsurvplot(
             fit, 
             data = mRNA_exprSet,
             surv.median.line = "hv",
             conf.int = TRUE,
             pval = TRUE,
             #facet.by="Grade",
             palette = c("pink","lightblue"),
             title="GSE43289")
p1

table(mRNA_exprSet$status,mRNA_exprSet$group)

library(meta)
library(Matrix)

dat=read.csv("data_forest.csv",header = T)

m1 <- metabin(ev.exp, 
              n.exp, 
              ev.cont, 
              n.cont,
              data = dat, 
              studlab = Study)

summary(m1)

forest(m1,
       col.study = "lightblue",
       col.square = "lightblue",
       col.inside = "white",
       col.diamond = "pink",
       col.diamond.lines = "pink",
       col.predict = "pink",
       col.fixed = "pink",
       col.random = "pink",
       lab.e = "High",
       lab.c = "Low"
       )
  
  
