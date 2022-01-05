DEG=data.table::fread("GlioVis - Visualization Tools for Glioma Datasets.csv",header = T)
library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
## 根据logfc降序排列基因
load("pathways_for_gsva_mets.Rdata")
list=org.Hs.egMsigdbC2REACTOME$gs
alldiff <- DEG
alldiff2 <- bitr(alldiff$V1,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = "org.Hs.eg.db",drop = T)
alldiff2=merge(alldiff2,alldiff,by.x="SYMBOL",by.y="V1")
## fgsea中输入的关键基因信息
alldiff2 <- alldiff2[order(alldiff2$logFC,decreasing = T),]
id <- alldiff2$logFC
names(id) <- alldiff2$ENTREZID
id

alldiff <- DEG[which(DEG$logFC>-1),]
alldiff <- alldiff[order(alldiff$logFC,decreasing = T),]
id <- alldiff$logFC
names(id) <- alldiff$V1
id

## fgsea中输入的关键通路信息
gmtfile <- "h.all.v7.4.symbols.gmt"
hallmark <- read.gmt(gmtfile)
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)

## Perform the fgsea analysis
fgseaRes <- fgsea(pathways = hallmark.list, 
                  stats = id,
                  minSize=1,
                  maxSize=10000,
                  nperm=1000)
sig <- fgseaRes[c(27,14,45,26,25,22),]
sig <- sig[order(sig$NES,decreasing = T),]
up=sig$pathway
sig <- sig[order(sig$NES,decreasing = F),]
down=sig[1:10]$pathway

dev.off()
plotGseaTable(hallmark.list[c(up)],id,fgseaRes,gseaParam = 0.3)
signature_collection
shiny::runGitHub("BioinformaticsFMRP/PanCanStem_Web/",subdir = "PanCanStem_Web")
SLC11A1_gsva=gsva(as.matrix(CGGA_exp1),
                           Pathways,
                           method='ssgsea',
                           kcdf='Gaussian',
                           abs.ranking=TRUE)
SLC11A1_gsva=as.data.frame(t(SLC11A1_gsva))
y=as.numeric(CGGA_exp$SLC11A1)
colnames <- colnames(SLC11A1_gsva)
data_s <- data.frame(colnames)
for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(SLC11A1_gsva[,i]),y,type="spearman")
  data_s[i,2] <- test$estimate                                            
  data_s[i,3] <- test$p.value
}
names(data_s) <- c("symbol","correlation","pvalue")
head(data_s)

data_s %>% 
  #filter(pvalue <0.05) %>% # 如果不想把p值大于0.05的放在图上，去掉最前面的#号
  ggplot(aes(correlation,
             forcats::fct_reorder(symbol,correlation))) +
  geom_segment(aes(xend=0,yend=symbol)) +
  geom_point(aes(col=pvalue,size=abs(correlation))) +
  scale_colour_gradientn(colours=c("#BC3C29","#0071B5")) +
  #scale_color_viridis_c(begin = 0.5, end = 1) +
  scale_size_continuous(range =c(2,8))+
  theme_bw() +
  theme(axis.title.y = element_text(angle = 90))+
  ylab(NULL)+
  xlab("Correlation")

library(ggpubr)
ggscatter(data = labels_tis, 
          x= "SLC11A1", 
          y = "V1",
          color = "#55A0CD",
          size = 1,
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "black",
                            fill = "#CAD9EC")
)+stat_cor(method = "pearson", label.x = 0,label.y = 1)+
  ylab("T-cell inflammatory signature")# Add correlation coefficient


