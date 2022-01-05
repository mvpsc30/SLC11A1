
#####################################################################################


y <- as.numeric(CGGA_exp[,"CD274"])
colnames <- colnames(CGGA_exp)
cor_data_df <- data.frame(colnames)
for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(CGGA_exp[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")
library(tidyverse)

data=cbind(CGGA_exp,CGGA_phe)


data=data.table::fread("TCGA.csv",header = T)

data=data %>% filter(Grade !="I")

data=na.omit(data)

summ_data <- data %>% 
  group_by(Grade) %>% 
  summarise(
    mean = mean(mRNA),
    sd = sd(mRNA),
    n = n()
  ) %>% 
  mutate(se = sd/sqrt(n),
         Grade = factor(Grade, levels = c('II', 'III', 'IV')))

data_plot <- data %>% 
  mutate(Grade = factor(Grade, levels = c('II', 'III', 'IV')))

# 绘图
library(gghalves)
library(ggpubr)
library(ggsignif)
library(ggsci)
p2=ggplot(data_plot, aes(x = Grade, y = mRNA, fill = Grade))+
  geom_half_violin(aes(fill = Grade),
                   position = position_nudge(x = .15, y = 0),
                   adjust=1.5, 
                   trim=T, 
                   colour=NA, 
                   side = 'r') +
  geom_point(aes(x = as.numeric(Grade) - 0.1,
                 y = mRNA,
                 color = Grade),
             position = position_jitter(width = .1),
             size = 1, shape = 20) +
  geom_boxplot(aes(x = Grade,y = mRNA, fill = Grade),
               outlier.shape = NA,
               width = .05,
               color = "black")+
  geom_point(data=summ_data,
             aes(x=Grade,y = mean,group = Grade, color = Grade),
             shape=18,
             size = 1.5,
             alpha=0.7,
             position = position_nudge(x = .1,y = 0)) +
  geom_errorbar(data =summ_data,
                aes(x = Grade, 
                    y = mean, 
                    group = Grade, 
                    colour = Grade,
                    ymin = mean-sd, 
                    ymax = mean+sd),
                width=.05,
                position=position_nudge(x = .1, y = 0)) +
  scale_colour_manual(values=c("#86BDFE", "#FFBC7E","#BC3C29")) +
  scale_fill_manual(values=c("#86BDFE", "#FFBC7E","#BC3C29")) +
  geom_signif(comparisons = list(c("II", "III"),
                                 c("II", "IV"),
                                 c("III", "IV")),
              y_position = c(12.5,13.5,12.9),
              map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05))+
  theme_classic()+
  ylab("SLC11A1")

p1+p2+p3+p4+p5+p6+
  plot_annotation(tag_levels = "A")

ggsave("plot2.pdf",height = 8.5,width = 12)


#####################################################################################

data_immune=data.table::fread("vital_sur.csv",header = T)
data_immune$cutoff_group

fit<- survfit(Surv(survival, status) ~ cutoff_group,
              data = data_immune)

ggsurvF <- ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = data_immune,            # data used to fit survival curves.
  risk.table = F,       # show risk table.
  pval = TRUE, 
  pval.method = F,# show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimates of survival curves.
  palette = "nejm",
  title="GSE43289",
  xlab = "Time in months",   # customize X axis label.
  break.time.by = 50,     # break X axis in time intervals by 500.
  ggtheme = theme_bw(), # customize plot and risk table with a theme.
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = c("High", "Low"),    # change legend labels.
)
library(patchwork)
library(cowplot)

(ggsurvA$plot)+(ggsurvB$plot)+(ggsurvC$plot)+
(ggsurvE$plot)+(ggsurvD$plot)+(ggsurvF$plot)+
  plot_layout(ncol = 3)+
  plot_annotation(tag_levels = 'A')

ggsave("plot.pdf",height = 9.5,width = 12)

#############################################################################################
library(Hmisc)
library(dplyr)
library(data.table)
immunemodulator=fread("immunomodulator.txt",header = F)
immunemodulator=immunemodulator[1:150,]

data=data[,c("SLC11A1",intersect(immunemodulator$V3,colnames(data)))]
library(ComplexHeatmap)
data$SLC11A1=ifelse(data$SLC11A1>median(data$SLC11A1),"High","Low")
data_high=data %>% filter(SLC11A1=="High")
data_low=data %>% filter(SLC11A1=="Low")

data=data[order(data$SLC11A1,decreasing = T),]
Heatmap(as.matrix(t(data[,2:ncol(data)])),
        show_row_dend = F,
        show_column_dend = F,
        show_row_names = T,
        show_column_names = F,
        cluster_columns = F)
library(pheatmap)
pheatmap(as.matrix(t(data[,2:ncol(data)])),
         cluster_cols = F,
         cluster_rows = F,
         show_colnames = F,
         scale = "none")
list=immunemodulator$V3[which(immunemodulator$V2=="MHC")]

data2=mRNA_exprSet[,c("type",immunemodulator$V3,"SLC11A1")]

pancor <- function(gene1,gene2,data2){
  data1 <- split(data2,data2$type)
  do.call(rbind,lapply(data1, function(x){
    dd  <- cor.test(as.numeric(x[,gene1]),as.numeric(x[,gene2]),type="pearson")
    data.frame(type=x$type[1],
               cor=dd$estimate,
               p.value=dd$p.value,
               gene=gene2)
  }))
}

plotdf=data.frame(type=NULL,
                  cor=NULL,
                  p.value=NULL,
                  gene=NULL)


for (i in list) {
  data3=pancor("SLC11A1",i,data2)
  plotdf=rbind(plotdf,data3)
}

colnames(plotdf)

library(ggplot2)
plotdf$pvalue=ifelse(plotdf$p.value>=0.05,"p≥0.05","p<0.05")

ggplot(plotdf,
       aes(type,
           gene,
           shape = pvalue,
           color = cor))+ 
  geom_point(size=3)+
  scale_shape_manual(values = c(15,7))+
  scale_color_gradient2(low = "#2b8cbe",
                        mid = "white",
                        high = "#e41a1c",limits = c(-1, 1))+
  theme_bw()+
  theme(axis.title.y= element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle =90,
                                   hjust = 0,
                                   vjust = 0
        ))

##############################################################################################

DEG=data.table::fread("DEG.csv",header = T)
library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

alldiff <- DEG
alldiff2 <- bitr(alldiff$V1,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = "org.Hs.eg.db",drop = T)
alldiff2=merge(alldiff2,alldiff,by.x="SYMBOL",by.y="V1")

alldiff2 <- alldiff2[order(alldiff2$logFC,decreasing = T),]
id <- alldiff2$logFC
names(id) <- alldiff2$ENTREZID
id

gmtfile <- "KEGG_metabolism.gmt"
hallmark <- read.gmt(gmtfile)
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)
list=org.Hs.egMsigdbC2REACTOME$gs
## Perform the fgsea analysis
fgseaRes <- fgsea(pathways = list, 
                  stats = id,
                  minSize=1,
                  maxSize=10000,
                  nperm=100)
sig <- fgseaRes[fgseaRes$pval<0.05,]
sig <- sig[order(sig$NES,decreasing = T),]
up=sig[1:10]$pathway
sig <- sig[order(sig$NES,decreasing = F),]
down=sig[1:10]$pathway

dev.off()
plotGseaTable(list[c(up,down)],id,fgseaRes,gseaParam = 0.5)

