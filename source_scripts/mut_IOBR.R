library(IOBR)
library(TCGAmutations)
library(maftools)
GBM=TCGAmutations::tcga_load(study = "GBM")
LGG=TCGAmutations::tcga_load(study = "LGG")
GBMLGG=merge_mafs(list(GBM,LGG))

mut=GBMLGG@data

names(mut_list)
#> [1] "all"        "snp"        "indel"      "frameshift"
# choose SNP mutation matrix as input
mut_TGCA<-mut_list$all
rownames(mut_TGCA)<-substring(rownames(mut_TGCA),1,12)
mut_TGCA<-mut_TGCA[!duplicated(rownames(mut_TGCA)),]





remotes::install_github("ropensci/UCSCXenaTools")

library("UCSCXenaTools")
var_gbm<-XenaGenerate(subset = XenaCohorts =="GDC TCGA Glioblastoma (GBM)") %>% 
  XenaFilter(filterDatasets    = "TCGA-GBM.mutect2_snv.tsv") %>% 
  XenaQuery() %>%
  XenaDownload() %>% 
  XenaPrepare()

var_lgg<-XenaGenerate(subset = XenaCohorts =="GDC TCGA Lower Grade Glioma (LGG)") %>% 
  XenaFilter(filterDatasets    = "TCGA-LGG.mutect2_snv.tsv") %>% 
  XenaQuery() %>%
  XenaDownload() %>% 
  XenaPrepare()

var_merge=rbind(var_gbm,var_lgg)
#> This will check url status, please be patient.
#> All downloaded files will under directory /tmp/Rtmp0RFeyC.
#> The 'trans_slash' option is FALSE, keep same directory structure as Xena.
#> Creating directories for datasets...
#> Downloading TCGA-STAD.mutect2_snv.tsv.gz
head(var_stad)
#> # A tibble: 6 x 11
#>   Sample_ID gene  chrom  start    end ref   alt   Amino_Acid_Chan… effect filter
#>   <chr>     <chr> <chr>  <dbl>  <dbl> <chr> <chr> <chr>            <chr>  <chr> 
#> 1 TCGA-CD-… C1or… chr1  2.19e6 2.19e6 G     -     p.P72Rfs*87      frame… PASS  
#> 2 TCGA-CD-… ERRF… chr1  8.01e6 8.01e6 C     T     p.P327P          synon… panel…
#> 3 TCGA-CD-… CLCN6 chr1  1.18e7 1.18e7 G     A     p.S486N          misse… PASS  
#> 4 TCGA-CD-… PRAM… chr1  1.28e7 1.28e7 G     A     p.G341R          misse… panel…
#> 5 TCGA-CD-… PRAM… chr1  1.32e7 1.32e7 G     T     p.P148H          misse… PASS  
#> 6 TCGA-CD-… CELA… chr1  1.55e7 1.55e7 G     A     p.G59R           misse… PASS  
#> # … with 1 more variable: dna_vaf <dbl>
# Then function `make_mut_matrix` can be used to transform data frame into mutation matrix
mut_list2<-make_mut_matrix(mut_data               = var_merge,
                           category               = "multi",
                           Tumor_Sample_Barcode   = "Sample_ID",
                           Hugo_Symbol            = "gene",
                           Variant_Classification = "effect",
                           Variant_Type           = "Variant_Type")

names(mut_list2)
#> [1] "all"        "snp"        "indel"      "frameshift"
# choose SNP mutation matrix as input
mut<-mut_list2$all
rownames(mut)<-substring(rownames(mut),1,12)
mut<-mut[!duplicated(rownames(mut)),]
dim(mut)

TCGA=readRDS(file="TCGA_GBMLGG.Rds")
TCGA=TCGA$expr
TCGA[1:3,1:4]
TCGA$Sample=gsub("\\.","-",TCGA$Sample)
TCGA=TCGA[which(TCGA$Sample %in% rownames(mut)),]
dim(TCGA)
TCGA=TCGA[,c("Sample","SLC11A1")]



res<-find_mutations(mutation_matrix     = mut, 
                    signature_matrix    = TCGA,
                    id_signature_matrix = "Sample",
                    signature           = "SLC11A1",
                    min_mut_freq        = 0.01,
                    plot                = TRUE,
                    method              = "wilcoxon",
                    save_path           = "mutations2",
                    palette             = "nrc",
                    show_plot           = T,
                    width               = 8, #oncoprint的宽度
                    height              = 4, #oncopritn的高度
                    oncoprint_group_by  = "mean", #使用平均值的分组进行分组
                    oncoprint_col       = "#224444",
                    gene_counts         = 10) #oncoprint展示多少个基因


dt <- GBMLGG@data


res2<-res$wilcoxon_test
dt2 <-dt[dt$Hugo_Symbol %in% rownames(res2)[1:10],]

dt2$ID<-substring(dt2$Tumor_Sample_Barcode,1,12)

Quantile<-function(x){
  ifelse(x>quantile(x,.75),"Q1",ifelse(x>quantile(x,.5),"Q2",ifelse(x>quantile(x,.25),"Q3","Q4")))
}#设置函数，四分位分组

data<-as.data.frame(TCGA[,2])
new_quantile<-apply(data ,2, Quantile)

TCGA$group<-new_quantile

dt2<-dt2[,c("Hugo_Symbol","ID")]

colnames(dt2)<-c("gene","ID")

colnames(TCGA)<-c("ID","SLC11A1","group")

dt2<-left_join(dt2,TCGA)

dt2<-na.omit(dt2)

x=table(dt2$gene,dt2$group)

write.csv(x,"mutfish.csv")

x<-read.csv("mutfish.csv",header = T,row.names = 1)


# library(devtools)
# install_github("chrisamiller/fishplot")

library(fishplot)

x<-x/10
frac.table <- as.matrix(x[1:4,])
timepoints <- c(0,75,110,120) #决定了鱼头、鱼身、鱼尾的长度
parents <- c(0,1,1,3) 

fish <- createFishObject(frac.table, parents, timepoints)

#设置每种clone的颜色
fish <- setCol(fish, c("#099D79", "#70C7EC", "#E8262D","#2C3789"))
fish <- layoutClones(fish)

#pdf('fish.pdf', width=10, height=6)
par(mar = par()$mar + c(0,0,3,0)) #在底部留出画图例的地方

fishPlot(fish, shape = "spline", #spline，圆滑的；或polygon，直的
         #title.btm = "Clonal architecture of tumors", #左下角可以写字
         cex.title = 1.2, #字号
         pad.left = 0.25, #鱼头的边的倾斜角度
         
         vlines = timepoints[c(1:2,4)], #画三条竖直实线
         col.vline = "white", #线的颜色
         vlab = c("PMF", "sAML", "sAML\nREM"), #竖线对应的文字，在画竖线的前提下才能写字
         
         bg.col = c("#F1F2F2","#F1F2F2","#F1F2F2"), #灰色背景
         border = 0.1 #每个clone的轮廓线的宽度
)

#如果三条竖线要画成虚线，就运行下面这行
#abline(v=timepoints[c(1:2,4)], col="white", lty=2, lwd=1)

#添加clone的图例
par(xpd = T)
legend("bottomright", 
       inset=c(.7,-.3), #把图例画到图外，根据自己的鱼调整
       pch=16, bty="n", 
       col=fish@col, text.col = fish@col,
       legend = paste0(row.names(frac)," ",frac$PMF,"%"))

legend("bottomright", 
       inset=c(.3,-.3),
       pch=16, bty="n", 
       col=fish@col, text.col = fish@col,
       legend = paste0(frac$sAML,"%"))

legend("bottomright", 
       inset=c(.0,-.3),
       pch=16, bty="n", 
       col=fish@col, text.col = fish@col,
       legend = paste0(frac$Saml_REM,"%"))

