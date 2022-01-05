library(ggcor)
library(ggplot2)
library(GSVA)

gene_set=fread("signature annotation.txt",header = T)
gene_set$step=paste0("Step",gene_set$Steps,"_",
                     gene_set$ImmuneCellType)

head(gene_set)
list<- split(as.matrix(gene_set)[,1], gene_set[,5])
gsva_matrix<- gsva(as.matrix(CGGA_exp1), 
                   list,
                   method='ssgsea',
                   kcdf='Gaussian',
                   abs.ranking=TRUE)
gsva_matrix=as.data.frame(t(gsva_matrix))
gsva_matrix$ID=rownames(gsva_matrix)
list2=list("Cancer Immunity Cycle"=unique(gene_set$step))

sig_res<-iobr_cor_plot(pdata_group       = CGGA_phe,
                   id1                   = "Sample",
                   feature_data          = gsva_matrix,
                   id2                   = "ID",
                   target                = NULL,
                   group                 = "group",
                   is_target_continuous  = FALSE,
                   padj_cutoff           = 1,
                   index                 = 1,
                   category              = "signature",
                   signature_group       = list2[c(1)],
                   ProjectID             = "CGGA",
                   palette_box           = "nrc",
                   palette_corplot       = "pheatmap",
                   palette_heatmap       = 4,
                   feature_limit         = 100,
                   character_limit       = 30,
                   show_heatmap_col_name = FALSE,
                   show_col              = FALSE,
                   show_plot             = TRUE,
                   path                  = "signatures")

gsva_matrix<- gsva(as.matrix(CGGA_exp1), 
                   list,
                   method='ssgsea',
                   kcdf='Gaussian',
                   abs.ranking=TRUE)

gsva_matrix=as.data.frame(t(gsva_matrix))
gsva_matrix$SLC11A1=as.numeric(CGGA_exp1["SLC11A1",])


y <- as.numeric(gsva_matrix[,"SLC11A1"])
colnames <- colnames(gsva_matrix)
cor_data_df <- data.frame(colnames)
for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(gsva_matrix[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("env","r","p.value")
cor_data_df$spec="SLC11A1"

mantel2 <- cor_data_df %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p.value, breaks = c(-Inf,0.01, 0.05, Inf),
                  labels = c("< 0.001",
                             "0.01 - 0.05",
                             ">= 0.05")))

mantel2=mantel2[,c("spec","env","r","p.value","rd","pd" )]

mantel2$line=ifelse(mantel2$r>0,
                    "Positive",
                    "Negative")

mantel2$line=factor(mantel2$line,levels = c("Positive",
                                            "Negative"))

quickcor(gsva_matrix, type = "upper") +
  geom_square() +
  anno_link(aes(colour = pd, size = rd,linetype= line), data = mantel2) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_fill_gradient2(low = "#0071B5",high = "#BB3B29")+
  scale_colour_manual(values = c("pink","lightblue","gray")) +
  guides(size = guide_legend(title = "Pearson's r",
                             override.aes = list(colour = "grey10"), 
                             order = 2),
         colour = guide_legend(title = "Pearson's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r",order = 4),
         
         linetype=guide_legend(title = "Relations",
                               override.aes = list(colour = "grey10"), 
                               order = 3)
  )

