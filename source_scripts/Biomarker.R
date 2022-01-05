library(IOBR)
library(data.table)
immunemodulator=fread("immunomodulator.txt",header = F)
immunemodulator=immunemodulator[1:150,]
marker=immunemodulator$V3
marker=list(marker=marker)

immune_marker=fread("immune_marker.txt",header = T)

marker2=signature_collection[c("CD_8_T_effector",
                               "NK_cells_Danaher_et_al",
                               "Macrophages_Danaher_et_al",
                               "Th1_cells_Bindea_et_al",
                               "DC_Bindea_et_al")]
marker2=as.data.frame(unlist(marker2))
marker2$group=substring(rownames(marker2),1,8)
colnames(marker2)=c("gene","Effector_gene")

marker=immune_marker$ID
marker=list(marker=marker)


res<-my_iobr_cor_plot(pdata_group         = CGGA_phe,
                      id1                   = "Sample",
                      feature_data          = CGGA_exp1,
                      id2                   = "ID",
                      target                = NULL,
                      group                 = "group",
                      is_target_continuous  = FALSE,
                      padj_cutoff           = 1,
                      index                 = 1,
                      category              = "gene",
                      signature_group       = marker,    
                      ProjectID             = "CGGA",
                      palette_box           = "nrc",
                      palette_corplot       = "pheatmap",
                      palette_heatmap       = 4,
                      feature_limit         = 200,
                      character_limit       = 200,
                      show_heatmap_col_name = FALSE,
                      show_col              = FALSE,
                      show_plot             = TRUE,
                      path                  = "marker_immune")
length(unique(res$variables))
immunemodulator=immunemodulator[,c("V2","V3")]
res2=merge(res,immune_marker,by.x = "variables",by.y = "ID")
res3=res2%>%filter(variables !="CXCR4")
res3=res3%>%filter(variables !="CXCL12")
colnames(res2)
library(tidyHeatmap)
res2 %>%
  group_by(`TILs type`,target_group) %>%
  tidyHeatmap::heatmap(
    .column = ID,
    .row = variables,
    .value = value,
    show_column_names=F,
    palette_value = c("#0072B5","white","#D30D17"),
    palette_grouping = list(c("#D30D17","#0072B5",
                              "#7FC97F" ,"#BEAED4","#FDC086", 
                              "#FFFF99", "#386CB0",
                              "#F0027F" ,"#BF5B17" ,"#666666"))
  )%>%
  add_tile(`TILs type`)%>%
  tidyHeatmap::save_pdf("sz2.pdf",
                            width = 6,height =6)

RColorBrewer::brewer.pal(8, "Accent")
table(res2$variables,res2$V2)
