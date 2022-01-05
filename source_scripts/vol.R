library(IOBR)
signature_collection
DEG$adj.P.Val
DEG$logFC
alldiff=DEG
alldiff$type <- ifelse(alldiff$adj.P.Val>0.05,'No-Sig',
                       ifelse(alldiff$logFC>=1.5,'Up',
                              ifelse(alldiff$logFC<(-1.5),
                                     'Down','No-Sig')))
table(alldiff$type)

ggplot(alldiff,aes(logFC,
                   -log10(adj.P.Val),
                   fill=type))+
  geom_point(shape=21,aes(size=-log10(adj.P.Val)))+
  scale_fill_manual(values=c('#0071B5','gray','#BB3B29'))+
  scale_color_manual(values=c('gray60','black'))+
  geom_vline(xintercept=c(-1.5,1.5),lty=2,col="gray30",lwd=0.6) +
  geom_hline(yintercept = -log10(0.05),lty=2,col="gray30",lwd=0.6)+
  theme_bw(base_rect_size = 1)+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.grid = element_blank(),
        plot.title = element_text(family = 'regular',hjust = 0.5),
        legend.position = c(0.9, 1),
        legend.justification = c(0.5, 1),
        legend.key.height = unit(0.5,'cm'),
        legend.background = element_rect(fill = NULL, colour = "black",
                                         size = 0.5))+
  xlim(-4,4)+
  guides(size=F,color=F)+
  ylab('-log10 (adj.P.Val)')+xlab('log2 (Fold Change)')
