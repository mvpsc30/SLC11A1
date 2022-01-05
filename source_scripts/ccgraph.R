library(ggraph)
library(tidygraph)

signature_collection$Immune_Checkpoint
EMT=unique(c(signature_collection$EMT1,signature_collection$EMT2,signature_collection$EMT3))
EMT=list(EMT=EMT)

table(immunemodulator$V2)
Immunoinhibitor=immunemodulator[which(immunemodulator$V2=="Immunoinhibitor"),]$V3


batch_cor <- function(gene){
  y <- as.numeric(CGGA_exp1[gene,])
  rownames <- intersect(signature_collection$Immune_Checkpoint,rownames(CGGA_exp1))
  do.call(rbind,lapply(rownames, function(x){
    dd  <- cor.test(as.numeric(CGGA_exp1[x,]),y,type="pearson")
    data.frame(term="Immune_Checkpoint",genes=x,correlation=dd$estimate)
  }))
}

EMT_cor=batch_cor("SLC11A1")

Immunoinhibitor_cor=batch_cor("SLC11A1")

Immune_Checkpoint_cor=batch_cor("SLC11A1")

df=rbind(EMT_cor,Immunoinhibitor_cor,Immune_Checkpoint_cor)

source(file = "gather_graph_node.R")
source(file = "gather_graph_edge.R")


nodes <- gather_graph_node(df, index = c("term", "genes"), value = "correlation", root="all")
edges <- gather_graph_edge(df, index = c("term", "genes"), root = "all")
nodes <- nodes %>% mutate_at(c("node.level","node.branch"),as.character)
head(nodes, 10)
head(edges, 10)

geneSpecial=df

geneCol <- ifelse(geneSpecial$correlation>0.4,"a","normal")
names(geneCol) <- geneSpecial$genes
geneCol



nodes$color <- "normal"
nodes[nodes$node.short_name %in% geneSpecial$genes,]$color <- geneCol[nodes[nodes$node.short_name %in% geneSpecial$genes,]$node.short_name]
nodes[nodes$node.short_name %in% geneSpecial$genes,]
nodes$color <- factor(nodes$color, levels = unique(nodes$color))


graph <- tbl_graph(nodes, edges)

gc <- ggraph(graph, layout = 'dendrogram', circular = TRUE) + 
 
  geom_edge_diagonal(aes(color = node1.node.branch,
                         filter=node1.node.level!="all"), 
                     alpha = 1/3,edge_width=1) + 
  geom_node_point(aes(size = node.size, 
                      color = node.branch,
                      filter=node.level!="all"), alpha = 1/3) + 
  scale_size(range = c(0.5,40)) + 
  theme(legend.position = "none") + 
  
 
  scale_edge_color_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  
  geom_node_text(
    aes(
      x = 1.048 * x,
      y = 1.048 * y,
      label = node.short_name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = leaf,
      color = node.branch
    ),
    size = 6, hjust = 'outward') +

  geom_node_text(
    aes(label=node.short_name,
        filter = !leaf & (node.level != "all"),
        color = node.branch),
    fontface="bold",
    size=6,
    family="sans"
  ) + 
  theme(panel.background = element_rect(fill = NA)) +
  coord_cartesian(xlim=c(-1.3,1.3),ylim = c(-1.3,1.3)) #扩大坐标系

gc


















