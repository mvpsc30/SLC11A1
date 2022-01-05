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
library(ggplot2)
library(ggplotify)
plot=as.ggplot(forest(m1,
                      col.study = "#0071B5",
                      col.square = "#0071B5",
                      col.inside = "white",
                      col.diamond = "#BB3B29",
                      col.diamond.fixed	= "#BB3B29",
                      col.diamond.random= "#BB3B29",
                      col.diamond.lines = "#BB3B29",
                      col.predict = "#BB3B29",
                      col.inside.random="#BB3B29",
                      col.fixed = "#BB3B29",
                      col.random = "#BB3B29",
                      lab.e = "SLC11A1  High",
                      lab.c = "Low",
                      layout = "RevMan5"
              ))
dev.off()
