cog<-read.delim("expression4.tsv")
View(cog)
row.names(cog)<-cog$Id
cog<-cog[,-1]
View(cog)
dff <- as.matrix(cog)
View(dff)
library(gplots)
library(superheat)
library(pheatmap)
library(RColorBrewer)
Colors=brewer.pal(11,"Spectral")
#Colors=brewer.pal(9,"YlGnBu")
#Colors=c("blue","yellow")
#Colors=colorRampPalette(Colors)(200)
pheatmap((dff), color = Colors, show_rownames = 7, show_colnames = 7, display_numbers =T, number_format = "%.2f",fontsize_number = 4, fontsize_row = 4, fontsize_col = 6, angle_col= 315, cluster_rows = FALSE, cluster_cols = FALSE)

superheat(X=dff, X.text = round(as.matrix(dff), 1),  X.text.size = 2.0, heat.pal = c("#b35806", "white", "#542788"),bottom.label.text.angle = 90,left.label.text.size = 2.5, bottom.label.text.size = 3, left.label.col = "white", bottom.label.col = "white",legend.vspace = 0.0)

#superheat(X=dff, X.text = round(as.matrix(dff), 1),  X.text.size = 2.1, heat.pal = c("#b35806", "white", "#542788"),bottom.label.text.angle = 90,left.label.text.size = 2.5, bottom.label.text.size = 3, left.label.col = "white", bottom.label.col = "white",legend.vspace = 0.0, row.dendrogram = TRUE, col.dendrogram = TRUE)
#heatmap.2((dff),cexCol =0.7,cexRow=1.0,col=Colors,trace = "none",margins = c(5,4),density.info = "none",scale = "none",lhei = c(2,8),lwid=c(2,8),srtCol=90,key = TRUE,key.xlab="COG copy number",key.title =NA)
#heatmap.2(t(dff),cexCol = 1.0,cexRow=1.2,col = heat.colors(256),trace = "none",margins = c(10,12),density.info = "none",scale = "none",lmat=rbind( c(0, 3, 4), c(2,1,1 ) ), lwid=c(1.5, 4, 2.6 ),srtCol=90,key = TRUE,key.xlab="Gene copy number",key.title =NA)
#superheat(X=dff, X.text = round(as.matrix(dff), 1),  X.text.size = 3, heat.pal = c("#b35806", "white", "#542788"), heat.pal.values = c(0, 0.5, 1))
