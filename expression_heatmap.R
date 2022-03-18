cog<-read.delim("expression.tsv")
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

