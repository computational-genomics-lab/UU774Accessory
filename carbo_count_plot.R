car<-read.delim("hetero_count_boxh.tsv")
View(car)
carr<-stack(car)
View(carr)
library(ggplot2)
g <- ggplot(carr, aes(values, ind))
g + geom_count(aes(x=ind, y=values),col="tomato3") + theme_classic()+theme(axis.text.x  = element_text( angle=45,face="bold",hjust = 1),axis.text.y  = element_text(face="bold")) + ylab("")  + xlab("")+ ggtitle("")
ggsave("carbohydrate_count.tiff",units= "mm", width=170, height= 165, dpi=300, compression = 'lzw')

##################################################################################################################################################
# boxplot
carr[carr== 0] <- 0.1
carr<-na.omit(carr)
carbo1<-log2(carr$values)
#carbo1<-(carbo$Count)
carbo2 <-cbind(carr,carbo1)
View(carbo2)
ggplot(carr, aes(x=ind, y=values,fill=ind)) + geom_boxplot() + theme_classic() + theme(axis.text.x  = element_text( angle=60,face="bold",hjust = 1),axis.text.y  = element_text(face="bold")) + ylab("") + xlab("") + ggtitle("")

ggsave("demo.tiff",units= "mm", width=170, height= 140, dpi=300, compression = 'lzw')
