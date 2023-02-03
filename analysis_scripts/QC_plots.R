par(mfrow=c(1,2))

GQfilter=read.csv("empirical_data/GQFilter_missingness.csv",header=T)
Can40filter=read.csv("empirical_data/Can40Filter_missingness.csv",header=T)

GQfilter$col=GQfilter$pop
GQfilter$col=gsub("Can40","#F3C127",GQfilter$col)
GQfilter$col=gsub("Lof07","#A62303",GQfilter$col)
GQfilter$col=gsub("Lof11","#D42D05",GQfilter$col)
GQfilter$col=gsub("Lof14","#FB8901",GQfilter$col)
GQfilter$col=gsub("Can13","#A2801A",GQfilter$col)

barplot(GQfilter$fmiss,col=GQfilter$col,ylab="Proportion Missing",ylim=c(0,1))
axis(1,labels=c("Can40","Lof07","Lof11","Lof14","Can13"),at=c(12,38,68,95,120))
#text(x=c(90,90),y=c(0.95,0.8),labels=c("All Loci","GQ Filter = 30, nloci = 346,290"))
text(x=1.5,y=0.95,labels="(a)",cex=2)

Can40filter$col=Can40filter$pop
Can40filter$col=gsub("Can40","red",Can40filter$col)
Can40filter$col=gsub("Lof07","orange",Can40filter$col)
Can40filter$col=gsub("Lof11","blue",Can40filter$col)
Can40filter$col=gsub("Lof14","purple",Can40filter$col)
Can40filter$col=gsub("Can13","gold",Can40filter$col)

barplot(Can40filter$fmiss,col=GQfilter$col,ylab="Proportion Missing",ylim=c(0,1))
axis(1,labels=c("Can40","Lof07","Lof11","Lof14","Can13"),at=c(12,38,68,95,120))
#text(x=c(90,90),y=c(0.95,0.8),labels=c("Can40 Low Missingness Loci","GQ Filter = 30, nloci = 112,082"))
text(x=1.5,y=0.95,labels="(b)",cex=2)
