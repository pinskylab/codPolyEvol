library(tidyverse)

Can40_afreqs=read.table("data/Can40Loci_afreqs/plink2.Can40.afreq",header=T)
Can13_afreqs=read.table("data/Can40Loci_afreqs/plink2.Can13.afreq",header=T)

Can40_afreqs$ContFreq=Can13_afreqs$ALT_FREQS
Can40_afreqs$FreqChange=Can40_afreqs$ContFreq-Can40_afreqs$ALT_FREQS

Can40_afreqs$AbsFreqChange=abs(Can40_afreqs$ContFreq-Can40_afreqs$ALT_FREQS)

Lof07_afreqs=read.table("data/Can40Loci_afreqs/plink2.Lof07.afreq",header=T)
Lof11_afreqs=read.table("data/Can40Loci_afreqs/plink2.Lof11.afreq",header=T)
Lof14_afreqs=read.table("data/Can40Loci_afreqs/plink2.Lof14.afreq",header=T)


Lof07_afreqs$Freq11=Lof11_afreqs$ALT_FREQS
Lof07_afreqs$Freq14=Lof14_afreqs$ALT_FREQS

Lof07_afreqs$FreqChange11=Lof07_afreqs$Freq11-Lof07_afreqs$ALT_FREQS
Lof07_afreqs$FreqChange14=Lof07_afreqs$Freq14-Lof07_afreqs$ALT_FREQS

### change covariance - Canada vs Norways!

Can40_afreqs$FreqChange
Lof07_afreqs$FreqChange11
Lof07_afreqs$FreqChange14

covdf=data.frame(Can40_afreqs$CHROM,Can40_afreqs$FreqChange,Lof07_afreqs$FreqChange11,Lof07_afreqs$FreqChange14)

covdf$FreqChange1114=Lof07_afreqs$Freq14-Lof07_afreqs$Freq11

covdf=na.omit(covdf)

cor_Can40_Lof11=cov(covdf$Can40_afreqs.FreqChange,covdf$Lof07_afreqs.FreqChange11)/sqrt(var(covdf$Can40_afreqs.FreqChange)*var(covdf$Lof07_afreqs.FreqChange11))
cor_Can40_Lof14=cov(covdf$Can40_afreqs.FreqChange,covdf$Lof07_afreqs.FreqChange14)/sqrt(var(covdf$Can40_afreqs.FreqChange)*var(covdf$Lof07_afreqs.FreqChange14))

cor_all=(cov(covdf$Can40_afreqs.FreqChange,covdf$Lof07_afreqs.FreqChange11)+cov(covdf$Can40_afreqs.FreqChange,covdf$Lof07_afreqs.FreqChange14))/(sqrt(var(covdf$Can40_afreqs.FreqChange)*var(covdf$Lof07_afreqs.FreqChange11))+sqrt(var(covdf$Can40_afreqs.FreqChange)*var(covdf$Lof07_afreqs.FreqChange14)))
cor_all=(cov(covdf$Can40_afreqs.FreqChange,covdf$Lof07_afreqs.FreqChange11)+cov(covdf$Can40_afreqs.FreqChange,covdf$Lof07_afreqs.FreqChange14))/(sqrt(var(covdf$Can40_afreqs.FreqChange)*var(covdf$Lof07_afreqs.FreqChange11))+sqrt(var(covdf$Can40_afreqs.FreqChange)*var(covdf$Lof07_afreqs.FreqChange14)))


### by LG
LGs=unique(covdf$Can40_afreqs.CHROM)

cordf=data.frame(matrix(ncol=3,nrow=23))
colnames(cordf)=c("LG","cor11","cor14")

for (i in 1:(length(LGs)-1)) {
  c=LGs[i]
  cor_4011_LG=covdf[which(covdf$Can40_afreqs.CHROM==c),]
  cordf$LG[i]=c
  cordf$cor11[i]=cov(cor_4011_LG$Can40_afreqs.FreqChange,cor_4011_LG$Lof07_afreqs.FreqChange11)/sqrt(var(cor_4011_LG$Can40_afreqs.FreqChange)*var(cor_4011_LG$Lof07_afreqs.FreqChange11))
  cordf$cor14[i]=cov(cor_4011_LG$Can40_afreqs.FreqChange,cor_4011_LG$Lof07_afreqs.FreqChange14)/sqrt(var(cor_4011_LG$Can40_afreqs.FreqChange)*var(cor_4011_LG$Lof07_afreqs.FreqChange14))
}

plot(x=NULL,y=NULL,xlim=c(0,24),ylim=c(-0.15,0.15),xlab="Linkage Group",ylab="Convergence Correlation")
points(x=seq(1,23)-0.1,y=cordf$cor11,col="blue",pch=19)
points(x=seq(1,23)+0.1,y=cordf$cor14,col="green",pch=19)
abline(v=seq(0,23)+0.5,lty=2,lwd=0.5)
abline(h=0,lty=2)
legend(x=8,y=-0.05,legend=c("Canada-Norway 2011","Canada-Norway 2014"),col=c("blue","green"),pch=c(19,19))

### change covariance - Norway07-11 vs Norway01-14!

cor_Lof11_Lof14=cov(covdf$Lof07_afreqs.FreqChange11,covdf$Lof07_afreqs.FreqChange14)/sqrt(var(covdf$Lof07_afreqs.FreqChange11)*var(covdf$Lof07_afreqs.FreqChange14))

LGs=unique(covdf$Can40_afreqs.CHROM)

cordf1114=data.frame(matrix(ncol=3,nrow=23))
colnames(cordf1114)=c("LG","cor1114")

for (i in 1:(length(LGs)-1)) {
  c=LGs[i]
  cor_4011_LG=covdf[which(covdf$Can40_afreqs.CHROM==c),]
  cordf1114$LG[i]=c
  cordf1114$cor1114[i]=cov(cor_4011_LG$Lof07_afreqs.FreqChange11,cor_4011_LG$Lof07_afreqs.FreqChange14)/sqrt(var(cor_4011_LG$Lof07_afreqs.FreqChange11)*var(cor_4011_LG$Lof07_afreqs.FreqChange14))
}

plot(x=NULL,y=NULL,xlim=c(0,24),ylim=c(-0.65,0.65),xlab="Linkage Group",ylab="Convergence Correlation")
points(x=seq(1,23)-0.1,y=cordf1114$cor1114,col="purple",pch=19)
abline(v=seq(0,23)+0.5,lty=2,lwd=0.5)
abline(h=0,lty=2)
legend(x=8,y=-0.05,legend=c("Norway 2011 - Norway 2014"),col=c("purple"),pch=c(19,19))

### change covariance - Norway11-14 vs Canada!

Lof07_afreqs$FreqChange1114=Lof07_afreqs$Freq14-Lof07_afreqs$Freq11

covdf1114=data.frame(Can40_afreqs$CHROM,Can40_afreqs$FreqChange,Lof07_afreqs$FreqChange1114)

covdf1114=na.omit(covdf1114)

cor_Lof1114_Can=cov(covdf1114$Lof07_afreqs.FreqChange1114,covdf1114$Can40_afreqs.FreqChange)/sqrt(var(covdf1114$Lof07_afreqs.FreqChange1114)*var(covdf1114$Can40_afreqs.FreqChange))

LGs=unique(covdf1114$Can40_afreqs.CHROM)

cordfCan1114=data.frame(matrix(ncol=3,nrow=23))
colnames(cordfCan1114)=c("LG","cor1114")

for (i in 1:(length(LGs)-1)) {
  c=LGs[i]
  cor_4011_LG=covdf1114[which(covdf1114$Can40_afreqs.CHROM==c),]
  cordfCan1114$LG[i]=c
  cordfCan1114$cor1114[i]=cov(cor_4011_LG$Lof07_afreqs.FreqChange1114,cor_4011_LG$Can40_afreqs.FreqChange)/sqrt(var(cor_4011_LG$Lof07_afreqs.FreqChange1114)*var(cor_4011_LG$Can40_afreqs.FreqChange))
}

plot(x=NULL,y=NULL,xlim=c(0,24),ylim=c(-0.65,0.65),xlab="Linkage Group",ylab="Convergence Correlation")
points(x=seq(1,23)-0.1,y=cordfCan1114$cor1114,col="purple",pch=19)
abline(v=seq(0,23)+0.5,lty=2,lwd=0.5)
abline(h=0,lty=2)
legend(x=8,y=-0.05,legend=c("Norway 2011 - Norway 2014"),col=c("purple"),pch=c(19,19))

#bootstrap genome-wide value
cor_all_boot=data.frame(matrix(ncol=4,nrow=100))
colnames(cor_all_boot)=c("convcor1","convcor2","convcor3","convcor4")

for (i in 1:100) {
  covdfboot=covdf[sample(nrow(covdf),nrow(covdf),replace=TRUE),]
  cor_all_boot$convcor1[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
  cor_all_boot$convcor2[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
  cor_all_boot$convcor3[i]=cov(covdfboot$Lof07_afreqs.FreqChange11,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Lof07_afreqs.FreqChange11)*var(covdfboot$Lof07_afreqs.FreqChange14))
  cor_all_boot$convcor4[i]=cov(covdfboot$FreqChange1114,covdfboot$Can40_afreqs.FreqChange)/sqrt(var(covdfboot$FreqChange1114)*var(covdfboot$Can40_afreqs.FreqChange))
}

#bootstrap LGs

LGboot=data.frame(matrix(ncol=5,nrow=(length(LGs)-1)*100))
colnames(LGboot)=c("LG","convcor1","convcor2","convcor3","convcor4")

for (i in 1:(length(LGs)-1)) {
  c=LGs[i]
  cor_4011_LG=covdf[which(covdf$Can40_afreqs.CHROM==c),]
  for (j in 1:100) {
    LGboot$LG[(100*(i-1))+j]=c
    covdfboot=cor_4011_LG[sample(nrow(cor_4011_LG),nrow(cor_4011_LG),replace=TRUE),]
    LGboot$convcor1[(100*(i-1))+j]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
    LGboot$convcor2[(100*(i-1))+j]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
    LGboot$convcor3[(100*(i-1))+j]=cov(covdfboot$Lof07_afreqs.FreqChange11,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Lof07_afreqs.FreqChange11)*var(covdfboot$Lof07_afreqs.FreqChange14))
    LGboot$convcor4[(100*(i-1))+j]=cov(covdfboot$FreqChange1114,covdfboot$Can40_afreqs.FreqChange)/sqrt(var(covdfboot$FreqChange1114)*var(covdfboot$Can40_afreqs.FreqChange))
  }
}

#### compare ####

library(RColorBrewer)
cols=brewer.pal(name="Set1",n=5)

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(x=NULL,y=NULL,xlim=c(1,23),ylim=c(-0.12,0.70),xlab="Linkage Group",ylab="ConvCor",xaxt="n")
points(x=seq(1,23),y=cordf1114$cor1114,col=cols[4],pch=19,cex=0.8)
points(x=seq(1,23)-0.15,y=cordf$cor11,col=cols[2],pch=19,cex=0.8)
points(x=seq(1,23),y=cordfCan1114$cor1114,col=cols[5],pch=19,cex=0.8)
points(x=seq(1,23)+0.15,y=cordf$cor14,col=cols[3],pch=19,cex=0.8)
for (i in 0:23){
  lines(x=c(i,i)+0.5,y=c(-0.12,0.70),lty=2,col="gray")
}
for (i in 1:23){
  lines(x=c(i-0.15,i-0.15),y=c(quantile(LGboot$convcor1[which(LGboot$LG==LGs[i])],0.025),quantile(LGboot$convcor1[which(LGboot$LG==LGs[i])],0.975)),col=cols[2],lwd=2)
  lines(x=c(i+0.15,i+0.15),y=c(quantile(LGboot$convcor2[which(LGboot$LG==LGs[i])],0.025),quantile(LGboot$convcor2[which(LGboot$LG==LGs[i])],0.975)),col=cols[3],lwd=2)
  lines(x=c(i,i),y=c(quantile(LGboot$convcor3[which(LGboot$LG==LGs[i])],0.025),quantile(LGboot$convcor3[which(LGboot$LG==LGs[i])],0.975)),col=cols[4],lwd=2)
  lines(x=c(i,i),y=c(quantile(LGboot$convcor4[which(LGboot$LG==LGs[i])],0.025),quantile(LGboot$convcor4[which(LGboot$LG==LGs[i])],0.975)),col=cols[5],lwd=2)
  }
lines(x=c(0.5,23.5),y=c(0,0),lty=2)
legend("topright",inset=c(-0.2,0),bty="n",legend=c(expression('ConvCor'[1]),expression('ConvCor'[2]),expression('ConvCor'[3]),expression('ConvCor'[4])),col=c(cols[2],cols[3],cols[4],cols[5]),pch=c(19,19,19,19),cex=1)
axis(side=1,at=seq(1,23),label=c(LGs[1:23]))

### se's over linkage groups
se <- function(x) sqrt(var(x) / length(x))
se(cordf1114$cor1114)
se(cordf$cor11)
se(cordf$cor14)
se(cordfCan1114$cor1114)

### within inversions? ###

quals=read.table("data/Can40Loci_afreqs/out.lqual",header=T)

Can40_afreqs$POS=quals$POS
Can40_LG1inv_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG01" & Can40_afreqs$POS>9114741 & Can40_afreqs$POS<26192489),]
Can40_LG2inv_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG02" & Can40_afreqs$POS>18489307 & Can40_afreqs$POS<24050282),]
Can40_LG7inv_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG07" & Can40_afreqs$POS>13606502 & Can40_afreqs$POS<23016726),]
Can40_LG12inv_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG12" & Can40_afreqs$POS>589105 & Can40_afreqs$POS<13631347),]

allinvind=c(rownames(Can40_LG1inv_afreqs),rownames(Can40_LG2inv_afreqs),rownames(Can40_LG7inv_afreqs),rownames(Can40_LG12inv_afreqs))

Can40_noinv_afreqs=Can40_afreqs[!(rownames(Can40_afreqs) %in% allinvind),]

Lof07_afreqs$POS=quals$POS
Lof07_LG1inv_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG01" & Lof07_afreqs$POS>9114741 & Lof07_afreqs$POS<26192489),]
Lof07_LG2inv_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG02" & Lof07_afreqs$POS>18489307 & Lof07_afreqs$POS<24050282),]
Lof07_LG7inv_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG07" & Lof07_afreqs$POS>13606502 & Lof07_afreqs$POS<23016726),]
Lof07_LG12inv_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG12" & Lof07_afreqs$POS>589105 & Lof07_afreqs$POS<13631347),]

allinvind=c(rownames(Lof07_LG1inv_afreqs),rownames(Lof07_LG2inv_afreqs),rownames(Lof07_LG7inv_afreqs),rownames(Lof07_LG12inv_afreqs))

Lof07_noinv_afreqs=Lof07_afreqs[!(rownames(Lof07_afreqs) %in% allinvind),]

covdf_noinv=na.omit(data.frame(Can40_noinv_afreqs$CHROM,Can40_noinv_afreqs$FreqChange,Lof07_noinv_afreqs$FreqChange11,Lof07_noinv_afreqs$FreqChange14))
cor11_noinv=cov(covdf_noinv$Can40_noinv_afreqs.FreqChange,covdf_noinv$Lof07_noinv_afreqs.FreqChange11)/sqrt(var(covdf_noinv$Can40_noinv_afreqs.FreqChange)*var(covdf_noinv$Lof07_noinv_afreqs.FreqChange11))
cor14_noinv=cov(covdf_noinv$Can40_noinv_afreqs.FreqChange,covdf_noinv$Lof07_noinv_afreqs.FreqChange14)/sqrt(var(covdf_noinv$Can40_noinv_afreqs.FreqChange)*var(covdf_noinv$Lof07_noinv_afreqs.FreqChange14))

cor_noinv_boot=data.frame(matrix(ncol=2,nrow=100))
colnames(cor_noinv_boot)=c("convcor1","convcor2")

for (i in 1:100) {
  covdfboot=covdf_noinv[sample(nrow(covdf_noinv),nrow(covdf_noinv),replace=TRUE),]
  colnames(covdfboot)=c("Can40_afreqs.CHROM","Can40_afreqs.FreqChange","Lof07_afreqs.FreqChange11","Lof07_afreqs.FreqChange14")
  cor_noinv_boot$convcor1[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
  cor_noinv_boot$convcor2[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
}

covdf_LG1inv=na.omit(data.frame(Can40_LG1inv_afreqs$CHROM,Can40_LG1inv_afreqs$FreqChange,Lof07_LG1inv_afreqs$FreqChange11,Lof07_LG1inv_afreqs$FreqChange14))
cor11_LG1inv=cov(covdf_LG1inv$Can40_LG1inv_afreqs.FreqChange,covdf_LG1inv$Lof07_LG1inv_afreqs.FreqChange11)/sqrt(var(covdf_LG1inv$Can40_LG1inv_afreqs.FreqChange)*var(covdf_LG1inv$Lof07_LG1inv_afreqs.FreqChange11))
cor14_LG1inv=cov(covdf_LG1inv$Can40_LG1inv_afreqs.FreqChange,covdf_LG1inv$Lof07_LG1inv_afreqs.FreqChange14)/sqrt(var(covdf_LG1inv$Can40_LG1inv_afreqs.FreqChange)*var(covdf_LG1inv$Lof07_LG1inv_afreqs.FreqChange14))

cor_LG1inv_boot=data.frame(matrix(ncol=2,nrow=100))
colnames(cor_LG1inv_boot)=c("convcor1","convcor2")

for (i in 1:100) {
  covdfboot=covdf_LG1inv[sample(nrow(covdf_LG1inv),nrow(covdf_LG1inv),replace=TRUE),]
  colnames(covdfboot)=c("Can40_afreqs.CHROM","Can40_afreqs.FreqChange","Lof07_afreqs.FreqChange11","Lof07_afreqs.FreqChange14")
  cor_LG1inv_boot$convcor1[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
  cor_LG1inv_boot$convcor2[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
}

covdf_LG2inv=na.omit(data.frame(Can40_LG2inv_afreqs$CHROM,Can40_LG2inv_afreqs$FreqChange,Lof07_LG2inv_afreqs$FreqChange11,Lof07_LG2inv_afreqs$FreqChange14))
cor11_LG2inv=cov(covdf_LG2inv$Can40_LG2inv_afreqs.FreqChange,covdf_LG2inv$Lof07_LG2inv_afreqs.FreqChange11)/sqrt(var(covdf_LG2inv$Can40_LG2inv_afreqs.FreqChange)*var(covdf_LG2inv$Lof07_LG2inv_afreqs.FreqChange11))
cor14_LG2inv=cov(covdf_LG2inv$Can40_LG2inv_afreqs.FreqChange,covdf_LG2inv$Lof07_LG2inv_afreqs.FreqChange14)/sqrt(var(covdf_LG2inv$Can40_LG2inv_afreqs.FreqChange)*var(covdf_LG2inv$Lof07_LG2inv_afreqs.FreqChange14))

cor_LG2inv_boot=data.frame(matrix(ncol=2,nrow=100))
colnames(cor_LG2inv_boot)=c("convcor1","convcor2")

for (i in 1:100) {
  covdfboot=covdf_LG2inv[sample(nrow(covdf_LG2inv),nrow(covdf_LG2inv),replace=TRUE),]
  colnames(covdfboot)=c("Can40_afreqs.CHROM","Can40_afreqs.FreqChange","Lof07_afreqs.FreqChange11","Lof07_afreqs.FreqChange14")
  cor_LG2inv_boot$convcor1[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
  cor_LG2inv_boot$convcor2[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
}

covdf_LG7inv=na.omit(data.frame(Can40_LG7inv_afreqs$CHROM,Can40_LG7inv_afreqs$FreqChange,Lof07_LG7inv_afreqs$FreqChange11,Lof07_LG7inv_afreqs$FreqChange14))
cor11_LG7inv=cov(covdf_LG7inv$Can40_LG7inv_afreqs.FreqChange,covdf_LG7inv$Lof07_LG7inv_afreqs.FreqChange11)/sqrt(var(covdf_LG7inv$Can40_LG7inv_afreqs.FreqChange)*var(covdf_LG7inv$Lof07_LG7inv_afreqs.FreqChange11))
cor14_LG7inv=cov(covdf_LG7inv$Can40_LG7inv_afreqs.FreqChange,covdf_LG7inv$Lof07_LG7inv_afreqs.FreqChange14)/sqrt(var(covdf_LG7inv$Can40_LG7inv_afreqs.FreqChange)*var(covdf_LG7inv$Lof07_LG7inv_afreqs.FreqChange14))

cor_LG7inv_boot=data.frame(matrix(ncol=2,nrow=100))
colnames(cor_LG7inv_boot)=c("convcor1","convcor2")

for (i in 1:100) {
  covdfboot=covdf_LG7inv[sample(nrow(covdf_LG7inv),nrow(covdf_LG7inv),replace=TRUE),]
  colnames(covdfboot)=c("Can40_afreqs.CHROM","Can40_afreqs.FreqChange","Lof07_afreqs.FreqChange11","Lof07_afreqs.FreqChange14")
  cor_LG7inv_boot$convcor1[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
  cor_LG7inv_boot$convcor2[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
}

covdf_LG12inv=na.omit(data.frame(Can40_LG12inv_afreqs$CHROM,Can40_LG12inv_afreqs$FreqChange,Lof07_LG12inv_afreqs$FreqChange11,Lof07_LG12inv_afreqs$FreqChange14))
cor11_LG12inv=cov(covdf_LG12inv$Can40_LG12inv_afreqs.FreqChange,covdf_LG12inv$Lof07_LG12inv_afreqs.FreqChange11)/sqrt(var(covdf_LG12inv$Can40_LG12inv_afreqs.FreqChange)*var(covdf_LG12inv$Lof07_LG12inv_afreqs.FreqChange11))
cor14_LG12inv=cov(covdf_LG12inv$Can40_LG12inv_afreqs.FreqChange,covdf_LG12inv$Lof07_LG12inv_afreqs.FreqChange14)/sqrt(var(covdf_LG12inv$Can40_LG12inv_afreqs.FreqChange)*var(covdf_LG12inv$Lof07_LG12inv_afreqs.FreqChange14))

cor_LG12inv_boot=data.frame(matrix(ncol=2,nrow=100))
colnames(cor_LG12inv_boot)=c("convcor1","convcor2")

for (i in 1:100) {
  covdfboot=covdf_LG12inv[sample(nrow(covdf_LG12inv),nrow(covdf_LG12inv),replace=TRUE),]
  colnames(covdfboot)=c("Can40_afreqs.CHROM","Can40_afreqs.FreqChange","Lof07_afreqs.FreqChange11","Lof07_afreqs.FreqChange14")
  cor_LG12inv_boot$convcor1[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
  cor_LG12inv_boot$convcor2[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
}

Can40_LG1noinv1_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG01" & Can40_afreqs$POS<9114741),]
Can40_LG1noinv2_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG01" & Can40_afreqs$POS>26192489),]
Can40_LG1noinv_afreqs=rbind(Can40_LG1noinv1_afreqs,Can40_LG1noinv2_afreqs)
Lof07_LG1noinv1_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG01" & Lof07_afreqs$POS<9114741),]
Lof07_LG1noinv2_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG01" & Lof07_afreqs$POS>26192489),]
Lof07_LG1noinv_afreqs=rbind(Lof07_LG1noinv1_afreqs,Lof07_LG1noinv2_afreqs)
covdf_LG1noinv=na.omit(data.frame(Can40_LG1noinv_afreqs$CHROM,Can40_LG1noinv_afreqs$FreqChange,Lof07_LG1noinv_afreqs$FreqChange11,Lof07_LG1noinv_afreqs$FreqChange14))
cor11_LG1noinv=cov(covdf_LG1noinv$Can40_LG1noinv_afreqs.FreqChange,covdf_LG1noinv$Lof07_LG1noinv_afreqs.FreqChange11)/sqrt(var(covdf_LG1noinv$Can40_LG1noinv_afreqs.FreqChange)*var(covdf_LG1noinv$Lof07_LG1noinv_afreqs.FreqChange11))
cor14_LG1noinv=cov(covdf_LG1noinv$Can40_LG1noinv_afreqs.FreqChange,covdf_LG1noinv$Lof07_LG1noinv_afreqs.FreqChange14)/sqrt(var(covdf_LG1noinv$Can40_LG1noinv_afreqs.FreqChange)*var(covdf_LG1noinv$Lof07_LG1noinv_afreqs.FreqChange14))

cor_LG1noinv_boot=data.frame(matrix(ncol=2,nrow=100))
colnames(cor_LG1noinv_boot)=c("convcor1","convcor2")

for (i in 1:100) {
  covdfboot=covdf_LG1noinv[sample(nrow(covdf_LG1noinv),nrow(covdf_LG1noinv),replace=TRUE),]
  colnames(covdfboot)=c("Can40_afreqs.CHROM","Can40_afreqs.FreqChange","Lof07_afreqs.FreqChange11","Lof07_afreqs.FreqChange14")
  cor_LG1noinv_boot$convcor1[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
  cor_LG1noinv_boot$convcor2[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
}

Can40_LG2noinv1_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG02" & Can40_afreqs$POS<18489307),]
Can40_LG2noinv2_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG02" & Can40_afreqs$POS>24050282),]
Can40_LG2noinv_afreqs=rbind(Can40_LG2noinv1_afreqs,Can40_LG2noinv2_afreqs)
Lof07_LG2noinv1_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG02" & Lof07_afreqs$POS<18489307),]
Lof07_LG2noinv2_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG02" & Lof07_afreqs$POS>24050282),]
Lof07_LG2noinv_afreqs=rbind(Lof07_LG2noinv1_afreqs,Lof07_LG2noinv2_afreqs)
covdf_LG2noinv=na.omit(data.frame(Can40_LG2noinv_afreqs$CHROM,Can40_LG2noinv_afreqs$FreqChange,Lof07_LG2noinv_afreqs$FreqChange11,Lof07_LG2noinv_afreqs$FreqChange14))
cor11_LG2noinv=cov(covdf_LG2noinv$Can40_LG2noinv_afreqs.FreqChange,covdf_LG2noinv$Lof07_LG2noinv_afreqs.FreqChange11)/sqrt(var(covdf_LG2noinv$Can40_LG2noinv_afreqs.FreqChange)*var(covdf_LG2noinv$Lof07_LG2noinv_afreqs.FreqChange11))
cor14_LG2noinv=cov(covdf_LG2noinv$Can40_LG2noinv_afreqs.FreqChange,covdf_LG2noinv$Lof07_LG2noinv_afreqs.FreqChange14)/sqrt(var(covdf_LG2noinv$Can40_LG2noinv_afreqs.FreqChange)*var(covdf_LG2noinv$Lof07_LG2noinv_afreqs.FreqChange14))

cor_LG2noinv_boot=data.frame(matrix(ncol=2,nrow=100))
colnames(cor_LG2noinv_boot)=c("convcor1","convcor2")

for (i in 1:100) {
  covdfboot=covdf_LG2noinv[sample(nrow(covdf_LG2noinv),nrow(covdf_LG2noinv),replace=TRUE),]
  colnames(covdfboot)=c("Can40_afreqs.CHROM","Can40_afreqs.FreqChange","Lof07_afreqs.FreqChange11","Lof07_afreqs.FreqChange14")
  cor_LG2noinv_boot$convcor1[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
  cor_LG2noinv_boot$convcor2[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
}

Can40_LG7noinv1_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG07" & Can40_afreqs$POS<18489307),]
Can40_LG7noinv2_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG07" & Can40_afreqs$POS>24050282),]
Can40_LG7noinv_afreqs=rbind(Can40_LG7noinv1_afreqs,Can40_LG7noinv2_afreqs)
Lof07_LG7noinv1_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG07" & Lof07_afreqs$POS<18489307),]
Lof07_LG7noinv2_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG07" & Lof07_afreqs$POS>24050282),]
Lof07_LG7noinv_afreqs=rbind(Lof07_LG7noinv1_afreqs,Lof07_LG7noinv2_afreqs)
covdf_LG7noinv=na.omit(data.frame(Can40_LG7noinv_afreqs$CHROM,Can40_LG7noinv_afreqs$FreqChange,Lof07_LG7noinv_afreqs$FreqChange11,Lof07_LG7noinv_afreqs$FreqChange14))
cor11_LG7noinv=cov(covdf_LG7noinv$Can40_LG7noinv_afreqs.FreqChange,covdf_LG7noinv$Lof07_LG7noinv_afreqs.FreqChange11)/sqrt(var(covdf_LG7noinv$Can40_LG7noinv_afreqs.FreqChange)*var(covdf_LG7noinv$Lof07_LG7noinv_afreqs.FreqChange11))
cor14_LG7noinv=cov(covdf_LG7noinv$Can40_LG7noinv_afreqs.FreqChange,covdf_LG7noinv$Lof07_LG7noinv_afreqs.FreqChange14)/sqrt(var(covdf_LG7noinv$Can40_LG7noinv_afreqs.FreqChange)*var(covdf_LG7noinv$Lof07_LG7noinv_afreqs.FreqChange14))

cor_LG7noinv_boot=data.frame(matrix(ncol=2,nrow=100))
colnames(cor_LG7noinv_boot)=c("convcor1","convcor2")

for (i in 1:100) {
  covdfboot=covdf_LG7noinv[sample(nrow(covdf_LG7noinv),nrow(covdf_LG7noinv),replace=TRUE),]
  colnames(covdfboot)=c("Can40_afreqs.CHROM","Can40_afreqs.FreqChange","Lof07_afreqs.FreqChange11","Lof07_afreqs.FreqChange14")
  cor_LG7noinv_boot$convcor1[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
  cor_LG7noinv_boot$convcor2[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
}

Can40_LG12noinv1_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG12" & Can40_afreqs$POS<18489307),]
Can40_LG12noinv2_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG12" & Can40_afreqs$POS>24050282),]
Can40_LG12noinv_afreqs=rbind(Can40_LG12noinv1_afreqs,Can40_LG12noinv2_afreqs)
Lof07_LG12noinv1_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG12" & Lof07_afreqs$POS<18489307),]
Lof07_LG12noinv2_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG12" & Lof07_afreqs$POS>24050282),]
Lof07_LG12noinv_afreqs=rbind(Lof07_LG12noinv1_afreqs,Lof07_LG12noinv2_afreqs)
covdf_LG12noinv=na.omit(data.frame(Can40_LG12noinv_afreqs$CHROM,Can40_LG12noinv_afreqs$FreqChange,Lof07_LG12noinv_afreqs$FreqChange11,Lof07_LG12noinv_afreqs$FreqChange14))
cor11_LG12noinv=cov(covdf_LG12noinv$Can40_LG12noinv_afreqs.FreqChange,covdf_LG12noinv$Lof07_LG12noinv_afreqs.FreqChange11)/sqrt(var(covdf_LG12noinv$Can40_LG12noinv_afreqs.FreqChange)*var(covdf_LG12noinv$Lof07_LG12noinv_afreqs.FreqChange11))
cor14_LG12noinv=cov(covdf_LG12noinv$Can40_LG12noinv_afreqs.FreqChange,covdf_LG12noinv$Lof07_LG12noinv_afreqs.FreqChange14)/sqrt(var(covdf_LG12noinv$Can40_LG12noinv_afreqs.FreqChange)*var(covdf_LG12noinv$Lof07_LG12noinv_afreqs.FreqChange14))

cor_LG12noinv_boot=data.frame(matrix(ncol=2,nrow=100))
colnames(cor_LG12noinv_boot)=c("convcor1","convcor2")

for (i in 1:100) {
  covdfboot=covdf_LG12noinv[sample(nrow(covdf_LG12noinv),nrow(covdf_LG12noinv),replace=TRUE),]
  colnames(covdfboot)=c("Can40_afreqs.CHROM","Can40_afreqs.FreqChange","Lof07_afreqs.FreqChange11","Lof07_afreqs.FreqChange14")
  cor_LG12noinv_boot$convcor1[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
  cor_LG12noinv_boot$convcor2[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
}

##### SNPs in CDS

Can40_CDS_afreqs=read.table("data/afreqs_Can40_gadmor2_CDS/plink2.Can40.afreq",header=T)
Can13_CDS_afreqs=read.table("data/afreqs_Can40_gadmor2_CDS/plink2.Can13.afreq",header=T)
Can40_CDS_afreqs$ContFreq=Can13_CDS_afreqs$ALT_FREQS
Can40_CDS_afreqs$FreqChange=Can40_CDS_afreqs$ContFreq-Can40_CDS_afreqs$ALT_FREQS
Can40_CDS_afreqs$AbsFreqChange=abs(Can40_CDS_afreqs$ContFreq-Can40_CDS_afreqs$ALT_FREQS)

Lof07_CDS_afreqs=read.table("data/afreqs_Can40_gadmor2_CDS/plink2.Lof07.afreq",header=T)
Lof11_CDS_afreqs=read.table("data/afreqs_Can40_gadmor2_CDS/plink2.Lof11.afreq",header=T)
Lof14_CDS_afreqs=read.table("data/afreqs_Can40_gadmor2_CDS/plink2.Lof14.afreq",header=T)
Lof07_CDS_afreqs$Freq11=Lof11_CDS_afreqs$ALT_FREQS
Lof07_CDS_afreqs$Freq14=Lof14_CDS_afreqs$ALT_FREQS
Lof07_CDS_afreqs$FreqChange11=Lof07_CDS_afreqs$Freq11-Lof07_CDS_afreqs$ALT_FREQS
Lof07_CDS_afreqs$FreqChange14=Lof07_CDS_afreqs$Freq14-Lof07_CDS_afreqs$ALT_FREQS

covdf_CDS=data.frame(Can40_CDS_afreqs$CHROM,Can40_CDS_afreqs$FreqChange,Lof07_CDS_afreqs$FreqChange11,Lof07_CDS_afreqs$FreqChange14)
covdf_CDS$FreqChange1114=Lof07_CDS_afreqs$Freq14-Lof07_CDS_afreqs$Freq11
covdf_CDS=na.omit(covdf_CDS)

cor_CDS_Can40_Lof11=cov(covdf_CDS$Can40_CDS_afreqs.FreqChange,covdf_CDS$Lof07_CDS_afreqs.FreqChange11)/sqrt(var(covdf_CDS$Can40_CDS_afreqs.FreqChange)*var(covdf_CDS$Lof07_CDS_afreqs.FreqChange11))
cor_CDS_Can40_Lof14=cov(covdf_CDS$Can40_CDS_afreqs.FreqChange,covdf_CDS$Lof07_CDS_afreqs.FreqChange14)/sqrt(var(covdf_CDS$Can40_CDS_afreqs.FreqChange)*var(covdf_CDS$Lof07_CDS_afreqs.FreqChange14))
cor_CDS_Lof11_Lof14=cov(covdf_CDS$Lof07_CDS_afreqs.FreqChange11,covdf_CDS$Lof07_CDS_afreqs.FreqChange14)/sqrt(var(covdf_CDS$Lof07_CDS_afreqs.FreqChange11)*var(covdf_CDS$Lof07_CDS_afreqs.FreqChange14))
cor_CDS_Lof1114_Can=cov(covdf_CDS$FreqChange1114,covdf_CDS$Can40_CDS_afreqs.FreqChange)/sqrt(var(covdf_CDS$FreqChange1114)*var(covdf_CDS$Can40_CDS_afreqs.FreqChange))

### bootstrap CDS

cor_CDS_boot=data.frame(matrix(ncol=4,nrow=100))
colnames(cor_CDS_boot)=c("convcor1","convcor2","convcor3","convcor4")

for (i in 1:100) {
  covdfboot=covdf_CDS[sample(nrow(covdf_CDS),nrow(covdf_CDS),replace=TRUE),]
  cor_CDS_boot$convcor1[i]=cov(covdfboot$Can40_CDS_afreqs.FreqChange,covdfboot$Lof07_CDS_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_CDS_afreqs.FreqChange)*var(covdfboot$Lof07_CDS_afreqs.FreqChange11))
  cor_CDS_boot$convcor2[i]=cov(covdfboot$Can40_CDS_afreqs.FreqChange,covdfboot$Lof07_CDS_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_CDS_afreqs.FreqChange)*var(covdfboot$Lof07_CDS_afreqs.FreqChange14))
  cor_CDS_boot$convcor3[i]=cov(covdfboot$Lof07_CDS_afreqs.FreqChange11,covdfboot$Lof07_CDS_afreqs.FreqChange14)/sqrt(var(covdfboot$Lof07_CDS_afreqs.FreqChange11)*var(covdfboot$Lof07_CDS_afreqs.FreqChange14))
  cor_CDS_boot$convcor4[i]=cov(covdfboot$FreqChange1114,covdfboot$Can40_CDS_afreqs.FreqChange)/sqrt(var(covdfboot$FreqChange1114)*var(covdfboot$Can40_CDS_afreqs.FreqChange))
}

#### plotting subsets

par(las=2,mai=c(2,1,1,1))

par(las=2,mar=c(10.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(x=NULL,y=NULL,xlim=c(1.3,8.2),ylim=c(-0.11,0.2),xlab="",ylab="ConvCor",xaxt="n")
points(x=1.5-0.1,y=cor11_LG1inv,col=cols[2],pch=19)
points(x=1.5+0.1,y=cor14_LG1inv,col=cols[3],pch=19)
points(x=2-0.1,y=cor11_LG1noinv,col=cols[2],pch=19)
points(x=2+0.1,y=cor14_LG1noinv,col=cols[3],pch=19)
points(x=3-0.1,y=cor11_LG2inv,col=cols[2],pch=19)
points(x=3+0.1,y=cor14_LG2inv,col=cols[3],pch=19)
points(x=3.5-0.1,y=cor11_LG2noinv,col=cols[2],pch=19)
points(x=3.5+0.1,y=cor14_LG2noinv,col=cols[3],pch=19)
points(x=4.5-0.1,y=cor11_LG7inv,col=cols[2],pch=19)
points(x=4.5+0.1,y=cor14_LG7inv,col=cols[3],pch=19)
points(x=5-0.1,y=cor11_LG7noinv,col=cols[2],pch=19)
points(x=5+0.1,y=cor14_LG7noinv,col=cols[3],pch=19)
points(x=6-0.1,y=cor11_LG12inv,col=cols[2],pch=19)
points(x=6+0.1,y=cor14_LG12inv,col=cols[3],pch=19)
points(x=6.5-0.1,y=cor11_LG12noinv,col=cols[2],pch=19)
points(x=6.5+0.1,y=cor14_LG12noinv,col=cols[3],pch=19)
points(x=7.5-0.1,y=cor_CDS_Can40_Lof11,col=cols[2],pch=19)
points(x=7.5+0.1,y=cor_CDS_Can40_Lof14,col=cols[3],pch=19)
points(x=8-0.1,y=cor_Can40_Lof11,col=cols[2],pch=19)
points(x=8+0.1,y=cor_Can40_Lof14,col=cols[3],pch=19)
lines(x=c(1.5-0.1,1.5-0.1),y=c(quantile(cor_LG1inv_boot$convcor1,0.025),quantile(cor_LG1inv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(1.5+0.1,1.5+0.1),y=c(quantile(cor_LG1inv_boot$convcor2,0.025),quantile(cor_LG1inv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(2-0.1,2-0.1),y=c(quantile(cor_LG1noinv_boot$convcor1,0.025),quantile(cor_LG1noinv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(2+0.1,2+0.1),y=c(quantile(cor_LG1noinv_boot$convcor2,0.025),quantile(cor_LG1noinv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(3-0.1,3-0.1),y=c(quantile(cor_LG2inv_boot$convcor1,0.025),quantile(cor_LG2inv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(3+0.1,3+0.1),y=c(quantile(cor_LG2inv_boot$convcor2,0.025),quantile(cor_LG2inv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(3.5-0.1,3.5-0.1),y=c(quantile(cor_LG2noinv_boot$convcor1,0.025),quantile(cor_LG2noinv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(3.5+0.1,3.5+0.1),y=c(quantile(cor_LG2noinv_boot$convcor2,0.025),quantile(cor_LG2noinv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(4.5-0.1,4.5-0.1),y=c(quantile(cor_LG7inv_boot$convcor1,0.025),quantile(cor_LG7inv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(4.5+0.1,4.5+0.1),y=c(quantile(cor_LG7inv_boot$convcor2,0.025),quantile(cor_LG7inv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(5-0.1,5-0.1),y=c(quantile(cor_LG7noinv_boot$convcor1,0.025),quantile(cor_LG7noinv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(5+0.1,5+0.1),y=c(quantile(cor_LG7noinv_boot$convcor2,0.025),quantile(cor_LG7noinv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(6-0.1,6-0.1),y=c(quantile(cor_LG12inv_boot$convcor1,0.025),quantile(cor_LG12inv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(6+0.1,6+0.1),y=c(quantile(cor_LG12inv_boot$convcor2,0.025),quantile(cor_LG12inv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(6.5-0.1,6.5-0.1),y=c(quantile(cor_LG12noinv_boot$convcor1,0.025),quantile(cor_LG12noinv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(6.5+0.1,6.5+0.1),y=c(quantile(cor_LG12noinv_boot$convcor2,0.025),quantile(cor_LG12noinv_boot$convcor2,0.975)),col=cols[3],lwd=2)

lines(x=c(7.5-0.1,7.5-0.1),y=c(quantile(cor_CDS_boot$convcor1,0.025),quantile(cor_CDS_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(7.5+0.1,7.5+0.1),y=c(quantile(cor_CDS_boot$convcor2,0.025),quantile(cor_CDS_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(8-0.1,8-0.1),y=c(quantile(cor_all_boot$convcor1,0.025),quantile(cor_LG7inv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(8+0.1,8+0.1),y=c(quantile(cor_all_boot$convcor2,0.025),quantile(cor_LG7inv_boot$convcor2,0.975)),col=cols[3],lwd=2)

lines(x=c(1.3,8.2),y=c(0,0),lty=2)
#legend("topright",bty="n",inset=c(-0.2,0),legend=c(expression('ConvCor'[1]),expression('ConvCor'[2])),col=c(cols[2],cols[3]),pch=c(19,19))
axis(side=1,at=c(1.5,2,3,3.5,4.5,5,6,6.5,7.5,8),label=c("LG01 inversion","LG01 outside inversion","LG02 inversion","LG02 outside inversion","LG07 inversion","LG07 outside inversion","LG12 inversion","LG12 outside inversion","Coding SNPs","All SNPs"))


#### Fig1 plots together

#par(las=2,mar=c(10.1, 0.1, 2.1, 0.1), xpd=TRUE,mfrow=c(2,1))

#library(RColorBrewer)
cols=brewer.pal(name="Set1",n=5)

par(mar=c(3.1, 4.1, 4.1, 0.1),xpd=TRUE,mfrow=c(2,1),las=2)

plot(x=NULL,y=NULL,xlim=c(1,23),ylim=c(-0.12,0.70),xlab="Linkage Group",ylab="ConvCor",xaxt="n")
points(x=seq(1,23),y=cordf1114$cor1114,col=cols[4],pch=19,cex=0.8)
points(x=seq(1,23)-0.15,y=cordf$cor11,col=cols[2],pch=19,cex=0.8)
points(x=seq(1,23),y=cordfCan1114$cor1114,col=cols[5],pch=19,cex=0.8)
points(x=seq(1,23)+0.15,y=cordf$cor14,col=cols[3],pch=19,cex=0.8)
for (i in 0:23){
  lines(x=c(i,i)+0.5,y=c(-0.12,0.70),lty=2,col="gray")
}
for (i in 1:23){
  lines(x=c(i-0.15,i-0.15),y=c(quantile(LGboot$convcor1[which(LGboot$LG==LGs[i])],0.025),quantile(LGboot$convcor1[which(LGboot$LG==LGs[i])],0.975)),col=cols[2],lwd=2)
  lines(x=c(i+0.15,i+0.15),y=c(quantile(LGboot$convcor2[which(LGboot$LG==LGs[i])],0.025),quantile(LGboot$convcor2[which(LGboot$LG==LGs[i])],0.975)),col=cols[3],lwd=2)
  lines(x=c(i,i),y=c(quantile(LGboot$convcor3[which(LGboot$LG==LGs[i])],0.025),quantile(LGboot$convcor3[which(LGboot$LG==LGs[i])],0.975)),col=cols[4],lwd=2)
  lines(x=c(i,i),y=c(quantile(LGboot$convcor4[which(LGboot$LG==LGs[i])],0.025),quantile(LGboot$convcor4[which(LGboot$LG==LGs[i])],0.975)),col=cols[5],lwd=2)
}
lines(x=c(0.5,23.5),y=c(0,0),lty=2)
#legend("topright",inset=c(-0.2,0),bty="n",legend=c(expression('ConvCor'[1]),expression('ConvCor'[2]),expression('ConvCor'[3]),expression('ConvCor'[4])),col=c(cols[2],cols[3],cols[4],cols[5]),pch=c(19,19,19,19),cex=1)
axis(side=1,at=seq(1,23),label=c(LGs[1:23]))

text(x=-1.5,y=0.68,labels="(a)",cex=1.5)
points(x=c(1,7,13,19),y=rep(1,4),col=c(cols[2],cols[3],cols[4],cols[5]),pch=c(19,19,19,19))
text(x=c(3,9,15,21),y=rep(1,4),labels=c(expression('ConvCor'[1]),expression('ConvCor'[2]),expression('ConvCor'[3]),expression('ConvCor'[4])))
text(x=c(3,9,15,21),y=rep(0.9,4),labels=c("Canada 1940-2013,","Canada 1940-2013,","Norway 1907-2011,","Canada 1940-2013,"))
text(x=c(3,9,15,21),y=rep(0.8,4),labels=c("Norway 1907-2011","Norway 1907-2014","Norway 1907-2014","Norway 2011-2014"))


par(mar=c(6.1, 4.1, 2.1, 0.1))
    
plot(x=NULL,y=NULL,xlim=c(1.3,8.2),ylim=c(-0.11,0.2),xlab="",ylab="ConvCor",xaxt="n")
points(x=1.5-0.1,y=cor11_LG1inv,col=cols[2],pch=19)
points(x=1.5+0.1,y=cor14_LG1inv,col=cols[3],pch=19)
points(x=2-0.1,y=cor11_LG1noinv,col=cols[2],pch=19)
points(x=2+0.1,y=cor14_LG1noinv,col=cols[3],pch=19)
points(x=3-0.1,y=cor11_LG2inv,col=cols[2],pch=19)
points(x=3+0.1,y=cor14_LG2inv,col=cols[3],pch=19)
points(x=3.5-0.1,y=cor11_LG2noinv,col=cols[2],pch=19)
points(x=3.5+0.1,y=cor14_LG2noinv,col=cols[3],pch=19)
points(x=4.5-0.1,y=cor11_LG7inv,col=cols[2],pch=19)
points(x=4.5+0.1,y=cor14_LG7inv,col=cols[3],pch=19)
points(x=5-0.1,y=cor11_LG7noinv,col=cols[2],pch=19)
points(x=5+0.1,y=cor14_LG7noinv,col=cols[3],pch=19)
points(x=6-0.1,y=cor11_LG12inv,col=cols[2],pch=19)
points(x=6+0.1,y=cor14_LG12inv,col=cols[3],pch=19)
points(x=6.5-0.1,y=cor11_LG12noinv,col=cols[2],pch=19)
points(x=6.5+0.1,y=cor14_LG12noinv,col=cols[3],pch=19)
points(x=7.5-0.1,y=cor_CDS_Can40_Lof11,col=cols[2],pch=19)
points(x=7.5+0.1,y=cor_CDS_Can40_Lof14,col=cols[3],pch=19)
points(x=8-0.1,y=cor_Can40_Lof11,col=cols[2],pch=19)
points(x=8+0.1,y=cor_Can40_Lof14,col=cols[3],pch=19)
lines(x=c(1.5-0.1,1.5-0.1),y=c(quantile(cor_LG1inv_boot$convcor1,0.025),quantile(cor_LG1inv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(1.5+0.1,1.5+0.1),y=c(quantile(cor_LG1inv_boot$convcor2,0.025),quantile(cor_LG1inv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(2-0.1,2-0.1),y=c(quantile(cor_LG1noinv_boot$convcor1,0.025),quantile(cor_LG1noinv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(2+0.1,2+0.1),y=c(quantile(cor_LG1noinv_boot$convcor2,0.025),quantile(cor_LG1noinv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(3-0.1,3-0.1),y=c(quantile(cor_LG2inv_boot$convcor1,0.025),quantile(cor_LG2inv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(3+0.1,3+0.1),y=c(quantile(cor_LG2inv_boot$convcor2,0.025),quantile(cor_LG2inv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(3.5-0.1,3.5-0.1),y=c(quantile(cor_LG2noinv_boot$convcor1,0.025),quantile(cor_LG2noinv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(3.5+0.1,3.5+0.1),y=c(quantile(cor_LG2noinv_boot$convcor2,0.025),quantile(cor_LG2noinv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(4.5-0.1,4.5-0.1),y=c(quantile(cor_LG7inv_boot$convcor1,0.025),quantile(cor_LG7inv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(4.5+0.1,4.5+0.1),y=c(quantile(cor_LG7inv_boot$convcor2,0.025),quantile(cor_LG7inv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(5-0.1,5-0.1),y=c(quantile(cor_LG7noinv_boot$convcor1,0.025),quantile(cor_LG7noinv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(5+0.1,5+0.1),y=c(quantile(cor_LG7noinv_boot$convcor2,0.025),quantile(cor_LG7noinv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(6-0.1,6-0.1),y=c(quantile(cor_LG12inv_boot$convcor1,0.025),quantile(cor_LG12inv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(6+0.1,6+0.1),y=c(quantile(cor_LG12inv_boot$convcor2,0.025),quantile(cor_LG12inv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(6.5-0.1,6.5-0.1),y=c(quantile(cor_LG12noinv_boot$convcor1,0.025),quantile(cor_LG12noinv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(6.5+0.1,6.5+0.1),y=c(quantile(cor_LG12noinv_boot$convcor2,0.025),quantile(cor_LG12noinv_boot$convcor2,0.975)),col=cols[3],lwd=2)

lines(x=c(7.5-0.1,7.5-0.1),y=c(quantile(cor_CDS_boot$convcor1,0.025),quantile(cor_CDS_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(7.5+0.1,7.5+0.1),y=c(quantile(cor_CDS_boot$convcor2,0.025),quantile(cor_CDS_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(8-0.1,8-0.1),y=c(quantile(cor_all_boot$convcor1,0.025),quantile(cor_LG7inv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(8+0.1,8+0.1),y=c(quantile(cor_all_boot$convcor2,0.025),quantile(cor_LG7inv_boot$convcor2,0.975)),col=cols[3],lwd=2)

lines(x=c(1.3,8.2),y=c(0,0),lty=2)
#legend("topright",bty="n",inset=c(-0.2,0),legend=c(expression('ConvCor'[1]),expression('ConvCor'[2])),col=c(cols[2],cols[3]),pch=c(19,19))
axis(side=1,at=c(1.5,2,3,3.5,4.5,5,6,6.5,7.5,8),label=c("Inversion","Outside","Inversion","Outside","Inversion","Outside","Inversion","Outside","Coding","All SNPs"))

text(x=c(1.8,3.3,4.8,6.3,7.8),y=rep(-0.35,5),labels=c("LG01","LG02","LG07","LG12","Genome"))

#points(x=c(1.2,3.2,5.2,7.2),y=rep(0.25,4),col=c(cols[2],cols[3],cols[4],cols[5]),pch=c(19,19,19,19))
#text(x=c(1.8,3.8,5.8,7.8),y=rep(0.25,4),labels=c(expression('ConvCor'[1]),expression('ConvCor'[2]),expression('ConvCor'[3]),expression('ConvCor'[4])))

text(x=0.5,y=0.25,labels="(b)",cex=1.5)

#########################################
###### original- no filter ##############
#########################################


Can40_afreqs=read.table("data/AllLoci_afreqs/plink2.Can40.afreq",header=T)
Can13_afreqs=read.table("data/AllLoci_afreqs/plink2.Can13.afreq",header=T)

Can40_afreqs$ContFreq=Can13_afreqs$ALT_FREQS
Can40_afreqs$FreqChange=Can40_afreqs$ContFreq-Can40_afreqs$ALT_FREQS

Lof07_afreqs=read.table("data/AllLoci_afreqs/plink2.Lof07.afreq",header=T)
Lof11_afreqs=read.table("data/AllLoci_afreqs/plink2.Lof11.afreq",header=T)
Lof14_afreqs=read.table("data/AllLoci_afreqs/plink2.Lof14.afreq",header=T)

Lof07_afreqs$Freq11=Lof11_afreqs$ALT_FREQS
Lof07_afreqs$Freq14=Lof14_afreqs$ALT_FREQS

Lof07_afreqs$FreqChange11=Lof07_afreqs$Freq11-Lof07_afreqs$ALT_FREQS
Lof07_afreqs$FreqChange14=Lof07_afreqs$Freq14-Lof07_afreqs$ALT_FREQS


Can40_afreqs$FreqChange
Lof07_afreqs$FreqChange11
Lof07_afreqs$FreqChange14

covdf=data.frame(Can40_afreqs$CHROM,Can40_afreqs$FreqChange,Lof07_afreqs$FreqChange11,Lof07_afreqs$FreqChange14)

covdf$FreqChange1114=Lof07_afreqs$Freq14-Lof07_afreqs$Freq11

covdf=na.omit(covdf)

cor_Can40_Lof11=cov(covdf$Can40_afreqs.FreqChange,covdf$Lof07_afreqs.FreqChange11)/sqrt(var(covdf$Can40_afreqs.FreqChange)*var(covdf$Lof07_afreqs.FreqChange11))
cor_Can40_Lof14=cov(covdf$Can40_afreqs.FreqChange,covdf$Lof07_afreqs.FreqChange14)/sqrt(var(covdf$Can40_afreqs.FreqChange)*var(covdf$Lof07_afreqs.FreqChange14))

cor_all=(cov(covdf$Can40_afreqs.FreqChange,covdf$Lof07_afreqs.FreqChange11)+cov(covdf$Can40_afreqs.FreqChange,covdf$Lof07_afreqs.FreqChange14))/(sqrt(var(covdf$Can40_afreqs.FreqChange)*var(covdf$Lof07_afreqs.FreqChange11))+sqrt(var(covdf$Can40_afreqs.FreqChange)*var(covdf$Lof07_afreqs.FreqChange14)))
cor_all=(cov(covdf$Can40_afreqs.FreqChange,covdf$Lof07_afreqs.FreqChange11)+cov(covdf$Can40_afreqs.FreqChange,covdf$Lof07_afreqs.FreqChange14))/(sqrt(var(covdf$Can40_afreqs.FreqChange)*var(covdf$Lof07_afreqs.FreqChange11))+sqrt(var(covdf$Can40_afreqs.FreqChange)*var(covdf$Lof07_afreqs.FreqChange14)))

### by LG
LGs=unique(covdf$Can40_afreqs.CHROM)

cordf=data.frame(matrix(ncol=3,nrow=23))
colnames(cordf)=c("LG","cor11","cor14")

for (i in 1:(length(LGs)-1)) {
  c=LGs[i]
  cor_4011_LG=covdf[which(covdf$Can40_afreqs.CHROM==c),]
  cordf$LG[i]=c
  cordf$cor11[i]=cov(cor_4011_LG$Can40_afreqs.FreqChange,cor_4011_LG$Lof07_afreqs.FreqChange11)/sqrt(var(cor_4011_LG$Can40_afreqs.FreqChange)*var(cor_4011_LG$Lof07_afreqs.FreqChange11))
  cordf$cor14[i]=cov(cor_4011_LG$Can40_afreqs.FreqChange,cor_4011_LG$Lof07_afreqs.FreqChange14)/sqrt(var(cor_4011_LG$Can40_afreqs.FreqChange)*var(cor_4011_LG$Lof07_afreqs.FreqChange14))
}

plot(x=NULL,y=NULL,xlim=c(0,24),ylim=c(-0.15,0.15),xlab="Linkage Group",ylab="Convergence Correlation")
points(x=seq(1,23)-0.1,y=cordf$cor11,col="blue",pch=19)
points(x=seq(1,23)+0.1,y=cordf$cor14,col="green",pch=19)
abline(v=seq(0,23)+0.5,lty=2,lwd=0.5)
abline(h=0,lty=2)
legend(x=8,y=-0.05,legend=c("Canada-Norway 2011","Canada-Norway 2014"),col=c("blue","green"),pch=c(19,19))

### change covariance - Norway07-11 vs Norway01-14!

cor_Lof11_Lof14=cov(covdf$Lof07_afreqs.FreqChange11,covdf$Lof07_afreqs.FreqChange14)/sqrt(var(covdf$Lof07_afreqs.FreqChange11)*var(covdf$Lof07_afreqs.FreqChange14))

LGs=unique(covdf$Can40_afreqs.CHROM)

cordf1114=data.frame(matrix(ncol=3,nrow=23))
colnames(cordf1114)=c("LG","cor1114")

for (i in 1:(length(LGs)-1)) {
  c=LGs[i]
  cor_4011_LG=covdf[which(covdf$Can40_afreqs.CHROM==c),]
  cordf1114$LG[i]=c
  cordf1114$cor1114[i]=cov(cor_4011_LG$Lof07_afreqs.FreqChange11,cor_4011_LG$Lof07_afreqs.FreqChange14)/sqrt(var(cor_4011_LG$Lof07_afreqs.FreqChange11)*var(cor_4011_LG$Lof07_afreqs.FreqChange14))
}

plot(x=NULL,y=NULL,xlim=c(0,24),ylim=c(-0.65,0.65),xlab="Linkage Group",ylab="Convergence Correlation")
points(x=seq(1,23)-0.1,y=cordf1114$cor1114,col="purple",pch=19)
abline(v=seq(0,23)+0.5,lty=2,lwd=0.5)
abline(h=0,lty=2)
legend(x=8,y=-0.05,legend=c("Norway 2011 - Norway 2014"),col=c("purple"),pch=c(19,19))

### change covariance - Norway11-14 vs Canada!

Lof07_afreqs$FreqChange1114=Lof07_afreqs$Freq14-Lof07_afreqs$Freq11

covdf1114=data.frame(Can40_afreqs$CHROM,Can40_afreqs$FreqChange,Lof07_afreqs$FreqChange1114)

covdf1114=na.omit(covdf1114)

cor_Lof1114_Can=cov(covdf1114$Lof07_afreqs.FreqChange1114,covdf1114$Can40_afreqs.FreqChange)/sqrt(var(covdf1114$Lof07_afreqs.FreqChange1114)*var(covdf1114$Can40_afreqs.FreqChange))

LGs=unique(covdf1114$Can40_afreqs.CHROM)

cordfCan1114=data.frame(matrix(ncol=3,nrow=23))
colnames(cordfCan1114)=c("LG","cor1114")

for (i in 1:(length(LGs)-1)) {
  c=LGs[i]
  cor_4011_LG=covdf1114[which(covdf1114$Can40_afreqs.CHROM==c),]
  cordfCan1114$LG[i]=c
  cordfCan1114$cor1114[i]=cov(cor_4011_LG$Lof07_afreqs.FreqChange1114,cor_4011_LG$Can40_afreqs.FreqChange)/sqrt(var(cor_4011_LG$Lof07_afreqs.FreqChange1114)*var(cor_4011_LG$Can40_afreqs.FreqChange))
}

plot(x=NULL,y=NULL,xlim=c(0,24),ylim=c(-0.65,0.65),xlab="Linkage Group",ylab="Convergence Correlation")
points(x=seq(1,23)-0.1,y=cordfCan1114$cor1114,col="purple",pch=19)
abline(v=seq(0,23)+0.5,lty=2,lwd=0.5)
abline(h=0,lty=2)
legend(x=8,y=-0.05,legend=c("Norway 2011 - Norway 2014"),col=c("purple"),pch=c(19,19))

#bootstrap genome-wide value
cor_all_boot=data.frame(matrix(ncol=4,nrow=100))
colnames(cor_all_boot)=c("convcor1","convcor2","convcor3","convcor4")

for (i in 1:100) {
  covdfboot=covdf[sample(nrow(covdf),nrow(covdf),replace=TRUE),]
  cor_all_boot$convcor1[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
  cor_all_boot$convcor2[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
  cor_all_boot$convcor3[i]=cov(covdfboot$Lof07_afreqs.FreqChange11,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Lof07_afreqs.FreqChange11)*var(covdfboot$Lof07_afreqs.FreqChange14))
  cor_all_boot$convcor4[i]=cov(covdfboot$FreqChange1114,covdfboot$Can40_afreqs.FreqChange)/sqrt(var(covdfboot$FreqChange1114)*var(covdfboot$Can40_afreqs.FreqChange))
}

#bootstrap LGs

LGboot=data.frame(matrix(ncol=5,nrow=(length(LGs)-1)*100))
colnames(LGboot)=c("LG","convcor1","convcor2","convcor3","convcor4")

for (i in 1:(length(LGs)-1)) {
  c=LGs[i]
  cor_4011_LG=covdf[which(covdf$Can40_afreqs.CHROM==c),]
  for (j in 1:100) {
    LGboot$LG[(100*(i-1))+j]=c
    covdfboot=cor_4011_LG[sample(nrow(cor_4011_LG),nrow(cor_4011_LG),replace=TRUE),]
    LGboot$convcor1[(100*(i-1))+j]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
    LGboot$convcor2[(100*(i-1))+j]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
    LGboot$convcor3[(100*(i-1))+j]=cov(covdfboot$Lof07_afreqs.FreqChange11,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Lof07_afreqs.FreqChange11)*var(covdfboot$Lof07_afreqs.FreqChange14))
    LGboot$convcor4[(100*(i-1))+j]=cov(covdfboot$FreqChange1114,covdfboot$Can40_afreqs.FreqChange)/sqrt(var(covdfboot$FreqChange1114)*var(covdfboot$Can40_afreqs.FreqChange))
  }
}

#### compare ####

library(RColorBrewer)
cols=brewer.pal(name="Set1",n=5)

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(x=NULL,y=NULL,xlim=c(1,23),ylim=c(-0.3,0.75),xlab="Linkage Group",ylab="ConvCor",xaxt="n")
points(x=seq(1,23),y=cordf1114$cor1114,col=cols[4],pch=19,cex=0.8)
points(x=seq(1,23)-0.15,y=cordf$cor11,col=cols[2],pch=19,cex=0.8)
points(x=seq(1,23),y=cordfCan1114$cor1114,col=cols[5],pch=19,cex=0.8)
points(x=seq(1,23)+0.15,y=cordf$cor14,col=cols[3],pch=19,cex=0.8)
for (i in 0:23){
  lines(x=c(i,i)+0.5,y=c(-0.3,0.75),lty=2,col="gray")
}
for (i in 1:23){
  lines(x=c(i-0.15,i-0.15),y=c(quantile(LGboot$convcor1[which(LGboot$LG==LGs[i])],0.025),quantile(LGboot$convcor1[which(LGboot$LG==LGs[i])],0.975)),col=cols[2],lwd=2)
  lines(x=c(i+0.15,i+0.15),y=c(quantile(LGboot$convcor2[which(LGboot$LG==LGs[i])],0.025),quantile(LGboot$convcor2[which(LGboot$LG==LGs[i])],0.975)),col=cols[3],lwd=2)
  lines(x=c(i,i),y=c(quantile(LGboot$convcor3[which(LGboot$LG==LGs[i])],0.025),quantile(LGboot$convcor3[which(LGboot$LG==LGs[i])],0.975)),col=cols[4],lwd=2)
  lines(x=c(i,i),y=c(quantile(LGboot$convcor4[which(LGboot$LG==LGs[i])],0.025),quantile(LGboot$convcor4[which(LGboot$LG==LGs[i])],0.975)),col=cols[5],lwd=2)
}
lines(x=c(0.5,23.5),y=c(0,0),lty=2)
legend("topright",inset=c(-0.2,0),bty="n",legend=c(expression('ConvCor'[1]),expression('ConvCor'[2]),expression('ConvCor'[3]),expression('ConvCor'[4])),col=c(cols[2],cols[3],cols[4],cols[5]),pch=c(19,19,19,19),cex=1)
axis(side=1,at=seq(1,23),label=c(LGs[1:23]))

### se's over linkage groups
se <- function(x) sqrt(var(x) / length(x))
se(cordf1114$cor1114)
se(cordf$cor11)
se(cordf$cor14)
se(cordfCan1114$cor1114)

### within inversions? ###

quals=read.table("data/AllLoci_afreqs/out.lqual",header=T)

Can40_afreqs$POS=quals$POS
Can40_LG1inv_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG01" & Can40_afreqs$POS>9114741 & Can40_afreqs$POS<26192489),]
Can40_LG2inv_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG02" & Can40_afreqs$POS>18489307 & Can40_afreqs$POS<24050282),]
Can40_LG7inv_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG07" & Can40_afreqs$POS>13606502 & Can40_afreqs$POS<23016726),]
Can40_LG12inv_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG12" & Can40_afreqs$POS>589105 & Can40_afreqs$POS<13631347),]

allinvind=c(rownames(Can40_LG1inv_afreqs),rownames(Can40_LG2inv_afreqs),rownames(Can40_LG7inv_afreqs),rownames(Can40_LG12inv_afreqs))

Can40_noinv_afreqs=Can40_afreqs[!(rownames(Can40_afreqs) %in% allinvind),]

Lof07_afreqs$POS=quals$POS
Lof07_LG1inv_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG01" & Lof07_afreqs$POS>9114741 & Lof07_afreqs$POS<26192489),]
Lof07_LG2inv_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG02" & Lof07_afreqs$POS>18489307 & Lof07_afreqs$POS<24050282),]
Lof07_LG7inv_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG07" & Lof07_afreqs$POS>13606502 & Lof07_afreqs$POS<23016726),]
Lof07_LG12inv_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG12" & Lof07_afreqs$POS>589105 & Lof07_afreqs$POS<13631347),]

allinvind=c(rownames(Lof07_LG1inv_afreqs),rownames(Lof07_LG2inv_afreqs),rownames(Lof07_LG7inv_afreqs),rownames(Lof07_LG12inv_afreqs))

Lof07_noinv_afreqs=Lof07_afreqs[!(rownames(Lof07_afreqs) %in% allinvind),]

covdf_noinv=na.omit(data.frame(Can40_noinv_afreqs$CHROM,Can40_noinv_afreqs$FreqChange,Lof07_noinv_afreqs$FreqChange11,Lof07_noinv_afreqs$FreqChange14))
cor11_noinv=cov(covdf_noinv$Can40_noinv_afreqs.FreqChange,covdf_noinv$Lof07_noinv_afreqs.FreqChange11)/sqrt(var(covdf_noinv$Can40_noinv_afreqs.FreqChange)*var(covdf_noinv$Lof07_noinv_afreqs.FreqChange11))
cor14_noinv=cov(covdf_noinv$Can40_noinv_afreqs.FreqChange,covdf_noinv$Lof07_noinv_afreqs.FreqChange14)/sqrt(var(covdf_noinv$Can40_noinv_afreqs.FreqChange)*var(covdf_noinv$Lof07_noinv_afreqs.FreqChange14))

cor_noinv_boot=data.frame(matrix(ncol=2,nrow=100))
colnames(cor_noinv_boot)=c("convcor1","convcor2")

for (i in 1:100) {
  covdfboot=covdf_noinv[sample(nrow(covdf_noinv),nrow(covdf_noinv),replace=TRUE),]
  colnames(covdfboot)=c("Can40_afreqs.CHROM","Can40_afreqs.FreqChange","Lof07_afreqs.FreqChange11","Lof07_afreqs.FreqChange14")
  cor_noinv_boot$convcor1[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
  cor_noinv_boot$convcor2[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
}

covdf_LG1inv=na.omit(data.frame(Can40_LG1inv_afreqs$CHROM,Can40_LG1inv_afreqs$FreqChange,Lof07_LG1inv_afreqs$FreqChange11,Lof07_LG1inv_afreqs$FreqChange14))
cor11_LG1inv=cov(covdf_LG1inv$Can40_LG1inv_afreqs.FreqChange,covdf_LG1inv$Lof07_LG1inv_afreqs.FreqChange11)/sqrt(var(covdf_LG1inv$Can40_LG1inv_afreqs.FreqChange)*var(covdf_LG1inv$Lof07_LG1inv_afreqs.FreqChange11))
cor14_LG1inv=cov(covdf_LG1inv$Can40_LG1inv_afreqs.FreqChange,covdf_LG1inv$Lof07_LG1inv_afreqs.FreqChange14)/sqrt(var(covdf_LG1inv$Can40_LG1inv_afreqs.FreqChange)*var(covdf_LG1inv$Lof07_LG1inv_afreqs.FreqChange14))

cor_LG1inv_boot=data.frame(matrix(ncol=2,nrow=100))
colnames(cor_LG1inv_boot)=c("convcor1","convcor2")

for (i in 1:100) {
  covdfboot=covdf_LG1inv[sample(nrow(covdf_LG1inv),nrow(covdf_LG1inv),replace=TRUE),]
  colnames(covdfboot)=c("Can40_afreqs.CHROM","Can40_afreqs.FreqChange","Lof07_afreqs.FreqChange11","Lof07_afreqs.FreqChange14")
  cor_LG1inv_boot$convcor1[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
  cor_LG1inv_boot$convcor2[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
}

covdf_LG2inv=na.omit(data.frame(Can40_LG2inv_afreqs$CHROM,Can40_LG2inv_afreqs$FreqChange,Lof07_LG2inv_afreqs$FreqChange11,Lof07_LG2inv_afreqs$FreqChange14))
cor11_LG2inv=cov(covdf_LG2inv$Can40_LG2inv_afreqs.FreqChange,covdf_LG2inv$Lof07_LG2inv_afreqs.FreqChange11)/sqrt(var(covdf_LG2inv$Can40_LG2inv_afreqs.FreqChange)*var(covdf_LG2inv$Lof07_LG2inv_afreqs.FreqChange11))
cor14_LG2inv=cov(covdf_LG2inv$Can40_LG2inv_afreqs.FreqChange,covdf_LG2inv$Lof07_LG2inv_afreqs.FreqChange14)/sqrt(var(covdf_LG2inv$Can40_LG2inv_afreqs.FreqChange)*var(covdf_LG2inv$Lof07_LG2inv_afreqs.FreqChange14))

cor_LG2inv_boot=data.frame(matrix(ncol=2,nrow=100))
colnames(cor_LG2inv_boot)=c("convcor1","convcor2")

for (i in 1:100) {
  covdfboot=covdf_LG2inv[sample(nrow(covdf_LG2inv),nrow(covdf_LG2inv),replace=TRUE),]
  colnames(covdfboot)=c("Can40_afreqs.CHROM","Can40_afreqs.FreqChange","Lof07_afreqs.FreqChange11","Lof07_afreqs.FreqChange14")
  cor_LG2inv_boot$convcor1[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
  cor_LG2inv_boot$convcor2[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
}

covdf_LG7inv=na.omit(data.frame(Can40_LG7inv_afreqs$CHROM,Can40_LG7inv_afreqs$FreqChange,Lof07_LG7inv_afreqs$FreqChange11,Lof07_LG7inv_afreqs$FreqChange14))
cor11_LG7inv=cov(covdf_LG7inv$Can40_LG7inv_afreqs.FreqChange,covdf_LG7inv$Lof07_LG7inv_afreqs.FreqChange11)/sqrt(var(covdf_LG7inv$Can40_LG7inv_afreqs.FreqChange)*var(covdf_LG7inv$Lof07_LG7inv_afreqs.FreqChange11))
cor14_LG7inv=cov(covdf_LG7inv$Can40_LG7inv_afreqs.FreqChange,covdf_LG7inv$Lof07_LG7inv_afreqs.FreqChange14)/sqrt(var(covdf_LG7inv$Can40_LG7inv_afreqs.FreqChange)*var(covdf_LG7inv$Lof07_LG7inv_afreqs.FreqChange14))

cor_LG7inv_boot=data.frame(matrix(ncol=2,nrow=100))
colnames(cor_LG7inv_boot)=c("convcor1","convcor2")

for (i in 1:100) {
  covdfboot=covdf_LG7inv[sample(nrow(covdf_LG7inv),nrow(covdf_LG7inv),replace=TRUE),]
  colnames(covdfboot)=c("Can40_afreqs.CHROM","Can40_afreqs.FreqChange","Lof07_afreqs.FreqChange11","Lof07_afreqs.FreqChange14")
  cor_LG7inv_boot$convcor1[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
  cor_LG7inv_boot$convcor2[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
}

covdf_LG12inv=na.omit(data.frame(Can40_LG12inv_afreqs$CHROM,Can40_LG12inv_afreqs$FreqChange,Lof07_LG12inv_afreqs$FreqChange11,Lof07_LG12inv_afreqs$FreqChange14))
cor11_LG12inv=cov(covdf_LG12inv$Can40_LG12inv_afreqs.FreqChange,covdf_LG12inv$Lof07_LG12inv_afreqs.FreqChange11)/sqrt(var(covdf_LG12inv$Can40_LG12inv_afreqs.FreqChange)*var(covdf_LG12inv$Lof07_LG12inv_afreqs.FreqChange11))
cor14_LG12inv=cov(covdf_LG12inv$Can40_LG12inv_afreqs.FreqChange,covdf_LG12inv$Lof07_LG12inv_afreqs.FreqChange14)/sqrt(var(covdf_LG12inv$Can40_LG12inv_afreqs.FreqChange)*var(covdf_LG12inv$Lof07_LG12inv_afreqs.FreqChange14))

cor_LG12inv_boot=data.frame(matrix(ncol=2,nrow=100))
colnames(cor_LG12inv_boot)=c("convcor1","convcor2")

for (i in 1:100) {
  covdfboot=covdf_LG12inv[sample(nrow(covdf_LG12inv),nrow(covdf_LG12inv),replace=TRUE),]
  colnames(covdfboot)=c("Can40_afreqs.CHROM","Can40_afreqs.FreqChange","Lof07_afreqs.FreqChange11","Lof07_afreqs.FreqChange14")
  cor_LG12inv_boot$convcor1[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
  cor_LG12inv_boot$convcor2[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
}

Can40_LG1noinv1_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG01" & Can40_afreqs$POS<9114741),]
Can40_LG1noinv2_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG01" & Can40_afreqs$POS>26192489),]
Can40_LG1noinv_afreqs=rbind(Can40_LG1noinv1_afreqs,Can40_LG1noinv2_afreqs)
Lof07_LG1noinv1_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG01" & Lof07_afreqs$POS<9114741),]
Lof07_LG1noinv2_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG01" & Lof07_afreqs$POS>26192489),]
Lof07_LG1noinv_afreqs=rbind(Lof07_LG1noinv1_afreqs,Lof07_LG1noinv2_afreqs)
covdf_LG1noinv=na.omit(data.frame(Can40_LG1noinv_afreqs$CHROM,Can40_LG1noinv_afreqs$FreqChange,Lof07_LG1noinv_afreqs$FreqChange11,Lof07_LG1noinv_afreqs$FreqChange14))
cor11_LG1noinv=cov(covdf_LG1noinv$Can40_LG1noinv_afreqs.FreqChange,covdf_LG1noinv$Lof07_LG1noinv_afreqs.FreqChange11)/sqrt(var(covdf_LG1noinv$Can40_LG1noinv_afreqs.FreqChange)*var(covdf_LG1noinv$Lof07_LG1noinv_afreqs.FreqChange11))
cor14_LG1noinv=cov(covdf_LG1noinv$Can40_LG1noinv_afreqs.FreqChange,covdf_LG1noinv$Lof07_LG1noinv_afreqs.FreqChange14)/sqrt(var(covdf_LG1noinv$Can40_LG1noinv_afreqs.FreqChange)*var(covdf_LG1noinv$Lof07_LG1noinv_afreqs.FreqChange14))

cor_LG1noinv_boot=data.frame(matrix(ncol=2,nrow=100))
colnames(cor_LG1noinv_boot)=c("convcor1","convcor2")

for (i in 1:100) {
  covdfboot=covdf_LG1noinv[sample(nrow(covdf_LG1noinv),nrow(covdf_LG1noinv),replace=TRUE),]
  colnames(covdfboot)=c("Can40_afreqs.CHROM","Can40_afreqs.FreqChange","Lof07_afreqs.FreqChange11","Lof07_afreqs.FreqChange14")
  cor_LG1noinv_boot$convcor1[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
  cor_LG1noinv_boot$convcor2[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
}

Can40_LG2noinv1_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG02" & Can40_afreqs$POS<18489307),]
Can40_LG2noinv2_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG02" & Can40_afreqs$POS>24050282),]
Can40_LG2noinv_afreqs=rbind(Can40_LG2noinv1_afreqs,Can40_LG2noinv2_afreqs)
Lof07_LG2noinv1_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG02" & Lof07_afreqs$POS<18489307),]
Lof07_LG2noinv2_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG02" & Lof07_afreqs$POS>24050282),]
Lof07_LG2noinv_afreqs=rbind(Lof07_LG2noinv1_afreqs,Lof07_LG2noinv2_afreqs)
covdf_LG2noinv=na.omit(data.frame(Can40_LG2noinv_afreqs$CHROM,Can40_LG2noinv_afreqs$FreqChange,Lof07_LG2noinv_afreqs$FreqChange11,Lof07_LG2noinv_afreqs$FreqChange14))
cor11_LG2noinv=cov(covdf_LG2noinv$Can40_LG2noinv_afreqs.FreqChange,covdf_LG2noinv$Lof07_LG2noinv_afreqs.FreqChange11)/sqrt(var(covdf_LG2noinv$Can40_LG2noinv_afreqs.FreqChange)*var(covdf_LG2noinv$Lof07_LG2noinv_afreqs.FreqChange11))
cor14_LG2noinv=cov(covdf_LG2noinv$Can40_LG2noinv_afreqs.FreqChange,covdf_LG2noinv$Lof07_LG2noinv_afreqs.FreqChange14)/sqrt(var(covdf_LG2noinv$Can40_LG2noinv_afreqs.FreqChange)*var(covdf_LG2noinv$Lof07_LG2noinv_afreqs.FreqChange14))

cor_LG2noinv_boot=data.frame(matrix(ncol=2,nrow=100))
colnames(cor_LG2noinv_boot)=c("convcor1","convcor2")

for (i in 1:100) {
  covdfboot=covdf_LG2noinv[sample(nrow(covdf_LG2noinv),nrow(covdf_LG2noinv),replace=TRUE),]
  colnames(covdfboot)=c("Can40_afreqs.CHROM","Can40_afreqs.FreqChange","Lof07_afreqs.FreqChange11","Lof07_afreqs.FreqChange14")
  cor_LG2noinv_boot$convcor1[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
  cor_LG2noinv_boot$convcor2[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
}

Can40_LG7noinv1_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG07" & Can40_afreqs$POS<18489307),]
Can40_LG7noinv2_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG07" & Can40_afreqs$POS>24050282),]
Can40_LG7noinv_afreqs=rbind(Can40_LG7noinv1_afreqs,Can40_LG7noinv2_afreqs)
Lof07_LG7noinv1_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG07" & Lof07_afreqs$POS<18489307),]
Lof07_LG7noinv2_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG07" & Lof07_afreqs$POS>24050282),]
Lof07_LG7noinv_afreqs=rbind(Lof07_LG7noinv1_afreqs,Lof07_LG7noinv2_afreqs)
covdf_LG7noinv=na.omit(data.frame(Can40_LG7noinv_afreqs$CHROM,Can40_LG7noinv_afreqs$FreqChange,Lof07_LG7noinv_afreqs$FreqChange11,Lof07_LG7noinv_afreqs$FreqChange14))
cor11_LG7noinv=cov(covdf_LG7noinv$Can40_LG7noinv_afreqs.FreqChange,covdf_LG7noinv$Lof07_LG7noinv_afreqs.FreqChange11)/sqrt(var(covdf_LG7noinv$Can40_LG7noinv_afreqs.FreqChange)*var(covdf_LG7noinv$Lof07_LG7noinv_afreqs.FreqChange11))
cor14_LG7noinv=cov(covdf_LG7noinv$Can40_LG7noinv_afreqs.FreqChange,covdf_LG7noinv$Lof07_LG7noinv_afreqs.FreqChange14)/sqrt(var(covdf_LG7noinv$Can40_LG7noinv_afreqs.FreqChange)*var(covdf_LG7noinv$Lof07_LG7noinv_afreqs.FreqChange14))

cor_LG7noinv_boot=data.frame(matrix(ncol=2,nrow=100))
colnames(cor_LG7noinv_boot)=c("convcor1","convcor2")

for (i in 1:100) {
  covdfboot=covdf_LG7noinv[sample(nrow(covdf_LG7noinv),nrow(covdf_LG7noinv),replace=TRUE),]
  colnames(covdfboot)=c("Can40_afreqs.CHROM","Can40_afreqs.FreqChange","Lof07_afreqs.FreqChange11","Lof07_afreqs.FreqChange14")
  cor_LG7noinv_boot$convcor1[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
  cor_LG7noinv_boot$convcor2[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
}

Can40_LG12noinv1_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG12" & Can40_afreqs$POS<18489307),]
Can40_LG12noinv2_afreqs=Can40_afreqs[which(Can40_afreqs$CHROM=="LG12" & Can40_afreqs$POS>24050282),]
Can40_LG12noinv_afreqs=rbind(Can40_LG12noinv1_afreqs,Can40_LG12noinv2_afreqs)
Lof07_LG12noinv1_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG12" & Lof07_afreqs$POS<18489307),]
Lof07_LG12noinv2_afreqs=Lof07_afreqs[which(Lof07_afreqs$CHROM=="LG12" & Lof07_afreqs$POS>24050282),]
Lof07_LG12noinv_afreqs=rbind(Lof07_LG12noinv1_afreqs,Lof07_LG12noinv2_afreqs)
covdf_LG12noinv=na.omit(data.frame(Can40_LG12noinv_afreqs$CHROM,Can40_LG12noinv_afreqs$FreqChange,Lof07_LG12noinv_afreqs$FreqChange11,Lof07_LG12noinv_afreqs$FreqChange14))
cor11_LG12noinv=cov(covdf_LG12noinv$Can40_LG12noinv_afreqs.FreqChange,covdf_LG12noinv$Lof07_LG12noinv_afreqs.FreqChange11)/sqrt(var(covdf_LG12noinv$Can40_LG12noinv_afreqs.FreqChange)*var(covdf_LG12noinv$Lof07_LG12noinv_afreqs.FreqChange11))
cor14_LG12noinv=cov(covdf_LG12noinv$Can40_LG12noinv_afreqs.FreqChange,covdf_LG12noinv$Lof07_LG12noinv_afreqs.FreqChange14)/sqrt(var(covdf_LG12noinv$Can40_LG12noinv_afreqs.FreqChange)*var(covdf_LG12noinv$Lof07_LG12noinv_afreqs.FreqChange14))

cor_LG12noinv_boot=data.frame(matrix(ncol=2,nrow=100))
colnames(cor_LG12noinv_boot)=c("convcor1","convcor2")

for (i in 1:100) {
  covdfboot=covdf_LG12noinv[sample(nrow(covdf_LG12noinv),nrow(covdf_LG12noinv),replace=TRUE),]
  colnames(covdfboot)=c("Can40_afreqs.CHROM","Can40_afreqs.FreqChange","Lof07_afreqs.FreqChange11","Lof07_afreqs.FreqChange14")
  cor_LG12noinv_boot$convcor1[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange11))
  cor_LG12noinv_boot$convcor2[i]=cov(covdfboot$Can40_afreqs.FreqChange,covdfboot$Lof07_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_afreqs.FreqChange)*var(covdfboot$Lof07_afreqs.FreqChange14))
}

##### SNPs in CDS

Can40_CDS_afreqs=read.table("data/afreqs_all_gadmor2_CDS/plink2.Can40.afreq",header=T)
Can13_CDS_afreqs=read.table("data/afreqs_all_gadmor2_CDS/plink2.Can13.afreq",header=T)
Can40_CDS_afreqs$ContFreq=Can13_CDS_afreqs$ALT_FREQS
Can40_CDS_afreqs$FreqChange=Can40_CDS_afreqs$ContFreq-Can40_CDS_afreqs$ALT_FREQS
Can40_CDS_afreqs$AbsFreqChange=abs(Can40_CDS_afreqs$ContFreq-Can40_CDS_afreqs$ALT_FREQS)

Lof07_CDS_afreqs=read.table("data/afreqs_all_gadmor2_CDS/plink2.Lof07.afreq",header=T)
Lof11_CDS_afreqs=read.table("data/afreqs_all_gadmor2_CDS/plink2.Lof11.afreq",header=T)
Lof14_CDS_afreqs=read.table("data/afreqs_all_gadmor2_CDS/plink2.Lof14.afreq",header=T)
Lof07_CDS_afreqs$Freq11=Lof11_CDS_afreqs$ALT_FREQS
Lof07_CDS_afreqs$Freq14=Lof14_CDS_afreqs$ALT_FREQS
Lof07_CDS_afreqs$FreqChange11=Lof07_CDS_afreqs$Freq11-Lof07_CDS_afreqs$ALT_FREQS
Lof07_CDS_afreqs$FreqChange14=Lof07_CDS_afreqs$Freq14-Lof07_CDS_afreqs$ALT_FREQS

covdf_CDS=data.frame(Can40_CDS_afreqs$CHROM,Can40_CDS_afreqs$FreqChange,Lof07_CDS_afreqs$FreqChange11,Lof07_CDS_afreqs$FreqChange14)
covdf_CDS$FreqChange1114=Lof07_CDS_afreqs$Freq14-Lof07_CDS_afreqs$Freq11
covdf_CDS=na.omit(covdf_CDS)

cor_CDS_Can40_Lof11=cov(covdf_CDS$Can40_CDS_afreqs.FreqChange,covdf_CDS$Lof07_CDS_afreqs.FreqChange11)/sqrt(var(covdf_CDS$Can40_CDS_afreqs.FreqChange)*var(covdf_CDS$Lof07_CDS_afreqs.FreqChange11))
cor_CDS_Can40_Lof14=cov(covdf_CDS$Can40_CDS_afreqs.FreqChange,covdf_CDS$Lof07_CDS_afreqs.FreqChange14)/sqrt(var(covdf_CDS$Can40_CDS_afreqs.FreqChange)*var(covdf_CDS$Lof07_CDS_afreqs.FreqChange14))
cor_CDS_Lof11_Lof14=cov(covdf_CDS$Lof07_CDS_afreqs.FreqChange11,covdf_CDS$Lof07_CDS_afreqs.FreqChange14)/sqrt(var(covdf_CDS$Lof07_CDS_afreqs.FreqChange11)*var(covdf_CDS$Lof07_CDS_afreqs.FreqChange14))
cor_CDS_Lof1114_Can=cov(covdf_CDS$FreqChange1114,covdf_CDS$Can40_CDS_afreqs.FreqChange)/sqrt(var(covdf_CDS$FreqChange1114)*var(covdf_CDS$Can40_CDS_afreqs.FreqChange))

### bootstrap CDS

cor_CDS_boot=data.frame(matrix(ncol=4,nrow=100))
colnames(cor_CDS_boot)=c("convcor1","convcor2","convcor3","convcor4")

for (i in 1:100) {
  covdfboot=covdf_CDS[sample(nrow(covdf_CDS),nrow(covdf_CDS),replace=TRUE),]
  cor_CDS_boot$convcor1[i]=cov(covdfboot$Can40_CDS_afreqs.FreqChange,covdfboot$Lof07_CDS_afreqs.FreqChange11)/sqrt(var(covdfboot$Can40_CDS_afreqs.FreqChange)*var(covdfboot$Lof07_CDS_afreqs.FreqChange11))
  cor_CDS_boot$convcor2[i]=cov(covdfboot$Can40_CDS_afreqs.FreqChange,covdfboot$Lof07_CDS_afreqs.FreqChange14)/sqrt(var(covdfboot$Can40_CDS_afreqs.FreqChange)*var(covdfboot$Lof07_CDS_afreqs.FreqChange14))
  cor_CDS_boot$convcor3[i]=cov(covdfboot$Lof07_CDS_afreqs.FreqChange11,covdfboot$Lof07_CDS_afreqs.FreqChange14)/sqrt(var(covdfboot$Lof07_CDS_afreqs.FreqChange11)*var(covdfboot$Lof07_CDS_afreqs.FreqChange14))
  cor_CDS_boot$convcor4[i]=cov(covdfboot$FreqChange1114,covdfboot$Can40_CDS_afreqs.FreqChange)/sqrt(var(covdfboot$FreqChange1114)*var(covdfboot$Can40_CDS_afreqs.FreqChange))
}

#### plotting subsets

par(las=2,mai=c(2,1,1,1))

par(las=2,mar=c(10.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(x=NULL,y=NULL,xlim=c(1.3,8.2),ylim=c(-0.03,0.32),xlab="",ylab="ConvCor",xaxt="n")
points(x=1.5-0.1,y=cor11_LG1inv,col=cols[2],pch=19)
points(x=1.5+0.1,y=cor14_LG1inv,col=cols[3],pch=19)
points(x=2-0.1,y=cor11_LG1noinv,col=cols[2],pch=19)
points(x=2+0.1,y=cor14_LG1noinv,col=cols[3],pch=19)
points(x=3-0.1,y=cor11_LG2inv,col=cols[2],pch=19)
points(x=3+0.1,y=cor14_LG2inv,col=cols[3],pch=19)
points(x=3.5-0.1,y=cor11_LG2noinv,col=cols[2],pch=19)
points(x=3.5+0.1,y=cor14_LG2noinv,col=cols[3],pch=19)
points(x=4.5-0.1,y=cor11_LG7inv,col=cols[2],pch=19)
points(x=4.5+0.1,y=cor14_LG7inv,col=cols[3],pch=19)
points(x=5-0.1,y=cor11_LG7noinv,col=cols[2],pch=19)
points(x=5+0.1,y=cor14_LG7noinv,col=cols[3],pch=19)
points(x=6-0.1,y=cor11_LG12inv,col=cols[2],pch=19)
points(x=6+0.1,y=cor14_LG12inv,col=cols[3],pch=19)
points(x=6.5-0.1,y=cor11_LG12noinv,col=cols[2],pch=19)
points(x=6.5+0.1,y=cor14_LG12noinv,col=cols[3],pch=19)
points(x=7.5-0.1,y=cor_CDS_Can40_Lof11,col=cols[2],pch=19)
points(x=7.5+0.1,y=cor_CDS_Can40_Lof14,col=cols[3],pch=19)
points(x=8-0.1,y=cor_Can40_Lof11,col=cols[2],pch=19)
points(x=8+0.1,y=cor_Can40_Lof14,col=cols[3],pch=19)
lines(x=c(1.5-0.1,1.5-0.1),y=c(quantile(cor_LG1inv_boot$convcor1,0.025),quantile(cor_LG1inv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(1.5+0.1,1.5+0.1),y=c(quantile(cor_LG1inv_boot$convcor2,0.025),quantile(cor_LG1inv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(2-0.1,2-0.1),y=c(quantile(cor_LG1noinv_boot$convcor1,0.025),quantile(cor_LG1noinv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(2+0.1,2+0.1),y=c(quantile(cor_LG1noinv_boot$convcor2,0.025),quantile(cor_LG1noinv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(3-0.1,3-0.1),y=c(quantile(cor_LG2inv_boot$convcor1,0.025),quantile(cor_LG2inv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(3+0.1,3+0.1),y=c(quantile(cor_LG2inv_boot$convcor2,0.025),quantile(cor_LG2inv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(3.5-0.1,3.5-0.1),y=c(quantile(cor_LG2noinv_boot$convcor1,0.025),quantile(cor_LG2noinv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(3.5+0.1,3.5+0.1),y=c(quantile(cor_LG2noinv_boot$convcor2,0.025),quantile(cor_LG2noinv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(4.5-0.1,4.5-0.1),y=c(quantile(cor_LG7inv_boot$convcor1,0.025),quantile(cor_LG7inv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(4.5+0.1,4.5+0.1),y=c(quantile(cor_LG7inv_boot$convcor2,0.025),quantile(cor_LG7inv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(5-0.1,5-0.1),y=c(quantile(cor_LG7noinv_boot$convcor1,0.025),quantile(cor_LG7noinv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(5+0.1,5+0.1),y=c(quantile(cor_LG7noinv_boot$convcor2,0.025),quantile(cor_LG7noinv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(6-0.1,6-0.1),y=c(quantile(cor_LG12inv_boot$convcor1,0.025),quantile(cor_LG12inv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(6+0.1,6+0.1),y=c(quantile(cor_LG12inv_boot$convcor2,0.025),quantile(cor_LG12inv_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(6.5-0.1,6.5-0.1),y=c(quantile(cor_LG12noinv_boot$convcor1,0.025),quantile(cor_LG12noinv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(6.5+0.1,6.5+0.1),y=c(quantile(cor_LG12noinv_boot$convcor2,0.025),quantile(cor_LG12noinv_boot$convcor2,0.975)),col=cols[3],lwd=2)

lines(x=c(7.5-0.1,7.5-0.1),y=c(quantile(cor_CDS_boot$convcor1,0.025),quantile(cor_CDS_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(7.5+0.1,7.5+0.1),y=c(quantile(cor_CDS_boot$convcor2,0.025),quantile(cor_CDS_boot$convcor2,0.975)),col=cols[3],lwd=2)
lines(x=c(8-0.1,8-0.1),y=c(quantile(cor_all_boot$convcor1,0.025),quantile(cor_LG7inv_boot$convcor1,0.975)),col=cols[2],lwd=2)
lines(x=c(8+0.1,8+0.1),y=c(quantile(cor_all_boot$convcor2,0.025),quantile(cor_LG7inv_boot$convcor2,0.975)),col=cols[3],lwd=2)

lines(x=c(1.3,8.2),y=c(0,0),lty=2)
legend("topright",bty="n",inset=c(-0.2,0),legend=c(expression('ConvCor'[1]),expression('ConvCor'[2])),col=c(cols[2],cols[3]),pch=c(19,19))
axis(side=1,at=c(1.5,2,3,3.5,4.5,5,6,6.5,7.5,8),label=c("LG01 inversion","LG01 outside inversion","LG02 inversion","LG02 outside inversion","LG07 inversion","LG07 outside inversion","LG12 inversion","LG12 outside inversion","Coding SNPs","All SNPs"))


### loci segregating across populations
Can40_afreqs=read.table("data/Can40Loci_afreqs/plink2.Can40.afreq",header=T)
Can13_afreqs=read.table("data/Can40Loci_afreqs/plink2.Can13.afreq",header=T)
Lof07_afreqs=read.table("data/Can40Loci_afreqs/plink2.Lof07.afreq",header=T)
Lof11_afreqs=read.table("data/Can40Loci_afreqs/plink2.Lof11.afreq",header=T)
Lof14_afreqs=read.table("data/Can40Loci_afreqs/plink2.Lof14.afreq",header=T)

head(Can40_afreqs)

which(Can40_afreqs==0)

length(which(Can40_afreqs$ALT_FREQS==0))
length(which(Can40_afreqs$ALT_FREQS>0 & Can40_afreqs$ALT_FREQS<1))
length(which(Can40_afreqs$ALT_FREQS==1))
length(which(is.na(Can40_afreqs$ALT_FREQS)))

segCan40=which(Can40_afreqs$ALT_FREQS>0 & Can40_afreqs$ALT_FREQS<1)
segCan13=which(Can13_afreqs$ALT_FREQS>0 & Can13_afreqs$ALT_FREQS<1)
segLof07=which(Lof07_afreqs$ALT_FREQS>0 & Lof07_afreqs$ALT_FREQS<1)
segLof11=which(Lof11_afreqs$ALT_FREQS>0 & Lof11_afreqs$ALT_FREQS<1)
segLof14=which(Lof14_afreqs$ALT_FREQS>0 & Lof14_afreqs$ALT_FREQS<1)

length(segCan40)
length(segCan13)
length(segLof07)
length(segLof11)
length(segLof14)

length(intersect(segCan40,segCan13))
length(intersect(segLof07,segLof11))
length(intersect(segLof07,segLof14))

length(intersect(segCan40,segLof07))/nrow(Can40_afreqs)
length(intersect(segCan13,segLof11))/nrow(Can40_afreqs)
length(intersect(segCan13,segLof14))/nrow(Can40_afreqs)

length(Reduce(intersect,list(segCan13,segLof11,segLof14)))/nrow(Can40_afreqs)


length(Reduce(intersect,list(segCan40,segCan13,segLof07,segLof11,segLof14)))
length(Reduce(intersect,list(segCan40,segCan13,segLof07,segLof11,segLof14)))/nrow(Can40_afreqs)
