library(tidyverse)
library(ggpubr)
library(RColorBrewer)
cols=brewer.pal(name="Set1",n=5)

##### load in the datasets for each migration scenario and calculate ConvCor statistics

covdf_m2n=data.frame(m=rep(NA,20),q=rep(NA,20),it=rep(NA,20),Can40Lof11=rep(NA,20),Can40Lof14=rep(NA,20),Lof11Lof14=rep(NA,20),Can40Lof1114=rep(NA,20))

for (i in 1:20){
  Can13file=paste("slim_outputs/slimout_afreqs/K7000_m1e-2_neutral_maf/i",i,"_maf.Can13.afreq",sep="")
  Can40file=paste("slim_outputs/slimout_afreqs/K7000_m1e-2_neutral_maf/i",i,"_maf.Can40.afreq",sep="")
  Lof07file=paste("slim_outputs/slimout_afreqs/K7000_m1e-2_neutral_maf/i",i,"_maf.Lof07.afreq",sep="")
  Lof11file=paste("slim_outputs/slimout_afreqs/K7000_m1e-2_neutral_maf/i",i,"_maf.Lof11.afreq",sep="")
  Lof14file=paste("slim_outputs/slimout_afreqs/K7000_m1e-2_neutral_maf/i",i,"_maf.Lof14.afreq",sep="")
  Can13_afreqs=read.table(Can13file,header=F)
  Can40_afreqs=read.table(Can40file,header=F)
  Lof07_afreqs=read.table(Lof07file,header=F)
  Lof11_afreqs=read.table(Lof11file,header=F)
  Lof14_afreqs=read.table(Lof14file,header=F)
  FreqChange=data.frame(Can40_afreqs$V5,Can13_afreqs$V5,Lof07_afreqs$V5,Lof11_afreqs$V5,Lof14_afreqs$V5)
  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(x, is.nan))
  FreqChange[is.nan(FreqChange)] <- 0
  FreqChange$Can=FreqChange$Can13_afreqs.V5-FreqChange$Can40_afreqs.V5
  FreqChange$Lof11=FreqChange$Lof11_afreqs.V5-FreqChange$Lof07_afreqs.V5
  FreqChange$Lof14=FreqChange$Lof14_afreqs.V5-FreqChange$Lof07_afreqs.V5
  FreqChange$Lof1114=FreqChange$Lof14_afreqs.V5-FreqChange$Lof11_afreqs.V5
  covdf_m2n$m[i]=0.01
  covdf_m2n$q[i]=0
  covdf_m2n$it[i]=i
  covdf_m2n$Can40Lof11[i]=cov(FreqChange$Can,FreqChange$Lof11)/sqrt(var(FreqChange$Can)*var(FreqChange$Lof11))
  covdf_m2n$Can40Lof14[i]=cov(FreqChange$Can,FreqChange$Lof14)/sqrt(var(FreqChange$Can)*var(FreqChange$Lof14))
  covdf_m2n$Lof11Lof14[i]=cov(FreqChange$Lof11,FreqChange$Lof14)/sqrt(var(FreqChange$Lof11)*var(FreqChange$Lof14))
  covdf_m2n$Can40Lof1114[i]=cov(FreqChange$Can,FreqChange$Lof1114)/sqrt(var(FreqChange$Can)*var(FreqChange$Lof1114))
}


covdf_m3n=data.frame(m=rep(NA,20),q=rep(NA,20),it=rep(NA,20),Can40Lof11=rep(NA,20),Can40Lof14=rep(NA,20),Lof11Lof14=rep(NA,20),Can40Lof1114=rep(NA,20))

for (i in 1:20){
  Can13file=paste("slim_outputs/slimout_afreqs/K7000_m1e-3_neutral_maf/i",i,"_maf.Can13.afreq",sep="")
  Can40file=paste("slim_outputs/slimout_afreqs/K7000_m1e-3_neutral_maf/i",i,"_maf.Can40.afreq",sep="")
  Lof07file=paste("slim_outputs/slimout_afreqs/K7000_m1e-3_neutral_maf/i",i,"_maf.Lof07.afreq",sep="")
  Lof11file=paste("slim_outputs/slimout_afreqs/K7000_m1e-3_neutral_maf/i",i,"_maf.Lof11.afreq",sep="")
  Lof14file=paste("slim_outputs/slimout_afreqs/K7000_m1e-3_neutral_maf/i",i,"_maf.Lof14.afreq",sep="")
  Can13_afreqs=read.table(Can13file,header=F)
  Can40_afreqs=read.table(Can40file,header=F)
  Lof07_afreqs=read.table(Lof07file,header=F)
  Lof11_afreqs=read.table(Lof11file,header=F)
  Lof14_afreqs=read.table(Lof14file,header=F)
  FreqChange=data.frame(Can40_afreqs$V5,Can13_afreqs$V5,Lof07_afreqs$V5,Lof11_afreqs$V5,Lof14_afreqs$V5)
  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(x, is.nan))
  FreqChange[is.nan(FreqChange)] <- 0
  FreqChange$Can=FreqChange$Can13_afreqs.V5-FreqChange$Can40_afreqs.V5
  FreqChange$Lof11=FreqChange$Lof11_afreqs.V5-FreqChange$Lof07_afreqs.V5
  FreqChange$Lof14=FreqChange$Lof14_afreqs.V5-FreqChange$Lof07_afreqs.V5
  FreqChange$Lof1114=FreqChange$Lof14_afreqs.V5-FreqChange$Lof11_afreqs.V5
  covdf_m3n$m[i]=0.001
  covdf_m3n$q[i]=0
  covdf_m3n$it[i]=i
  covdf_m3n$Can40Lof11[i]=cov(FreqChange$Can,FreqChange$Lof11)/sqrt(var(FreqChange$Can)*var(FreqChange$Lof11))
  covdf_m3n$Can40Lof14[i]=cov(FreqChange$Can,FreqChange$Lof14)/sqrt(var(FreqChange$Can)*var(FreqChange$Lof14))
  covdf_m3n$Lof11Lof14[i]=cov(FreqChange$Lof11,FreqChange$Lof14)/sqrt(var(FreqChange$Lof11)*var(FreqChange$Lof14))
  covdf_m3n$Can40Lof1114[i]=cov(FreqChange$Can,FreqChange$Lof1114)/sqrt(var(FreqChange$Can)*var(FreqChange$Lof1114))
}


covdf_m4n=data.frame(m=rep(NA,20),q=rep(NA,20),it=rep(NA,20),Can40Lof11=rep(NA,20),Can40Lof14=rep(NA,20),Lof11Lof14=rep(NA,20),Can40Lof1114=rep(NA,20))

for (i in 1:20){
  Can13file=paste("slim_outputs/slimout_afreqs/K7000_m1e-4_neutral_maf/i",i,"_maf.Can13.afreq",sep="")
  Can40file=paste("slim_outputs/slimout_afreqs/K7000_m1e-4_neutral_maf/i",i,"_maf.Can40.afreq",sep="")
  Lof07file=paste("slim_outputs/slimout_afreqs/K7000_m1e-4_neutral_maf/i",i,"_maf.Lof07.afreq",sep="")
  Lof11file=paste("slim_outputs/slimout_afreqs/K7000_m1e-4_neutral_maf/i",i,"_maf.Lof11.afreq",sep="")
  Lof14file=paste("slim_outputs/slimout_afreqs/K7000_m1e-4_neutral_maf/i",i,"_maf.Lof14.afreq",sep="")
  Can13_afreqs=read.table(Can13file,header=F)
  Can40_afreqs=read.table(Can40file,header=F)
  Lof07_afreqs=read.table(Lof07file,header=F)
  Lof11_afreqs=read.table(Lof11file,header=F)
  Lof14_afreqs=read.table(Lof14file,header=F)
  FreqChange=data.frame(Can40_afreqs$V5,Can13_afreqs$V5,Lof07_afreqs$V5,Lof11_afreqs$V5,Lof14_afreqs$V5)
  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(x, is.nan))
  FreqChange[is.nan(FreqChange)] <- 0
  FreqChange$Can=FreqChange$Can13_afreqs.V5-FreqChange$Can40_afreqs.V5
  FreqChange$Lof11=FreqChange$Lof11_afreqs.V5-FreqChange$Lof07_afreqs.V5
  FreqChange$Lof14=FreqChange$Lof14_afreqs.V5-FreqChange$Lof07_afreqs.V5
  FreqChange$Lof1114=FreqChange$Lof14_afreqs.V5-FreqChange$Lof11_afreqs.V5
  covdf_m4n$m[i]=0.0001
  covdf_m4n$q[i]=0
  covdf_m4n$it[i]=i
  covdf_m4n$Can40Lof11[i]=cov(FreqChange$Can,FreqChange$Lof11)/sqrt(var(FreqChange$Can)*var(FreqChange$Lof11))
  covdf_m4n$Can40Lof14[i]=cov(FreqChange$Can,FreqChange$Lof14)/sqrt(var(FreqChange$Can)*var(FreqChange$Lof14))
  covdf_m4n$Lof11Lof14[i]=cov(FreqChange$Lof11,FreqChange$Lof14)/sqrt(var(FreqChange$Lof11)*var(FreqChange$Lof14))
  covdf_m4n$Can40Lof1114[i]=cov(FreqChange$Can,FreqChange$Lof1114)/sqrt(var(FreqChange$Can)*var(FreqChange$Lof1114))
}


covdf_m5n=data.frame(m=rep(NA,20),q=rep(NA,20),it=rep(NA,20),Can40Lof11=rep(NA,20),Can40Lof14=rep(NA,20),Lof11Lof14=rep(NA,20),Can40Lof1114=rep(NA,20))

for (i in 1:20){
  Can13file=paste("slim_outputs/slimout_afreqs/K7000_m1e-5_neutral_maf/i",i,"_maf.Can13.afreq",sep="")
  Can40file=paste("slim_outputs/slimout_afreqs/K7000_m1e-5_neutral_maf/i",i,"_maf.Can40.afreq",sep="")
  Lof07file=paste("slim_outputs/slimout_afreqs/K7000_m1e-5_neutral_maf/i",i,"_maf.Lof07.afreq",sep="")
  Lof11file=paste("slim_outputs/slimout_afreqs/K7000_m1e-5_neutral_maf/i",i,"_maf.Lof11.afreq",sep="")
  Lof14file=paste("slim_outputs/slimout_afreqs/K7000_m1e-5_neutral_maf/i",i,"_maf.Lof14.afreq",sep="")
  Can13_afreqs=read.table(Can13file,header=F)
  Can40_afreqs=read.table(Can40file,header=F)
  Lof07_afreqs=read.table(Lof07file,header=F)
  Lof11_afreqs=read.table(Lof11file,header=F)
  Lof14_afreqs=read.table(Lof14file,header=F)
  FreqChange=data.frame(Can40_afreqs$V5,Can13_afreqs$V5,Lof07_afreqs$V5,Lof11_afreqs$V5,Lof14_afreqs$V5)
  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(x, is.nan))
  FreqChange[is.nan(FreqChange)] <- 0
  FreqChange$Can=FreqChange$Can13_afreqs.V5-FreqChange$Can40_afreqs.V5
  FreqChange$Lof11=FreqChange$Lof11_afreqs.V5-FreqChange$Lof07_afreqs.V5
  FreqChange$Lof14=FreqChange$Lof14_afreqs.V5-FreqChange$Lof07_afreqs.V5
  FreqChange$Lof1114=FreqChange$Lof14_afreqs.V5-FreqChange$Lof11_afreqs.V5
  covdf_m5n$m[i]=0.00001
  covdf_m5n$q[i]=0
  covdf_m5n$it[i]=i
  covdf_m5n$Can40Lof11[i]=cov(FreqChange$Can,FreqChange$Lof11)/sqrt(var(FreqChange$Can)*var(FreqChange$Lof11))
  covdf_m5n$Can40Lof14[i]=cov(FreqChange$Can,FreqChange$Lof14)/sqrt(var(FreqChange$Can)*var(FreqChange$Lof14))
  covdf_m5n$Lof11Lof14[i]=cov(FreqChange$Lof11,FreqChange$Lof14)/sqrt(var(FreqChange$Lof11)*var(FreqChange$Lof14))
  covdf_m5n$Can40Lof1114[i]=cov(FreqChange$Can,FreqChange$Lof1114)/sqrt(var(FreqChange$Can)*var(FreqChange$Lof1114))
}


covdf_m6n=data.frame(m=rep(NA,20),q=rep(NA,20),it=rep(NA,20),Can40Lof11=rep(NA,20),Can40Lof14=rep(NA,20),Lof11Lof14=rep(NA,20),Can40Lof1114=rep(NA,20))

for (i in 1:20){
  Can13file=paste("slim_outputs/slimout_afreqs/K7000_m1e-6_neutral_maf/i",i,"_maf.Can13.afreq",sep="")
  Can40file=paste("slim_outputs/slimout_afreqs/K7000_m1e-6_neutral_maf/i",i,"_maf.Can40.afreq",sep="")
  Lof07file=paste("slim_outputs/slimout_afreqs/K7000_m1e-6_neutral_maf/i",i,"_maf.Lof07.afreq",sep="")
  Lof11file=paste("slim_outputs/slimout_afreqs/K7000_m1e-6_neutral_maf/i",i,"_maf.Lof11.afreq",sep="")
  Lof14file=paste("slim_outputs/slimout_afreqs/K7000_m1e-6_neutral_maf/i",i,"_maf.Lof14.afreq",sep="")
  Can13_afreqs=read.table(Can13file,header=F)
  Can40_afreqs=read.table(Can40file,header=F)
  Lof07_afreqs=read.table(Lof07file,header=F)
  Lof11_afreqs=read.table(Lof11file,header=F)
  Lof14_afreqs=read.table(Lof14file,header=F)
  FreqChange=data.frame(Can40_afreqs$V5,Can13_afreqs$V5,Lof07_afreqs$V5,Lof11_afreqs$V5,Lof14_afreqs$V5)
  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(x, is.nan))
  FreqChange[is.nan(FreqChange)] <- 0
  FreqChange$Can=FreqChange$Can13_afreqs.V5-FreqChange$Can40_afreqs.V5
  FreqChange$Lof11=FreqChange$Lof11_afreqs.V5-FreqChange$Lof07_afreqs.V5
  FreqChange$Lof14=FreqChange$Lof14_afreqs.V5-FreqChange$Lof07_afreqs.V5
  FreqChange$Lof1114=FreqChange$Lof14_afreqs.V5-FreqChange$Lof11_afreqs.V5
  covdf_m6n$m[i]=0.000001
  covdf_m6n$q[i]=0
  covdf_m6n$it[i]=i
  covdf_m6n$Can40Lof11[i]=cov(FreqChange$Can,FreqChange$Lof11)/sqrt(var(FreqChange$Can)*var(FreqChange$Lof11))
  covdf_m6n$Can40Lof14[i]=cov(FreqChange$Can,FreqChange$Lof14)/sqrt(var(FreqChange$Can)*var(FreqChange$Lof14))
  covdf_m6n$Lof11Lof14[i]=cov(FreqChange$Lof11,FreqChange$Lof14)/sqrt(var(FreqChange$Lof11)*var(FreqChange$Lof14))
  covdf_m6n$Can40Lof1114[i]=cov(FreqChange$Can,FreqChange$Lof1114)/sqrt(var(FreqChange$Can)*var(FreqChange$Lof1114))
}

### combine migration scenarios

covdf_n=rbind(covdf_m2n,covdf_m3n,covdf_m4n,covdf_m5n,covdf_m6n)

plot(log10(covdf_n$m),covdf_n$Can40Lof11,xlab="Log10(Migration Rate)",ylab="Convergence Correlation",main="Canada 1940-2013 vs Norway 1907-2011")
plot(log10(covdf_n$m),covdf_n$Can40Lof14,xlab="Log10(Migration Rate)",ylab="Convergence Correlation",main="Canada 1940-2013 vs Norway 1907-2014")
plot(log10(covdf_n$m),covdf_n$Can40Lof1114,xlab="Log10(Migration Rate)",ylab="Convergence Correlation",main="Canada 1940-2013 vs Norway 2011-2014")
plot(log10(covdf_n$m),covdf_n$Lof11Lof14,xlab="Log10(Migration Rate)",ylab="Convergence Correlation",main="Norway 1907-2011 vs Norway 1907-2014")

mean(covdf_n$Can40Lof11)
mean(covdf_n$Can40Lof14)

## migration on log scale

covdf_n$logm=as.factor(log10(covdf_n$m))

## plot convcor

convcor1plot=ggplot(data=covdf_n,aes(x=logm,y=Can40Lof11)) + 
  geom_boxplot(color=cols[2]) +
  ylim(-0.08,0.15) +
  xlab(expression('Log'[10]*' Migration Rate')) +
  ylab(expression('ConvCor'[1])) +
  geom_hline(yintercept=0.085,linetype="dashed") +
  annotate("rect",xmin=0,xmax=6,ymin=quantile(cor_all_boot$convcor1,0.025),ymax=quantile(cor_all_boot$convcor1,0.975),alpha=0.4,fill=cols[2]) +
  theme_bw()

convcor2plot=ggplot(data=covdf_n,aes(x=logm,y=Can40Lof14)) + 
  geom_boxplot(color=cols[3]) +
  ylim(-0.08,0.15) +
  xlab(expression('Log'[10]*' Migration Rate')) +
  ylab(expression('ConvCor'[2])) +
  geom_hline(yintercept=0.082,linetype="dashed") +
  annotate("rect",xmin=0,xmax=6,ymin=quantile(cor_all_boot$convcor2,0.025),ymax=quantile(cor_all_boot$convcor2,0.975),alpha=0.4,fill=cols[3]) +
  theme_bw()

convcor3plot=ggplot(data=covdf_n,aes(x=logm,y=Lof11Lof14)) + 
  geom_boxplot(color=cols[4]) +
  ylim(0.38,0.62) +
  xlab(expression('Log'[10]*' Migration Rate')) +
  ylab(expression('ConvCor'[3])) +
  geom_hline(yintercept=0.519,linetype="dashed") +
  annotate("rect",xmin=0,xmax=6,ymin=quantile(cor_all_boot$convcor3,0.025),ymax=quantile(cor_all_boot$convcor3,0.975),alpha=0.4,fill=cols[4]) +
  theme_bw()

convcor4plot=ggplot(data=covdf_n,aes(x=logm,y=Can40Lof1114)) + 
  geom_boxplot(color=cols[5]) +
  ylim(-0.1,0.1) +
  xlab(expression('Log'[10]*' Migration Rate')) +
  ylab(expression('ConvCor'[4])) +
  geom_hline(yintercept=-0.026,linetype="dashed") +
  annotate("rect",xmin=0,xmax=6,ymin=quantile(cor_all_boot$convcor4,0.025),ymax=quantile(cor_all_boot$convcor4,0.975),alpha=0.4,fill=cols[5]) +
  theme_bw()

ggarrange(convcor1plot,convcor2plot,convcor3plot,convcor4plot,ncol=2,nrow=2)


### just CIs ####

convcor1plot=ggplot(data=covdf_n,aes(x=logm,y=Can40Lof11)) + 
#  geom_boxplot(color=cols[2]) +
  ylim(-0.08,0.15) +
  scale_x_discrete(labels=c("-6","-5","-4","-3","-2")) +
  xlab(expression('Log'[10]*' Migration Rate')) +
  ylab(expression('ConvCor'[1])) +
  geom_hline(yintercept=0.085,linetype="dashed") +
  annotate("rect",xmin=0,xmax=6,ymin=quantile(cor_all_boot$convcor1,0.025),ymax=quantile(cor_all_boot$convcor1,0.975),alpha=0.4,fill=cols[2]) +
  theme_bw()

convcor2plot=ggplot(data=covdf_n,aes(x=logm,y=Can40Lof14)) + 
#  geom_boxplot(color=cols[3]) +
  ylim(-0.08,0.15) +
  scale_x_discrete(labels=c("-6","-5","-4","-3","-2")) +
  xlab(expression('Log'[10]*' Migration Rate')) +
  ylab(expression('ConvCor'[2])) +
  geom_hline(yintercept=0.082,linetype="dashed") +
  annotate("rect",xmin=0,xmax=6,ymin=quantile(cor_all_boot$convcor2,0.025),ymax=quantile(cor_all_boot$convcor2,0.975),alpha=0.4,fill=cols[3]) +
  theme_bw()

convcor3plot=ggplot(data=covdf_n,aes(x=logm,y=Lof11Lof14)) + 
#  geom_boxplot(color=cols[4]) +
  ylim(0.38,0.62) +
  scale_x_discrete(labels=c("-6","-5","-4","-3","-2")) +
  xlab(expression('Log'[10]*' Migration Rate')) +
  ylab(expression('ConvCor'[3])) +
  geom_hline(yintercept=0.519,linetype="dashed") +
  annotate("rect",xmin=0,xmax=6,ymin=quantile(cor_all_boot$convcor3,0.025),ymax=quantile(cor_all_boot$convcor3,0.975),alpha=0.4,fill=cols[4]) +
  theme_bw()

convcor4plot=ggplot(data=covdf_n,aes(x=logm,y=Can40Lof1114)) + 
#  geom_boxplot(color=cols[5]) +
  ylim(-0.1,0.1) +
  scale_x_discrete(labels=c("-6","-5","-4","-3","-2")) +
  xlab(expression('Log'[10]*' Migration Rate')) +
  ylab(expression('ConvCor'[4])) +
  geom_hline(yintercept=-0.026,linetype="dashed") +
  annotate("rect",xmin=0,xmax=6,ymin=quantile(cor_all_boot$convcor4,0.025),ymax=quantile(cor_all_boot$convcor4,0.975),alpha=0.4,fill=cols[5]) +
  theme_bw()

ggarrange(convcor1plot,convcor2plot,convcor3plot,convcor4plot,ncol=2,nrow=2)
