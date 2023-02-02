library(ggplot2)
library(ggpubr)

## read simulated FST datasets for each migration scenario

mig2_n<-as.matrix(read.table("slim_outputs/slimout_fsts/K7000_m1e-2_neutral_maf_fsts"))
mig2_n_d<-data.frame(t(mig2_n))
names(mig2_n_d)<-mig2_n_d[1,]
mig2_n_d$m=0.01
mig2_n_d$qtl=0
mig2_n_d$I=seq(0,20)
mig2_n_d=mig2_n_d[-1,]

mig3_n<-as.matrix(read.table("slim_outputs/slimout_fsts/K7000_m1e-3_neutral_maf_fsts"))
mig3_n_d<-data.frame(t(mig3_n))
names(mig3_n_d)<-mig3_n_d[1,]
mig3_n_d$m=0.001
mig3_n_d$qtl=0
mig3_n_d$I=seq(0,20)
mig3_n_d=mig3_n_d[-1,]

mig4_n<-as.matrix(read.table("slim_outputs/slimout_fsts/K7000_m1e-4_neutral_maf_fsts"))
mig4_n_d<-data.frame(t(mig4_n))
names(mig4_n_d)<-mig4_n_d[1,]
mig4_n_d$m=0.0001
mig4_n_d$qtl=0
mig4_n_d$I=seq(0,20)
mig4_n_d=mig4_n_d[-1,]

mig5_n<-as.matrix(read.table("slim_outputs/slimout_fsts/K7000_m1e-5_neutral_maf_fsts"))
mig5_n_d<-data.frame(t(mig5_n))
names(mig5_n_d)<-mig5_n_d[1,]
mig5_n_d$m=0.00001
mig5_n_d$qtl=0
mig5_n_d$I=seq(0,20)
mig5_n_d=mig5_n_d[-1,]

mig6_n<-as.matrix(read.table("slim_outputs/slimout_fsts/K7000_m1e-6_neutral_maf_fsts"))
mig6_n_d<-data.frame(t(mig6_n))
names(mig6_n_d)<-mig6_n_d[1,]
mig6_n_d$m=0.000001
mig6_n_d$qtl=0
mig6_n_d$I=seq(0,20)
mig6_n_d=mig6_n_d[-1,]

## combine migration scenarios and create isolate temporal / spatial comparisons

mig_n<-rbind(mig2_n_d,mig3_n_d,mig4_n_d,mig5_n_d,mig6_n_d)

tmig=data.frame(fst=c(mig_n$ctfst,mig_n$ltfst),mig=c(mig_n$m,mig_n$m))
tmig$logm=log10(tmig$mig)
tmig$cm=as.factor(tmig$logm)
tmig$fst=as.numeric(tmig$fst)

smig=data.frame(fst=c(mig_n$clcfst,mig_n$clhfst),mig=c(mig_n$m,mig_n$m))
smig$logm=log10(smig$mig)
smig$cm=as.factor(smig$logm)
smig$fst=as.numeric(smig$fst)

## plot simulated FSTs

tmpfstplot=ggplot(data=tmig,aes(x=cm,y=fst)) + 
  geom_boxplot() +
  ylim(-0.005,0.014) +
  xlab(expression('Log'[10]*' Migration Rate')) +
  ylab(expression('Temporal F'['ST'])) +
  geom_hline(yintercept=0.012,linetype="dashed") +
  theme_bw()

spfstplot=ggplot(data=smig,aes(x=cm,y=fst)) + 
  geom_boxplot() +
  ylim(-0.005,0.14) +
  xlab(expression('Log'[10]*' Migration Rate')) +
  xlab("") +
  ylab(expression('Spatial F'['ST'])) +
  geom_hline(yintercept=0.11,linetype="dashed") +
  theme_bw()

ggarrange(spfstplot,tmpfstplot,ncol=2,nrow=1)
