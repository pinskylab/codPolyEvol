mutages=read.table("slimulation/SLiM_MutationAges.txt",header=T)

mutages_nonfixed=mutages[which(mutages$Prev<129830),]
length(which(mutages_nonfixed$Tick>57400))/nrow(mutages_nonfixed)

mutages_nonfixed_m1=mutages[which(mutages$Prev<129830 & mutages$Type=="m1"),]
length(which(mutages_nonfixed_m1$Tick>57400))/nrow(mutages_nonfixed_m1)
mutages_nonfixed_m1_shared=mutages[which(mutages$Prev<129830 & mutages$Type=="m1" & mutages$Tick<57400),]
mutages_nonfixed_m1_postsplit=mutages[which(mutages$Prev<129830 & mutages$Type=="m1" & mutages$Tick>57400),]
plot(density(mutages_nonfixed_m1_shared$Prev/129830))
plot(density(mutages_nonfixed_m1_postsplit$Prev/129830))
mean(mutages_nonfixed_m1_shared$Prev/129830)
mean(mutages_nonfixed_m1_postsplit$Prev/129830)

mutages_nonfixed_m2=mutages[which(mutages$Prev<129830 & mutages$Type=="m2"),]
length(which(mutages_nonfixed_m2$Tick>57400))/nrow(mutages_nonfixed_m2)
mutages_nonfixed_m2_shared=mutages[which(mutages$Prev<129830 & mutages$Type=="m2" & mutages$Tick<57400),]
mutages_nonfixed_m2_postsplit=mutages[which(mutages$Prev<129830 & mutages$Type=="m2" & mutages$Tick>57400),]
plot(density(mutages_nonfixed_m2_shared$Prev/129830))
plot(density(mutages_nonfixed_m2_postsplit$Prev/129830))
mean(mutages_nonfixed_m2_shared$Prev/129830)
mean(mutages_nonfixed_m2_postsplit$Prev/129830)

mutages_nonfixed_m3=mutages[which(mutages$Prev<129830 & mutages$Type=="m3"),]
length(which(mutages_nonfixed_m3$Tick>57400))/nrow(mutages_nonfixed_m3)
mean(mutages_nonfixed_m3$Prev/129830)

#par(mar = c(bottom, left, top, right))
#The default value for mar is c(5.1, 4.1, 4.1, 2.1)
par(mfrow=c(3,1),mar = c(4.2, 4.5, 1.2, 0.1))
options(scipen=999)

plot((63940-mutages_nonfixed_m2$Tick)*10,mutages_nonfixed_m2$Prev/129830,xlim=c(0,639400),ylab="",xlab="",col="gray",cex=0.5,cex.axis=1.2)
abline(v=65300,lty=2)
plot((63940-mutages_nonfixed_m1$Tick)*10,mutages_nonfixed_m1$Prev/129830,xlim=c(0,639400),ylab="Frequency",xlab="",col="purple",cex=0.5,cex.lab=1.5,cex.axis=1.2)
abline(v=65300,lty=2)
plot((63940-mutages_nonfixed_m3$Tick)*10,mutages_nonfixed_m3$Prev/129830,xlim=c(0,639400),ylab="",xlab="Mutation Age (years)",col="red",cex=0.5,cex.lab=1.5,cex.axis=1.2)
abline(v=65300,lty=2)


par(mfrow=c(3,1),mar = c(4.2, 4.5, 1.2, 0.1))
options(scipen=999)

plot((63940-mutages_nonfixed_m2$Tick),mutages_nonfixed_m2$Prev/129830,xlim=c(0,639400),ylab="",xlab="",col="gray",cex=0.5,cex.axis=1.2)
abline(v=65300,lty=2)
plot((63940-mutages_nonfixed_m1$Tick),mutages_nonfixed_m1$Prev/129830,xlim=c(0,639400),ylab="Frequency",xlab="",col="purple",cex=0.5,cex.lab=1.5,cex.axis=1.2)
abline(v=65300,lty=2)
plot((63940-mutages_nonfixed_m3$Tick),mutages_nonfixed_m3$Prev/129830,xlim=c(0,639400),ylab="",xlab="Mutation Age (generations)",col="red",cex=0.5,cex.lab=1.5,cex.axis=1.2)
abline(v=65300,lty=2)


## ggtime!

library(ggplot2)

cbPalette <- c("#E69F00", "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

mutages_nonfixed$Age=639400-(mutages_nonfixed$Tick*10)
mutages_nonfixed$Frequency=mutages_nonfixed$Prev/129830
mutages_nonfixed$Type=gsub("m1","QTL",mutages_nonfixed$Type)
mutages_nonfixed$Type=gsub("m2","Neutral",mutages_nonfixed$Type)
mutages_nonfixed$Type=gsub("m3","Deleterious",mutages_nonfixed$Type)

options(scipen=999)

ggplot(aes(x=Age,y=Frequency,color=Type),data=mutages_nonfixed)+
  geom_point(show.legend=FALSE)+
  scale_colour_manual(values=cbPalette)+
  facet_grid(Type~.)+
  geom_vline(xintercept=65300,linetype="dashed")+
  labs(x="Age (Years)")+
  theme_bw()+
  theme(text=element_text(size=15))

## nodata

ggplot(aes(x=Age,y=Frequency,color=Type),data=mutages_nonfixed)+
# geom_point(show.legend=FALSE)+
  xlim(0,639400)+
  ylim(0,1)+
  scale_colour_manual(values=cbPalette)+
  facet_grid(Type~.)+
  geom_vline(xintercept=65300,linetype="dashed")+
  labs(x="Age (Years)")+
  theme_bw()+
  theme(text=element_text(size=15))




##################### 4 types! #####################

mutages=read.table("slimulation/SLiM_MutationAges_4type.txt",header=T)

mutages_nonfixed=mutages[which(mutages$Prev<132990),]
length(which(mutages_nonfixed$Tick>57400))/nrow(mutages_nonfixed)

mutages_nonfixed_m1=mutages[which(mutages$Prev<132990 & mutages$Type=="m1"),]
length(which(mutages_nonfixed_m1$Tick>57400))/nrow(mutages_nonfixed_m1)
mutages_nonfixed_m1_shared=mutages[which(mutages$Prev<132990 & mutages$Type=="m1" & mutages$Tick<57400),]
mutages_nonfixed_m1_postsplit=mutages[which(mutages$Prev<132990 & mutages$Type=="m1" & mutages$Tick>57400),]
plot(density(mutages_nonfixed_m1_shared$Prev/132990))
plot(density(mutages_nonfixed_m1_postsplit$Prev/132990))
mean(mutages_nonfixed_m1_shared$Prev/132990)
mean(mutages_nonfixed_m1_postsplit$Prev/132990)

plot(density(mutages_nonfixed_m1_shared$Prev/132990))
plot(density(mutages_nonfixed_m1_postsplit$Prev/132990))
plot(density(mutages_nonfixed_m1_shared$SelCoef))
plot(density(mutages_nonfixed_m1_postsplit$SelCoef))


mutages_nonfixed_m2=mutages[which(mutages$Prev<132990 & mutages$Type=="m2"),]
length(which(mutages_nonfixed_m2$Tick>57400))/nrow(mutages_nonfixed_m2)
mutages_nonfixed_m2_shared=mutages[which(mutages$Prev<132990 & mutages$Type=="m2" & mutages$Tick<57400),]
mutages_nonfixed_m2_postsplit=mutages[which(mutages$Prev<132990 & mutages$Type=="m2" & mutages$Tick>57400),]
plot(density(mutages_nonfixed_m2_shared$Prev/132990))
plot(density(mutages_nonfixed_m2_postsplit$Prev/132990))
mean(mutages_nonfixed_m2_shared$Prev/132990)
mean(mutages_nonfixed_m2_postsplit$Prev/132990)


mutages_nonfixed_m3=mutages[which(mutages$Prev<132990 & mutages$Type=="m3"),]
length(which(mutages_nonfixed_m3$Tick>57400))/nrow(mutages_nonfixed_m3)
mutages_nonfixed_m3_shared=mutages[which(mutages$Prev<132990 & mutages$Type=="m3" & mutages$Tick<57400),]
mutages_nonfixed_m3_postsplit=mutages[which(mutages$Prev<132990 & mutages$Type=="m3" & mutages$Tick>57400),]
plot(density(mutages_nonfixed_m3_shared$Prev/132990))
plot(density(mutages_nonfixed_m3_postsplit$Prev/132990))
mean(mutages_nonfixed_m3_shared$Prev/132990)
mean(mutages_nonfixed_m3_postsplit$Prev/132990)

mean(mutages_nonfixed_m3$Prev/132990)


plot(density(mutages_nonfixed_m3_shared$Prev/132990))
plot(density(mutages_nonfixed_m3_postsplit$Prev/132990))
plot(density(mutages_nonfixed_m3_shared$SelCoef))
plot(density(mutages_nonfixed_m3_postsplit$SelCoef))


mutages_nonfixed_m4=mutages[which(mutages$Prev<132990 & mutages$Type=="m4"),]
length(which(mutages_nonfixed_m4$Tick>57400))/nrow(mutages_nonfixed_m4)
mutages_nonfixed_m4_shared=mutages[which(mutages$Prev<132990 & mutages$Type=="m4" & mutages$Tick<57400),]
mutages_nonfixed_m4_postsplit=mutages[which(mutages$Prev<132990 & mutages$Type=="m4" & mutages$Tick>57400),]
plot(density(mutages_nonfixed_m4_shared$Prev/132990))
plot(density(mutages_nonfixed_m4_postsplit$Prev/132990))
mean(mutages_nonfixed_m4_shared$Prev/132990)
mean(mutages_nonfixed_m4_postsplit$Prev/132990)

plot(density(mutages_nonfixed_m4_shared$Prev/132990))
plot(density(mutages_nonfixed_m4_postsplit$Prev/132990))
plot(density(mutages_nonfixed_m4_shared$SelCoef))
plot(density(mutages_nonfixed_m4_postsplit$SelCoef))


## ggtime!

library(ggplot2)

cbPalette <- c("#E69F00", "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

mutages_nonfixed$Age=639300-(mutages_nonfixed$Tick*10)
mutages_nonfixed$Frequency=mutages_nonfixed$Prev/132990
mutages_nonfixed$Type=gsub("m1","QTL",mutages_nonfixed$Type)
mutages_nonfixed$Type=gsub("m2","Neutral",mutages_nonfixed$Type)
mutages_nonfixed$Type=gsub("m3","Deleterious1",mutages_nonfixed$Type)
mutages_nonfixed$Type=gsub("m4","Deleterious0",mutages_nonfixed$Type)

mutages_nonfixed %>% filter(Type!="Deleterious1") -> mutages_nonfixed_m3

options(scipen=999)

ggplot(aes(x=Age,y=Frequency,color=Type),data=mutages_nonfixed_m3)+
  geom_point(show.legend=FALSE)+
  scale_colour_manual(values=cbPalette)+
  facet_grid(Type~.)+
  geom_vline(xintercept=65300,linetype="dashed")+
  labs(x="Age (Years)")+
  theme_bw()+
  theme(text=element_text(size=15))

## nodata

ggplot(aes(x=Age,y=Frequency,color=Type),data=mutages_nonfixed)+
  # geom_point(show.legend=FALSE)+
  xlim(0,639400)+
  ylim(0,1)+
  scale_colour_manual(values=cbPalette)+
  facet_grid(Type~.)+
  geom_vline(xintercept=65300,linetype="dashed")+
  labs(x="Age (Years)")+
  theme_bw()+
  theme(text=element_text(size=15))


x <- rgamma(n=1000, shape=5, rate=3)


######## some densities

library(dplyr)

mutages_nonfixed %>%
  mutate(Split = ifelse(Age>65300,"Before Split","After Split")) %>%
  filter(Type=="Deleterious1") %>% 
  ggplot(aes(x=SelCoef,fill=Split)) + geom_density(alpha = 0.5)

mutages_nonfixed %>%
  mutate(Split = ifelse(Age>65300,"Before Split","After Split")) %>%
  filter(Type=="Deleterious0") %>% 
  ggplot(aes(x=SelCoef,fill=Split)) + geom_density(alpha = 0.5) + labs(title="Deleterious Mutations")

mutages_nonfixed %>%
  mutate(Split = ifelse(Age>65300,"Before Split","After Split")) %>%
  filter(Type=="QTL") %>% 
  ggplot(aes(x=SelCoef,fill=Split)) + geom_density(alpha = 0.5) + labs(title="QTL Mutations")









##################### 3 types! #####################

mutages=read.table("slimulation/SLiM_MutationAges_3type.txt",header=T)

mutages_nonfixed=mutages[which(mutages$Prev<131458),]
length(which(mutages_nonfixed$Tick>57400))/nrow(mutages_nonfixed)

mutages_nonfixed_m1=mutages[which(mutages$Prev<131458 & mutages$Type=="m1"),]
length(which(mutages_nonfixed_m1$Tick>57400))/nrow(mutages_nonfixed_m1)
mutages_nonfixed_m1_shared=mutages[which(mutages$Prev<131458 & mutages$Type=="m1" & mutages$Tick<57400),]
mutages_nonfixed_m1_postsplit=mutages[which(mutages$Prev<131458 & mutages$Type=="m1" & mutages$Tick>57400),]
plot(density(mutages_nonfixed_m1_shared$Prev/131458))
plot(density(mutages_nonfixed_m1_postsplit$Prev/131458))
mean(mutages_nonfixed_m1_shared$Prev/131458)
mean(mutages_nonfixed_m1_postsplit$Prev/131458)

plot(density(mutages_nonfixed_m1_shared$Prev/131458))
plot(density(mutages_nonfixed_m1_postsplit$Prev/131458))
plot(density(mutages_nonfixed_m1_shared$SelCoef))
plot(density(mutages_nonfixed_m1_postsplit$SelCoef))


mutages_nonfixed_m2=mutages[which(mutages$Prev<131458 & mutages$Type=="m2"),]
length(which(mutages_nonfixed_m2$Tick>57400))/nrow(mutages_nonfixed_m2)
mutages_nonfixed_m2_shared=mutages[which(mutages$Prev<131458 & mutages$Type=="m2" & mutages$Tick<57400),]
mutages_nonfixed_m2_postsplit=mutages[which(mutages$Prev<131458 & mutages$Type=="m2" & mutages$Tick>57400),]
plot(density(mutages_nonfixed_m2_shared$Prev/131458))
plot(density(mutages_nonfixed_m2_postsplit$Prev/131458))
mean(mutages_nonfixed_m2_shared$Prev/131458)
mean(mutages_nonfixed_m2_postsplit$Prev/131458)


mutages_nonfixed_m3=mutages[which(mutages$Prev<131458 & mutages$Type=="m3"),]
length(which(mutages_nonfixed_m3$Tick>57400))/nrow(mutages_nonfixed_m3)
mutages_nonfixed_m3_shared=mutages[which(mutages$Prev<131458 & mutages$Type=="m3" & mutages$Tick<57400),]
mutages_nonfixed_m3_postsplit=mutages[which(mutages$Prev<131458 & mutages$Type=="m3" & mutages$Tick>57400),]
plot(density(mutages_nonfixed_m3_shared$Prev/131458))
plot(density(mutages_nonfixed_m3_postsplit$Prev/131458))
mean(mutages_nonfixed_m3_shared$Prev/131458)
mean(mutages_nonfixed_m3_postsplit$Prev/131458)

mean(mutages_nonfixed_m3$Prev/131458)


plot(density(mutages_nonfixed_m3_shared$Prev/131458))
plot(density(mutages_nonfixed_m3_postsplit$Prev/131458))
plot(density(mutages_nonfixed_m3_shared$SelCoef))
plot(density(mutages_nonfixed_m3_postsplit$SelCoef))




## ggtime!

library(ggplot2)

cbPalette <- c("#E69F00", "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

mutages_nonfixed$Age=639300-(mutages_nonfixed$Tick*10)
mutages_nonfixed$Frequency=mutages_nonfixed$Prev/131458
mutages_nonfixed$Type=gsub("m1","QTL",mutages_nonfixed$Type)
mutages_nonfixed$Type=gsub("m2","Neutral",mutages_nonfixed$Type)
mutages_nonfixed$Type=gsub("m3","Deleterious",mutages_nonfixed$Type)

mutages_nonfixed %>% filter(Type!="Deleterious") -> mutages_nonfixed_m3

options(scipen=999)

ggplot(aes(x=Age,y=Frequency,color=Type),data=mutages_nonfixed_m3)+
  geom_point(show.legend=FALSE)+
  scale_colour_manual(values=cbPalette)+
  facet_grid(Type~.)+
  geom_vline(xintercept=65300,linetype="dashed")+
  labs(x="Age (Years)")+
  theme_bw()+
  theme(text=element_text(size=15))

ggplot(aes(x=Age,y=Frequency,color=Type),data=mutages_nonfixed)+
  geom_point(show.legend=FALSE)+
  scale_colour_manual(values=cbPalette)+
  facet_grid(Type~.)+
  geom_vline(xintercept=65300,linetype="dashed")+
  labs(x="Age (Years)")+
  theme_bw()+
  theme(text=element_text(size=15))

ggplot(aes(x=Age,y=Frequency,color=Type),data=mutages_nonfixed)+
  # geom_point(show.legend=FALSE)+
  xlim(0,639400)+
  ylim(0,1)+
  scale_colour_manual(values=cbPalette)+
  facet_grid(Type~.)+
  geom_vline(xintercept=65300,linetype="dashed")+
  labs(x="Age (Years)")+
  theme_bw()+
  theme(text=element_text(size=15))


x <- rgamma(n=1000, shape=5, rate=3)


######## some densities

library(dplyr)

mutages_nonfixed %>%
  mutate(Age = ifelse(Age>65300,"Before Split","After Split")) %>%
  filter(Type=="Deleterious") %>% 
  ggplot(aes(x=SelCoef,fill=Age)) + geom_density(alpha = 0.5) + labs(title="Deleterious Mutations") +
  theme_bw()+
  labs(x="Fitness Effect",y="Density")+
  theme(text=element_text(size=15))


mutages_nonfixed %>%
  mutate(Age = ifelse(Age>65300,"Before Split","After Split")) %>%
  filter(Type=="QTL") %>% 
  ggplot(aes(x=SelCoef,fill=Age)) + geom_density(alpha = 0.5) + labs(title="QTL Mutations") +
  theme_bw()+
  labs(x="Phenotypic Effect",y="Density")+
  theme(text=element_text(size=15))


#### applying maf

mutages=read.table("slimulation/SLiM_MutationAges_3type.txt",header=T)

mutages_nonfixed=mutages[which(mutages$Prev<129830),]
length(which(mutages_nonfixed$Tick>57400))/nrow(mutages_nonfixed)

mutages_maf=mutages[which(mutages$Prev<129830 & mutages$Prev>129830*0.05),]

nrow(mutages_nonfixed)
nrow(mutages_maf)

length(which(mutages_maf$Tick>57400))/nrow(mutages_maf)

mutages_maf_m1=mutages[which(mutages$Prev<129830 & mutages$Prev>129830*0.05 & mutages$Type=="m1"),]
length(which(mutages_maf_m1$Tick>57400))/nrow(mutages_maf_m1)

mutages_maf_m2=mutages[which(mutages$Prev<129830 & mutages$Prev>129830*0.05 & mutages$Type=="m2"),]
length(which(mutages_maf_m2$Tick>57400))/nrow(mutages_maf_m2)

mutages_maf_m3=mutages[which(mutages$Prev<129830 & mutages$Prev>129830*0.05 & mutages$Type=="m3"),]
length(which(mutages_maf_m1$Tick>57400))/nrow(mutages_maf_m3)
