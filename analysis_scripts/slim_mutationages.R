library(ggplot2)
library(tidyr)

## read simulated dataset with 3 mutation types

mutages=read.table("slim_outputs/SLiM_MutationAges_3type.txt",header=T)

mutages_nonfixed=mutages[which(mutages$Prev<131458),]
length(which(mutages_nonfixed$Tick>57400))/nrow(mutages_nonfixed)

## calculate frequencies and ages for QTL mutations before and after population split

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

## calculate frequencies and ages for neutral mutations before and after population split

mutages_nonfixed_m2=mutages[which(mutages$Prev<131458 & mutages$Type=="m2"),]
length(which(mutages_nonfixed_m2$Tick>57400))/nrow(mutages_nonfixed_m2)
mutages_nonfixed_m2_shared=mutages[which(mutages$Prev<131458 & mutages$Type=="m2" & mutages$Tick<57400),]
mutages_nonfixed_m2_postsplit=mutages[which(mutages$Prev<131458 & mutages$Type=="m2" & mutages$Tick>57400),]
plot(density(mutages_nonfixed_m2_shared$Prev/131458))
plot(density(mutages_nonfixed_m2_postsplit$Prev/131458))
mean(mutages_nonfixed_m2_shared$Prev/131458)
mean(mutages_nonfixed_m2_postsplit$Prev/131458)

## calculate frequencies and ages for deleterious mutations before and after population split

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




## plot age + frequency by type

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


######## density plots of effect sizes for QTL and deleterious mutations

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


#### applying minor allele filter and recalculating prevalences for different mutation types 

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
