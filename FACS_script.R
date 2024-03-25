library(flowCore)
library(ggcyto)
library(tidyverse)
library(knitr)

fs <- read.flowSet(path = "C:/Users/ccape/Downloads/2020-07-08-FACS-data/2020-07-08-FACS-data",pattern = ".fcs",alter.names = T)
#acces phenotype information
pData(fs)[1:3,]
pData(fs)$well <- gsub(".*_.*_(.*)_.*.fcs","\\1",sampleNames(fs)) # extract well from name and add new 'well' column
pData(fs)[1:3,]

colnames(fs)
#labling your wells by florescent protien.
colnames(fs)[colnames(fs)=="FITC.A"] <- "GFP"
colnames(fs)[colnames(fs)=="Pacific.Blue.A"] <- "BFP"

#gating 
gs <- GatingSet(fs)

g.singlets <- polygonGate(filterId = "Singlets","FSC.A"=c(2e4,25e4,25e4,2e4),"FSC.H"=c(0e4,12e4,18e4,6e4)) # define gate
ggcyto(gs[[1]],aes(x=FSC.A,y=FSC.H),subset="root")+geom_hex(bins = 200)+geom_gate(g.singlets)+ggcyto_par_set(limits = "instrument") # check gate

add(gs,g.singlets) # add gate to GatingSet
recompute(gs) # recompute 

ggcyto(gs[[1]],aes(x=FSC.A,y=SSC.A),subset="root")+geom_hex(bins = 200)+ggcyto_par_set(limits = "instrument")

#  checking where the singlets you filtered out show up in FSC vs SSC.

ggcyto(gs[[1]],aes(x=FSC.A,y=SSC.A),subset="root")+geom_hex(bins = 200)+ggcyto_par_set(limits = "instrument")

ggcyto(gs[[1]],aes(x=FSC.A,y=SSC.A),subset="Singlets")+geom_hex(bins = 200)+ggcyto_par_set(limits = "instrument")

#plotting all samples
ggcyto(gs,aes(x=FSC.A,y=FSC.H),subset="root")+geom_hex(bins = 100)+geom_gate("Singlets")+
  geom_stats(adjust = 0.8)+ggcyto_par_set(limits = "instrument")+
  facet_wrap(~well,ncol = 10)

#separating live from dead cells.

g.live <- polygonGate(filterId = "Live","FSC.A"=c(5e4,24e4,24e4,8e4),"SSC.A"=c(0,2e4,15e4,5e4)) # define gate
ggcyto(gs[[1]],aes(x=FSC.A,y=SSC.A),subset="Singlets")+geom_hex(bins = 200)+geom_gate(g.live)+ggcyto_par_set(limits = "instrument") # check gate
add(gs,g.live,parent="Singlets") # add gate to GatingSet
recompute(gs) # recompute GatingSet


ggcyto(gs,aes(x=FSC.A,y=SSC.A),subset="Singlets")+geom_hex(bins = 100)+geom_gate("Live")+
  geom_stats(adjust = 0.8)+ggcyto_par_set(limits = "instrument")+
  facet_wrap(~well,ncol = 10)
#focusing on GFP positve cells

g.gfp <- rectangleGate(filterId="GFP positive","GFP"=c(1000, Inf)) # set gate
ggcyto(gs[[1]],aes(x=GFP),subset="Live")+geom_density(fill="forestgreen")+geom_gate(g.gfp)+ggcyto_par_set(limits = "instrument")+scale_x_flowJo_biexp() # check gate


add(gs,g.gfp,parent="Live") # add gate to GatingSet
recompute(gs) # recalculate Gatingset


ggcyto(gs,aes(x=GFP),subset="Live",)+geom_density(fill="forestgreen")+geom_gate("GFP positive")+
  geom_stats(adjust = 0.1,y=0.002,digits = 1)+
  ggcyto_par_set(limits = "instrument")+scale_x_flowJo_biexp()+
  facet_wrap(~well,ncol = 10)

ps <- gs_pop_get_count_with_meta(gs)
ps <- ps %>% mutate(percent_of_parent=Count/ParentCount)
ps_df<-ps %>% select(sampleName,well,Population,Count,ParentCount,percent_of_parent) %>% head() %>% data.frame() 

write.csv(ps_df,"C:/Users/ccape/Downloads/FACS_results")


# orginal code from https://jchellmuth.com/posts/FACS-with-R/ added/changed lines to add with reporting.

