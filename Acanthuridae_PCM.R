#Ecomorphology of coral reef fishes
#written by Michalis Mihalitsis
#-------------------------------------------------------------------------------------
rm(list=ls()) #clears R workspace (helps prevent errors from old data)
#-------------------------------------------------------------------------------------
setwd("")
library(readxl)
#load packages
library(stats)
library(phytools)
library(ape)
library(ggtree)
library(tidyverse)
library(treeio)
library(geiger)
##### FIRST WE IMPORT THE PHYLOGENETIC TREE  ####
tr = read.tree('actinopt_12k_treePL.tre')
#then the dataset with thte traits
svl <- read.csv('Acanthuridae_frequency_traits.csv',header = T, stringsAsFactors=T)
svl2 = svl[,2:7]
rownames(svl2) = svl$Species
species = c('Acanthurus_japonicus','Acanthurus_leucosternon','Naso_lituratus','Naso_elegans',
            'Prionurus_laticlavius', 'Zebrasoma_desjardinii','Zebrasoma_velifer',
            'Acanthurus_bahianus','Acanthurus_coeruleus','Acanthurus_chirurgus',
            'Acanthurus_nigrofuscus','Acanthurus_xanthopterus',
            'Ctenochaetus_strigosus','Acanthurus_lineatus','Paracanthurus_hepatus','Zanclus_cornutus',
            'Kyphosus_vaigiensis','Apolemichthys_xanthurus','Siganus_magnificus','Chaetodon_kleinii')
tr = drop.tip(tr,tr$tip.label[-match(species, tr$tip.label)])
p=ggtree(tr, color="black", size=1,ladderize=T)%<+% svl2+
  geom_tiplab()+
  #geom_tippoint(aes(shape = trait, color = trait,size=3))+
  coord_cartesian(clip = 'off')+
  #geom_treescale(x=4, y=3.5)
  theme_tree2(plot.margin=margin(20, 10, 20, 6))
gheatmap(p, svl2, offset=30, width=0.3,colnames=T,colnames_angle=-90)+
  scale_fill_viridis_c(option="D", name="continuous\nvalue")
library(deeptime)
p <- ggtree(tr,ladderize = T) +
  coord_geo(xlim = c(-110, 30), ylim = c(-2, Ntip(tr)+2), neg = TRUE, abbrv = FALSE) +
  geom_tiplab()+
  scale_x_continuous(breaks = seq(-110, 0, 20), labels = abs(seq(-110, 0, 20))) +
  theme_tree2()
#revts(p)
#######################################################################################
###################################################################################
toothdatall = read_xlsx('surgeonfish_teeth.xlsx',sheet = 'Sheet3')
toothdatall = toothdatall[,c(2:5)]
toothdatall = toothdatall %>% group_by(Genus,species) %>% 
  summarise_at(c('uj_ratio','lj_ratio'),mean,na.rm = TRUE) %>% ungroup()
toothdatall$sp = paste(toothdatall$Genus, toothdatall$species, sep="_")
#
surtree = read.tree('Acanthuridae_nwk.tree')
species = c('Acanthurus_bahianus','Acanthurus_chirurgus','Acanthurus_japonicus','Acanthurus_leucosternon',
            'Acanthurus_lineatus','Acanthurus_nigrofuscus','Acanthurus_xanthopterus','Naso_lituratus',
            'Naso_elegans','Paracanthurus_hepatus','Prionurus_laticlavius',
            'Zebrasoma_desjardinii','Zebrasoma_velifer','Acanthurus_coeruleus',
            'Acanthurus_achilles','Acanthurus_blochii','Acanthurus_guttatus',
            'Acanthurus_gahhm','Acanthurus_maculiceps','Acanthurus_mata','Acanthurus_nigricans',
            'Acanthurus_nigroris','Acanthurus_nigricauda','Acanthurus_olivaceus','Acanthurus_pyroferus','Acanthurus_sohal',
            'Acanthurus_tennentii','Acanthurus_thompsoni','Acanthurus_triostegus','Naso_brevirostris','Naso_vlamingii',
            'Prionurus_punctatus','Zebrasoma_flavescens','Zebrasoma_scopas','Zebrasoma_xanthurum','Ctenochaetus_binotatus',
            'Ctenochaetus_striatus','Ctenochaetus_strigosus')
surtree = drop.tip(surtree,surtree$tip.label[-match(species, surtree$tip.label)])
#
svl2=toothdatall[,c(3)] %>% as.data.frame()
rownames(svl2) = toothdatall$sp
p=ggtree(surtree, color="black", size=0.1,ladderize=T)+ #%<+% svl2+
  geom_tiplab()+
  scale_x_continuous(labels=abs)+
  coord_cartesian(clip = 'off')
gheatmap(p, svl2, offset=30, width=0.09,colnames=T,colnames_angle=-90,colnames_offset_y = .25)+
  scale_fill_viridis_c(limits = c(0.01,1),option="B", name="Tooth symmetry",direction = -1)
#
library(deeptime)
p <- ggtree(surtree,ladderize = T) +
  coord_geo(xlim = c(-80, 30), ylim = c(-2, Ntip(surtree)+2), neg = TRUE, abbrv = FALSE) +
  geom_tiplab()+
  scale_x_continuous(breaks = seq(-80, 0, 20), labels = abs(seq(-80, 0, 20))) +
  theme_tree2()
revts(p)
##############################################################################
preds = read.csv('surgeonfish_predictions.csv',row.names = 1,stringsAsFactors=T)
preds$species2 <- paste(preds$Genus, preds$species, sep="_")
rownames(preds) = preds$species2
name.check(surtree,preds)
#### first, lets plot the character on the tree  #############################
svl = preds[,5] %>% as.data.frame()
rownames(svl) <- rownames(preds)
p=ggtree(surtree, color="black", size=1,ladderize=T) %<+% svl+
  geom_tiplab()+
  scale_x_continuous(labels=abs)+
  coord_cartesian(clip = 'off')
gheatmap(p, svl, offset=30, width=0.1,colnames=T,colnames_angle=-90,colnames_offset_y = .25)+
  scale_fill_viridis_d(option="D", name="Bite adaptation")
####
svl3 = data.frame(Genus_sp=rownames(svl),trait=as.numeric(svl$.))
rownames(svl3) = rownames(svl)
fmode<-setNames(svl3$trait,rownames(svl))
name.check(surtree,svl3)
fmode<-setNames(svl4$trait,rownames(svl4))
name.check(surtree,svl4)
traitmap = make.simmap(surtree,fmode,model='ER',nsim=100,Q='mcmc',vQ=0.01,
                       prior=list(use.empirical=T),samplefreq=10)
#write.simmap(traitmap,'traitmap.nex',format="phylip", version=1.5)
cols=set_names(c('#899DA4',"#C93311",'#EAD9AD'),c('1','2','3'))
cols=set_names(c('#446455',"#FDD262",'#C93311'),c('1','2','3'))
densityTree(traitmap,method="plotSimmap",alpha=0.1,lwd=4,nodes="intermediate",
            colors=cols,compute.consensus=F)
test=summary(traitmap) # to see the time spent in each state
test2=describe.simmap(traitmap)
transitions=countSimmap(traitmap)
#
fiter = fitDiscrete(surtree,fmode,model='ER',transform = 'none')
fitsym = fitDiscrete(surtree,fmode,model='SYM')
fitard = fitDiscrete(surtree,fmode,model='ARD')
AIC(fiter,fitsym,fitard)
plot(fiter,color=TRUE,lwd=3)
###############################################################################
############################ NOW WE MAKE A SIMMAP TO USE FOR OU MODELS  ########
preds = read.csv('for_OU.csv')
preds$species2 <- paste(preds$Genus, preds$species, sep="_")
rownames(preds) = preds$species2
name.check(surtree,preds)
fmode<-setNames(preds$uj_ratio,rownames(preds))
traitmap = make.simmap(surtree,fmode,model='ER',nsim=100,Q='mcmc',vQ=0.01,
                       prior=list(use.empirical=T),samplefreq=10)
plotSimmap(traitmap)
fiter = fitMk(surtree,fmode,model='ER')
fitard = fitMk(surtree,fmode,model='ARD')
fitsym = fitMk(surtree,fmode,model='SYM')
AIC(fiter,fitard,fitsym)
#
write.simmap(traitmap,'traitmap.nex',format="phylip", version=1.5)
#traitmap<-read.simmap("traitmap.nex",format="phylip",version=1.5)
##############################################################################
### NEXT WE RUN THE EVOLUTION MODELS
##############################################################################
uj_ratio=preds$uj_ratio
names(uj_ratio) = rownames(preds)
phylosig(surtree,uj_ratio,method='lambda',test = T)
evmod1 = fitContinuous(surtree,uj_ratio,model='BM') #brownian motion
evmod2 = fitContinuous(surtree,uj_ratio,model='EB') #early burst
evmod3 = fitContinuous(surtree,uj_ratio,model='OU')
evmod4 = fitContinuous(surtree,uj_ratio,model='rate_trend')
evmod5 = fitContinuous(surtree,uj_ratio,model='lambda')
evmod6 = fitContinuous(surtree,uj_ratio,model='kappa')
evmod7 = fitContinuous(surtree,uj_ratio,model='delta')
evmod8 = fitContinuous(surtree,uj_ratio,model='mean_trend')
evmod9 = fitContinuous(surtree,uj_ratio,model='white')
AIC(evmod1,evmod2,evmod3,evmod4,evmod5,evmod6,evmod7,evmod8,evmod9)
###
library(OUwie)
fmode<-setNames(preds$pred_uj,rownames(preds))
ecomorph.tree<-read.simmap('traitmap.nex',format = 'phylip',version=1.5)
ouwiedata = data.frame(Genus_species=rownames(preds),Reg=preds$pred_uj,
                       X=preds$uj_ratio)
fitbm=OUwie(ecomorph.tree,ouwiedata,model="BM1",simmap.tree = T)
fitbms=OUwie(ecomorph.tree,ouwiedata,model="BMS",simmap.tree = T,root.station = F)
fitou1=OUwie(ecomorph.tree,ouwiedata,model="OU1",simmap.tree = T,root.station = F)
fitoum=OUwie(ecomorph.tree,ouwiedata,model="OUM",simmap.tree = T,root.station = F)
fitoumv=OUwie(ecomorph.tree,ouwiedata,model="OUMV",simmap.tree = T,root.station = F)
fitouma=OUwie(ecomorph.tree,ouwiedata,model="OUMA",simmap.tree = T,root.station = F)
fitoumva=OUwie(ecomorph.tree,ouwiedata,model="OUMVA",simmap.tree = T,root.station = F)
aic= setNames(c(fitbm$AIC,fitbms$AIC,fitou1$AIC,fitoum$AIC,fitoumv$AIC,fitouma$AIC,fitoumva$AIC),
              c('BM1','BMS','OU1','OUM','OUMV','OUMA','OUMVA'))
aic
aicc= setNames(c(fitbm$AICc,fitbms$AICc,fitou1$AICc,fitoum$AICc,fitoumv$AICc,fitouma$AICc,fitoumva$AICc),
              c('BM1','BMS','OU1','OUM','OUMV','OUMA','OUMVA'))
aicc
#save.image(fitoumv,file='OUmodels.rdata')
#
phenogram(ecomorph.tree,uj_ratio,ftype="i",colors=setNames(c('#446455',"#FDD262",'#C93311'),c('1','2','3')),
          fsize=0.5,xlab="time (ma)",ylab="Tooth symmetry",spread.labels=TRUE,spread.cost=c(2.5,1.5))
ggplot(ouwiedata,aes(x=X))+
  geom_density(aes(fill=Reg),alpha=0.7)+
  geom_vline(xintercept = 0.908)+
  geom_vline(xintercept = 0.91206439+0.01890357,linetype='dashed')+
  geom_vline(xintercept = 0.628)+
  geom_vline(xintercept = 0.64499779+0.02506114,linetype='dashed')+
  geom_vline(xintercept = 0.129)+
  geom_vline(xintercept = 0.14660255-0.05120231,linetype='dashed')+
  theme_classic()


