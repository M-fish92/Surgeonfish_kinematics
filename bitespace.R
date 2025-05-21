#Surgeonfishes
#written by Michalis Mihalitsis
#-------------------------------------------------------------------------------------
rm(list=ls()) #clears R workspace (helps prevent errors from old data)
#-------------------------------------------------------------------------------------
setwd("")
library(readxl)
#load packages
library(stats)
library(ggeffects) #for effects plots in ggplot
library(tidyverse) #for data wrangling
##############################################################################################
library(vegan)
library(ape)
bites = read.csv('Acanthuridae_discrete_bites.csv',header = T, stringsAsFactors=T)
#
test = prcomp(bites[2:7],scale. = F,center = T)
biplot(test,var.axes=T,cex=0.8,xlabs=1:nrow(bites),
                pch = 16,xlab = "PC 1", ylab = "PC 2")
test2 = test$x[,c(1:2)] %>% as.data.frame()
test2 = cbind(test2,bites)
d=ggplot(test2,aes(x=PC1,y=PC2))+
  geom_point(aes(size=1))+
  theme_classic()
d + geom_density_2d_filled(contour_var = "ndensity",alpha=0.6)
#
library(rayshader)
library(reshape2)
if(run_documentation()) {
  plot_gg(test, multicore = TRUE, raytrace = TRUE, width = 7, height = 6, 
          scale = 500, windowsize = c(1200, 860), zoom = 0.6, phi = 30, theta = 30)
  render_snapshot()
}
###################################################################################
#################### NOW WE CALCULATE THE AVERAGE POSITION OF EACH SPECIES IN PC SPACE
ax1 = test2 %>% group_by(Species) %>% summarise(avg_ax1 = mean(PC1), sd_ax1 = sd(PC1))
ax2 = test2 %>% group_by(Species) %>% summarise(avg_ax2 = mean(PC2), sd_ax2 = sd(PC2))
survar = cbind(ax1,ax2) %>% as.data.frame()
survar = survar[,c(1,2,5)]
#
ggplot(survar,aes(x=avg_ax1,y=avg_ax2,label=Species))+
  geom_point()+
  geom_text()+
  theme_classic()
####### NOW WE CALCULATE KINEMATIC VARIANCE  #############
## FIRST BY SPECIES
var1 = test2 %>% group_by(Species) %>% summarise(var_ax1 = var(PC1))
var2 = test2 %>% group_by(Species) %>% summarise(var_ax2 = var(PC2))
var_sum = cbind(var1,var_ax2=var2$var_ax2)
var_sum_kin = var_sum %>% group_by(Species) %>% mutate(sumvar=sum(var_ax1,var_ax2))
ggplot(var_sum_kin,aes(x=reorder(Species,-sumvar),y=sumvar))+
  geom_point(size=3)+
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  theme_classic()








