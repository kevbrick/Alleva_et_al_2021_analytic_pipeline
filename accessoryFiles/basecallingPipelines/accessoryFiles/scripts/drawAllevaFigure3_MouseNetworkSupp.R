library(plyr)
library(ggplot2)
library(data.table)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(grid)
library(png)
library(igraph)
library(ggraph)


#############################
### theme7point
### KB June 29 2018
### Set the default theme to 7-point font
### suitable for publication
### ARGS: none
### OUTPUTS: Sets the theme
theme7point <- function() {
  theme_base_size <- 7
  theme_set(theme_bw(base_size = theme_base_size) %+replace%
              theme(axis.text = element_text(size=theme_base_size,
                                             color='black'),
                    panel.grid = element_blank(),
                    panel.border=element_blank(),
                    axis.line=element_line(size=.2),
                    axis.ticks = element_line(size=.2)))
}

#############################
### plotGraph
### KB October 29 2020
plotMouseGraph <- function(dfIn, dfStrains, seed = 102){
  
  edgeDF         <- dfIn %>% filter(N>0,zfa != parent) %>% select(from=parent,to=zfa,weight=N) %>% distinct()
  
  nodesAsParents <- dfIn %>% filter(N>0,zfa != parent) %>% select(from=parent) %>% 
    group_by(from) %>% count(name='nKids')
  
  nodesAsKids    <- dfIn %>% filter(N>0,zfa != parent) %>% select(from=zfa) %>% 
    group_by(from) %>% count(name='nParents')
  
  nodesDF <- nodesAsKids %>% full_join(nodesAsParents) %>% 
    mutate(numKids    = ifelse(is.na(nKids),0,nKids),
           numParents = ifelse(is.na(nParents),0,nParents)) %>%
    mutate(nEdges = numKids+numParents)
  
  graph  <- graph_from_data_frame(d=edgeDF, vertices=nodesDF, directed=TRUE) 
  
  dfN <- data.frame(name = V(graph)$name) %>% inner_join(dfStrains)
  
  labelNames  <- dfN$PrZFA
  labelNames[is.na(labelNames)] <- ''
  #labelStrain <- dfN$strain
  
  labelPubNames  <- dfN$pub
  labelPubNames[is.na(labelPubNames)] <- ''
  
  V(graph)$lbl     <- labelNames
  V(graph)$pub     <- labelPubNames
  V(graph)$strain  <- dfN$strain
  V(graph)$degree  <- degree(graph)
  V(graph)$hasKids <- V(graph)$numKids >0
  
  E(graph)$weight <- max(as.numeric(E(graph)$weight))-as.numeric(E(graph)$weight)+1
  
  theme7point()
  
  #r=sample(1:2000,1)
  #seed <- 102
  set.seed(seed)
  
  xy <- layout_with_fr(graph,niter = 2000,start.temp = 1) 

  gGraphMain <- ggraph(graph,layout = "manual",x = xy[,1],y = xy[,2])+
    geom_edge_fan(aes(color=as.numeric(weight)),alpha=0.6, 
                  arrow = arrow(angle=10,
                                length=unit(.15,'inches'),
                                ends = "last",
                                type = "closed")) +
    geom_node_point(aes(size = numParents,fill=hasKids),shape=21,color='black',alpha=1) + 
    scale_size_continuous('# Parents') +
    scale_fill_manual('Parent allele',values=c('firebrick','forestgreen')) +
    scale_edge_color_gradient('# Switches',low='grey20',high=alpha('dodgerblue1',.5)) +
    theme(legend.position='right',
          legend.key.size = unit(0.3,'cm'),
          plot.margin = unit(c(0,0,0,0),'cm')) 
  #geom_node_label(aes(label=labelNames),repel=TRUE) 
  #scale_edge_color_manual(values=c('dodgerblue1','grey70','grey50','grey20'))+
  #scale_edge_alpha_manual(values=c(.3,.5,.7,.9)) +   
  
  gGraphStrain <- ggraph(graph,layout = "manual",x = xy[,1],y = xy[,2])+
    geom_edge_fan(aes(color=as.numeric(weight)),alpha=.6, 
                  arrow = arrow(angle=10,
                                length=unit(.15,'inches'),
                                ends = "last",
                                type = "closed")) +
    geom_node_point(aes(size = numParents,fill=strain),shape=21,color='black',alpha=1) + 
    scale_size_continuous('# Parents') +
    scale_fill_manual('Strain',values=c('pink','yellow','green','dodgerblue1','grey50')) +
    scale_edge_color_gradient('# Switches',low='grey20',high=alpha('dodgerblue1',.5)) +
    theme(legend.position='right',
          legend.key.size = unit(0.3,'cm'),
          plot.margin = unit(c(0,0,0,0),'cm')) 
  
  gPubStrains <-   ggraph(graph,layout = "manual",x = xy[,1],y = xy[,2])+
    geom_edge_fan(aes(color=as.numeric(weight)),alpha=.6, 
                  arrow = arrow(angle=10,
                                length=unit(.15,'inches'),
                                ends = "last",
                                type = "closed")) +
    geom_node_point(aes(size = numParents,fill=strain),shape=21,color='black',alpha=1) + 
    scale_size_continuous('# Parents') +
    scale_fill_manual('Strain',values=c('pink','yellow','green','dodgerblue1','grey50')) +
    scale_edge_color_gradient('# Switches',low='grey20',high=alpha('dodgerblue1',.5)) +
    theme(legend.position='right',
          legend.key.size = unit(0.3,'cm'),
          plot.margin = unit(c(0,0,0,0),'cm')) +
    geom_node_label(aes(label=labelPubNames),color='magenta',repel=TRUE) 
  
  #geom_node_label(aes(label=labelNames),repel=TRUE) +
  #scale_edge_color_manual(values=c('dodgerblue1','grey70','grey50','grey20'))+
  #scale_edge_alpha_manual(values=c(.3,.5,.7,.9)) +   
  
  gX2 <- ggarrange(gGraphMain,gPubStrains,ncol=1,nrow=2)
  
  return(list(gMain = gGraphMain,
              gStrain = gPubStrains,
              gX2 = gX2))
  
}

theme7point()

dfStrain      <- read.table('PrZFA_allele_details.mouse.txt',header=FALSE)
names(dfStrain) <- c('strain','name','id','pub','code')

x <- fread('PrZFA_relatedness.mouse.tab',header=FALSE)
names(x) <- c('zfa','parent','event','N')

allelesOrder <- as.character(sort(unique(sort(rbind(x$zfa,x$parent)))))

x <- x[!(x$zfa %in% c('PrZFA')),]

x$possible_parent                                          <- grepl('END',x$event)
x$nSwitches                                                <- x$N

x$zfa    <- gsub('^mm','',x$zfa)
x$parent <- gsub('^mm','',x$parent)

x <- x %>% 
  inner_join(dfStrain, by = c('zfa' = 'name')) %>% rename(child_strain = strain) %>% 
  inner_join(dfStrain, by = c('parent' = 'name')) %>% rename(parent_strain = strain) 

lMouse <- plotMouseGraph(x,dfStrain,seed=102)

ggarrange(lMouse$gMain,lMouse$gStrain,labels=c('A','B'),ncol=1,nrow=2,
          font.label = list(size=8,face='bold'))

ggsave('Alleva_recomb_Fig3_MouseNetworkSupplement.png',dpi=500,height=9,width=7)
ggsave('Alleva_recomb_Fig3_MouseNetworkSupplement.pdf',height=9,width=7)