library(ggrepel)
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

source('genericFunctions.R')

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
### drawParentagePie
### KB October 29 2020
drawParentagePie <- function(dfPie,sNames=c('Blood:A-A','Sperm:A-A','Known alleles','Novel alleles')){
  
  dfParentagePies <- dfPie %>%
    group_by(name) %>%
    mutate(name=factor(name,levels=sNames),
           type=factor(type,levels=c('Yes: 1 switch',
                                     'Yes: 2 switches', 
                                     'Yes: 3 switches',  
                                     'Yes: 4 switches', 
                                     'No')),
           cumsumpos=cumsum(n/total*100),
           pc=n/total*100,
           ypos=cumsumpos-(pc/2),
           lbl=paste0(round(n/total*100,0),'% (N=',n,')'),
           end = 2 * pi * cumsum(pc)/100,
           start = lag(end, default = 0),
           middle = 0.5 * (start + end),
           hjust = ifelse(middle > pi, 1, 0),
           vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))
  
  gPie <- ggplot(dfParentagePies, aes(x = "", y = pc, fill = type)) + 
    geom_bar(width = 1,stat="identity",color='black',lwd=.3) + 
    geom_label_repel(aes(x=1.5,label = lbl), 
                     size=7*5/14, 
                     position = position_stack(vjust = 0.5),
                     segment.size = 0.4,
                     show.legend=FALSE)  + 
    coord_polar(theta = "y") + 
    facet_wrap(~name,nrow=1,strip.position = 'top')+
    scale_fill_manual('',values=c('Yes: 1 switch'   = '#44AA99',
                                  'Yes: 2 switches' = '#6699CC', 
                                  'Yes: 3 switches' = '#88CCEE',  
                                  'Yes: 4 switches' = '#117733', 
                                  'No'              = '#999933'),
                      guide=guide_legend(reverse = FALSE)) +
    scale_color_manual('',values=c('Yes: 1 switch'   = '#44AA99',
                                   'Yes: 2 switches' = '#6699CC', 
                                   'Yes: 3 switches' = '#88CCEE',  
                                   'Yes: 4 switches' = '#117733', 
                                   'No'              = '#999933'),
                       guide=guide_legend(reverse = FALSE)) +
    theme_void() + 
    theme(legend.position='top',
          legend.key.size=unit(c(0.3),'cm',),
          legend.text = element_text(size=7),
          strip.text=element_text(size=7,face='bold'),
          plot.margin=unit(c(0,0,0,0),'cm'))
  
  return(gPie)
}

#############################
### plotByNSwitches
### KB October 29 2020
plotByNSwitches <- function(dfSwitch){
  
  dfSwitch$name <- factor(dfSwitch$name,levels = pieOrder)
  g <- ggplot(dfSwitch, aes(x=name,y=n/total*100,fill=type)) + 
    geom_bar(color='black',lwd=.2,stat='identity') + 
    geom_text(aes(label=paste0(round(n/total*100,0),'% [N = ',n,']')),
              position = position_stack(vjust = 0.5),
              size = 7*5/14) +
    scale_fill_manual('',
                      values=c('firebrick1','grey90','grey70','grey60','grey50')) + 
    xlab('PrZFA source') + 
    ylab('Alleles with putative parent (%)') + 
    coord_cartesian(ylim=c(0,100),expand = FALSE) + 
    theme(legend.position = c(0.075,.05), 
          legend.justification = c(0,0),
          legend.title=element_blank(),
          legend.key.size = unit(0.2,'cm'))
  
  gPie <- drawParentagePie(dfSwitch)
  
  return(list(gBar=g,
              gPie=gPie,
              data=dfSwitch))
}
###############################

#############################
### plotGraph
### KB October 29 2020
plotGraph <- function(dfIn, whatFig='main', seed = 102){
  
  edgeDF         <- dfIn %>% filter(N>0,zfa != parent, inPop,parent_inPop) %>% select(from=parent,to=zfa,weight=N,type=newAllele) %>% distinct()
  
  nodesAsParents <- dfIn %>% filter(N>0,zfa != parent,inPop,parent_inPop) %>% select(from=parent) %>% 
    group_by(from) %>% count(name='nKids')
  
  nodesAsKids    <- dfIn %>% filter(N>0,zfa != parent,inPop,parent_inPop) %>% select(from=zfa) %>% 
    group_by(from) %>% count(name='nParents')
  
  nodesDF <- nodesAsKids %>% full_join(nodesAsParents) %>% 
    mutate(numKids    = ifelse(is.na(nKids),0,nKids),
           numParents = ifelse(is.na(nParents),0,nParents)) %>%
    mutate(nEdges = numKids+numParents)
  
  graph  <- graph_from_data_frame(d=edgeDF, vertices=nodesDF, directed=TRUE) 
  
  labelNames <- V(graph)$name
  
  if (whatFig == 'main'){
    labelNames[!(labelNames %in% c('A','C','F','E','L5','L20','M18','L19','L9','B','L4','N'))] <- ''
  }
  
  if (whatFig == 'suppNew'){
    labelNames[!(grepl('^M',labelNames))] <- ''
  }  
  
  V(graph)$lbl    <- labelNames
  V(graph)$degree <- degree(graph)
  V(graph)$hasKids <- V(graph)$numKids >0
  
  E(graph)$weight <- 5-as.numeric(E(graph)$weight)
  
  E(graph)$weightN <- 0.01*(5-as.numeric(E(graph)$weight))
  theme7point()
  
  #r=sample(1:2000,1)
  #seed <- 102
  set.seed(seed)
  
  xy <- layout_with_mds(graph) 
  xy[which(V(graph)$name %in% c('L9','L20','M18')),1] <- xy[which(V(graph)$name %in% c('L9','L20','M18')),1] + 9 
  xy[which(V(graph)$name %in% c('L20')),2]            <- xy[which(V(graph)$name %in% c('L20')),2] + 2 
  
  gGraph <- ggraph(graph,layout = "manual",x = xy[,1],y = xy[,2])+
    geom_edge_fan(aes(color=weight,alpha=weight), 
                  arrow = arrow(angle=10,
                                length=unit(.15,'inches'),
                                ends = "last",
                                type = "closed")) +
    scale_edge_color_manual(values=c('dodgerblue1','grey70','grey50','grey20'))+
    scale_edge_alpha_manual(values=c(.3,.5,.7,.9)) + 
    geom_node_point(aes(size = numParents,fill=hasKids),shape=21,color='black',alpha=1) + 
    geom_node_label(aes(label=lbl),size=8*5/14,repel=TRUE) +
    scale_fill_manual(values=c('#CC6677','#117733')) +
    scale_alpha_continuous(range=c(.8,1)) + 
    theme(legend.position='none',
          plot.margin = unit(c(0,0,0,0),'cm')) 
  
  return(gGraph)
  
  
}

###############################

theme7point()

popAlleles <- read.table('PrZFA_alleles.foundInPops.txt')
names(popAlleles) <- c('allele','code')

dfRelatedness <- fread('PrZFA_relatedness.human.tab',header=FALSE)
names(dfRelatedness) <- c('zfa','parent','event','N')

allelesAll  <- as.character(sort(unique(sort(rbind(dfRelatedness$zfa,dfRelatedness$parent)))))
alleleOrder <- c(LETTERS,
                 str_sort(allelesAll[grep('^L[0-9]+$',perl=TRUE,allelesAll)],numeric=TRUE),
                 str_sort(allelesAll[grep('^(A|C|L[0-9]+)v',perl=TRUE,allelesAll)],numeric=TRUE),
                 str_sort(allelesAll[grep('^M[0-9]+$',perl=TRUE,allelesAll)],numeric=TRUE),
                 str_sort(allelesAll[grep('^baudat_[A-Z]$',perl=TRUE,allelesAll)],numeric=TRUE),
                 'Unk','noZFA','Known','Novel')


dfRelatedness <- dfRelatedness[!(dfRelatedness$zfa %in% c('PrZFA')),]

dfRelatedness$simple                                                   <- grepl(':s:',dfRelatedness$zfa)
dfRelatedness$complex                                                  <- grepl(':c:',dfRelatedness$zfa)
dfRelatedness$possible_parent                                          <- grepl('END',dfRelatedness$event)
dfRelatedness$inPop                                                    <- !(dfRelatedness$complex | dfRelatedness$simple)
dfRelatedness$inPop[dfRelatedness$zfa %in% popAlleles$allele & !dfRelatedness$inPop]           <- TRUE

dfRelatedness$parent_inPop                                             <- !(grepl('v:',dfRelatedness$parent))
dfRelatedness$parent_inPop[dfRelatedness$parent %in% popAlleles$allele & !dfRelatedness$inPop] <- TRUE

dfRelatedness$newAllele                                                <- grepl('^M',dfRelatedness$zfa)
dfRelatedness$nSwitches                                                <- dfRelatedness$N

dfRelatedness$zfa[dfRelatedness$zfa == 'baudat_I'] <- 'I'
dfRelatedness$zfa[dfRelatedness$zfa == 'baudat_F'] <- 'F'
dfRelatedness$zfa[dfRelatedness$zfa == 'baudat_G'] <- 'G'
dfRelatedness$zfa[dfRelatedness$zfa == 'baudat_H'] <- 'H'
dfRelatedness$zfa[dfRelatedness$zfa == 'Av:s:0053:M1S:A-A'] <- 'N'

dfRelatedness$parent[dfRelatedness$parent == 'baudat_I'] <- 'I'
dfRelatedness$parent[dfRelatedness$parent == 'baudat_F'] <- 'F'
dfRelatedness$parent[dfRelatedness$parent == 'baudat_G'] <- 'G'
dfRelatedness$parent[dfRelatedness$parent == 'baudat_H'] <- 'H'
dfRelatedness$parent[dfRelatedness$parent == 'Av:s:0053:M1S:A-A'] <- 'N'

dfRelatedness$inPop[dfRelatedness$zfa == 'N']           <- TRUE
dfRelatedness$parent_inPop[dfRelatedness$parent == 'N'] <- TRUE

plotDF <- dfRelatedness
plotDF$nSwitches[plotDF$nSwitches == 0] <- NA

sampleOrder <- c('Known','Novel','Blood:A-A',
                 paste0('Sperm:',c('A-A',
                                   'C-C',
                                   'C-A',
                                   'D-A',
                                   'A-L12',
                                   'L14-A',
                                   'A-L15',
                                   'A-L16',
                                   'A-L20',
                                   'A-L24',
                                   'A-L27',
                                   'C-L12',
                                   'C-L14',
                                   'L14-L16')))

pieOrder <- c('Known alleles','Novel alleles','Blood:A-A',
                 paste0('Sperm:',c('A-A',
                                   'C-C',
                                   'C-A',
                                   'D-A',
                                   'A-L12',
                                   'L14-A',
                                   'A-L15',
                                   'A-L16',
                                   'A-L20',
                                   'A-L24',
                                   'A-L27',
                                   'C-L12',
                                   'C-L14',
                                   'L14-L16')))

## Plot A: Pie charts of inferred parentage ------------------------------------
dfAgg <- plotDF %>% filter(grepl('A-A',zfa) & parent == 'A') %>%
  mutate(src = ifelse (grepl('B:',zfa),'Blood','Sperm')) %>%
  group_by(zfa,src) %>% summarise(minS = min(nSwitches)) %>%
  mutate(hasParent = minS > 0) %>%
  group_by(src) %>%
  add_count(name = 'total') %>%
  group_by(src,minS,hasParent,total) %>%
  count() %>%
  mutate(type = ifelse(minS == 0, 
                       'No', 
                       ifelse(minS == 1, 
                              paste0('Yes: ',minS,' switch'),
                              paste0('Yes: ',minS,' switches'))))

dfInit <- plotDF %>% filter(grepl("-",zfa)) %>% 
  mutate(p1       = gsub(perl = TRUE, pattern = '^.+:(.+)-.+', replacement ='\\1',zfa),
         p2       = gsub(perl = TRUE, pattern = '^.+:.+-(.+)', replacement ='\\1',zfa),
         pP1      = gsub(perl = TRUE, pattern = '^.+:(.+-.+)', replacement ='\\1',parent),
         parental = gsub(perl = TRUE, pattern = '^.+:(.+-.+)', replacement ='\\1',zfa)) %>%
  mutate(ptype = ifelse(p1 == parent | p2 == parent, 
                        'Parental',
                        ifelse (parental == pP1, 'P-Var','Other')),
         src = ifelse (grepl('B:',zfa),'Blood','Sperm'),
         homhet = ifelse(p1==p2,'Hom','Het'))

## PLOT A: Sperm/Blood variants ------------------------------------------------
dfAABloodSperm <- dfInit %>%
  filter(ptype == 'Parental' & dfInit$parental == 'A-A') %>%
  group_by(zfa,parental,src) %>% 
  summarise(minS = min(nSwitches,na.rm=TRUE)) %>%
  replace_na(list(minS=0)) %>%
  mutate(hasParent = minS > 0) %>%
  group_by(parental,src) %>%
  add_count(name = 'total') %>%
  group_by(src,parental,minS,hasParent,total) %>%
  count() %>%
  mutate(type = ifelse(minS == 0, 
                       'No', 
                       ifelse(minS == 1, 
                              paste0('Yes: ',minS,' switch'),
                              paste0('Yes: ',minS,' switches'))),
         name = paste0(src,":",parental))

lst_AA_Blood_and_Sperm <- plotByNSwitches(dfAABloodSperm)

## PLOT B: Known vs Novel  -----------------------------------------------------
dfKnownVNew <- plotDF %>% 
  filter(inPop) %>%
  mutate(src = ifelse (grepl('^M',zfa),'Novel','Known')) %>%
  group_by(zfa,src) %>% 
  summarise(minS = min(nSwitches,na.rm=TRUE)) %>%
  mutate(hasParent = minS > 0) %>%
  group_by(src) %>%
  add_count(name = 'total') %>%
  group_by(src,minS,hasParent,total) %>%
  count() %>%
  mutate(type = ifelse(is.na(minS), 
                       'No', 
                       ifelse(minS == 1, 
                              paste0('Yes: ',minS,' switch'),
                              paste0('Yes: ',minS,' switches')))) %>%
  mutate(name=paste0(src,' alleles')) %>% 
  select(name,n,total,type) 

lst_Known_V_New <- plotByNSwitches(dfKnownVNew)

## Plot AB : Pies  --------------------------------------------------------------
gPieLegend <- getGGLegend(lst_Known_V_New$gPie)

gPies <- ggarrange(lst_AA_Blood_and_Sperm$gPie + theme(legend.position='none'),
          lst_Known_V_New$gPie        +  theme(legend.position='none'),
          ncol=2,nrow=1,
          labels=c('A','B'),
          hjust=0,vjust=1,
          font.label = list(size=8,face='bold'))

gAB <- ggarrange(gPies,gPieLegend,ncol=1,nrow=2,heights=c(8,1))

# PLOT C: Network --------------------------------------------------------------
gC           <- plotGraph(dfRelatedness)
gNetworkSupp <- plotGraph(dfRelatedness,whatFig='allLabels')

gABC <- ggarrange(gAB,gC,
                 ncol=1,nrow=2,
                 labels=c('','C'),
                 hjust=0,vjust=1,
                 font.label = list(size=8,face='bold'),
                 heights=c(1,1.8))

ggsave('Alleva_et_al_Figure3.png',plot = gABC, dpi=500,height=5,width=6)
ggsave('Alleva_et_al_Figure3.png',plot = gABC, dpi=500,height=5,width=6)

ggsave('Alleva_et_al_Supplementary_Human_Network_Figure_FullNetwork.png',
       plot = gNetworkSupp, dpi=500,height=8.5,width=11)
ggsave('Alleva_et_al_Supplementary_Human_Network_Figure_FullNetwork.pdf',
       plot = gNetworkSupp ,height=8.5,width=11)


