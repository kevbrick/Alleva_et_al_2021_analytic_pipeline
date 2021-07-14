library(ggplot2)
library(plyr)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(RColorBrewer)

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
###############################

theme7point()

## Plotting function
plotMe <- function(df,numCols){
  g <- ggplot(data=df,aes(x=type,y=cover)) + 
    geom_boxplot(data=df %>% dplyr:::filter(coltype != 'ZF unique to one allele'),
                 color='grey60',varwidth=TRUE,outlier.alpha=0,lwd=.2) + 
    geom_hline(data=df %>% group_by(idGT,type) %>% dplyr:::filter(coltype=='ZF in neither allele') %>%summarize(m=quantile(cover,.98)),
               aes(yintercept=m),color='forestgreen',lty='dashed',lwd=.2) + 
    geom_hline(data=df %>% group_by(idGT,type) %>% dplyr:::filter(coltype=='ZF in neither allele') %>%summarize(m=quantile(cover,.95)),
               aes(yintercept=m),color='forestgreen',lty='dotted',lwd=.2) + 
    geom_point(data=df %>% dplyr:::filter(coltype != 'ZF unique to one allele'),
               position=position_jitter(h=0,width=.2),
               aes(color=coltype),
               size=.5) +
    geom_point(data=df %>% dplyr:::filter(coltype == 'ZF unique to one allele'),
               aes(color=coltype),
               size=.5) + 
    xlab('ZF type') + 
    ylab('Coverage\n(short-reads)') + 
    facet_wrap(~idGT,scales='free',ncol=numCols) + 
    scale_color_manual('',values=c('firebrick2','darkolivegreen3','dodgerblue2')) +
    theme(strip.background=element_blank(),
          strip.text=element_text(size=7,face='bold',hjust=.5),
          panel.border=element_rect(fill=NA,size=.2),
          legend.position='top',
          legend.key.size=unit(0.3,'cm'),
          legend.text = element_text(size=7))
  return(g)
}

getMorder <- function(gt){
  nret <- 999999999
  gts <- str_split(gt,"/")[[1]]
  for (g in gts){
    if (grepl("^M(\\d+)$",g)){
      g <- gsub("M(\\d+)","\\1",g)
      nret <- min(nret,as.numeric(g))
    }
  }
  return(nret)
}

###############################
dfAll <- read.table('allShortReadValidationData.tab') %>% 
  dplyr:::rename(id=V1,
         genotype=V2,
         zf=V3,
         type=V4,
         coverraw=V5,
         nZF=V6,
         cover=V7) %>%
  mutate(idGT=paste0(id," [",genotype,"]"),
         coltype=ifelse(type == 'Both',
                        'ZF in both alleles',
                        ifelse(type == 'None',
                               'ZF in neither allele',
                               'ZF unique to one allele'))) %>%
  rowwise %>%
  mutate(ord=getMorder(genotype)) %>% 
  dplyr:::filter(!grepl('Unk',genotype))

dfAll$idGT <- factor(dfAll$idGT,levels=unique(dfAll$idGT[order(dfAll$ord)]))

## Sort by type
uniqueTypes <- unique(dfAll$type)[!(unique(dfAll$type) %in% c('Both','None'))]
typeOrder <- uniqueTypes[order(uniqueTypes)]

dfAll$type <- factor(dfAll$type,levels=c("Both",typeOrder,"None"))

## Get unique IDs & count them
uIDs <- levels(dfAll$idGT)
nIDs <- length(uIDs)

## Set up plot parameters across multiple pages
plotSize <- 25
plotCols <- 5

froms <- seq(1       ,(nIDs-plotSize+1),plotSize)
tos   <- seq(plotSize,nIDs           ,plotSize)

nLast <- tos[length(tos)]
if (nIDs > nLast){
  froms <- c(froms,nLast+1)
  tos   <- c(tos,nIDs)
}

## Make each plot & save
for (i in 1:length(tos)){
  lstIDs <- uIDs[froms[i]:tos[i]]
  dfPlot <- dfAll[dfAll$idGT %in% lstIDs,]
  print(paste0("Plot ",i," : ",length(lstIDs)))
  gP <- plotMe(dfPlot,plotCols)
  print(gP)
  ggsave (paste0('Alleva_et_al_newAlleleValidation_Page',i,'.png'),
          width=6.5,height=1.3*ceiling(length(lstIDs)/plotCols),
          dpi=400)
  ggsave (paste0('Alleva_et_al_newAlleleValidation_Page',i,'.pdf'),
          width=6.5,height=1.3*ceiling(length(lstIDs)/plotCols))
}
