library(ggplot2)
library(plyr)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(RColorBrewer)
library(data.table)

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
## Set order of Alleles / ZFAs
setAlleleOrder <- function(df,col='id',rev=FALSE){
  allelesAll  <- sort(unique(sort(df[[col]])))
  alleleOrder <- c(str_sort(allelesAll[grep('^[A-Z]$',perl=TRUE,allelesAll)],numeric=TRUE),
                   str_sort(allelesAll[grep('^L[0-9]+$',perl=TRUE,allelesAll)],numeric=TRUE),
                   str_sort(allelesAll[grep('^(A|C|L[0-9]+)v',perl=TRUE,allelesAll)],numeric=TRUE),
                   str_sort(allelesAll[grep('^M[0-9]+$',perl=TRUE,allelesAll)],numeric=TRUE),
                   'Unk','noZFA')
  
  if (rev){
    df[[col]] <- factor(df[[col]],levels=rev(alleleOrder))
  }else{
    df[[col]] <- factor(df[[col]],levels=alleleOrder)
  }
  
  return(df)
}

## Part A
drawParentageForJeffreys2013 <- function(){
  dfParentage <- fread('parentage_Jeffreys_2013_all.txt',header=TRUE) %>%
    rowwise() %>%
    mutate(shortid=gsub("^(.+v:).+:(\\d+):M\\d+[SB].+","\\1\\2",id,perl=TRUE),
           type = ifelse(type=='Monoallelic','Uni-Parental','Bi-Parental'),
           has_recombinant=ifelse(code=='NORECOMBINANT',FALSE,TRUE),
           mannum=as.numeric(gsub("^M(\\d+)$","\\1",manid,perl=TRUE)),
           mansid=paste0(gsub("^M(\\d+)$","\\1",manid,perl=TRUE),sb),
           mangt1=gsub("^(\\S+)\\/\\S+$","\\1",mangt,perl=TRUE),
           mangt2=gsub("^\\S+\\/(\\S)+$","\\1",mangt,perl=TRUE),
           man=paste0("Man ", mannum, ifelse(sb=="Sperm","S","B"),": ", mangt),
           switchkey=ifelse(p1 == mangt1,switchkey,gsub(":4",":2",gsub(":2",":1",gsub(":1",":4",paste0(":",switchkey))))),
           switchkey=gsub("^:","",switchkey),
           nt=paste0(as.numeric(strsplit(switchkey,":")[[1]]),collapse = ""))
  
  dfUniORBi <- dfParentage %>%
    group_by(id,man,shortid,nswitches) %>% 
    summarise(nd=n_distinct(type)) 
  
  dfBestRecombinant <- dfParentage %>%
    group_by(man,shortid) %>% 
    slice_min(order_by = nswitches) %>% 
    slice_head()
  
  dfParentageProcessed <- inner_join(dfBestRecombinant,
                                     dfUniORBi,
                                     by = c('id','man','nswitches'))%>%
    mutate(type=ifelse(nd==2,"Uni/Bi",
                       ifelse(mangt1 == mangt2,
                              'Uni/Bi',
                              type)))
  
  dfParentagePlot <- reshape2:::melt.data.frame(dfParentageProcessed %>% 
                                                  mutate(type=ifelse(nswitches==99,"None",type)) %>%
                                                  group_by(man) %>% 
                                                  add_count(name='tot') %>% 
                                                  group_by(man,type,tot) %>% 
                                                  dplyr:::count(),
                                                measure.vars = c('n'),
                                                id.vars=c('type','man','tot')) 
  
  dfParentagePlot$type <- factor(dfParentagePlot$type,levels=c('None','Uni-Parental','Bi-Parental','Uni/Bi'))
  
  # Order men by # of alleles with no parentage
  dfMen <- dfParentagePlot %>% 
   mutate(pc=100*value/tot) %>% 
   dplyr:::filter(type!='None') %>% 
   group_by(man) %>% 
   tally(pc) %>% 
  mutate(mannum=gsub("Man\\s(\\d+)[SB].+$","\\1",man))
   
  dfParentagePlot$man <- factor(dfParentagePlot$man,
                                levels=dfMen$man[rev(order(dfMen$mannum))])
  
  ## Colorblind palette
  cbPalette <- c("#bbbbbb", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  gParentage <- ggplot(dfParentagePlot,
                       aes(y=man,x=value/tot*100,fill=as.factor(type))) +
    geom_bar(stat='identity',color='black',lwd=.2,width=1) + 
    geom_text(hjust=1,
              aes(label=paste0(round(value/tot*100),""),
                  color=as.factor(type)),
              position=position_stack(vjust = 0.9),size=7*5/14) + 
    geom_label(hjust=0,
               x=2,
               size=7*5/14,
               color='black',
               fill='#FFE5D0',
               label.size=NA,
               label.padding = unit(0.04,'cm'),
              aes(y=man,label=man)) + 
    theme(legend.position='top',
          legend.key.size=unit(0.2,'cm'),
          legend.text=element_text(size=7),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.margin = margin(c(0,0.2,0.2,0),unit='cm'),
          axis.title=element_text(size=7))+
    xlab('Parentage inferred (%)') + 
    ylab('') + 
    guides(fill = guide_legend(reverse = TRUE),
           color=FALSE) + 
    scale_fill_manual('',values=c('Uni/Bi' = cbPalette[7],
                                  'Bi-Parental' = cbPalette[5],
                                  'Uni-Parental' = cbPalette[4],
                                  'None' = cbPalette[1])) + 
    scale_color_manual('',values=c('Uni/Bi' = 'white',
                                   'Bi-Parental' = 'black',
                                   'Uni-Parental' = 'black',
                                   'None' = 'black')) + 
    coord_cartesian(expand=FALSE,clip='off')
  
  return(gParentage)
}

## Supplement
drawSupplementForJeffreys2013 <- function(){
  dfParentage <- fread('parentage_Jeffreys_2013_all.txt',header=TRUE) %>%
    rowwise() %>%
    mutate(shortid=gsub("^(.+v:).+:(\\d+):M\\d+[SB].+","\\1\\2",id,perl=TRUE),
           type = ifelse(type=='Monoallelic','Uni-Parental','Bi-Parental'),
           has_recombinant=ifelse(code=='NORECOMBINANT',FALSE,TRUE),
           mannum=as.numeric(gsub("^M(\\d+)$","\\1",manid,perl=TRUE)),
           mansid=paste0(gsub("^M(\\d+)$","\\1",manid,perl=TRUE),sb),
           mangt1=gsub("^(\\S+)\\/\\S+$","\\1",mangt,perl=TRUE),
           mangt2=gsub("^\\S+\\/(\\S)+$","\\1",mangt,perl=TRUE),
           man=paste0("Man ", mannum, ifelse(sb=="Sperm","S","B"),": ", mangt),
           switchkey=ifelse(p1 == mangt1,switchkey,gsub(":4",":2",gsub(":2",":1",gsub(":1",":4",paste0(":",switchkey))))),
           switchkey=gsub("^:","",switchkey),
           nt=paste0(as.numeric(strsplit(switchkey,":")[[1]]),collapse = ""),
           randn=sample(seq(0,1,0.001),1))
  
  dfUniORBi <- dfParentage %>%
    group_by(id,man,shortid,nswitches) %>% 
    summarise(nd=n_distinct(type)) 
  
  dfBestRecombinant <- dfParentage %>%
    group_by(man,shortid) %>% 
    slice_min(order_by = nswitches+randn) %>% 
    slice_head()
  
  dfParentageProcessed <- inner_join(dfBestRecombinant,
                                     dfUniORBi,
                                     by = c('id','shortid','man','nswitches'))%>%
    mutate(type=ifelse(nd==2,"Uni/Bi",
                       ifelse(mangt1 == mangt2,
                              'Uni/Bi',
                              type)),
           type=ifelse(nswitches==99,'None',type))
  
  for (i in 1:dim(dfParentageProcessed)[1]){
    pid <- as.numeric(strsplit(dfParentageProcessed$switchkey[i],":")[[1]])
    d <- data.frame(p=pid,
                    p1=dfParentageProcessed$p1[i],
                    p2=dfParentageProcessed$p2[i],
                    nid=dfParentageProcessed$shortid[i],
                    man=dfParentageProcessed$man[i],
                    mansid=dfParentageProcessed$mansid[i],
                    mannum=dfParentageProcessed$mannum[i],
                    n=1:length(pid),
                    y=dfParentageProcessed$shortid[i],
                    type=dfParentageProcessed$type[i]) %>% 
      rowwise() %>% 
      mutate(pid = p)
    
    if (dfParentageProcessed$shortid[i] == 'Rv:0075'){
      print(d)
    }
    dType <- d[1,]
    dType$n <- -1

    d <- rbind(dType,d)
    
    if (i == 1){
      dfAll <- d 
    }else{
      dfAll <- rbind(dfAll,d)
    }
  }
  
  dfAll$y   <- factor(dfAll$y,levels=rev(unique(dfAll$y[order(dfAll$nid)])))
  dfAll$man <- factor(dfAll$man,levels=unique(dfAll$man[order(dfAll$mannum)]))
  
  dfAll$pid[dfAll$pid==2 & (dfAll$p1 == dfAll$p2)] <- 1
  
  dSet1 <- dfAll %>% rowwise() %>% dplyr:::filter(mansid %in% c("1Blood","2Sperm","3Sperm","4Sperm"))
  dSet2 <- dfAll %>% rowwise() %>% dplyr:::filter(mansid %in% c("5Sperm","8Sperm","9Sperm","12Sperm"))
  dSet3 <- dfAll %>% rowwise() %>% dplyr:::filter(mansid %in% c("1Sperm","7Sperm","11Sperm","13Sperm"))
  dSet4 <- dfAll %>% rowwise() %>% dplyr:::filter(mannum %in% c(6,10,14,16))
  
  plotMe <- function(dfX) { 
    ## Colorblind palette
    cbPalette <- c("#bbbbbb", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
    g <- ggplot() + 
      geom_tile(data = dfX[dfX$n>0,], 
                aes(x=n,fill=factor(pid),y=y),
                color='black',lwd=.2) + 
      geom_point(data = dfX[dfX$n<0,], 
                aes(x=n,color=factor(type),y=y),
                shape=19,
                size=1.5) +
      scale_fill_manual("Parent",values = c('0'='grey70',
                                            '3'='green',
                                            '2'='dodgerblue1',
                                            '1'='yellow')) + 
      scale_color_manual("Type",values = c('Uni/Bi' = cbPalette[7],
                                           'Bi-Parental' = cbPalette[5],
                                           'Uni-Parental' = cbPalette[4],
                                           'None' = cbPalette[1])) + 
      theme(axis.ticks.y=element_blank(),
            strip.text=element_text(size=7)) + 
      ylab('') + xlab('ZF #') + 
      facet_wrap(~man,scales='free_y',nrow=1)
    
    return(g)
  }
  
  pList  <- list();
  szList <- list();
  for (m in unique(dfAll$man)){
    d <- dfAll %>% rowwise() %>% dplyr:::filter(man == m)  
    pList[[m]] <- plotMe(d)
    szList[[m]] <- dim(d)[1]
  }
  
  gSet1 <- plotMe(dSet1)
  gSet2 <- plotMe(dSet2)
  gSet3 <- plotMe(dSet3)
  gSet4 <- plotMe(dSet4)
  
  noX    <- theme(axis.text.x=element_blank(),plot.margin = unit(c(0,0,0,0),'cm'))
  noLeg  <- theme(legend.position='none')
  topLeg <- theme(legend.position='top', 
                  legend.key.size=unit(0.4,'cm'),
                  legend.margin = margin(c(0,0,0,0),'cm'))
  small <-  theme(axis.text.y=element_text(size=6),
                  plot.margin = unit(c(0,0,0,0),'cm'))
  
  ggarrange(gSet1 + noX + xlab('') + topLeg + small,
            gSet2 + noX + xlab('') + noLeg + small,
            gSet3 + noX + xlab('') + noLeg + small,
            gSet4 + noLeg + small,
            ncol=1,
            heights=c(3,5,10,13),
            align='v')
  
  g14 <- ggarrange(gSet1 + noX + xlab('') + topLeg + small,
                   gSet4 + noLeg + small,
                   ncol=1,
                   heights=c(1.7,5),
                   align='v')
  
  g23 <- ggarrange(gSet2 + noX + xlab('') + topLeg + small,
                   gSet3 + noLeg + small,
                   ncol=1,
                   heights=c(5,8.45),
                   align='v')

  ggsave('Alleva_et_al_RelatednessPlot_Supp_Page1.png',plot = g23, height=6.5,width=6.5,dpi=400)
  ggsave('Alleva_et_al_RelatednessPlot_Supp_Page1.pdf',plot = g23, height=6.5,width=6.5)

  ggsave('Alleva_et_al_RelatednessPlot_Supp_Page2.png',plot = g14, height=8.5,width=6.5,dpi=400)
  ggsave('Alleva_et_al_RelatednessPlot_Supp_Page2.pdf',plot = g14, height=8.5,width=6.5)
  
  return(list(gPage1=g14,
              gPage2=g23))
}

## Parts B & C
drawParentagePlotForPopulationAlleles <- function(){
  #dfInitialData <- fread('PrZFA_relatedness_popAllelesONLY_v3.tab',header=FALSE) %>%
  dfInitialData <- fread('parentage_Pop_alleles_ONLY.txt',header=TRUE) %>%
    dplyr:::rename(child=id) %>%
    mutate(nswitches=ifelse(nswitches==0,NA,nswitches),
           shortid=gsub("^(.+v:).+:(\\d+):M\\d+[SB].+","\\1\\2",child,perl=TRUE),
           pAC=ifelse(p1AC==p2AC,paste0(p1AC,"/",p1AC,"-type"),'A/C-type'),
           parentalAC=pAC,
           childAC=paste0(childAC,"-type"))
  
  dfPopAlleleParentageInit <- dfInitialData %>%
    select(shortid,childAC,parentalAC,nswitches) %>%
    group_by(shortid,childAC,parentalAC,nswitches) %>% 
    dplyr:::count() 
  
  dfPopAlleleParentage <- setAlleleOrder(dfPopAlleleParentageInit,'shortid',rev=TRUE)
  
  gAllAtype <- ggplot(dfPopAlleleParentage %>% 
                        dplyr:::filter(childAC=='A-type'),
                      aes(x=nswitches,y=shortid,fill=log10(n))) + 
    geom_tile(color='black') + 
    facet_grid(paste0("Child = ",childAC)~parentalAC) + 
    scale_fill_viridis_c() +
    theme(panel.grid.major.y = element_line(size=.2,color='grey80'),
          strip.text=element_text(size=7),
          panel.border = element_rect(fill=NA)) + 
    ylab('') +
    scale_x_continuous(breaks=c(1:12),
                       labels=c(1,'','','','',6,'','','','',11,'')) 
  
  gAllCtype <- ggplot(dfPopAlleleParentage %>% 
                        dplyr:::filter(childAC=='C-type'),
                      aes(x=nswitches,y=shortid,fill=log10(n))) + 
    geom_tile(color='black') + 
    facet_grid(paste0("Child = ",childAC)~parentalAC) + 
    scale_fill_viridis_c() +
    theme(panel.grid.major.y = element_line(size=.2,color='grey80'),
          strip.text=element_text(size=7),
          panel.border = element_rect(fill=NA)) + 
    ylab('') +
    scale_x_continuous(breaks=c(1:12),
                       labels=c(1,'','','','',6,'','','','',11,''))
  
  noX      <- theme(axis.text.x=element_blank())
  noM      <- theme(plot.margin = unit(c(0,0,0,0),'cm'))
  noXStrip <- theme(strip.background.x = element_blank(),strip.text.x=element_blank())
  noYStrip <- theme(strip.background.y = element_blank(),strip.text.y=element_blank())
  noLeg    <- theme(legend.position='none')
  custLeg  <- theme(legend.position=c(1,0.5),
                    legend.justification = c(1.05,.5),
                    legend.key.size=unit(0.4,'cm'),
                    legend.key.width = unit(0.2,'cm'),
                    legend.background = element_rect(fill='white'))
  
  custLeg  <- theme(legend.position='top',
                    legend.direction='horizontal',
                    legend.key.height = unit(0.1,'cm'),
                    legend.key.width = unit(0.5,'cm'),
                    legend.background = element_rect(fill='white'))
  
  dfAggregateBarPlots <- dfInitialData %>%
    select(childAC,parentalAC,nswitches) %>%
    group_by(childAC,parentalAC,nswitches) %>% 
    dplyr:::count() 
  
  gQAtype <- ggplot(dfAggregateBarPlots %>% dplyr:::filter(childAC == "A-type"),
                    aes(x=nswitches,y=n)) + 
    geom_bar(color='black',stat='identity',lwd=.2) + 
    facet_grid(paste0("Child = ",childAC)~parentalAC,scales='free_y') +
    theme(panel.border = element_rect(fill=NA),
          strip.text=element_text(size=7)) + 
    ylab('n') +  
    noXStrip + 
    noM +
    xlab('') + 
    coord_cartesian(xlim=c(1,12)) + 
    scale_x_continuous(breaks=c(1:12),
                       labels=c(1,'','','','',6,'','','','',11,''))
  
  gQCtype <- ggplot(dfAggregateBarPlots %>% dplyr:::filter(childAC == "C-type"),
                    aes(x=nswitches,y=n)) + 
    geom_bar(color='black',stat='identity',lwd=.2) + 
    facet_grid(paste0("Child = ",childAC)~parentalAC,scales='free_y') +
    theme(panel.border = element_rect(fill=NA),
          strip.text=element_text(size=7)) + 
    ylab('n') +  
    noXStrip + 
    noM + 
    xlab('') + 
    coord_cartesian(xlim=c(1,12)) + 
    scale_x_continuous(breaks=c(1:12),
                       labels=c(1,'','','','',6,'','','','',11,''))
  
  gAtypeMerge <- ggarrange(gAllAtype + noYStrip + noM + noX + custLeg + xlab('') + coord_cartesian(xlim=c(1,12)),
                           gQAtype + noYStrip + xlab('Template switches required (#)'),
                           ncol=1,nrow=2,
                           heights=c(8,2),
                           align='v',
                           labels=c(''),
                           hjust=0,vjust=1,
                           font.label = list(size=8,face='bold'))
  
  gCtypeMerge <- ggarrange(gAllCtype + noYStrip + noM + noX + noLeg + xlab('') + coord_cartesian(xlim=c(1,12)),
                           gQCtype + noYStrip + xlab('Template switches required (#)') ,
                           ncol=1,nrow=2,
                           heights=c(4,2),
                           align='v',
                           labels=c(''),
                           hjust=0,vjust=1,
                           font.label = list(size=8,face='bold'))
  
  return(list(gA=gAtypeMerge,
              gC=gCtypeMerge))
}

###############################

theme7point()

###############################

gJeffreys <- drawParentageForJeffreys2013()
gJeffSupp <- drawSupplementForJeffreys2013()
lstPop    <- drawParentagePlotForPopulationAlleles()

#################################################################################

gLeft <- ggarrange(gJeffreys,
                   lstPop$gC,
                   ncol=1,nrow=2,
                   heights=c(1,1.4),
                   labels=c('A','C'),
                   hjust=0,vjust=1,
                   font.label = list(size=8,face='bold'))

gAll <- ggarrange(gLeft,
                  lstPop$gA,
                  ncol=2,nrow=1,
                  labels=c('','B'),
                  hjust=0,vjust=1,
                  font.label = list(size=8,face='bold'))

ggsave('Alleva_et_al_RelatednessPlot.png',height=7,width=6.5,dpi=400)
ggsave('Alleva_et_al_RelatednessPlot.pdf',height=7,width=6.5)

