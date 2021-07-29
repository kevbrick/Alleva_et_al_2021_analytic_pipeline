library(ggplot2)
library(plyr)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(RColorBrewer)
library(ggtext)
library(data.table)

source('accessoryFiles/scripts/genericFunctions.R')

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

prdm9_allele_details <- read.table('atype_ctype.txt',header=TRUE,stringsAsFactors = FALSE,skip = 11) %>% 
  mutate(allele=short_ID)

atypes <- unique(prdm9_allele_details$short_ID[prdm9_allele_details$ACtype == 'A'])
ctypes <- unique(prdm9_allele_details$short_ID[prdm9_allele_details$ACtype == 'C'])

prdm9_DF <- read.table('prdm9_haplotypes.final.tab',header=TRUE,stringsAsFactors = FALSE)

## Modify Jeffreys Sperm/Blood alleles name
prdm9_DF$allele_1 <- gsub(pattern = ':[sc]:([0-9]+):.+$',':\\1',prdm9_DF$allele_1)
prdm9_DF$allele_2 <- gsub(pattern = ':[sc]:([0-9]+):.+$',':\\1',prdm9_DF$allele_2)
prdm9_DF$diploid  <- gsub(pattern = ':[sc]:([0-9]+):.+/(.+)$',':\\1/\\2',prdm9_DF$diploid)
prdm9_DF$diploid  <- gsub(pattern = '^(.+)/(.+):[sc]:([0-9]+):.+$','\\1/\\2:\\3',prdm9_DF$diploid)

## Set order of Alleles / ZFAs
allelesAll  <- as.character(sort(unique(sort(rbind(prdm9_DF$allele_1,prdm9_DF$allele_2)))))
alleleOrder <- c(str_sort(allelesAll[grep('^[A-Z]$',perl=TRUE,allelesAll)],numeric=TRUE),
                 str_sort(allelesAll[grep('^L[0-9]+$',perl=TRUE,allelesAll)],numeric=TRUE),
                 str_sort(allelesAll[grep('^(A|C|L[0-9]+)v',perl=TRUE,allelesAll)],numeric=TRUE),
                 str_sort(allelesAll[grep('^M[0-9]+$',perl=TRUE,allelesAll)],numeric=TRUE),
                 'Unk','noZFA')

allGenotypes  <- as.character(sort(unique(sort(prdm9_DF$diploid))))
genotypeOrder <- c(str_sort(allGenotypes[grep('^[A-Z][0-9]*/[A-Z][0-9]*$',perl=TRUE,allGenotypes)],numeric=TRUE),
                   str_sort(allGenotypes[grep('^[A-Z][0-9]*/(A|C|L[0-9]+)v.+$',perl=TRUE,allGenotypes)],numeric=TRUE),
                   str_sort(allGenotypes[grep('^(A|C|L[0-9]+)v.+/',perl=TRUE,allGenotypes)],numeric=TRUE))

### FOR FULL ALT ALLELE NAMES: 
# allelesAll  <- as.character(sort(unique(sort(rbind(prdm9_DF$allele_1,prdm9_DF$allele_2)))))
# alleleOrder <- c(str_sort(allelesAll[grep('^[A-Z]$',perl=TRUE,allelesAll)],numeric=TRUE),
#                  str_sort(allelesAll[grep('^L[0-9]+$',perl=TRUE,allelesAll)],numeric=TRUE),
#                  str_sort(allelesAll[grep('^(A|C|L[0-9]+)v',perl=TRUE,allelesAll)],numeric=TRUE),
#                  str_sort(allelesAll[grep('^M[0-9]+$',perl=TRUE,allelesAll)],numeric=TRUE),
#                  'Unk','noZFA')
# 
# allGenotypes  <- as.character(sort(unique(sort(prdm9_DF$diploid))))
# genotypeOrder <- c(str_sort(allGenotypes[grep('^[A-Z][0-9]*/[A-Z][0-9]*$',perl=TRUE,allGenotypes)],numeric=TRUE),
#                    str_sort(allGenotypes[grep('^[A-Z][0-9]*/(A|C|L[0-9]+)v.+$',perl=TRUE,allGenotypes)],numeric=TRUE),
#                    str_sort(allGenotypes[grep('^(A|C|L[0-9]+)v.+/',perl=TRUE,allGenotypes)],numeric=TRUE))

prdm9_DF$ok <- TRUE
prdm9_DF$ok[prdm9_DF$diploid == '/'] <- FALSE
prdm9_DF$ok[grepl('Unk|noZF|NA',prdm9_DF$diploid)] <- FALSE

prdm9_DF$typeHH  <- 'het'
prdm9_DF$typeHH[prdm9_DF$size_1 == prdm9_DF$size_2 & prdm9_DF$allele_1 == prdm9_DF$allele_2]  <- 'hom'
prdm9_DF$typeHH <- factor(prdm9_DF$type,levels=c('het','hom'))

prdm9_DF$type  <- 'unk'
prdm9_DF$type[prdm9_DF$size_1 != prdm9_DF$size_2]  <- 'Het (unequal #ZFs)'
prdm9_DF$type[prdm9_DF$size_1 == prdm9_DF$size_2 & prdm9_DF$allele_1 != prdm9_DF$allele_2]  <- 'Het (equal #ZFs)'
prdm9_DF$type[prdm9_DF$size_1 == prdm9_DF$size_2 & prdm9_DF$allele_1 == prdm9_DF$allele_2]  <- 'Hom'

prdm9_DF$type <- factor(prdm9_DF$type,levels=c('unk','Het (unequal #ZFs)','Het (equal #ZFs)','Hom'))

prdm9_DF       <- prdm9_DF %>% mutate(pop=ifelse(is_child,paste0(pop,"(T)"),pop))

## SET ORDER
popOrder <- c('LWK','YRI','YRI(T)','CHB','PJL','PEL','TSI','FIN')

prdm9_DF$pop <- factor(prdm9_DF$pop,popOrder)

### Pops DF
popsDF <- data.frame(pop = c('CHB','FIN','LWK','PEL','PJL','TSI','YRI'),
                     tot = c(120,100,120,70,108,114,179))

popsDF <- data.frame(pop = c('CHB','FIN','LWK','PEL','PJL','TSI','YRI','YRI(T)'),
                     tot = c(120,100,120,70,108,114,120,59))

popsDF$pop <- factor(popsDF$pop,popOrder)

## Make parent DF ##############################################################
plotDF <- prdm9_DF %>% filter(ok & pop != 'OTH') %>%
  select(pop,diploid,allele_1,allele_2,type,typeHH) 

plotDF$pop <- factor(plotDF$pop,popOrder)

## Figure drawing functions ####################################################
drawGenotypingSuccessFigs <- function(popOrder = c('LWK','YRI',"YRI(T)",'CHB','PJL','PEL','TSI','FIN'), drawFigs = FALSE){

  plotDF <- prdm9_DF %>% filter(ok & pop != 'OTH') %>%
    select(pop,diploid,allele_1,allele_2,type,typeHH,is_child)
  
  plotDF$pop <- factor(plotDF$pop,popOrder)
  
  ## Fig 1 #####################################################################
  plot1DF <- reshape2:::melt.data.frame(plotDF[,c('diploid','pop')] %>% 
                    group_by(pop) %>% 
                    count() %>% 
                    inner_join(popsDF) %>% 
                    mutate(yes=n,no=tot-n) %>% 
                    select(pop,yes,no))
  
  plot1DF$pop <- factor(plot1DF$pop,popOrder)
  
  plot1DF$variable <- factor(plot1DF$variable,levels=c('no','yes'))
  gGTCounts <- ggplot(plot1DF,aes(y=pop,x=value,fill=variable)) + 
    geom_bar(stat='identity',width=.75,lwd=.2) + 
    geom_text(aes(label=value),
              color='white',
              size=8*5/14,
              position = position_stack(vjust = 0.5)) +
    scale_fill_manual('PRDM9 diploid genotype inferred ?',
                      values=c('darkorange','grey40'),
                      guide = guide_legend(reverse = TRUE)) + 
    xlab('Individuals (#)') + 
    ylab('') + 
    theme(legend.position='top')
  
  ## Fig 2 #####################################################################
  gtsDF <- prdm9_DF %>% filter(pop != 'OTH') %>%
    select(pop,diploid,allele_1,allele_2,type,typeHH,ok,nseqs,nzfseqs) 
  
  gtsDF$pop <- factor(gtsDF$pop,popOrder)
  gZFvsStatus <- ggplot(gtsDF,aes(x=pop,y=nzfseqs,fill=ok)) + 
    geom_boxplot(position=position_dodge()) + 
    scale_y_log10()+
    scale_fill_manual('PRDM9 diploid genotype inferred ?',
                      values=c('grey40','darkorange')) + 
    xlab('') + 
    ylab('ZF array reads (#)') + 
    theme(legend.position='top')
  
  gtsDF$yesno <- 'no'
  gtsDF$yesno[gtsDF$ok] <- 'yes'
  gtsDF$yesno <- factor(gtsDF$yesno,levels=c('no','yes'))
  
  gZFvsStatusH <- ggplot(gtsDF,aes(x=pop,y=nzfseqs,fill=yesno)) + 
    geom_boxplot(position=position_dodge(),alpha=.2,lwd=.2,outlier.alpha=0) + 
    geom_point(color=alpha('white',0),size = 0.7, shape = 21, position = position_jitterdodge(),lwd=.2) +
    scale_y_log10(labels=fancy_scientific) +
    scale_fill_manual('PRDM9 diploid genotype inferred ?',
                      values=c('yes'='grey40','no'='darkorange')) +
    xlab('') + 
    ylab('Reads with valid ZF array (#)') + 
    theme(legend.position='top') + 
    annotation_logticks(sides='l',
                        size=.2,
                        short = unit(0.1,'cm'),
                        mid = unit(0.1,'cm'),
                        long = unit(0.2,'cm'))
  
  gZFvsStatusV <- ggplot(gtsDF,aes(y=pop,x=nzfseqs,fill=yesno)) + 
    geom_boxplot(position=position_dodge(),alpha=.2,lwd=.2,outlier.alpha=0) + 
    geom_point(color=alpha('white',0),size = 0.7, shape = 21, position = position_jitterdodge(),lwd=.2) +
    scale_x_log10(labels=fancy_scientific) +
    scale_fill_manual('PRDM9 diploid genotype inferred ?',
                      values=c('yes'='grey40','no'='darkorange')) +
    ylab('') + 
    xlab('Reads with valid ZF array (#)') + 
    theme(legend.position='top') + 
    annotation_logticks(sides='b',
                        size=.2,
                        short = unit(0.1,'cm'),
                        mid = unit(0.1,'cm'),
                        long = unit(0.2,'cm'))
  
  ## Fig 3 #####################################################################
  gLowCoverageH <- ggplot(gtsDF,aes(x=pop,y=nseqs,fill=yesno)) + 
    geom_boxplot(position=position_dodge(),alpha=.2,lwd=.2,outlier.alpha=0) + 
    geom_point(color=alpha('white',0),size = 0.7, shape = 21, position = position_jitterdodge(),lwd=.2) +
    scale_y_log10(labels=fancy_scientific) +
    scale_fill_manual('PRDM9 diploid genotype inferred ?',
                      values=c('yes'='grey40','no'='darkorange')) +
    xlab('') + 
    ylab('Raw reads (#)') + 
    theme(legend.position='top') + 
    annotation_logticks(sides='l',
                        size=.2,
                        short = unit(0.1,'cm'),
                        mid = unit(0.1,'cm'),
                        long = unit(0.2,'cm'))
  
  gLowCoverageV <- ggplot(gtsDF,aes(y=pop,x=nseqs,fill=yesno)) + 
    geom_boxplot(position=position_dodge(),alpha=.2,lwd=.2,outlier.alpha=0) + 
    geom_point(color=alpha('white',0),size = 0.7, shape = 21, position = position_jitterdodge(),lwd=.2) +
    scale_x_log10(labels=fancy_scientific) +
    scale_fill_manual('PRDM9 diploid genotype inferred ?',
                      values=c('yes'='grey40','no'='darkorange')) +
    ylab('') + 
    xlab('Raw reads (#)') + 
    theme(legend.position='top') + 
    annotation_logticks(sides='b',
                        size=.2,
                        short = unit(0.1,'cm'),
                        mid = unit(0.1,'cm'),
                        long = unit(0.2,'cm'))
  
  ## Supp Fig 1 ################################################################
  gLeg    <- get_legend(gGTCounts + theme(legend.background=element_blank(),
                                          legend.key.size=unit(0.3,'cm'))) 
  
  gLabels <- ggarrange(ggplot+theme_void() + theme(plot.margin=unit(c(0,0,.2,0),'cm')),
                       ggplot+theme_void() + theme(plot.margin=unit(c(0,0,.2,0),'cm')),
                       ggplot+theme_void() + theme(plot.margin=unit(c(0,0,0,0),'cm')),
                       ncol=3,nrow=1,
                       labels=c('E','F','G'),
                       font.label = list(size=8,face='bold',hjust=0,vjust=1))
  
  gSuppData <- ggarrange(gGTCounts     + theme(plot.margin=unit(c(0,0,0,0),'cm'), legend.position='none'),
                         gLowCoverageV + theme(plot.margin=unit(c(0,0,0,0),'cm'), legend.position='none', axis.text.y=element_blank()),
                         gZFvsStatusV  + theme(plot.margin=unit(c(0,0,0,0),'cm'), legend.position='none', axis.text.y=element_blank()),
                       ncol=3,nrow=1)

  gSuppGT <- ggarrange(gLabels,gLeg,gSuppData,nrow = 3,ncol=1,heights=c(1,1,20))  
  if (drawFigs){
    ggsave(paste0('PRDM9_diploid_GT_Counts.png'),gGTCounts,width=4,height=3)
    ggsave(paste0('PRDM9_diploid_GT_Counts.pdf'),gGTCounts,width=4,height=3)
    
    ggsave(paste0('PRDM9_ZFs_per_person.png'),gZFvsStatus,width=4,height=3)
    ggsave(paste0('PRDM9_ZFs_per_person.pdf'),gZFvsStatus,width=4,height=3)

    ggsave(paste0('PRDM9_raw_reads_per_person.png'),gLowCoverageH,width=4,height=3)
    ggsave(paste0('PRDM9_raw_reads_per_person.pdf'),gLowCoverageH,width=4,height=3)    
    
    ggsave(paste0('Alleva_et_al_Genotyping_By_Pop_VS_Reads.png'),gSuppGT,width=7,height=2)
    ggsave(paste0('Alleva_et_al_Genotyping_By_Pop_VS_Reads.pdf'),gSuppGT,width=7,height=2)
  }
  
  return(list(gSupp = gSuppGT, 
              gFig1 = gGTCounts,
              gFig2 = gZFvsStatus,
              gFig3 = gLowCoverageH))
}

drawDiploidGenotypes <- function(drawFigs=FALSE){
  
  ## Omit children here as it will bias pop estimates
  ## (Only relevant for YRI)
  plotDF <- prdm9_DF %>% 
    filter(ok & pop != 'OTH' & !is_child) %>%
    select(pop,diploid,allele_1,allele_2,type,typeHH) 
  
  ## Make DF ###################################################################
  plot2DFi <- as.data.frame(plotDF[,c('diploid','pop')] %>% 
                              group_by(pop) %>% 
                              add_count(name='popTot') %>% 
                              group_by(pop,diploid,popTot) %>% 
                              count(name='gtPopTot') %>% 
                              group_by(pop,diploid,popTot,gtPopTot) %>% 
                              summarize(n=gtPopTot,
                                        popTot=popTot,
                                        pc=gtPopTot/popTot*100,
                                        pcRnd=round(gtPopTot/popTot*100/10)*10))
  
  dfGenotypes      <- data.frame(diploid=genotypeOrder,grp=0)
  dfGenotypes$grp  <- rep(1:3,each=ceiling(length(dfGenotypes$grp)/3))[1:length(dfGenotypes$diploid)]
  
  plot2DF <- join(plot2DFi,dfGenotypes,by='diploid')
  plot2DF$pop <- factor(plot2DF$pop,popOrder)
  
  ord <- genotypeOrder[genotypeOrder %in% unique(sort(plot2DF$diploid))]
  
  plot2DF$diploid <- factor(plot2DF$diploid,rev(ord))
  
  ## Make Diploid genotypes figs #######################################################
  ## This bit is to add color to y-axis labels
  #  Uses ggtext to replace axis text with HTML span elements
  plot2DF <- plot2DF %>%
    mutate(
      allele1=gsub("^(\\S+)\\/.+$","\\1",diploid),
      allele2=gsub("^\\S+\\/(.+)$","\\1",diploid),
      y.label = paste("<span style = 'color: ",
                      ifelse(allele1 %in% atypes, "magenta", "forestgreen"),
                      ";'>",
                      allele1,
                      "</span>", 
                      "<span style = 'color: black;'>/</span>",
                      "<span style = 'color: ",
                      ifelse(allele2 %in% atypes, "magenta", "forestgreen"),
                      ";'>",
                      allele2,
                      "</span>", 
                      sep = ""))
  
  ## Order span levels
  lvl <- levels(plot2DF$diploid)
  n = 1
  
  for (x in levels(plot2DF$diploid)){
    lvl[n] <- plot2DF$y.label[plot2DF$diploid == x][1]
    n = n+1
  }
  
  plot2DF$y.label <- factor(plot2DF$y.label,lvl)
  
  ggplot(plot2DF,aes(x=pop,y=y.label,fill=pcRnd)) + 
    geom_tile(lwd=.05,color='black') + 
    geom_text(aes(label=round(pc,0)),
              size=7*5/14) +
    scale_fill_gradient('Freq(%)',high = 'firebrick1', low='white', na.value = 'white') + 
    xlab('Population') + ylab('Diploid PRDM9 genotype') +
    ggtitle('Genotype frequency (%)') +
    theme(strip.text=element_blank(),
          strip.background=element_blank(),
          legend.position = 'none',
          plot.title = element_text(size=7,hjust=0),
          axis.text.x=element_text(angle=90,hjust=1),
          axis.text.y = element_markdown()) + 
    facet_wrap(~grp,ncol=3,scales='free_y')

  gGenotypePCNoWrap <- ggplot(plot2DF,aes(x=pop,y=y.label,fill=pcRnd)) + 
    geom_tile(lwd=.05,color='black') + 
    geom_text(aes(label=round(pc,0)),
              size=7*5/14) +
    scale_fill_gradient('Freq(%)',high = 'firebrick1', low='white', na.value = 'white') + 
    xlab('Population') + ylab('Diploid PRDM9 genotype') +
    ggtitle('Genotype frequency (%)') +
    theme(strip.text=element_blank(),
          strip.background=element_blank(),
          legend.position = 'none',
          plot.title = element_text(size=7,hjust=0),
          axis.text.x=element_text(angle=90,hjust=1),
          axis.text.y = element_markdown())
  
  ## Original ... uncolored labels
  # gGenotypePCNoWrap <- ggplot(plot2DF,aes(x=pop,y=diploid,fill=pcRnd)) + 
  #   geom_tile(lwd=.05,color='black') + 
  #   geom_text(aes(label=round(pc,0)),
  #             size=7*5/14) +
  #   scale_fill_gradient('Freq(%)',high = 'firebrick1', low='white', na.value = 'white') + 
  #   xlab('Population') + ylab('Diploid PRDM9 genotype') +
  #   ggtitle('Genotype frequency (%)') +
  #   theme(strip.text=element_blank(),
  #         strip.background=element_blank(),
  #         legend.position = 'none',
  #         plot.title = element_text(size=7,hjust=0),
  #         axis.text.x=element_text(angle=90,hjust=1))
  
  gGenotypePC <- gGenotypePCNoWrap + facet_wrap(~grp,ncol=3,scales='free_y')
  
  gGenotypeNNoWrap <- ggplot(plot2DF,aes(x=pop,y=y.label,fill=n)) + 
    geom_tile(lwd=.05,color='black') + 
    geom_text(aes(label=round(n,0)),
              size=7*5/14) +
    scale_fill_gradient('Freq(%)',high = 'dodgerblue2', low='white', na.value = 'white') + 
    xlab('Population') + ylab('Diploid PRDM9 genotype') +
    ggtitle('Genotype frequency (Count)') +
    theme(strip.text=element_blank(),
          strip.background=element_blank(),
          legend.position = 'none',
          axis.text.x=element_text(angle=90,hjust=1),
          axis.text.y = element_markdown())
  
  gGenotypeN <- gGenotypeNNoWrap + facet_wrap(~grp,ncol=3,scales='free_y')
  
  ## Make HomHet fig ###########################################################
  # plot7DF <- plotDF %>%  select(pop,type) %>% group_by(pop,type) %>% count() %>%  
  #   inner_join(plotDF %>% group_by(pop) %>% count(name = "tot")) %>% 
  #   summarise(N=tot, pc=n/tot*100)
  
  plot7DF <- plotDF %>%  select(pop,type) %>% mutate(type=ifelse(type=='Hom','Hom','Het')) %>%
    group_by(pop,type) %>% count() %>%  
    inner_join(plotDF %>% group_by(pop) %>% count(name = "tot")) %>% 
    summarise(N=tot, pc=n/tot*100)
  
  plot7DF$pop <- factor(plot7DF$pop,levels=popOrder)
  
  gHomHet <- ggplot(plot7DF[plot7DF$type == 'Het',],aes(x=pc,y=pop)) + 
    geom_bar(stat='identity',width=.8,color='grey50',lwd=.1,fill='grey50') + 
    xlab('PRDM9 heterozygosity (%)') + ylab('') + 
    theme(legend.position='none',
          legend.key.size = unit(0.1,'cm'),
          legend.text = element_text(size=7))
  
  ## Draw figs #################################################################
  if (drawFigs){
    ggsave(paste0('PRDM9_diploid_GT_by_Pop_Freqs.png'),gLong,width=14,height=4)
    ggsave(paste0('PRDM9_diploid_GT_by_Pop_Freqs.pdf'),gLong,width=14,height=4)
    
    ggsave(paste0('PRDM9_diploid_GT_by_Pop_Freqs_wide.png'),gDiploidGT_By_Pop,width=7,height=6,dpi=500)
    ggsave(paste0('PRDM9_diploid_GT_by_Pop_Freqs_wide.pdf'),gDiploidGT_By_Pop,width=7,height=6,dpi=500)
  }
  
  return(list(gFig        = gGenotypePC,
              gFigN       = gGenotypeN,
              gFigNoWrap  = gGenotypePCNoWrap,
              gFigNNoWrap = gGenotypeNNoWrap,
              gHomHet     = gHomHet))
}

drawAlleleFrequencies <- function(drawFigs=FALSE){
  
  popOrder <- popOrder[popOrder != 'YRI(T)']
  
  ## Omit children here as it will bias pop estimates
  ## (Only relevant for YRI)
  plotDF <- prdm9_DF %>% 
    filter(ok & pop != 'OTH' & !is_child) %>%
    select(pop,diploid,allele_1,allele_2,type,typeHH) 
  
  plotDF$pop <- factor(plotDF$pop,levels=popOrder)
  
  ## Make DF ###################################################################
  plot4DF <- reshape2:::melt.data.frame(plotDF %>% 
                    select(pop,allele_1,allele_2),#[,c('pop','allele_1','allele_2')],
                  id.vars='pop',
                  value.name='allele') %>%
    group_by(pop) %>% 
    add_count(name='tot') %>% 
    group_by(pop,allele,tot) %>% 
    count() %>% 
    group_by(pop,allele,tot) %>% 
    left_join(plotDF %>% count(pop, name = "totPop")) %>%
    summarize(N=n,
              pc=n/tot*100,
              pcRnd=round(n/tot*100/10)*10,
              popPC=n/totPop*100) %>%
    inner_join(prdm9_allele_details %>% select(allele,ACtype), by='allele') %>%
    dplyr:::filter(!(allele %in% c('Unk','noZFA','Av:0053','L13')))
  
  plot4DF$allele <- factor(plot4DF$allele,levels=rev(alleleOrder))
  plot4DF$pop    <- factor(plot4DF$pop,popOrder)
  
  dfAlleles      <- data.frame(allele=rev(levels(plot4DF$allele)),grp=0)
  dfAlleles$grp  <- rep(1:3,each=ceiling(length(dfAlleles$grp)/3))[1:length(dfAlleles$allele)]
  
  p              <- join(plot4DF,dfAlleles,by='allele')
  plot4DF        <- p
  
  ##############################################################################
  ## This bit is to add color to y-axis labels
  #  Uses ggtext to replace axis text with HTML span elements
  plot4DF <- plot4DF %>%
    mutate(y.label = paste("<span style = 'color: ",
                           ifelse(allele %in% atypes, "magenta", "forestgreen"),
                           ";'>", allele, "</span>", sep = ""))
  
  ## Order span levels
  fAllele <- levels(plot4DF$allele)
  fAllele <- fAllele[!(fAllele %in% c('Unk','noZFA','Av:0053','L13'))]
  lvl <- vector(length=length(fAllele))
  n = 1
  
  for (x in fAllele){
    lvl[n] <- plot4DF$y.label[plot4DF$allele == x][1]
    n = n+1
  }
  
  plot4DF$y.label <- factor(plot4DF$y.label,lvl)
  ### Done adding label colors
  
  ## Fig 1 : Allele freqs (%) ##################################################
  gAlleleFreqPCNoWrap <- 
    ggplot(plot4DF,aes(x=pop,y=y.label,fill=pc)) + 
    geom_tile(lwd=.05,color='black') + 
    geom_text(aes(label=round(pc,0)),
              size=7*5/14) +
    scale_fill_gradient('%',high = 'firebrick1', 
                        low='white', 
                        na.value = 'white') + 
    xlab('Population') + ylab('PRDM9 allele') + 
    ggtitle('Allele frequency (%)') + 
    theme(strip.text=element_blank(),
          strip.background=element_blank(),
          legend.position = 'none',
          axis.text.x=element_text(angle=90,hjust=1),
          axis.text.y=element_markdown())
  
  gAlleleFreqPC <- gAlleleFreqPCNoWrap + facet_wrap(~grp,ncol=3,scales='free_y')
  
  ## Fig 2 : Allele counts #####################################################
  gAlleleFreqCntNoWrap <- ggplot(plot4DF,aes(x=pop,y=y.label,fill=N)) + 
    geom_tile(lwd=.05,color='black') + 
    geom_text(aes(label=round(N,0)),
              size=7*5/14) +
    scale_fill_gradient('%',high = 'dodgerblue2', 
                        low='white', 
                        na.value = 'white') + 
    xlab('Population') + ylab('PRDM9 allele') + 
    ggtitle('Allele frequency (Count)') + 
    theme(strip.text=element_blank(),
          strip.background=element_blank(),
          legend.position = 'none',
          axis.text.x=element_text(angle=90,hjust=1),
          axis.text.y=element_markdown())
  
  gAlleleFreqCnt <- gAlleleFreqCntNoWrap + 
    facet_wrap(~grp,ncol=3,scales='free_y') 
  
  ## Supp Fig. Pop specific alleles
  dfPopSpecific <- as.data.frame(plot4DF %>% mutate(n=N) %>% select(pop,tot,allele,n))
  
  lstSupp <- drawSamplingErrorFig(dfPopSpecific,nIterations = 1000)
  
  ggsave('Alleva_et_al_SuppFig_alleleFrequency_with_CIs_from_bootstrapped_sampling.png',plot=lstSupp$fig,dpi=500,width=6.5,height=9)
  ggsave('Alleva_et_al_SuppFig_alleleFrequency_with_CIs_from_bootstrapped_sampling.pdf',plot=lstSupp$fig,width=6.5,height=9)
  
  ## Fig 3 : Pop with one copy (%) #############################################
  dfA1  <- plotDF[,c('allele_1','pop')]; names(dfA1) <- c('allele','pop')
  dfA2  <- plotDF[,c('allele_2','pop')]; names(dfA2) <- c('allele','pop')
  dfA3 <- dfA2[dfA2$allele != dfA1$allele,]
  
  plot6DFi <- rbind(dfA1,dfA3)  %>% 
    group_by(pop,allele) %>% add_count() %>% 
    inner_join(popsDF) %>%
    summarise(N=n, 
              popPC = n/tot*100)
  
  plot6DF <- join(plot6DFi,dfAlleles,by='allele')
  plot6DF$allele <- factor(plot6DF$allele,levels=rev(alleleOrder))
  plot6DF$pop    <- factor(plot6DF$pop,popOrder)
  
  ##############################################################################
  ## This bit is to add color to y-axis labels
  #  Uses ggtext to replace axis text with HTML span elements
  plot6DF <- plot6DF %>%
    mutate(y.label = paste("<span style = 'color: ",
                           ifelse(allele %in% atypes, "magenta", "forestgreen"),
                           ";'>", allele, "</span>", sep = ""))
  
  ## Order span levels
  fAllele <- levels(plot6DF$allele)
  fAllele <- fAllele[!(fAllele %in% c('Unk','noZFA','Av:0053','L13'))]
  lvl <- vector(length=length(fAllele))
  n = 1
  
  for (x in fAllele){
    lvl[n] <- plot6DF$y.label[plot6DF$allele == x][1]
    n = n+1
  }
  
  plot6DF$y.label <- factor(plot6DF$y.label,lvl)
  ### Done adding label colors
  gOneCopyNoWrap <- ggplot(plot6DF,aes(x=pop,y=y.label,fill=popPC)) + 
    geom_tile(lwd=.05,color='black') + 
    geom_text(aes(label=round(popPC,0)),
              size=7*5/14,check_overlap = TRUE) +
    scale_fill_gradient(high = 'firebrick1', low='white', na.value = 'white') + 
    xlab('Population') + ylab('PRDM9 allele') + 
    ggtitle('Carriers (at least one allele) (%)') + 
    theme(strip.text=element_blank(),
          strip.background=element_blank(),
          legend.position = 'none',
          axis.text.x=element_text(angle=90,hjust=1),
          axis.text.y = element_markdown())
  
  gOneCopy <- gOneCopyNoWrap + facet_wrap(~grp,ncol=3,scales='free_y') 
  
  gOneCopyCntNoWrap <- ggplot(plot6DF,aes(x=pop,y=y.label,fill=N)) + 
    geom_tile(lwd=.05,color='black') + 
    geom_text(aes(label=round(N,0)),
              size=7*5/14,check_overlap = TRUE) +
    scale_fill_gradient(high = 'dodgerblue2', low='white', na.value = 'white') + 
    xlab('Population') + ylab('PRDM9 allele') + 
    ggtitle('Carriers (at least one allele) (Count)') + 
    theme(strip.text=element_blank(),
          strip.background=element_blank(),
          legend.position = 'none',
          axis.text.x=element_text(angle=90,hjust=1),
          axis.text.y = element_markdown())
  
  gOneCopyCnt <- gOneCopyCntNoWrap + facet_wrap(~grp,ncol=3,scales='free_y')  
    
  ## Pie Charts ################################################################
  n <- length(alleleOrder)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  col_DF <- data.frame(allele=alleleOrder,col=col_vector[1:n])
  
  lstPies <- list()
  
  for (pop in popOrder[popOrder != 'YRI(T)']){
    
    popAlleles <- unique(plot4DF$allele[plot4DF$pop == pop])
    colz <- col_DF$col[col_DF$allele %in% popAlleles]
    
    gPie <- ggplot(plot4DF[plot4DF$pop == pop,],aes(x="",y=N,fill=allele)) + 
      geom_bar(stat='identity',color='grey50',lwd=0.05) + 
      coord_polar("y", start=0) +
      scale_fill_manual('',values=as.character(rev(colz)),guide=guide_legend(reverse = TRUE)) +
      theme_void() + 
      theme(legend.position='none',
            plot.margin = unit(c(0,0,0,0),'cm')) 
    
    lstPies[[pop]] <- gPie
  }
  
  ## Stacked Bar Charts ########################################################
  plot4DF$simpleAllele <- paste0("Oth (",plot4DF$ACtype,'-type)')
  plot4DF$simpleAllele[plot4DF$allele %in% c('A','B','C','L14')] <- as.character(plot4DF$allele[plot4DF$allele %in% c('A','B','C','L14')])
  
  plot4DF$simpleAllele <- factor(plot4DF$simpleAllele,levels=rev(c('A','B','C','L14','Other','Oth (A-type)','Oth (C-type)')))
  
  gAllelesByPop_barchart <- ggplot(plot4DF,aes(x=N*100,y=pop,fill=simpleAllele)) + 
    geom_bar(color='black',
             position='fill',stat='identity',
             lwd=.1) + 
    scale_fill_manual('',values=c('#009E73','#CC79A7','#0072B2','#56B4E9','#D55E00','#E69F00'),
                      guide=guide_legend(reverse = TRUE)) + 
    theme(legend.position='top',
          legend.direction='horizontal',
          legend.key.size = unit(0.1,'cm'),
          legend.key.height = unit(0.2,'cm'),
          legend.key.width = unit(0.2,'cm'),
          legend.box.margin = margin(-0.2,0,-0.4,0,'cm'),
          legend.background=element_blank(),
          legend.text = element_text(size=7)) + 
    ylab('Population') + 
    xlab('PRDM9 alleles (%)') + 
    scale_x_continuous(breaks=c(0,.5,1),labels=c(0,50,100))

  dfPopAlleles <- rbind(plot4DF %>% mutate(plottype='ABC'   ,colorBy=simpleAllele),
                        plot4DF %>% mutate(plottype='ACtype',colorBy=paste0(ACtype,'-type'))) %>%
    group_by(plottype,pop,colorBy) %>%
    summarize(pc=sum(N)/mean(tot)*100,N=sum(N),tot=mean(tot))
  dfPopAlleles$colorBy <- factor(dfPopAlleles$colorBy,levels=rev(c('A','B','C','L14','Other','A-type','C-type')))
  
  gAllelesByPop_barchart_AC <- ggplot(dfPopAlleles,aes(x=pc,y=pop,fill=colorBy)) + 
    geom_bar(color='black',
             position='fill',
             stat='identity',
             lwd=.2) + 
    scale_fill_manual('',values=c('A'      = 'yellow',
                                  'B'      = 'darkolivegreen2',
                                  'C'      = 'firebrick',
                                  'L14'    = 'violet',
                                  'Other'  = 'grey70',
                                  'A-type' = 'gold',
                                  'C-type' = 'red'),
                      guide=guide_legend(reverse = TRUE)) + 
    theme(legend.position='top',
          legend.direction='horizontal',
          legend.key.size = unit(0.1,'cm'),
          legend.key.height = unit(0.3,'cm'),
          legend.key.width = unit(1,'cm'),
          legend.text = element_text(size=7)) + 
    ylab('Population') + 
    xlab('PRDM9 alleles (%)') + 
    scale_x_continuous(breaks=c(0,.5,1),labels=c(0,50,100)) + 
    facet_wrap(~plottype)
  
  ## Make pies merge figure
  gPieLast <- ggarrange(lstPies[[popOrder[7]]] + ggtitle(popOrder[7]),
                        as_ggplot(get_legend(gAllelesByPop_barchart)),
                        ncol=1,nrow=2,heights=c(1,2))
  
  gPie1to6 <- ggarrange(lstPies[[popOrder[1]]] + ggtitle(popOrder[1]),
                        lstPies[[popOrder[2]]] + ggtitle(popOrder[2]),
                        lstPies[[popOrder[3]]] + ggtitle(popOrder[3]),
                        lstPies[[popOrder[4]]] + ggtitle(popOrder[4]),
                        lstPies[[popOrder[5]]] + ggtitle(popOrder[5]),
                        lstPies[[popOrder[6]]] + ggtitle(popOrder[6]),
                        nrow=3,ncol=2)
  
  gPies <- ggarrange(gPie1to6,gPieLast,ncol=2,nrow=1,
                     widths=c(2,1))
  ## Lollipops !! 
  dfA1  <- plotDF[,c('allele_1','pop')]; names(dfA1) <- c('allele','pop')
  dfA2  <- plotDF[,c('allele_2','pop')]; names(dfA2) <- c('allele','pop')
  
  dfM <- rbind(
    prdm9_DF %>% mutate(allele=allele_1,code=allele_1_code) %>% select(allele,code) %>% distinct() %>% dplyr:::filter(!grepl("(UU|NA)",code,perl=TRUE)),
    prdm9_DF %>% mutate(allele=allele_2,code=allele_2_code) %>% select(allele,code) %>% distinct() %>% dplyr:::filter(!grepl("(UU|NA)",code,perl=TRUE))) %>%
    distinct() %>%
    mutate(newZF = ifelse(grepl("!",code),TRUE,FALSE))
  
  dfPlotLollipop <- rbind(dfA1,dfA2)  %>% 
    group_by(allele) %>% count() %>%
    mutate(type='Pop',col='Pop') %>% 
    inner_join(dfM, by='allele') %>% 
    mutate(type = ifelse(grepl(':',allele) ,'Blood/Sperm',type),
           type = ifelse(grepl('^M',allele),'Novel: Pub ZFs',type),
           type = ifelse(grepl('!',code),'Novel: New ZF',type)) %>%
    mutate(col = ifelse(grepl(':',allele)                 ,'Sperm',col),
           col = ifelse(grepl(':00(0[0-9]|1[1-8])',allele),'Blood',col),
           col = ifelse(grepl('^M',allele)                ,'Novel: Pub ZFs',col),
           col = ifelse(grepl('!',code)                   ,'Novel: New ZF',col)) 
  
  dfPlotLollipop$allele <- factor(dfPlotLollipop$allele,levels=rev(alleleOrder))
  dfPlotLollipop$type   <- factor(dfPlotLollipop$type,levels=c('Pop','Blood/Sperm','Novel: Pub ZFs','Novel: New ZF'))
  dfPlotLollipop$col    <- factor(dfPlotLollipop$col,levels=c('Pop','Blood','Sperm','Novel: Pub ZFs','Novel: New ZF'))
  
  gLollipop <- ggplot(dfPlotLollipop) + 
    geom_vline(xintercept=1,lty='dashed',lwd=.2) +
    geom_segment(aes(x=0.5,xend=n,y=allele,yend=allele),lwd=.2) + 
    geom_point(aes(x=n,y=allele,fill=col),shape=21,size=1.6,lwd=.2) +
    scale_x_log10() + 
    scale_fill_manual(values=c('salmon','yellow','dodgerblue2','magenta','grey50')) +
    coord_cartesian(xlim=c(.8,960),expand=TRUE) + 
    geom_text(data=dfPlotLollipop[dfPlotLollipop$n == 1,],
              aes(x=n+1,y=allele),size=9*5/14,label='*') + 
    facet_wrap(~type,ncol=4,scales='free_y') + 
    theme(strip.text=element_text(size=7,hjust=0,face='bold'),
          strip.background=element_blank(),
          legend.position = 'none',
          plot.margin = unit(c(0,0.2,0,0),'cm')) + 
    annotation_logticks(sides = 'b',
                        size  = 0.2,
                        long  = unit(0.1,'cm'),
                        mid   = unit(0.05,'cm'),
                        short = unit(0.05,'cm')) +
    xlab('Individuals with allele (N)') +
    ylab('PRDM9 allele')
  
  ## COUNTS OF RARE ALLELES
  dfRareALL <- rbind(dfA1,dfA2)  %>% 
    mutate(is_rare = ifelse(allele %in% c('A','B','C','L14'),'common','rare')) %>% 
    group_by(is_rare) %>% count() %>% mutate(pop="ALL") %>% 
    reshape2:::dcast(pop~is_rare,value.var='n')
  
  dfRarePOP <- rbind(dfA1,dfA2)  %>% 
    mutate(is_rare = ifelse(allele %in% c('A','B','C','L14'),'common','rare')) %>% 
    group_by(is_rare,pop) %>% count() %>%
    reshape2:::dcast(pop~is_rare)
  
  dfRareByPop <- rbind(dfRareALL,dfRarePOP) %>%
    mutate(pc=rare/(common+rare)*100)
  
  fwrite(x=dfRareByPop, file = 'Alleva_et_al_Rare_and_Common_allele_Freqs_by_pop.txt',quote = FALSE,row.names = FALSE)
  
  ## Draw all figs if you want #################################################
  if (drawFigs){
    ggsave(paste0('PRDM9_ZFallele_Freq_Percent.png'),gAlleleFreqPC,width=14,height=4)
    ggsave(paste0('PRDM9_ZFallele_Freq_Percent.pdf'),gAlleleFreqPC,width=14,height=4)
    
    ggsave(paste0('PRDM9_ZFallele_Freq_Counts.png'),gAlleleFreqCnt,width=14,height=4)
    ggsave(paste0('PRDM9_ZFallele_Freq_Counts.pdf'),gAlleleFreqCnt,width=14,height=4)
    
    ggsave(paste0('PRDM9_ZFallele_Freq_Counts_and_PC.png'),gAlleleFreqCntPC,width=14,height=4)
    ggsave(paste0('PRDM9_ZFallele_Freq_Counts_and_PC.pdf'),gAlleleFreqCntPC,width=14,height=4)
    
    ggsave(paste0('PRDM9_ZFalleles_by_Pop_OneCopy.png'),gOneCopy,width=14,height=4)
    ggsave(paste0('PRDM9_ZFalleles_by_Pop_OneCopy.pdf'),gOneCopy,width=14,height=4)
    
    ggsave(paste0('PRDM9_ZFalleles_by_Pop_Freqs_Count.png'),gOneCopyCnt,width=14,height=4)
    ggsave(paste0('PRDM9_ZFalleles_by_Pop_Freqs_Count.pdf'),gOneCopyCnt,width=14,height=4)
  }
  
  return(list(gFreqPC          = gAlleleFreqPC,
              gFreqCnt         = gAlleleFreqCnt,
              gOneCopyPC       = gOneCopy,
              gOneCopyN        = gOneCopyCnt,
              gFreqPCNoWrap    = gAlleleFreqPCNoWrap,
              gFreqCntNoWrap   = gAlleleFreqCntNoWrap,
              gOneCopyPCNoWrap = gOneCopyNoWrap,
              gOneCopyNNoWrap  = gOneCopyCntNoWrap,
              gPies            = gPies,
              allPies          = lstPies,
              gLollipop        = gLollipop,
              gBars            = gAllelesByPop_barchart,
              supp             = lstSupp))
}

drawSamplingErrorFig <- function(dfBegin,nIterations = 10000){
  
  for (currentPop in popOrder[popOrder != 'YRI(T)']){
    print (paste0("Checking : ",currentPop))
    
    dfRealPopCounts <- dfBegin[dfBegin$pop == currentPop,] %>% 
      mutate(popTot=tot,type='real',randNum=0) %>%
      select(pop,allele,n,popTot,type,randNum) 
    
    #nPopSize <- dfRealPopCounts$popTot[1]
    nPopSize <- sum(dfRealPopCounts$n)
    
    dfPlotFreqVRand <- dfRealPopCounts
    
    ## DF to join so that zeros are counted
    dfPopAllZeros <- data.frame(allele=unique(dfRealPopCounts$allele),
                                pop=currentPop,
                                popTot=nPopSize,
                                nZero = 0)
    
    for (i in 1:nIterations){
      ## Randomly choose a set of alleles
      vRandomAlleles <- sample(x       = dfRealPopCounts$allele,
                               size    = nPopSize,
                               replace = TRUE,
                               prob    = dfRealPopCounts$n)
      
      ## Make DF from randomly chosen alleles
      dfTmp <- data.frame(allele=as.character(vRandomAlleles)) %>% 
        group_by(allele) %>% 
        count() %>% 
        mutate(pop=currentPop,popTot=nPopSize,type='random',randNum=i) %>% 
        select(pop,popTot,allele,n,type,randNum)
      
      ## Add alleles with ZERO counts
      dfRandomAlleles <- dfPopAllZeros %>% join(dfTmp,by=c('allele','pop','popTot'))
      dfRandomAlleles$type <- 'random'
      dfRandomAlleles$randNum <- i
      dfRandomAlleles$n[is.na(dfRandomAlleles$n)] <- 0
      
      dfRandomAlleles <- dfRandomAlleles %>% select(pop,popTot,allele,n,type,randNum)
      
      ## Add this round of sampling to main dataframe
      dfPlotFreqVRand <- rbind(dfPlotFreqVRand,dfRandomAlleles) %>% 
        select(pop,allele,n,popTot,type,randNum)
    }
    
    ## Get summary statistics for random iterations
    dfError <- dfPlotFreqVRand[dfPlotFreqVRand$type=='random',] %>% 
      group_by(pop,allele) %>% 
      summarize(p1=quantile(n,.01),
                p5=quantile(n,.05),
                p95=quantile(n,.95),
                p99=quantile(n,.99),
                mean=mean(n),
                .groups='drop_last')
    
    ## Add this pop to main DF
    if (currentPop == popOrder[1]){
      dfAllPops <- dfRealPopCounts %>% inner_join(dfError,by=c('allele','pop'))
    }else{
      dfAllPops <- rbind(dfAllPops,dfRealPopCounts %>% inner_join(dfError,by=c('allele','pop')))
    }
  }
  
  dfAllPops$pop <- factor(dfAllPops$pop,levels=popOrder)
  
  ### Select Alleles
  ## Select only alleles where 1st %ile > 0
  selectAlleles        <- dfAllPops %>% group_by(allele) %>% 
    summarize(maxP1=max(p1)) %>% 
    filter(maxP1>0) %>% select(allele)
  
  selectAlleles$allele <- factor(selectAlleles$allele,levels=alleleOrder)
  
  numAlleles <- length(unique(selectAlleles$allele))
  
  dfAllCombos <- data.frame(allele=rep(unique(selectAlleles$allele),
                                       each = length(popOrder)),
                            pop=rep(popOrder,numAlleles),
                            nAlt=0)
  
  dfAllPopPlot <- dfAllCombos %>% 
    inner_join(dfAllPops,by=c('pop','allele')) %>% 
    mutate(n=ifelse(n>0,n,0))
  
  dfAllPopPlot$allele <- factor(dfAllPopPlot$allele,levels=alleleOrder)
  
  dfToPlot <- dfAllPopPlot[dfAllPopPlot$allele %in% selectAlleles$allele,]
  dfToPlot$pop <- factor(dfToPlot$pop,levels=popOrder)
  
  g <- ggplot(data=dfToPlot) +  
    geom_vline(xintercept=0,lwd=.2) + 
    geom_bar(aes(x=n/popTot*100,y=pop,fill=pop),stat='identity',position=position_dodge(),alpha=.1) +
    geom_errorbar(aes(fill=pop,y=pop,xmin=p1/popTot*100,xmax=p99/popTot*100),width=.5,lwd=.2,position=position_dodge())  + 
    geom_errorbar(aes(color=pop,fill=pop,y=pop,xmin=mean/popTot*100,xmax=mean/popTot*100),width=.8,lwd=.7,position=position_dodge())  + 
    facet_wrap(~allele,ncol=4,scales='free_x') + 
    scale_fill_brewer(palette='Dark2',guide = guide_legend(reverse = TRUE) ) +
    scale_color_brewer(palette='Dark2',guide = guide_legend(reverse = TRUE) ) +
    theme(legend.position=c(1,0),
          legend.justification=c(1,0),
          legend.title=element_blank(),
          legend.key.size=unit(c(0.3),'cm'),
          strip.text = element_text(size=7),
          legend.text = element_text(size=7),
          axis.line.y=element_blank(),axis.ticks.y=element_blank()) + 
    xlab('Allele frequency (%)') + 
    ylab('') + 
    coord_cartesian(ylim=c(.3,7.7),expand=FALSE,clip=FALSE)
  
  return(list(fig=g,
              data=dfAllPopPlot,
              plotData=dfToPlot))
}

## Now draw the figures ########################################################
lstSupp       <- drawGenotypingSuccessFigs()
lstDiploidGTs <- drawDiploidGenotypes()
lstAlleleFreq <- drawAlleleFrequencies()

## Supp fig first ##############################################################
ggsave(paste0('Alleva_et_al_SuppFig_Genotyping_By_Pop_VS_Reads.png'),lstSupp$gSupp,width=7,height=3)
ggsave(paste0('Alleva_et_al_SuppFig_Genotyping_By_Pop_VS_Reads.pdf'),lstSupp$gSupp,width=7,height=3)

## Now main Fig ################################################################
gTop <- ggarrange(lstAlleleFreq$gBars, lstDiploidGTs$gHomHet,
                  ncol=2,nrow=1,
                  widths=c(4,2.4),
                  labels=c('A','B'),
                  font.label=c('size'=9))

ggarrange(gTop, lstDiploidGTs$gFig, 
          ncol=1, nrow=2,
          heights=c(2,5),
          labels=c('','C'),
          font.label=c('size'=9))  

gBC <- ggarrange(lstAlleleFreq$gBars, lstDiploidGTs$gHomHet,
                  ncol=1,nrow=2,
                  heights=c(6,5),
                  labels=c('B','C'),
                  font.label=c('size'=8))

gABC <- ggarrange(lstAlleleFreq$gLollipop, gBC,
                 ncol=2,nrow=1,
                 widths=c(2,1),
                 labels=c('A',''),
                 font.label=c('size'=8))

gABCD <- ggarrange(gABC, lstDiploidGTs$gFig,
                  ncol=1,nrow=2,
                  heights=c(2,3),
                  labels=c('','D'),
                  font.label=c('size'=8))

ggsave('Alleva_et_al_Figure2.png',plot = gABCD, height=7,width=6.5,dpi=500)
ggsave('Alleva_et_al_Figure2.pdf',plot = gABCD, height=7,width=7)

ggsave('Alleva_et_al_SuppFig_PopPies.png',lstAlleleFreq$gPies,height=5,width=5,dpi=500)
ggsave('Alleva_et_al_SuppFig_PopPies.pdf',lstAlleleFreq$gPies,height=5,width=5)

ggsave('Alleva_et_al_SuppFig_PrZFA_alleleFreq_N_and_PC.png',lstAlleleFreq$gFreqBoth,height=10,width=7,dpi=500)
ggsave('Alleva_et_al_SuppFig_PrZFA_alleleFreq_N_and_PC.pdf',lstAlleleFreq$gFreqBoth,height=10,width=7)
  
ggsave('Alleva_et_al_SuppFig_PrZFA_alleleCarriers_PC.png',lstAlleleFreq$gOneCopyPC,height=10,width=7,dpi=500)
ggsave('Alleva_et_al_SuppFig_PrZFA_alleleCarriers_PC.pdf',lstAlleleFreq$gOneCopyPC,height=10,width=7)

ggsave('Alleva_et_al_SuppFig_alleleFrequency_with_CIs_from_bootstrapped_sampling.png',plot=lstAlleleFreq$supp$fig,dpi=500,width=5,height=6)
ggsave('Alleva_et_al_SuppFig_alleleFrequency_with_CIs_from_bootstrapped_sampling.pdf',plot=lstAlleleFreq$supp$fig,width=5,height=6)

gAlleles <- ggarrange(lstAlleleFreq$gFreqPCNoWrap    + ggtitle('Frequency (%)') + ylab('') + theme(plot.margin=unit(c(0,0,0,0),'cm')),
                      lstAlleleFreq$gFreqCntNoWrap   + ggtitle('(Count)')       + ylab('') + theme(axis.text.y = element_blank(),plot.margin=unit(c(0,0,0,0),'cm')),
                      lstAlleleFreq$gOneCopyPCNoWrap + ggtitle('Carriers (%)')  + ylab('') + theme(axis.text.y = element_blank(),plot.margin=unit(c(0,0,0,0),'cm')),
                      lstAlleleFreq$gOneCopyNNoWrap  + ggtitle('(Count)')       + ylab('') + theme(axis.text.y = element_blank(),plot.margin=unit(c(0,0.2,0,0),'cm')),
                      widths=c(3.3, 3, 2.3, 3),
                      ncol=4,
                      labels=c('A','','B',''),
                      font.label = list(size=8,face='bold',hjust=0,vjust=1))

ggarrange(gAlleles,
          lstDiploidGTs$gFigN + xlab('') + theme(plot.margin=unit(c(0,0,0,0),'cm')),
          ncol=1,nrow=2,
          heights=c(1.8,1),
          labels=c('','C'),
          font.label = list(size=8,face='bold',hjust=0,vjust=1))

ggsave('Alleva_et_al_SuppFig_PrZFA_diploinFreq_alleleFreq_and_carrierFreq.png',bg = 'white',height=10,width=6.5,dpi=500)
ggsave('Alleva_et_al_SuppFig_PrZFA_diploinFreq_alleleFreq_and_carrierFreq.pdf',bg = 'white',height=10,width=6.5)

