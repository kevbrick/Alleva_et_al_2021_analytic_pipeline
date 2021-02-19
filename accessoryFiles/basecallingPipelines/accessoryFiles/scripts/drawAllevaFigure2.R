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

## SET ORDER
popOrder <- c('LWK','YRI','CHB','PJL','PEL','TSI','FIN')

### Pops DF
popsDF <- data.frame(pop = c('CHB','FIN','LWK','PEL','PJL','TSI','YRI'),
                     tot = c(120,100,120,70,108,114,120))

popsDF$pop <- factor(popsDF$pop,popOrder)

## Make parent DF ##############################################################
plotDF <- prdm9_DF %>% filter(ok & pop != 'OTH') %>%
  select(pop,diploid,allele_1,allele_2,type,typeHH) 

plotDF$pop <- factor(plotDF$pop,popOrder)

## Figure drawing functions ####################################################
drawGenotypingSuccessFigs <- function(popOrder = c('LWK','YRI','CHB','PJL','PEL','TSI','FIN'), drawFigs = FALSE){
  
  ## Fig 1 #####################################################################
  plot1DF <- melt(plotDF[,c('diploid','pop')] %>% group_by(pop) %>% count() %>% 
                    inner_join(popsDF) %>% mutate(yes=n,no=tot-n) %>% select(pop,yes,no))
  
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
  
  ## Fig 3 #####################################################################
  gtsDF$yesno <- 'no'
  gtsDF$yesno[gtsDF$ok] <- 'yes'
  gtsDF$yesno <- factor(gtsDF$yesno,levels=c('yes','no'))
  
  gLowCoverage <- ggplot(gtsDF,aes(x=pop,y=nseqs,fill=yesno)) + 
    geom_boxplot(position=position_dodge(),alpha=.2,lwd=.2,outlier.alpha=0) + 
    geom_point(color=alpha('white',0),size = 0.7, shape = 21, position = position_jitterdodge(),lwd=.2) +
    scale_y_log10()+
    scale_fill_manual('PRDM9 diploid genotype inferred ?',
                      values=c('grey40','darkorange')) +
    xlab('') + 
    ylab('Raw reads (#)') + 
    theme(legend.position='top')
  
  ## Supp Fig 1 ################################################################
  gSuppGT <- ggarrange(gGTCounts,gLowCoverage,ncol=2,nrow=1,labels=c('A','B'),font.label = c(size=9))
  
  if (drawFigs){
    ggsave(paste0('PRDM9_diploid_GT_Counts.png'),gGTCounts,width=4,height=3)
    ggsave(paste0('PRDM9_diploid_GT_Counts.pdf'),gGTCounts,width=4,height=3)
    
    ggsave(paste0('PRDM9_ZFs_per_person.png'),gZFvsStatus,width=4,height=3)
    ggsave(paste0('PRDM9_ZFs_per_person.pdf'),gZFvsStatus,width=4,height=3)

    ggsave(paste0('PRDM9_raw_reads_per_person.png'),gLowCoverage,width=4,height=3)
    ggsave(paste0('PRDM9_raw_reads_per_person.pdf'),gLowCoverage,width=4,height=3)    
    
    ggsave(paste0('Alleva_et_al_Genotyping_By_Pop_VS_Reads.png'),gSuppGT,width=7,height=3)
    ggsave(paste0('Alleva_et_al_Genotyping_By_Pop_VS_Reads.pdf'),gSuppGT,width=7,height=3)
  }
  
  return(list(gSupp = gSuppGT, 
              gFig1 = gGTCounts,
              gFig2 = gZFvsStatus,
              gFig3 = gLowCoverage))
}

drawDiploidGenotypes <- function(drawFigs=FALSE){
  
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
  
  ## Make Long-form figs #######################################################
  gGenotypePC <- ggplot(plot2DF,aes(x=pop,y=diploid,fill=pcRnd)) + 
    geom_tile(lwd=.05,color='black') + 
    geom_text(aes(label=round(pc,0)),
              size=7*5/14) +
    scale_fill_gradient('Freq(%)',high = 'firebrick1', low='white', na.value = 'white') + 
    xlab('Population') + ylab('Diploid PRDM9 genotype') +
    ggtitle('Genotype frequency (%)') +
    facet_wrap(~grp,ncol=3,scales='free_y') + 
    theme(strip.text=element_blank(),
          strip.background=element_blank(),
          legend.position = 'none',
          plot.title = element_text(size=7,hjust=0),
          axis.text.x=element_text(angle=90,hjust=1))
  
  gGenotypeN <- ggplot(plot2DF,aes(x=pop,y=diploid,fill=n)) + 
    geom_tile(lwd=.05,color='black') + 
    geom_text(aes(label=round(n,0)),
              size=7*5/14) +
    scale_fill_gradient('Freq(%)',high = 'dodgerblue2', low='white', na.value = 'white') + 
    xlab('Population') + ylab('Diploid PRDM9 genotype') +
    ggtitle('Genotype frequency (Count)') +
    facet_wrap(~grp,ncol=3,scales='free_y') + 
    theme(strip.text=element_blank(),
          strip.background=element_blank(),
          legend.position = 'none',
          axis.text.x=element_text(angle=90,hjust=1))
  
  ## Make HomHet fig ###########################################################
  # plot7DF <- plotDF %>%  select(pop,type) %>% group_by(pop,type) %>% count() %>%  
  #   inner_join(plotDF %>% group_by(pop) %>% count(name = "tot")) %>% 
  #   summarise(N=tot, pc=n/tot*100)
  
  plot7DF <- plotDF %>%  select(pop,type) %>% mutate(type=ifelse(type=='Hom','Hom','Het')) %>%
    group_by(pop,type) %>% count() %>%  
    inner_join(plotDF %>% group_by(pop) %>% count(name = "tot")) %>% 
    summarise(N=tot, pc=n/tot*100)
  
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
  
  return(list(gFig    = gGenotypePC,
              gFigN   = gGenotypeN,
              gHomHet = gHomHet))
}

drawSamplingErrorFig <- function(dfBegin,nIterations = 1000){
  
  for (currentPop in popOrder){
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
  
  g <- ggplot(data=dfAllPopPlot[dfAllPopPlot$allele %in% selectAlleles$allele,]) +  
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
              data=dfAllPopPlot))
}

drawAlleleFrequencies <- function(drawFigs=FALSE){
  
  ## Make DF ###################################################################
  plot4DF <- melt(plotDF[,c('pop','allele_1','allele_2')],
                  id.vars='pop',value.name='allele') %>%
    group_by(pop) %>% 
    add_count(name='tot') %>% 
    group_by(pop,allele,tot) %>% 
    count() %>% 
    group_by(pop,allele,tot) %>% 
    left_join(plotDF %>% count(pop, name = "totPop")) %>%
    summarize(N=n,
              pc=n/tot*100,
              pcRnd=round(n/tot*100/10)*10,
              popPC=n/totPop*100)
  
  plot4DF$allele <- factor(plot4DF$allele,levels=rev(alleleOrder))
  plot4DF$pop    <- factor(plot4DF$pop,popOrder)
  
  dfAlleles      <- data.frame(allele=rev(levels(plot4DF$allele)),grp=0)
  dfAlleles$grp  <- rep(1:3,each=ceiling(length(dfAlleles$grp)/3))[1:length(dfAlleles$allele)]
  
  p              <- join(plot4DF,dfAlleles,by='allele')
  plot4DF        <- p
  
  ## Fig 1 : Allele freqs (%) ##################################################
  gAlleleFreqPC <- ggplot(plot4DF,aes(x=pop,y=allele,fill=pc)) + 
    geom_tile(lwd=.05,color='black') + 
    geom_text(aes(label=round(pc,0)),
              size=7*5/14) +
    scale_fill_gradient('%',high = 'firebrick1', 
                        low='white', 
                        na.value = 'white') + 
    xlab('Population') + ylab('PRDM9 allele') + 
    ggtitle('Allele frequency (%)') + 
    facet_wrap(~grp,ncol=3,scales='free_y') + 
    theme(strip.text=element_blank(),
          strip.background=element_blank(),
          legend.position = 'none',
          axis.text.x=element_text(angle=90,hjust=1))
  
  gAlleleFreqCnt <- ggplot(plot4DF,aes(x=pop,y=allele,fill=N)) + 
    geom_tile(lwd=.05,color='black') + 
    geom_text(aes(label=round(N,0)),
              size=7*5/14) +
    scale_fill_gradient('%',high = 'dodgerblue2', 
                        low='white', 
                        na.value = 'white') + 
    xlab('Population') + ylab('PRDM9 allele') + 
    ggtitle('Allele frequency (Count)') + 
    facet_wrap(~grp,ncol=3,scales='free_y') + 
    theme(strip.text=element_blank(),
          strip.background=element_blank(),
          legend.position = 'none',
          axis.text.x=element_text(angle=90,hjust=1))
  
  ## Fig 2 : Allele counts #####################################################
  gAlleleFreqCnt <- ggplot(plot4DF,aes(x=pop,y=allele,fill=N)) + 
    geom_tile(lwd=.05,color='black') + 
    geom_text(aes(label=round(N,0)),
              size=7*5/14) +
    scale_fill_gradient('%',high = 'dodgerblue2', 
                        low='white', 
                        na.value = 'white') + 
    xlab('Population') + ylab('PRDM9 allele') + 
    ggtitle('Allele frequency (Count)') + 
    facet_wrap(~grp,ncol=3,scales='free_y') + 
    theme(strip.text=element_blank(),
          strip.background=element_blank(),
          legend.position = 'none',
          axis.text.x=element_text(angle=90,hjust=1))

  ## Supp Fig. Pop specific alleles
  dfPopSpecific <- as.data.frame(plot4DF %>% mutate(n=N) %>% select(pop,tot,allele,n))
  
  lstSupp <- drawSamplingErrorFig(dfPopSpecific,nIterations = 1000)
  
  ggsave('Alleva_et_al_SuppFig_alleleFrequency_with_CIs_from_bootstrapped_sampling.png',plot=lstSupp$fig,dpi=500,width=5,height=6)
  ggsave('Alleva_et_al_SuppFig_alleleFrequency_with_CIs_from_bootstrapped_sampling.pdf',plot=lstSupp$fig,width=5,height=6)
  
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
  
  gOneCopy <- ggplot(plot6DF,aes(x=pop,y=allele,fill=popPC)) + 
    geom_tile(lwd=.05,color='black') + 
    geom_text(aes(label=round(popPC,0)),
              size=7*5/14,check_overlap = TRUE) +
    scale_fill_gradient(high = 'firebrick1', low='white', na.value = 'white') + 
    xlab('Population') + ylab('PRDM9 allele') + 
    ggtitle('Carriers (at least one allele) (%)') + 
    facet_wrap(~grp,ncol=3,scales='free_y') + 
    theme(strip.text=element_blank(),
          strip.background=element_blank(),
          legend.position = 'none',
          axis.text.x=element_text(angle=90,hjust=1))
  
  gOneCopyCnt <- ggplot(plot6DF,aes(x=pop,y=allele,fill=N)) + 
    geom_tile(lwd=.05,color='black') + 
    geom_text(aes(label=round(N,0)),
              size=7*5/14,check_overlap = TRUE) +
    scale_fill_gradient(high = 'dodgerblue2', low='white', na.value = 'white') + 
    xlab('Population') + ylab('PRDM9 allele') + 
    ggtitle('Carriers (at least one allele) (Count)') + 
    facet_wrap(~grp,ncol=3,scales='free_y') + 
    theme(strip.text=element_blank(),
          strip.background=element_blank(),
          legend.position = 'none',
          axis.text.x=element_text(angle=90,hjust=1))
  
  ## Pie Charts ################################################################
  n <- length(alleleOrder)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  col_DF <- data.frame(allele=alleleOrder,col=col_vector[1:n])
  
  lstPies <- list()
  
  for (pop in popOrder){
    
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
  plot4DF$simpleAllele <- 'Other'
  plot4DF$simpleAllele[plot4DF$allele %in% c('A','B','C')] <- as.character(plot4DF$allele[plot4DF$allele %in% c('A','B','C')])
  
  plot4DF$simpleAllele <- factor(plot4DF$simpleAllele,levels=rev(c('A','B','C','Other')))
  
  gAllelesByPop_barchart <- ggplot(plot4DF,aes(x=N*100,y=pop,fill=simpleAllele)) + 
    geom_bar(color='grey50',
             position='fill',stat='identity',
             lwd=.1) + 
    scale_fill_manual('',values=c('yellow','darkolivegreen2','violet','grey70'),
                      guide=guide_legend(reverse = TRUE)) + 
    theme(legend.position='top',
          legend.direction='horizontal',
          legend.key.size = unit(0.1,'cm'),
          legend.key.height = unit(0.3,'cm'),
          legend.text = element_text(size=7)) + 
    ylab('Population') + 
    xlab('PRDM9 alleles (%)') + 
    scale_x_continuous(breaks=c(0,.5,1),labels=c(0,50,100))
  
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

  dfPlotLollipop <- rbind(dfA1,dfA2)  %>% 
    group_by(allele) %>% count() %>%
    mutate(type='Pop',col='Pop')
  
  dfPlotLollipop$allele <- factor(dfPlotLollipop$allele,levels=rev(alleleOrder))
  dfPlotLollipop$type[grepl(':',dfPlotLollipop$allele)] <- 'Blood/Sperm'
  dfPlotLollipop$type[grepl('^M',dfPlotLollipop$allele)] <- 'Novel: New ZF'
  dfPlotLollipop$type[dfPlotLollipop$allele %in% paste0('M',c(2,4,6,13:16,18:20,22:27,29,31))] <- 'Novel: Pub ZFs'
  
  dfPlotLollipop$col[grepl(':',dfPlotLollipop$allele)] <- 'Sperm'
  dfPlotLollipop$col[grepl('000[29]',dfPlotLollipop$allele)] <- 'Blood'
  dfPlotLollipop$col[grepl('^M',dfPlotLollipop$allele)] <- 'Novel: New ZF'
  dfPlotLollipop$col[dfPlotLollipop$allele %in% paste0('M',c(2,4,6,13:16,18:20,22:27,29,31))] <- 'Novel: Pub ZFs'
  
  dfPlotLollipop$type <- factor(dfPlotLollipop$type,levels=c('Pop','Blood/Sperm','Novel: Pub ZFs','Novel: New ZF'))
  dfPlotLollipop$col  <- factor(dfPlotLollipop$col,levels=c('Pop','Blood','Sperm','Novel: Pub ZFs','Novel: New ZF'))
  
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
  
  return(list(gFreqPC    = gAlleleFreqPC,
              gFreqCnt   = gAlleleFreqCnt,
              gOneCopyPC = gOneCopy,
              gOneCopyN  = gOneCopyCnt,
              gPies      = gPies,
              allPies    = lstPies,
              gLollipop  = gLollipop,
              gBars      = gAllelesByPop_barchart))
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