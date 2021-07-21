source('accessoryFiles/scripts/genericFunctions.R')

library(data.table)
library(ggplot2)
library(plyr)
library(tidyverse)

theme7point()

## Fills out sparse data frame with 0s where necessary
## built in R functions don't do this quite right !!
fillOut <- function(x){
  for (i in unique(x$phenotype)){
    snpID <- x$ID[x$phenotype == i][1]
    for (j in unique(x$SNP[x$phenotype == i])){
      for (k in unique(x$prdm9GT)){
        xx <- x %>% dplyr:::filter(phenotype == i,
                                   SNP==j,
                                   prdm9GT == k)
        if (dim(xx)[1] == 0){
          xAdd <- data.frame(ID=snpID,
                             SNP=j,
                             prdm9GT=k,
                             tot=99999,
                             phenotypeID=i,
                             cnt=0,
                             phenotype=i)
          x <- rbind(x,xAdd)
        }
        #print(paste0(i,":",j,":",k, "-->",dim(xx)[1]))
        
      } 
    } 
  }
  return(x)
}

fillOutPop <- function(x){
  pops <- unique(x$pop)
  for (i in unique(x$phenotype)){
    snpID <- x$ID[x$phenotype == i][1]
    GT    <- unique(x$prdm9GT[x$phenotype == i])
    SNPs  <- unique(x$SNP[x$phenotype == i])
    
    for (j in SNPs){
      for (k in GT){
        for (p in pops){
          xx <- x %>% dplyr:::filter(phenotype == i,
                                     SNP==j,
                                     prdm9GT == k,
                                     pop==p)
          if (dim(xx)[1] == 0){
            xAdd <- data.frame(ID=snpID,
                               SNP=j,
                               prdm9GT=k,
                               tot=99999,
                               pop=p,
                               phenotypeID=i,
                               cnt=0,
                               phenotype=i)
            x <- rbind(x,xAdd)
          }
        }
      }
    }
  }
   
  return(x)
}

getSNPpenetrance <- function(phenotypeID = 'C',snpid=NULL){
  
  dfCMH <- read.table(paste0('prdm9AS.',phenotypeID,'.mod.cmh'),header=TRUE)
  
  dfCMH$col="Others"
  dfCMH$col[dfCMH$selected>0] <- 'Selected'
  
  ## Suck in phenotypes and individual genotypes again ... just a bit more convenient
  dfPhenotypes   <- read.table('phenotypes.txt',header=TRUE)
  dfGT           <- read.table('genotypes.txt',header=FALSE)
  names(dfGT)    <- c('FID','IID','diploid')
  dfClust        <- read.table('clusters.txt',header=FALSE)
  names(dfClust) <- c('FID','IID','pop')
  
  dfPRDM9 <- join(dfGT,dfClust,by="FID") %>% select(FID,pop,diploid)
  
  ## Get SNP IDs
  if (!is.null(snpid)){
    dfRSID     <- read.table(paste0('prdm9AS_keepers.select',phenotypeID,'.map'), header=FALSE) %>% 
      select(ID=V2) %>%
      inner_join(dfCMH, by=c("ID" = "SNP")) %>%
      dplyr:::filter(ID == snpid) %>%
      mutate(SNP=paste0("SNP",row_number()))
  }else{
    dfRSID     <- read.table(paste0('prdm9AS_keepers.select',phenotypeID,'.map'), header=FALSE) %>% 
      select(ID=V2) %>%
      inner_join(dfCMH, by=c("ID" = "SNP")) %>%
      dplyr:::filter(ID != 'rs77023486' & ID != 'rs141586808' & ID != 'rs138354146') %>%
      arrange(P) %>%
      mutate(SNP=paste0("SNP",row_number()))
    
    if (dfRSID$P[1] == dfRSID$P[dfRSID$ID == 'rs6889665']){
      dfRSID <- dfRSID[dfRSID$ID == 'rs6889665',]
    }else{
      dfRSID <- dfRSID[1,]
    }
  }
  
  SNPID <- dfRSID$ID[1]
  print(phenotypeID)
  print(SNPID)
  
  #dfGenotypes <- read.table(paste0('prdm9AS_keepers.select',phenotypeID,'.ped'),header=FALSE)
  dfGenotypes <- read.table(paste0('prdm9AS_keepers.select',phenotypeID,'.SNPsequences.tab'), header=TRUE)
  
  names(dfGenotypes)[1] <- 'FID'
  
  dfGenotypes  <- join(dfGenotypes, dfPhenotypes, by='FID')
  dfGenotypes  <- join(dfGenotypes, dfPRDM9     , by='FID')
  
  dfGenotypes  <- dfGenotypes[!is.na(dfGenotypes$pop),]
  
  if (phenotypeID %in% c('Alike','Clike')){
    dfGenotypes <- dfGenotypes %>% rowwise() %>% 
      mutate(diploid = ifelse(Alike == 2, ifelse(Clike == 2,'Alike/Clike','Alike/Alike'),'Clike/Clike'))
  }
  
  dfMGenotypes <- reshape2:::melt.data.frame(dfGenotypes %>% 
                                               dplyr:::rename(gt1=paste0(SNPID,"_A"), 
                                                      gt2=paste0(SNPID,"_B")) %>%
                                               select(FID,A,B,C,pop,diploid,gt1,gt2),
                                             id.vars=c("FID","A","B","C","pop",'diploid'),
                                             measure_vars=paste("gt",1:2)) 
  
  dfMGenotypes$ID <- 'NA'
  
  dfMGenotypes$ID <- dfRSID$ID[1]
  
  # dfMGenotypes <- reshape2:::melt.data.frame(dfGenotypes %>% 
  #                                              select(FID,A,B,C,pop,diploid,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16),
  #                                            id.vars=c("FID","A","B","C","pop",'diploid'),
  #                                            measure_vars=paste("V",7:8)) 
  # dfMGenotypes$SNP <- 'NA'
  # 
  # dfMGenotypes$SNP[dfMGenotypes$variable %in% paste0("V",7:8)]   <- 'SNP1'
  # #dfMGenotypes$SNP[dfMGenotypes$variable %in% paste0("V",9:10)]  <- 'SNP2'
  # #dfMGenotypes$SNP[dfMGenotypes$variable %in% paste0("V",11:12)] <- 'SNP3'
  # #dfMGenotypes$SNP[dfMGenotypes$variable %in% paste0("V",13:14)] <- 'SNP4'
  # #dfMGenotypes$SNP[dfMGenotypes$variable %in% paste0("V",15:16)] <- 'SNP5'
  
  dfMGenotypes <- join(dfMGenotypes,dfRSID,by="ID") %>% 
    select('FID', 'pop', 'diploid', SNP='value', 'ID')
  
  ### For coarse resolution
  dfMGenotypes$prdm9GTsimple <- NA
  dfMGenotypes$prdm9GTsimple[dfMGenotypes$diploid == paste0(phenotypeID,'/',phenotypeID)]  <- "HOM"
  dfMGenotypes$prdm9GTsimple[  grepl(paste0('^',phenotypeID,'|',phenotypeID,'$'),dfMGenotypes$diploid)  & is.na(dfMGenotypes$prdm9GTsimple)] <- "HET"
  dfMGenotypes$prdm9GTsimple[!(grepl(paste0('^',phenotypeID,'|',phenotypeID,'$'),dfMGenotypes$diploid)) & is.na(dfMGenotypes$prdm9GTsimple)] <- "NONE"
  dfMGenotypes$prdm9GTsimple <- factor(dfMGenotypes$prdm9GTsimple,levels=c('NONE','HET','HOM'))
  
  ### For detailed resolution
  dfMGenotypes$prdm9GTdetail <- NA
  dfMGenotypes$prdm9GTdetail[dfMGenotypes$diploid == paste0(phenotypeID,'/',phenotypeID)]  
  dfMGenotypes$prdm9GTdetail[dfMGenotypes$diploid == paste0(phenotypeID,'/',phenotypeID)]  <- paste0(phenotypeID,'/',phenotypeID)
  dfMGenotypes$prdm9GTdetail[dfMGenotypes$diploid == 'B/B']     <- 'B/B'
  dfMGenotypes$prdm9GTdetail[dfMGenotypes$diploid == 'A/A']     <- 'A/A'
  dfMGenotypes$prdm9GTdetail[dfMGenotypes$diploid == 'C/C']     <- 'C/C'
  dfMGenotypes$prdm9GTdetail[dfMGenotypes$diploid == 'L14/L14'] <- 'L14/L14'
  dfMGenotypes$prdm9GTdetail[dfMGenotypes$diploid == 'L4/L4']   <- 'L4/L4'
  dfMGenotypes$prdm9GTdetail[dfMGenotypes$diploid == 'A/L4']    <- 'A/L4'
  dfMGenotypes$prdm9GTdetail[dfMGenotypes$diploid == 'A/L14']   <- 'A/L14'
  dfMGenotypes$prdm9GTdetail[dfMGenotypes$diploid == 'A/L19']   <- 'A/L19'
  dfMGenotypes$prdm9GTdetail[  grepl(paste0('^',phenotypeID,'|',phenotypeID,'$'),dfMGenotypes$diploid)  & is.na(dfMGenotypes$prdm9GTdetail)] <- paste0(phenotypeID,'-carrier')
  dfMGenotypes$prdm9GTdetail[!(grepl(paste0('^',phenotypeID,'|',phenotypeID,'$'),dfMGenotypes$diploid)) & is.na(dfMGenotypes$prdm9GTdetail)] <- 'Other alleles'
  
  dfMGenotypes$phenotypeID <- phenotypeID
  
  dfMGenotypes <- dfMGenotypes[!is.na(dfMGenotypes$pop),]
  
  return(dfMGenotypes)
}

plotSingleSNPFig <- function(df2Plot,lbl='numpc'){
  g <- ggplot(df2Plot,aes(fill=cnt/tot*100,x=SNP,y=prdm9GT)) + 
    geom_tile(color='black',size=0.25) + 
    theme(panel.border = element_rect(size=.2,fill=NA),
          plot.title = element_text(size=7,hjust=0.5), 
          legend.position = 'top',
          legend.key.width = unit(1,'cm'),
          legend.key.height=unit(0.2,'cm'),
          legend.text=element_text(size=7)) + 
    ylab('PRDM9 genotype') + 
    xlab('SNP haplotype') + 
    scale_fill_gradient2('Allele frequency (%)',low='white',mid='yellow',high='red',midpoint=50) + 
    ggtitle (paste0("Phenotype: ",df2Plot$phenotypeID[1]))
  
  if ("pop" %in% names(df2Plot)){
    g <- g + facet_grid(pop~ID,scales='free_x') 
  }else{
    g <- g + facet_grid(~ID,scales='free_x') 
  }
  
  if (lbl == 'numpc'){
    g <- g + geom_text(aes(label=paste0(round(cnt/tot*100),"%\n(n = ",cnt,")")),size=7*5/14,check_overlap = TRUE)
  }
  
  if (lbl == 'num'){
    g <- g + geom_text(aes(label=cnt),size=7*5/14,check_overlap = TRUE)
  }
  
  if (lbl == 'pc'){
    g <- g + geom_text(aes(label=paste0(round(cnt/tot*100),"%")),size=7*5/14,check_overlap = TRUE)
  }
  
  return(g)
}

getPlotData <- function(dfGT){
  
  #-----------------------------------------------------------
  # Make DFs for plotting:
  
  ## Simple: All
  dfPlotSimple <- dfGT %>% 
    group_by(prdm9GTsimple,ID,phenotypeID) %>% 
    add_count(name='tot') %>% 
    group_by(ID,SNP,prdm9GT = prdm9GTsimple,tot,phenotypeID) %>% 
    count(name = 'cnt')
  
  ## Simple: By Pop 
  dfPlotSimpleByPop <-dfGT %>% 
    group_by(prdm9GTsimple,ID,pop,phenotypeID) %>% 
    add_count(name='tot') %>% 
    group_by(ID,SNP,prdm9GT = prdm9GTsimple,tot,pop,phenotypeID) %>% 
    count(name = 'cnt')
  
  ## Detail: All
  dfPlotDetail <- dfGT %>% 
    group_by(prdm9GTdetail,ID,phenotypeID) %>% 
    add_count(name='tot') %>% 
    group_by(ID,SNP,prdm9GT = prdm9GTdetail,tot,phenotypeID) %>% 
    count(name = 'cnt')
  
  ## Simple: By Pop 
  dfPlotDetailByPop <- dfGT %>% 
    group_by(prdm9GTdetail,ID,pop,phenotypeID) %>% 
    add_count(name='tot') %>% 
    group_by(ID,SNP,prdm9GT = prdm9GTdetail,tot,pop,phenotypeID) %>% 
    count(name = 'cnt')
  
  return(list(simple_all = dfPlotSimple,
              simple_pop = dfPlotSimpleByPop,
              detail_all = dfPlotDetail,
              detail_pop = dfPlotDetailByPop))
}

plinkplot <- function(x, name, chosenSNPs=NULL){
  df <- read.table(x,header=TRUE)
  
  df$col="Others"
  
  if (!is.null(chosenSNPs)){
    df$selected <- 0
    df$selected[df$SNP %in% chosenSNPs] <- 1
  }
  
  df$col[df$selected>0] <- 'Selected'

  df$col[df$SNP == "5:23532534" | df$SNP == "rs6889665"] <- "rs6889665"
  df$col[df$SNP == "5:23492211" | df$SNP == "rs10043097"] <- "rs10043097"
  
  df$col <- factor(df$col, levels=c('Selected',
                                    'rs6889665',
                                    'rs10043097',
                                    'Others'))
  df$test <- name
  g <- ggplot(df,aes(x=BP/1000000,y=-10*log10(P),color=factor(col))) + 
    geom_point(size=2) + 
    geom_point(data=df[df$col=="rs6889665",],size=3) + 
    coord_cartesian(xlim=c(23.3,23.7)) + 
    scale_color_manual(values=c('red','dodgerblue2',"pink",'grey50'))
  
  return(list(data=df,
              plot=g))
}

drawManhattanPlot <- function(dfX, d=NULL, midpoint=23527218){
  
  if (is.null(d)){
    dfChk <- dfX %>%
      mutate(pos=BP/1000000) 
    lowLimit <- 0
    highLimit <- max(dfChk$BP)/1000000
    d <- max(dfChk$BP)
  }else{
    lowLimit  <- (midpoint-d)/1000000
    highLimit <- (midpoint+d)/1000000
  }
  
  dfPlot <- dfX %>%
    mutate(pos=BP/1000000) %>%
    dplyr:::filter(pos > lowLimit, pos < highLimit) %>%
    mutate(test=gsub("like","-type",test))
  
  dfPrdm9 <- read.table('PRDM9.transcript.bed',header=FALSE) %>% 
    mutate(cs=V1, from=V2/1000000, to=V3/1000000)
  dfExons <- read.table('PRDM9.exons.bed'     ,header=FALSE) %>% 
    mutate(cs=V1, from=V2/1000000, to=V3/1000000)
  dfPrZFA <- read.table('PRDM9.ZFDomains.bed' ,header=FALSE) %>% 
    mutate(cs=V1, from=V2/1000000, to=V3/1000000)
  
  gPRDM9 <- ggplot(dfPlot) + 
    geom_rect(data=dfPrdm9,aes(xmin=from,xmax=to),
              ymin=-3,ymax =70,
              fill='chartreuse',alpha=.05)
  
  if (d <= 500000){
    gPRDM9 <- gPRDM9 + 
      geom_rect(data=dfPrdm9,aes(xmin=from,xmax=to),ymin=-3,ymax =70,fill='forestgreen',alpha=.1) + 
      geom_segment(y=-3,yend=-3,data=dfPrdm9,aes(x=from,xend=to),lwd=.5,color='forestgreen') +
      geom_segment(y=-3,yend=-3,data=dfExons,aes(x=from,xend=to),lwd=3,color='forestgreen') +
      geom_segment(y=-3,yend=-3,data=dfPrZFA,aes(x=from,xend=to),lwd=3,color='grey40') 
  }
  
  if (d > 1000000){
    gPRDM9 <- gPRDM9 + 
      geom_vline(data=dfPrdm9,aes(xintercept=from),
                 color='forestgreen',alpha=.2)
  }
  
  g <- gPRDM9 + 
    geom_point(aes(x=pos, y=-1*log10(P), color=factor(col)), size=.1) + 
    geom_point(data=dfPlot[dfPlot$col=="Selected",],aes(x=pos, y=-1*log10(P), color=factor(col)), size=.5) + 
    geom_point(data=dfPlot[dfPlot$col=="rs6889665",],aes(x=pos, y=-1*log10(P), color=factor(col)), size=.4) + 
    geom_point(data=dfPlot[dfPlot$col=="rs10043097",],aes(x=pos, y=-1*log10(P), color=factor(col)), size=.4) + 
    facet_wrap(~test,nrow=1) +
    scale_color_manual('SNP type',values=c('dodgerblue2','red',"magenta",'grey50')) + 
    xlab('Position on hg38 chr5 (Mb)') + 
    ylab('Association score\n(-log10(P))') + 
    theme(legend.position=c(0,1),
          legend.justification=c(0,1),
          legend.background=element_blank(),
          legend.key.size=unit(0.3,'cm')) 
  
  if (d <= 500000){
    g <- g + coord_cartesian(xlim=c(lowLimit,highLimit), ylim=c(-5,75),expand=FALSE)   
  }else{
    g <- g + coord_cartesian(xlim=c(lowLimit,highLimit), ylim=c(-1,75),expand=FALSE)   
  }
  
  return(g)
}

drawPenetrancePlots <- function(SNP=NULL){
  for (p in c('Alike','Clike','A','B','C','L14')){
    dS <- getSNPpenetrance(phenotypeID = p, snpid=SNP) 
    df <- getPlotData(dS)  
    dS$phenotype <- p
    df$simple_all$phenotype <- p
    df$simple_pop$phenotype <- p
    df$detail_all$phenotype <- p
    df$detail_pop$phenotype <- p
    
    if (p == 'Alike'){
      dfI  <- dS
      dfSA <- df$simple_all
      dfSP <- df$simple_pop
      dfDA <- df$detail_all
      dfDP <- df$detail_pop
    }else{
      dfI  <- rbind(dfI,dS)
      dfSA <- rbind(dfSA,df$simple_all)
      dfSP <- rbind(dfSP,df$simple_pop)
      dfDA <- rbind(dfDA,df$detail_all)
      dfDP <- rbind(dfDP,df$detail_pop)
    }
  }
  
  dfSA                <- fillOut(dfSA) %>% 
                         mutate(phenotype = gsub("like","-type",phenotype),
                                facetID   = paste0(phenotype," : ",ID),
                                pc       = ifelse(cnt/tot>0,cnt/tot*100,NA))

  dfSA$phenotype      <- factor(dfSA$phenotype,levels=c('A-type','C-type','A','B','C','L14'))
  dfSP$pop            <- factor(dfSP$pop,levels=c('FIN','TSI','PEL','PJL','CHB','YRI','LWK'))
  dfSA$prdm9GT        <- factor(dfSA$prdm9GT,levels=c('NONE','HET','HOM'))

  dfSA$facetID        <- factor(dfSA$facetID,levels=c(dfSA$facetID[dfSA$phenotype == 'A-type'][1],
                                                      dfSA$facetID[dfSA$phenotype == 'C-type'][1],
                                                      dfSA$facetID[dfSA$phenotype == 'A'][1],
                                                      dfSA$facetID[dfSA$phenotype == 'B'][1],
                                                      dfSA$facetID[dfSA$phenotype == 'C'][1],
                                                      dfSA$facetID[dfSA$phenotype == 'L14'][1]))

  gSA <- ggplot(dfSA,aes(x=SNP,y=prdm9GT, fill=pc)) + 
    geom_tile(color='black',lwd=.2) + 
    facet_wrap(~facetID,scales='free_x',ncol=6) + 
    geom_text(aes(label=paste0(round(cnt/tot*100),"%\nn=",cnt)),size=7*5/14) + 
    scale_fill_gradient2(low='white',mid='pink',high='red',midpoint=50,na.value = 'grey90') + 
    ylab('') + 
    theme(legend.position='none',strip.text=element_text(size=7)) + 
    xlab('Nucleotide')
  
  gSABubblePlot <- ggplot(dfSA,aes(x=SNP,y=prdm9GT)) + 
    geom_point(shape=21,color='white',aes(size=pc,fill=pc))+ 
    scale_size_continuous(range=c(1,15)) + 
    facet_wrap(~facetID,scales='free_x',ncol=6) + 
    geom_text(aes(label=paste0(round(cnt/tot*100),"%\nn=",cnt)),size=7*5/14) + 
    scale_fill_gradient2(low='grey70',mid='pink',high='firebrick1',midpoint=50,na.value = 'grey90') + 
    ylab('') + 
    theme(legend.position='none',
          panel.border=element_rect(fill=NA),
          strip.text=element_text(size=7)) + 
    xlab('Nucleotide')
  
  gBar <- ggplot(dfSA ,aes(x=SNP,fill=prdm9GT, y=cnt/tot*100))  + 
    geom_bar(stat='identity',position=position_dodge()) + 
    facet_grid(prdm9GT~paste0(phenotype," : ",ID),scales='free_x') +
    theme(panel.border=element_rect(fill=NA),
          panel.grid.major.y = element_line(color='grey70',size=.3),
          legend.position='none')+ geom_hline(yintercept=0,lwd=.3) +
    geom_hline(yintercept=50,color='magenta',lty='dashed',lwd=.3) + 
    scale_y_continuous(breaks=c(0,25,50,75,100),
                       labels=c('0','',"50",'',100)) + 
    ylab('Pentrance (%)') 
  
  dfSP                <- fillOutPop(dfSP) %>% 
                         mutate(phenotype = gsub("like","-type",phenotype),
                                facetID   = paste0(phenotype," : ",ID),
                                pc       = ifelse(cnt/tot>0,cnt/tot*100,NA))
  
  dfSP$phenotype      <- factor(dfSP$phenotype,levels=c('A-type','C-type','A','B','C','L14'))
  dfSP$prdm9GT        <- factor(dfSP$prdm9GT,levels=c('NONE','HET','HOM'))
  dfSP$pop            <- factor(dfSP$pop,levels=c('FIN','TSI','PEL','PJL','CHB','YRI','LWK'))
  dfSP$facetID        <- factor(dfSP$facetID,levels=c(dfSP$facetID[dfSP$phenotype == 'A-type'][1],
                                                      dfSP$facetID[dfSP$phenotype == 'C-type'][1],
                                                      dfSP$facetID[dfSP$phenotype == 'A'][1],
                                                      dfSP$facetID[dfSP$phenotype == 'B'][1],
                                                      dfSP$facetID[dfSP$phenotype == 'C'][1],
                                                      dfSP$facetID[dfSP$phenotype == 'L14'][1]))
  
  gSP <- ggplot(dfSP,aes(x=SNP,y=prdm9GT, fill=pc)) + 
    geom_tile(color='black',lwd=.2) + 
    facet_grid(pop~facetID,scales='free_x') + 
    geom_text(aes(label=cnt),size=7*5/14) + 
    scale_fill_gradient2(low='dodgerblue1',mid='green',high='yellow',midpoint=50,na.value = 'grey90') + 
    ylab('') + 
    theme(legend.position='none',strip.text=element_text(size=7),panel.border=element_rect(fill=NA)) + 
    xlab('Nucleotide')
  #geom_text(aes(label=paste0(round(cnt/tot*100),"%\nn=",cnt)),size=7*5/14) + 
  
  gSPBubblePlot <- ggplot(dfSP,aes(x=SNP,y=prdm9GT)) + 
    geom_point(shape=21,color='white',aes(size=pc, fill=pc))+ 
    facet_grid(pop~facetID,scales='free_x') + 
    geom_text(aes(label=cnt),size=7*5/14) + 
    scale_fill_gradient2(low='grey70',mid='pink',high='firebrick1',midpoint=50,na.value = 'grey90') + 
    scale_size_continuous(range=c(1,12)) + 
    ylab('') + 
    theme(legend.position='none',strip.text=element_text(size=7),panel.border=element_rect(fill=NA)) + 
    xlab('Nucleotide')
  
  return(list(gSimple=gSA,
              gSimpleBar=gBar,
              gSimpleBubb=gSABubblePlot,
              gPop=gSP,
              gPopBubb=gSPBubblePlot,
              dataSimple=dfSA,
              dataPop=dfSP))
}

lAll   <- drawPenetrancePlots(NULL)
lHinch <- drawPenetrancePlots('rs6889665')

SNPsUsed <- lAll$dataSimple[,c('ID','phenotype')] %>% distinct()

#-----------------------------------------------------------
# Manhattan Plots

gAl <- plinkplot('prdm9AS.Alike.mod.cmh',"Alike carriers",SNPsUsed %>% dplyr:::filter(phenotype == 'A-type') %>% select(ID))
gCl <- plinkplot('prdm9AS.Clike.mod.cmh',"Clike carriers",SNPsUsed %>% dplyr:::filter(phenotype == 'C-type') %>% select(ID))

dfLike <- rbind(gAl$data,
                gCl$data)

dfLike$test <- factor(dfLike$test,levels=paste0(c('A-type','C-type','A','C','B','L14'),' carriers'))

gA <- plinkplot('prdm9AS.A.mod.cmh',"A carriers",chosenSNPs = SNPsUsed %>% dplyr:::filter(phenotype == 'A') %>% select(ID))
gB <- plinkplot('prdm9AS.B.mod.cmh',"B carriers",SNPsUsed %>% dplyr:::filter(phenotype == 'B') %>% select(ID))
gC <- plinkplot('prdm9AS.C.mod.cmh',"C carriers",SNPsUsed %>% dplyr:::filter(phenotype == 'C') %>% select(ID))
gL14 <- plinkplot('prdm9AS.L14.mod.cmh',"L14 carriers",SNPsUsed %>% dplyr:::filter(phenotype == 'L14') %>% select(ID))

df4X <- rbind(gA$data,
              gB$data,
              gC$data,
              gL14$data)

df4X$test <- factor(df4X$test,levels=paste0(c('Alike','Clike','A','B','C','L14'),' carriers'))

gLikeCS   <- drawManhattanPlot(dfLike)
g4XCS     <- drawManhattanPlot(df4X)

gLike40M  <- drawManhattanPlot(dfLike, d = 20000000)
g4X40M    <- drawManhattanPlot(df4X  , d = 20000000)
gLike150k <- drawManhattanPlot(dfLike, d = 50000)
g4X150k   <- drawManhattanPlot(df4X  , d = 50000)

ltop <- theme(legend.position='top',
              legend.title=element_blank(),
              legend.key.size=unit(0.2,'cm'),
              legend.justification = 'center')

maxM    <- theme(plot.margin = unit(c(0,0,0,0),'cm'))

noStrip <- theme(strip.background=element_blank(),
                 strip.text=element_blank())

normStripTxt <- theme(strip.text=element_text(size=7))

box <- theme(panel.border=element_rect(fill=NA))
                 
nl <- theme(legend.position='none')

ggLeg  <- get_legend(gLike150k + ltop + maxM + theme(legend.text=element_text(size=7)))
xLAll  <- scale_x_continuous(breaks=c(0,10,20,30,40),
                             labels=c('',10,'',30,''))

xL150  <- scale_x_continuous(breaks=c(23.5,23.55),
                             labels=c('23.50',23.55))

ggLike <- ggarrange(gLikeCS  + box + xLAll + nl + maxM + normStripTxt + xlab(''),
                    gLike150k + box + xL150 + nl + maxM + noStrip,
                    nrow=2,ncol=1)

gg4X   <- ggarrange(g4XCS  + box + xLAll + nl + maxM  + ylab('') + normStripTxt + xlab(''),
                    g4X150k + box + xL150 + nl + maxM  + ylab('') + noStrip,
                    nrow=2,ncol=1)

ggAll  <- ggarrange(ggLike,
                    gg4X,
                    nrow=1,ncol=2,
                    widths=c(2,4))

ggFinal <- ggarrange(ggLeg,ggAll,
                     ncol=1,nrow=2,
                     heights=c(1,10))

# ggsave('Alleva_et_all_ManhattanPlots.png',
#        plot=ggFinal,dpi=400,height=6*.5,width=14*.5)

gHits   <- ggarrange(lHinch$gSimpleBubb + xlab('') + maxM,
                     lAll$gSimpleBubb + maxM,
                     ncol=1,
                     labels=c('B','C'),
                     font.label = list(size=8,face='bold'),
                     hjust=0,vjust=1)

gOut    <- ggarrange(ggFinal,
                     gHits,
                     ncol=1,
                     labels=c('A',''),
                     font.label = list(size=8,face='bold'),
                     hjust=0,vjust=1,
                     heights=c(1,1.4))

ggsave('Alleva_et_al_Figure_SNP_Associations.png',gOut,width=8,height=7)
ggsave('Alleva_et_al_Figure_SNP_Associations.pdf',gOut,width=8,height=7)

gSupp   <- lAll$gPopBubb
ggsave('Alleva_et_al_Supplement_Associations_Plot.png',gSupp,width=8,height=9)
ggsave('Alleva_et_al_Supplement_Associations_Plot.pdf',gSupp,width=8,height=9)
