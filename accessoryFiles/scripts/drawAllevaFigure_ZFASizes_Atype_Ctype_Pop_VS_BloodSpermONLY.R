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

###############################

prdm9_allele_details <- read.table('atype_ctype.txt',
                                   header=TRUE,
                                   stringsAsFactors = FALSE,
                                   skip = 11)

allelesInPop_Jeffreys2013 <- c(paste0("L",c(1,2,3,5,6,10,11,15,19,21,22)),'C','E')

formatPval <- function(p,pName=NULL,type,altP=0.05){
  
  if (is.null(pName)){
    pNm <- 'P'
  }else{
    pNm <- paste0("P(",pName,")")
  }
  
  ## For P == 0, show P < altP
  if (p == 0){
    if (p > 0.001){
      return(paste0(pNm,' < ',round(altP,3),'; ',type))    
    }
    
    p <- round(log10(altP))
    
    return(paste0(pNm,' < 10',altP,'; ',type)) 
  }else{
    if (p > 0.01){
      return(paste0('n.s. [',pNm,' = ',round(p,2),']; ',type))
    }
    
    if (p > 0.001){
      return(paste0(pNm,' = ',round(p,3),'; ',type))    
    }
    
    p<-round(log10(p))
    
    return(paste0(pNm,' = 10',p,'; ',type)) 
  }
}

compare_pop_with_BS_alleles <- function(df){
  
  set.seed(42)
  
  dfStat <- prdm9_allele_details %>% 
    rowwise() %>% 
    mutate(zf_len = str_length(zf_code)/2,
           foundWhere=ifelse(in_pop,"POP","BS")) %>% 
    dplyr:::filter(grepl("(A|C|L14)v:",short_ID) | short_ID %in% allelesInPop_Jeffreys2013) %>% 
    select(ID,short_ID,zf_len,ACtype,foundWhere) %>%
    dplyr:::filter(grepl(":(A|C|L14)-(A|C|L14)$",ID))
  
  dA_All <- dfStat$zf_len[dfStat$ACtype == 'A']
  dC_All <- dfStat$zf_len[dfStat$ACtype == 'C']

  dA_Pop <- dfStat$zf_len[dfStat$foundWhere == 'POP' & dfStat$ACtype == 'A']
  dC_Pop <- dfStat$zf_len[dfStat$foundWhere == 'POP' & dfStat$ACtype == 'C']
  dA_BS  <- dfStat$zf_len[dfStat$foundWhere == 'BS'  & dfStat$ACtype == 'A']
  dC_BS  <- dfStat$zf_len[dfStat$foundWhere == 'BS'  & dfStat$ACtype == 'C']

  statWxAtype <- formatPval(wilcox.test(dA_BS,dA_Pop,exact=FALSE)$p.value,type="Wilcoxon")
  statWxCtype <- formatPval(wilcox.test(dC_BS,dC_Pop,exact=FALSE)$p.value,type="Wilcoxon")
  statWxBS    <- formatPval(wilcox.test(dA_BS,dC_BS,exact=FALSE)$p.value,type="Wilcoxon")
  statWxPop   <- formatPval(wilcox.test(dA_Pop,dC_Pop,exact=FALSE)$p.value,type="Wilcoxon")
  
  vAtop <- as.numeric(names(sort(summary(as.factor(dA_Pop)),decreasing=TRUE))[1:2])
  vCtop <- as.numeric(names(sort(summary(as.factor(dC_Pop)),decreasing=TRUE))[1:2])
  nA    <- length(dA_Pop)
  nC    <- length(dC_Pop)
  
  nObs  <- c(sum(dA_Pop == vAtop[1]),
            sum(dA_Pop == vAtop[2]),
            sum(dC_Pop == vCtop[1]),
            sum(dC_Pop == vCtop[2]))
  
  nPerm <- 100000
  v <- matrix(nrow = 4, ncol = nPerm)
  for (i in 1:nPerm){
    v[1,i] <- sum(sample(dA_All,nA,TRUE) == vAtop[1])
    v[2,i] <- sum(sample(dA_All,nA,TRUE) == vAtop[2])
    v[3,i] <- sum(sample(dC_All,nC,TRUE) == vCtop[1])
    v[4,i] <- sum(sample(dC_All,nC,TRUE) == vCtop[2])
  }
  
  nP <- c(formatPval(sum(v[1,] >= nObs[1])/nPerm,vAtop[1],'Permutation',altP=1/nPerm),
          formatPval(sum(v[2,] >= nObs[2])/nPerm,vAtop[2],'Permutation',altP=1/nPerm),
          formatPval(sum(v[3,] >= nObs[3])/nPerm,vCtop[1],'Permutation',altP=1/nPerm),
          formatPval(sum(v[4,] >= nObs[4])/nPerm,vCtop[2],'Permutation',altP=1/nPerm))
  
  dfPlotMe <- prdm9_allele_details %>% 
    rowwise() %>% 
    mutate(zf_len     = str_length(zf_code)/2,
           foundWhere = ifelse(in_pop,"In population","Blood/Sperm ONLY"),
           pAC        = ifelse(ACtype == 'A',
                                paste0("A-type\n",statWxAtype,"\n",nP[1],"\n",nP[2]),
                                paste0("C-type\n",statWxCtype,"\n",nP[3],"\n",nP[4])),
           pPopBS     = ifelse(in_pop, statWxPop, statWxBS),
           ACtype     = paste0(ACtype,'-type')) %>% 
    dplyr:::filter(grepl("(A|C|L14)v:",short_ID) | short_ID %in% namedinSB) %>% 
    dplyr:::filter(grepl(":(A|C|L14)-(A|C|L14)$",ID)) %>%
    select(short_ID,zf_len,ACtype,foundWhere,pAC,pPopBS) %>% 
    group_by(foundWhere,ACtype,zf_len,pAC,pPopBS) %>% 
    count() %>% 
    group_by(ACtype,foundWhere) %>% 
    add_tally(n,name="tot")
  
  g <- ggplot(dfPlotMe) + 
    geom_bar(aes(x=zf_len,y=n/tot*100,fill=ACtype),
             color='black',lwd=.3,
             position=position_dodge(),
             width=.75,stat="identity") + 
    facet_grid(foundWhere+pPopBS~pAC) + 
    geom_text(x=Inf,y=Inf,
              hjust=1,vjust=1,
              aes(label=paste0("\nN = ",tot,"  ")),
              size=7*5/14,check_overlap=TRUE) + 
    coord_cartesian(xlim=c(6,23),ylim=c(0,70),expand=FALSE) + 
    scale_y_continuous(breaks=c(0,25,50)) + 
    scale_x_continuous(minor_breaks=seq(0,40,1)) + 
    theme(legend.position='none',
          strip.text=element_text(size=7),
          panel.grid = element_line(size=.2,color='grey80')) +
    xlab('ZF array size (# ZFs)') + 
    ylab('Alleles (%)') +
    scale_fill_manual(values = c('A-type' = 'purple',
                                 'C-type' = 'chartreuse3'))

  return (g)
}

gLengthComp <- compare_pop_with_BS_alleles(prdm9_allele_details)

ggsave('Alleva_et_al_AlleleSizesAnalysis_BS_v_Pop_Alike_V_Clike_ONLY_ACL14.png',
       plot = gLengthComp,
       height=3.5,width=4,dpi=400)

ggsave('Alleva_et_al_AlleleSizesAnalysis_BS_v_Pop_Alike_V_Clike_ONLY_ACL14.pdf',
       plot = gLengthComp,
       height=3.5,width=4)
