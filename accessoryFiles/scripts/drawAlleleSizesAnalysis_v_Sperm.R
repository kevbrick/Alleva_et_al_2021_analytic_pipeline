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
checkSpermVariantDets <- function(){
  dfAll <- fread('PrZFA_alleles.details.txt',skip=11,header=TRUE)

  dfBS <- fread('humanPRDM9alleles.BloodandSpermVariants.txt',header=FALSE) %>%
    rename(V1='man',V2='zf_code',V3='ID') %>%
    inner_join(dfAll, 
               by=c('ID','zf_code')) %>%
    select(man, ID, short_ID, zf_code, 
           published_allele, in_pop, ABDdist, 
           CBDdist, ACtype) %>% 
    rowwise() %>% 
    mutate(zf_len     = str_length(zf_code)/2,
           foundWhere = ifelse(in_pop,"In population","Blood/Sperm ONLY"),
           ACtype     = paste0(ACtype,'-type'),
           parentGT   = gsub('^.+:(\\S+)$','\\1',man,perl=TRUE),
           pAllele1   = gsub('-.+$','',parentGT,perl=TRUE),
           pAllele2   = gsub('^.+-','',parentGT,perl=TRUE)) %>% 
    inner_join(dfAll %>% rowwise %>% mutate(p1len=str_length(zf_code)/2, p1AC=ACtype) %>%  select(ID,p1len,p1AC), 
               by=c('pAllele1'='ID')) %>%
    inner_join(dfAll %>% rowwise %>% mutate(p2len=str_length(zf_code)/2, p2AC=ACtype) %>% select(ID,p2len,p2AC), 
               by=c('pAllele2'='ID')) %>%
    select(ID,short_ID,zf_len,ACtype,man,parentGT,pAllele1,pAllele2,ABDdist,CBDdist,p1len,p2len,p1AC,p2AC) %>%
    dplyr:::filter(p1len==p2len, p1AC==p2AC) %>% 
    mutate(zfdiff = zf_len-p1len)

  dfPlot <- dfBS %>% 
           group_by(p1AC) %>% 
           add_count(name='tot') %>%
           group_by(p1AC,zfdiff,tot) %>%
           count()
  
  gLen <- ggplot(dfPlot)+ 
    geom_bar(aes(x=zfdiff,y=n/tot*100, fill=p1AC), 
             stat='identity',
             color='black',lwd=.2,width=.7) + 
    geom_vline(aes(xintercept=0),color='red') + 
    facet_grid(p1AC~.) + 
    scale_fill_manual('', values = c('A' = 'purple',
                                     'C' = 'chartreuse3')) +
    xlab('Child - Parent allele\nlength (# ZFs)') + 
    ylab('Sperm variants (%)') + 
    theme(legend.position='none',
          strip.text=element_blank(),
          strip.background=element_blank()) 
  
  dfPlot <- dfBS %>% 
    mutate(score=ifelse(p1AC=='A',ABDdist,CBDdist)) %>%
    group_by(p1AC) %>% 
    add_count(name='tot') %>%
    group_by(p1AC,score,tot) %>%
    count()
  
  gScore <- ggplot(dfPlot)+ 
    geom_bar(aes(x=score,y=n/tot*100, fill=p1AC), 
             stat='identity',
             color='black',lwd=.2,width=.7) + 
    facet_grid(p1AC~.) + 
    scale_fill_manual('', values = c('A' = 'purple',
                                     'C' = 'chartreuse3')) +
    xlab('Binding site distance (# AAs)') +
  ylab('Sperm variants (%)') + 
    theme(legend.position='none',
          strip.text=element_blank(),
          strip.background=element_blank()) 
  
  return(list(gLen=gLen,
              gSc=gScore))
}

getPopAlleleSizeDistByType <- function(){
  dfStatMe <- prdm9_allele_details %>% 
    rowwise() %>% 
    mutate(zf_len     = str_length(zf_code)/2,
           ACtype     = paste0(ACtype,'-type')) %>% 
    dplyr:::filter(in_pop) 
  
  # p <- formatPval(wilcox.test(dfStatMe$zf_len[dfStatMe$ACtype == 'A-type'],dfStatMe$zf_len[dfStatMe$ACtype == 'C-type'])$p.value,
  #                 type='Wilcoxon')
  p <- wilcox.test(dfStatMe$zf_len[dfStatMe$ACtype == 'A-type'],
                   dfStatMe$zf_len[dfStatMe$ACtype == 'C-type'])$p.value
  
  dfPlotMe <- prdm9_allele_details %>% 
    rowwise() %>% 
    mutate(zf_len     = str_length(zf_code)/2,
           ACtype     = paste0(ACtype,'-type')) %>% 
    dplyr:::filter(in_pop) %>%
    group_by(ACtype,zf_len) %>% 
    count() 
  
  p  <- paste("'P ='~","~10^-5~", ";Wilx")
  
  gTop  <- ggplot(dfStatMe,aes(x=zf_len,fill=ACtype)) + 
    geom_density(alpha=.2,adjust=2,lwd=.2) +
    scale_x_continuous(minor_breaks=seq(8,19,1),breaks=seq(8,19,5)) + 
    theme(legend.position=c(0.05,1),
          legend.justification=c(0,1),
          legend.title=element_blank(),
          legend.background=element_blank(),
          legend.key.height=unit(0.3,'cm'),
          legend.key.width=unit(0.3,'cm'),
          strip.text=element_text(size=7),
          panel.grid = element_line(size=.2,color='grey80')) +
    ylab('Frequency') +
    xlab('') + 
    scale_fill_manual('', values = c('A-type' = 'purple',
                                     'C-type' = 'chartreuse3'))+
    coord_cartesian(xlim=c(7.3,19.7),expand=FALSE) +
    annotate(geom='text',label=paste("P == 10 ^ -5 "),parse=TRUE,x=Inf,y=Inf,hjust=1,vjust=1,size=7*5/14) + 
    annotate(geom='text',label=paste0("  Median\n  (A-type)\n  = ",median(dfStatMe$zf_len[dfStatMe$ACtype == 'A-type']), ' ZFs'),
             x=-Inf,y=.1,hjust=0,vjust=0,size=7*5/14,color="purple") + 
    annotate(geom='text',label=paste0("  Median\n  (C-type)\n  = ",median(dfStatMe$zf_len[dfStatMe$ACtype == 'C-type']), ' ZFs'),
             x=Inf,y=.1,hjust=1,vjust=0,size=7*5/14,color="chartreuse3")
  
  gBott <- ggplot(dfPlotMe,aes(x=zf_len,y=n,fill=ACtype)) + 
    geom_bar(stat='identity',color='black',lwd=.2,
             position=position_dodge2(preserve = "single"),width=.6) + 
    scale_x_continuous(minor_breaks=seq(8,19,1),breaks=seq(8,19,5)) + 
    theme(legend.position='none',
          strip.text=element_text(size=7),
          panel.grid = element_line(size=.2,color='grey80')) +
    ylab('Alleles (#)') +
    xlab('Allele length (# ZFs)') +
    scale_fill_manual(values = c('A-type' = 'purple',
                                 'C-type' = 'chartreuse3'))
  coord_cartesian(xlim=c(7.2,19.8),expand=FALSE)  
  
  gRet <- ggarrange(gTop,gBott,ncol=1,nrow=2,align='v')
  
  return(gRet)
}

gPop <- getPopAlleleSizeDistByType()
lstSperm <- checkSpermVariantDets()

ggarrange(gPop,
          lstSperm$gLen,
          lstSperm$gSc,
          ncol=3,nrow=1,
          labels=c('D','E','F'),
          font.label = list(size=8,face='bold'),
          hjust=0,vjust=1)

ggsave('Alleva_et_al_AlleleSizesAnalysis_v_Sperm.png',
       height=2.6,width=8,dpi=400)

ggsave('Alleva_et_al_AlleleSizesAnalysis_v_Sperm.pdf',
       height=2.6,width=8)