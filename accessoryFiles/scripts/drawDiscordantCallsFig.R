library(plyr,include.only = "join")
library(dplyr)
library(stringr)
library(ggplot2)
library(maps)
library(RColorBrewer)
library(ggpubr)
library(tidyverse)
library(reshape2)
library(ggrepel)
library(png)
library(grid)
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

theme7point()
options(scipen = 999) ## To disable scientific notation
###############################
getGT <- function(retGT=1,gt1,gt2,s11,s12,s21,s22){
  gt1 <- gsub("M\\d+($|\\/)","New\\1",gt1,perl=TRUE)
  gt2 <- gsub("M\\d+($|\\/)","New\\1",gt2,perl=TRUE)
  splitGT1 <- strsplit(gt1,"/")
  splitGT2 <- strsplit(gt2,"/")
  if(s11 == s21){
    if (s12 == s22){
      gt1 = paste0("-/-")
      gt2 = paste0("-/-")
      if (retGT==1){
        return(gt1)
      }else{
        return(gt2)
      }
    }else{
      gt1 = paste0("-/",splitGT1[[1]][2])
      gt2 = paste0("-/",splitGT2[[1]][2])
      if (retGT==1){
        return(gt1)
      }else{
        return(gt2)
      }
    }
  }
  
  if(s11 == s22){
    if (s12 == s21){
      gt1 = paste0("-/-")
      gt2 = paste0("-/-")             
      if (retGT==1){
        return(gt1)
      }else{
        return(gt2)
      }
    }else{
      gt1 = paste0("-/",splitGT1[[1]][2])
      gt2 = paste0(splitGT2[[1]][1],'/-')
      if (retGT==1){
        return(gt1)
      }else{
        return(gt2)
      }
    }
  }
  
  if(s12 == s21){
    if (s11 == s22){
      gt1 = paste0("-/-")
      gt2 = paste0("-/-")             
      if (retGT==1){
        return(gt1)
      }else{
        return(gt2)
      }
    }else{
      gt1 = paste0(splitGT1[[1]][1],'/-')
      gt2 = paste0("-/",splitGT2[[1]][2])
      if (retGT==1){
        return(gt1)
      }else{
        return(gt2)
      }
    }
  }
  
  if(s12 == s22){
    if (s11 == s21){
      gt1 = paste0("-/-")
      gt2 = paste0("-/-")             
      if (retGT==1){
        return(gt1)
      }else{
        return(gt2)
      }
    }else{
      gt1 = paste0(splitGT1[[1]][1],'/-')
      gt2 = paste0(splitGT2[[1]][1],'/-')
      if (retGT==1){
        return(gt1)
      }else{
        return(gt2)
      }
    }
  }
  
  if (retGT==1){
    return(gt1)
  }else{
    return(gt2)
  }
  
}

drawZFsizesComparison <- function(dfNanopore, dfPacbio, sameLen=TRUE, nagree=0, labelAlleles=TRUE){
  if (nagree > 0){
    dfDiscordant <- dfNanopore %>%
      inner_join(dfPacbio, by="id") %>%
      filter(!grepl('^A[AN]-',id)) %>%
      mutate(agree  = ifelse((bonSeq1 == pbSeq1 | bonSeq1 == pbSeq2) & 
                               (bonSeq2 == pbSeq1 | bonSeq2 == pbSeq2), TRUE, FALSE),
             homhet = ifelse(bonHH == pbHH, pbHH, paste0(bonHH, ":", pbHH))) %>% 
      select(id,bonGT,bonHH,bonN,pbGT,pbHH,pbN,agree,homhet,pbSz1,bonSz1,pbSz2,bonSz2,pbSeq1,pbSeq2,bonSeq1,bonSeq2) %>%
      dplyr:::filter(agree)
    
    dfDiscordant <- dfDiscordant[sample(1:length(dfDiscordant$id),nagree),]
  }else{
    dfDiscordant <- dfNanopore %>%
      inner_join(dfPacbio, by="id") %>%
      filter(!grepl('^A[AN]-',id)) %>%
      mutate(agree  = ifelse((bonSeq1 == pbSeq1 | bonSeq1 == pbSeq2) & 
                               (bonSeq2 == pbSeq1 | bonSeq2 == pbSeq2), TRUE, FALSE),
             homhet = ifelse(bonHH == pbHH, pbHH, paste0(bonHH, ":", pbHH))) %>% 
      select(id,bonGT,bonHH,bonN,pbGT,pbHH,pbN,agree,homhet,pbSz1,bonSz1,pbSz2,bonSz2,pbSeq1,pbSeq2,bonSeq1,bonSeq2) %>%
      dplyr:::filter(!agree)
  }
  
  dfZFsizes_nanopore <- fread('nanopore_ZFAsizes.txt', header=TRUE) %>% 
    mutate(id=gsub('\\.\\S+$','',id,perl=TRUE),type='Nanopore') %>% 
    select(id,type,size)
  
  dfZFsizes_pacbio   <- fread('pacbio_ZFAsizes.txt'  , header=TRUE) %>%
    mutate(id=gsub('\\.\\S+$','',id,perl=TRUE),type='Pacbio'  ) %>%
    select(id,type,size)
  
  if (sameLen){
    dfSelect <- dfDiscordant %>%
      mutate(pbC=paste0(pbSz1,":",pbSz2),
             bonC1=paste0(bonSz1,":",bonSz2),
             bonC2=paste0(bonSz2,":",bonSz1)) %>%
      dplyr:::filter((pbC == bonC1 | pbC == bonC2))
  }else{
    dfSelect <- dfDiscordant %>%
      mutate(pbC=paste0(pbSz1,":",pbSz2),
             bonC1=paste0(bonSz1,":",bonSz2),
             bonC2=paste0(bonSz2,":",bonSz1)) %>%
      dplyr:::filter(!(pbC == bonC1 | pbC == bonC2))
  }
  dfPB  <- join(dfSelect, dfZFsizes_pacbio  , by='id')
  dfONT <- join(dfSelect, dfZFsizes_nanopore, by='id')
  
  dfPBONTplot <- rbind(dfONT,dfPB) %>% 
    rowwise() %>% 
    mutate(lbl    = id,
           gtPacbio   = getGT(1,pbGT,bonGT,pbSeq1,pbSeq2,bonSeq1,bonSeq2),
           gtNanopore = getGT(2,pbGT,bonGT,pbSeq1,pbSeq2,bonSeq1,bonSeq2),
           lblPB  = paste0("\n\n",gtPacbio," : ",paste0(pbSz1,"/",pbSz2)),
           lblONT = paste0("\n",gtNanopore," : ",paste0(bonSz1,":",bonSz2)))

  g <- ggplot(dfPBONTplot,aes(x=size,fill=type)) +
    geom_bar(position=position_dodge(),
             color='black',lwd=.2,width=.6) + 
    xlab('# ZFs') + ylab('# reads') + 
    theme(legend.position=c(1,0),
          legend.title=element_blank(),
          legend.justification=c(1,0),
          strip.background=element_rect(fill='grey90'),
          strip.text=element_blank(),
          legend.key.size=unit(0.5,'cm'),
          panel.grid.major.x = element_line(size=.2,linetype='dotted',color='grey70')) + 
    facet_wrap(~id) +
    scale_fill_manual(values=c('#0072B2','#E69F00')) +
    geom_text(x=-Inf,y=Inf,hjust=0,vjust=1,
              check_overlap = TRUE, size=7*5/14, 
              fontface='bold',
              color='black',
              aes(label=id)) + 
    geom_text(x=-Inf,y=Inf,hjust=0,vjust=1,
              check_overlap = TRUE, size=7*5/14, 
              color='dodgerblue4',
              aes(label=lblONT)) + 
    geom_text(x=-Inf,y=Inf,hjust=0,vjust=1,
              check_overlap = TRUE, size=7*5/14, 
              color='#E69F00',
              aes(label=lblPB)) + 
    coord_cartesian(xlim=c(8,20)) + 
    scale_x_continuous(breaks=seq(9,19,2))
  
  g
  
  return(g)
}

plotLengthBiases <- function(dfNanopore, dfPacbio){
  dfOK <- dfNanopore %>%
    inner_join(dfPacbio, by="id") %>%
    filter(!grepl('^A[AN]-',id)) %>%
    mutate(agree  = ifelse((bonSeq1 == pbSeq1 | bonSeq1 == pbSeq2) & 
                             (bonSeq2 == pbSeq1 | bonSeq2 == pbSeq2), TRUE, FALSE),
           homhet = ifelse(bonHH == pbHH, pbHH, paste0(bonHH, ":", pbHH))) %>% 
    select(id,bonGT,bonHH,bonN,pbGT,pbHH,pbN,agree,homhet,pbSz1,bonSz1,pbSz2,bonSz2) %>%
    dplyr:::filter(agree)

  dfZFsizes_nanopore <- fread('nanopore_ZFAsizes.txt', header=TRUE) %>% mutate(id=gsub('\\.\\S+$','',id,perl=TRUE),type='Nanopore') %>% select(id,type,size)
  dfZFsizes_pacbio   <- fread('pacbio_ZFAsizes.txt'  , header=TRUE) %>% mutate(id=gsub('\\.\\S+$','',id,perl=TRUE),type='Pacbio'  ) %>% select(id,type,size)
  
  dfONT <- join(dfOK,dfZFsizes_nanopore,by='id')
  dfPB  <- join(dfOK,dfZFsizes_pacbio,by='id')
  
  ## Remove homozygotes
  dfPlotX2 <- rbind(dfONT,dfPB) %>%
    group_by(id,type) %>%
    add_count(name='tot') %>%
    dplyr:::filter(size == pbSz1 | size == pbSz2) %>%
    dplyr:::filter(pbSz1 != pbSz2) %>%
    group_by(id,type,size,tot) %>%
    count() %>% 
    mutate(pc=n/tot*100)
  
  gAmpBias <- ggplot(dfPlotX2,aes(x=size,y=pc,group=factor(size))) + 
    geom_hline(yintercept=50,color='magenta',lwd=.2,lty='dotted') + 
    geom_boxplot(fill='salmon',lwd=.4,varwidth=TRUE,outlier.size=.1) + 
    scale_y_continuous(breaks=c(0,25,50,100)) + 
    scale_x_continuous(breaks=seq(7,19,2)) +
    facet_wrap(~type,ncol=1) +
    theme(legend.position='none',
          strip.text=element_text(face='bold',size=7)) + 
    xlab('# ZFs') + 
    ylab('Percent of reads in heterozygote (%)')
  
  return(gAmpBias)
}
################################################################################

prdm9_Bonito <- read.table('prdm9_haplotypes.bonito1d.tab', header=TRUE, stringsAsFactors = FALSE)
prdm9_Pacbio <- read.table('prdm9_haplotypes.pacbio.tab'  , header=TRUE, stringsAsFactors = FALSE)

dfBon <- prdm9_Bonito %>% 
  filter(!grepl('(Unk|noZFA|NA)',perl=TRUE,diploid)) %>%
  mutate(homhet = ifelse(homhet == 'het' & size_1 == size_2,'Het\n(equal)',
                         ifelse(homhet == 'het' & size_1 != size_2, 'Het\n(unequal)', 'Hom'))) %>%
  select(id,bonGT = diploid, 
         bonHH = homhet, 
         bonSeq1 = seq_1, 
         bonSeq2 = seq_2, 
         bonN = nzfseqs,
         bonSz1 = size_1,
         bonSz2 = size_2)  

dfPB <- prdm9_Pacbio %>% 
  filter(!grepl('(Unk|noZFA|NA)',perl=TRUE,diploid)) %>%
  mutate(homhet = ifelse(homhet == 'het' & size_1 == size_2,'Het\n(equal)',
                         ifelse(homhet == 'het' & size_1 != size_2, 'Het\n(unequal)', 'Hom'))) %>%
  select(id,pbGT = diploid, 
         pbHH = homhet, 
         pbSeq1 = seq_1, 
         pbSeq2 = seq_2, 
         pbN = nzfseqs,
         pbSz1 = size_1,
         pbSz2 = size_2) 

gDiffLen   <- drawZFsizesComparison(dfBon,dfPB, FALSE, FALSE) + theme(legend.position='none')
gSameLen   <- drawZFsizesComparison(dfBon,dfPB, TRUE)
gSameLenOK <- drawZFsizesComparison(dfBon,dfPB, TRUE, nagree=6)

gLenBiases <- plotLengthBiases(dfBon, dfPB)
g1 <- ggarrange(gSameLenOK + theme(legend.position='none'),gSameLen,gDiffLen,
                heights = c(2,2,3),
                ncol=1,nrow=3,
                labels =c('A','B','C'),
                font.label = c(size=9,hjust=0,vjust=1))

gFig <- ggarrange(g1,gLenBiases,
          ncol=2,nrow=1,
          widths=c(3,1),
          labels=c('','D'),
          font.label = c(size=9,hjust=0,vjust=1))

ggsave('Alleva_et_al_discordant_calls_PB_V_ONT.png',gFig,height=7,width=7,dpi=300)
ggsave('Alleva_et_al_discordant_calls_PB_V_ONT.pdf',gFig,height=7,width=7)

strcmpK <- function(a,b){
sa      <- strsplit(a,"")
sb      <- strsplit(b,"")
diffs   <- (sa[[1]] != sb[[1]])
diffPos <- which(diffs)
return(list(n=diffPos,
seqA=sa[[1]][diffPos],
seqB=sb[[1]][diffPos]))
}
