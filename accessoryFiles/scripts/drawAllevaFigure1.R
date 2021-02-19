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

###########################################
drawMap <- function(bw=TRUE){
  world <- map_data("world")
  worldplot <- ggplot(fill='blue') +
    geom_polygon(data = world, aes(x=long, y = lat, group = group)) + 
    coord_fixed(1.3)
  
  dfSet <- data.frame(country   =c('Italy'  , 'Peru'    , 'China'   , 'Kenya' ,'Finland' ,'Nigeria','Pakistan'),
                      city      =c('Tuscany', 'Lima'    , 'Beijing' , 'Webuye','Helsinki','Ibadan' ,'Lahore'),
                      code      =c('TSI'    , 'PEL'     , 'CHB'     , "LWK"   , "FIN"    , "YRI"   , "PJL"),
                      N         =c(      114,         70,        120,      120,       100,      120,   108),
                      latitude  =c(43.355099, -12.045969,  39.90647 ,  0.6166 , 60.166628, 7.378606, 31.562192),
                      longtitude=c(11.029560, -77.030581, 116.391195, 34.7666 , 24.943508, 3.896993, 74.322852))
  
  ## Let's ditch many of the unnecessary elements
  plain_theme <- theme(
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5),
    legend.position='none'
  )
  
  if (bw){
    popColz <- rep('black',7)
  }else{
    popColz <- c('f95955','e77700','ffff00','3abe40','00afff','505050','7564d5')
  }
  
  if (bw){
  gMap <- ggplot() + 
    geom_polygon(data = world, 
                 aes(x=long, y = lat, group = group), 
                 fill='grey20',lwd=.1) + 
    geom_polygon(data = world[world$region %in% c('Italy'  , 'Peru'    , 'China'   , 'Kenya' ,'Finland' ,'Nigeria','Pakistan'),], 
                 aes(x=long, y = lat, group = group), 
                 fill='grey20',color=NA,lwd=.1) + 
    geom_point(data=dfSet,size=3,aes(y=latitude,x=longtitude),color='red') + 
    coord_cartesian(ylim=c(-52,70),
                    xlim=c(-165,175)) + 
    geom_label_repel(data=dfSet,
               aes(y=latitude,
                   x=longtitude,
                   label=paste0(code,": N = ",N),
                   color=code),
               size=7*5/14,
               lwd=0) +
    scale_color_manual(values=popColz)+
    plain_theme
  }else{
    gMap <- ggplot() + 
      geom_polygon(data = world, 
                   aes(x=long, y = lat, group = group), 
                   fill='grey20',lwd=.1) + 
      geom_polygon(data = world[world$region %in% c('Italy'  , 'Peru'    , 'China'   , 'Kenya' ,'Finland' ,'Nigeria','Pakistan'),], 
                   aes(x=long, y = lat, group = group), 
                   fill='grey20',color=NA,lwd=.1) + 
      geom_point(data=dfSet,size=3,aes(y=latitude,x=longtitude,color=code)) + 
      coord_cartesian(ylim=c(-52,70),
                      xlim=c(-165,175)) + 
      geom_label_repel(data=dfSet,
                       aes(y=latitude,
                           x=longtitude,
                           label=paste0(code,": N = ",N),
                           color=code),
                       size=7*5/14,
                       lwd=0) +
      scale_color_manual(values=popColz)+
      plain_theme
  }

  return(gMap)
}

###########################################
drawPRDM9alleles <- function(){
  ################################
  col_vector<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8',
                '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                '#bcf60c', '#fabebe', '#222222', '#e6beff',
                '#9a6324', '#999999', '#800000', '#aaffc3', 
                '#808000', '#ffd8b1', '#000075', '#808080')
  
  pr <- read.table('accessoryFiles/otherdata/prdm9alleles.txt',header=TRUE)
  
  pr$col <- 'black'
  pr$col[pr$domain == 'ZF'] <- pr$name[pr$domain == 'ZF']
  
  pr$allele <- factor(pr$allele, levels=c('L4','C','N','B','A'))
  
  prlines <- pr %>% group_by(allele) %>% filter(to == max(to)) %>% select(allele, to) %>% mutate(from = 1, to = to+4)
  
  dfAnnot <- pr %>% group_by(domain) %>% filter(from == min(from)) %>% 
    select(domain,from,to) %>% 
    mutate(to = ifelse (domain == 'ZF', 1026, to)) %>% distinct
  
  dfAnnot$domain <- factor(dfAnnot$domain,levels=c('KRAB','SET','ZF','ZF-array'))
  dfAnnot$domain[dfAnnot$domain == 'ZF'] <- 'ZF-array'
  
  gPRDM9 <- ggplot() + 
    geom_segment(data=prlines,aes(x=1,xend=to,y=allele,yend=allele)) + 
    geom_rect(data=pr,aes(xmin=from,
                          xmax=to,
                          ymin = as.numeric(allele) - 0.25,
                          ymax = as.numeric(allele) + 0.25,
                          fill=col),
              color='black',lwd=.1) + 
    scale_fill_manual(values=col_vector) + 
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=7),
          legend.position='none') + 
    xlab('') + ylab('') + 
    geom_text(data=dfAnnot, aes(x=(from+to)/2,label=domain),y=5.5,hjust=0.5,size=7*5/14) + 
    coord_cartesian(ylim=c(0.7,5.7),expand=FALSE,clip='off')

  return(gPRDM9)
}

###########################################
drawONTvPB <- function(){
  
  prdm9_Bonito <- read.table('prdm9_haplotypes.bonito.tab',header=TRUE,stringsAsFactors = FALSE)
  prdm9_Pacbio <- read.table('prdm9_haplotypes.pacbio.tab',header=TRUE,stringsAsFactors = FALSE)
  
  dfBon <- prdm9_Bonito %>% 
    filter(!grepl('(Unk|noZFA|NA)',perl=TRUE,diploid)) %>%
    mutate(homhet = ifelse(homhet == 'het' & size_1 == size_2,'Het\n(equal)',
                           ifelse(homhet == 'het' & size_1 != size_2, 'Het\n(unequal)', 'Hom'))) %>%
    select(id,bonGT = diploid, bonHH = homhet, bonSeq1 = seq_1, bonSeq2 = seq_2, bonN = nzfseqs) 
  
  dfPB <- prdm9_Pacbio %>% 
    filter(!grepl('(Unk|noZFA|NA)',perl=TRUE,diploid)) %>%
    mutate(homhet = ifelse(homhet == 'het' & size_1 == size_2,'Het\n(equal)',
                           ifelse(homhet == 'het' & size_1 != size_2, 'Het\n(unequal)', 'Hom'))) %>%
    select(id,pbGT = diploid, pbHH = homhet, pbSeq1 = seq_1, pbSeq2 = seq_2, pbN = nzfseqs) 
  
  df <- dfBon %>%
    inner_join(dfPB, by="id") %>%
    filter(!grepl('^A[AN]-',id)) %>%
    mutate(agree  = ifelse((bonSeq1 == pbSeq1 | bonSeq1 == pbSeq2) & 
                             (bonSeq2 == pbSeq1 | bonSeq2 == pbSeq2), TRUE, FALSE),
           homhet = ifelse(bonHH == pbHH, pbHH, paste0(bonHH, ":", pbHH))) %>% 
    group_by(pbHH) %>%
    add_count() %>% 
    select(homhet,agree,n) %>% 
    group_by(pbHH, agree, n) %>%
    count()%>%
    mutate(pc=nn/n*100) 
  
  dfOvFALSE <- data.frame(pbHH='Overall',agree=FALSE,n=sum(df$nn),nn=sum(df$nn[!df$agree]),pc=sum(df$nn[!df$agree])/sum(df$nn)*100) 
  dfOvTRUE <- data.frame(pbHH='Overall',agree=TRUE ,n=sum(df$nn),nn=sum(df$nn[df$agree]) ,pc=sum(df$nn[df$agree])/sum(df$nn)*100)

  df <- rbind(df,dfOvFALSE,dfOvTRUE)
  
  df$pbHH <- factor(df$pbHH, levels=rev(c('Overall','Hom','Het\n(unequal)','Het\n(equal)')))
  
  gPC <- ggplot(df[df$agree,],aes(y=pbHH,x=pc,fill=pbHH)) + geom_bar(width = .75, stat='identity') + 
    scale_fill_manual(values=c('grey50','grey50','grey50','firebrick')) + ylab('') + xlab('Concordant (%)') +
    theme(legend.position='none',
          axis.text.y = element_text(hjust=0.5)) + 
    geom_hline(yintercept=3.5,lty='dashed',lwd=.2) +
    coord_cartesian(xlim=c(0,100),ylim=c(0.5,4.5),expand=FALSE,clip='off') +
    geom_text(x=5,aes(label=paste0(round(pc,1),"%; N=",nn)),size=7*5/14,color='white',hjust=0)
  
  gN <- ggplot(df[df$agree,],aes(x=pbHH,y=nn,fill=agree)) + geom_bar(stat='identity') +  
    scale_fill_manual(values=c('grey50','firebrick')) + xlab('') + ylab('Concordant (#)') + 
    theme(legend.position='none')+ 
    geom_label(y=25,aes(label=paste0("N = ",nn)),size=8*5/14,fill='white')
  
  gPBvONT <- ggarrange(gN,gPC,ncol=1,labels =c('E',''),font.label = c(size=9))
  
  return(list(gPC = gPC,
              gN = gN,
              gX2 = gPBvONT))
}

## START ##################################
gMap <- drawMap()
gAlleles <- drawPRDM9alleles()
gONTvPB <- drawONTvPB()

gACD <- ggarrange(gMap,gAlleles,gONTvPB$gPC + theme(plot.margin=unit(c(0,0.2,0,0),'cm')),
                 ncol=3,nrow=1,
                 widths=c(2,1,1),
                 labels =c('A','C','D'),
                 font.label = c(size=9))

img <- readPNG("accessoryFiles/otherdata/Alleva_Schematic.png")
gB <- rasterGrob(img, interpolate=TRUE)

gAll <- ggarrange(gACD,gB,
                  ncol=1,nrow=2,
                  heights=c(1,1.8),
                  labels =c('','B'),
                  font.label = c(size=9))

ggsave('Alleva_et_al_Figure1.png',gAll,height=5.1, width=7.2,dpi = 500)
ggsave('Alleva_et_al_Figure1.pdf',gAll,height=5.7, width=7.2)