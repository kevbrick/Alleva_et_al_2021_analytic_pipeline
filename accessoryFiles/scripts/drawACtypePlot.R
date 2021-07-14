source('accessoryFiles/scripts/genericFunctions.R')

library(ggplot2)
library(data.table)
library(plyr)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(png)
library(motifStack)
library(Biostrings)
library('universalmotif')

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

theme7point()

###############################
checkSpermVariantDets <- function(){
  dfAll <- fread('PrZFA_alleles.details.txt',skip=11,header=TRUE)
  
  dfBS <- fread('humanPRDM9alleles.BloodandSpermVariants.txt',header=FALSE) %>%
    dplyr:::rename(man=V1,
                   zf_code=V2,
                   ID=V3) %>%
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
  prdm9_allele_details <- read.table('atype_ctype.txt',header=TRUE,stringsAsFactors = FALSE,skip = 11) %>%
    mutate(allele=short_ID)
  
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
  
  gRet <- ggarrange(gTop + theme(legend.position='top',plot.margin=unit(c(0,0,0,0),'cm')),
                    gBott+theme(plot.margin=unit(c(0,0,0,0),'cm')),
                    ncol=1,nrow=2,align='v')
  
  return(gRet)
}

makeMotifsPlot <- function(){
  
  ## Functions -------------------------------------------------------------------
  parsePWM <- function(pwmFile = 'selectPRDM9.AA.pwm',reqID='') {
    ## Read in PWM IDs from pwm_predict file; skip 9 heasder rows
    namesTab <- read.table(pipe(paste0('grep [a-z] ',pwmFile)),skip = 6)
    pwmIDs   <- gsub('>','',namesTab$V1)
    
    nid <- which(pwmIDs == reqID)
    
    pwm <- read.table(pwmFile,skip = 9+((nid-1)*5),nrows=4)
    
    write.table(x = as.matrix(pwm), file='tmp.pfm',col.names = FALSE,row.names = FALSE,sep="\t",quote = FALSE)
    return('tmp.pfm')
  }
  
  getPredMotif <- function(pfmFile,rc=TRUE,name='motif'){
    pfm <- read.table(pfmFile)
    
    if(pfm$V1[1] %in% c('A','C','T','G')){
      rownames(pfm) <- pfm$V1
      pfm <- pfm[,2:dim(pfm)[2]]
    }
    
    if (rc){
      pfm <- rev(pfm)
      rownames(pfm) <- c('T','G','C','A')
    }else{
      rownames(pfm) <- c('A','C','G','T')
    }
    
    pfm <- as.data.frame(pfm %>% add_rownames() %>% arrange(rowname) %>% select(-rowname))
    rownames(pfm) <- c('A','C','G','T')
    
    motif <- new("pfm", mat=as.matrix(pfm), name=name)
    
    #plot(motif)
    return(list(pfm=pfm,
                motif=motif,
                rc=rc,
                name=name,
                file=pfmFile,
                len=dim(pfm)[2]))
  }
  ## END Functions -------------------------------------------------------------------
  
  dfDets   <- fread('PrZFA_alleles.details.txt',header=TRUE,skip = 11) 
  namesTab <- read.table(pipe(paste0('grep [a-z] ','PrZFA_Alleles.AA.pwm')),skip = 6)
  pwmIDs   <- gsub('>','',namesTab$V1)
  
  motifs <- list()
  lstAllMotifs <- list()
  nMotif <- 1
  for (m in pwmIDs){
    if (dfDets$in_pop[dfDets$ID == m] || dfDets$in_pop[dfDets$short_ID == m]){
      zz <- getPredMotif(parsePWM("PrZFA_Alleles.AA.pwm",m),rc=TRUE,m)
      m <- gsub("^(.+):[a-z]:([0-9]+):.+$","\\1:\\2",m,perl=TRUE)
      motifs[[m]]      <- zz$motif
      motifs[[m]]@name <- m
    }
    df <- data.frame(xmin=1,xmax=(dim(zz$pfm)[2])+1,ymin=0,ymax=1)
    
    df$mot <- list(zz$motif)
    
    gMotif <- ggplot(df) + 
      geom_vline(xintercept=seq(1,
                                ceiling(df$xmax/3)*3,3),
                 lwd=.2,color='grey80') +
      geom_motif(ic.scale=TRUE,
                 aes(xmin=xmin,
                     xmax=xmax,
                     ymin=ymin,
                     ymax=ymax,
                     motif=mot)) + 
      coord_cartesian(xlim=c(.8,df$xmax+.2),ylim=c(0,round(max(df$mot[[1]]$mat))),expand=FALSE) +
      scale_y_continuous(breaks=seq(0,round(max(df$mot[[1]]$mat)))) + 
      scale_x_continuous(breaks=seq(1,(max(df$xmax)-3),3)+1.5,
                         labels=1:(max(df$xmax)/3)) + 
      xlab(paste0('Zinc finger # (',m,")")) + 
      ylab("Score")

    lstAllMotifs[[nMotif]] <- gMotif
    ggsave(paste0('Alleva_et_al_Motifs_',m,'.png'),gMotif,height=.8,width=4)
    nMotif <- nMotif+1
  }
  
  hc <- clusterMotifs(motifs)
  phylog <- ade4::hclust2phylog(hc)
  
  pdf('motifAlignment.pdf',height=6,width=6)
  motifPiles(phylog,motifs,clabel.leaves=.5,r.tree=.1)
  dev.off()
}

################################################################################
dfAC <- fread('PrZFA_alleles.details.txt',header=TRUE,skip = 11) %>% 
  dplyr:::filter(in_pop == TRUE) %>%
  mutate(dA = ABDdist,
         dC = CBDdist) %>%
  group_by(dA,dC,ACtype) %>% 
  summarize(alleles = toString(short_ID),
            n=n()) %>%
  select(alleles,dA,dC,ACtype,n) %>%
  mutate(l1    = gsub("^(.{25}.*?),","\\1\n",alleles),
         l2    = gsub("^(.{52}.*?),","\\1\n",l1),
         l3    = gsub("^(.{79}.*?),","\\1\n",l2),
         label = gsub("^(.{111}.*?),","\\1\n",l3))

g <- ggplot(dfAC,aes(x=dA,y=dC)) + 
  geom_line(data=data.frame(dA=c(0,23),dC=c(0,23))) +
  geom_hline(lty='dotted',yintercept=0,lwd=.2,color='grey70') + 
  geom_vline(lty='dotted',xintercept=0,lwd=.2,color='grey70') + 
  geom_label_repel(aes(label=label,color=factor(ACtype)),size=7*5/14) +
  geom_point(aes(size=n,
                 fill=factor(ACtype)),
             shape=21,
             position=position_jitterdodge()) + 
  xlab(bquote("Distance to "*PRDM9[A]*" binding site")) + 
  ylab(bquote("Distance to "*PRDM9[C]*" binding site")) + 
  coord_cartesian(xlim=c(-0.6,20),ylim=c(-0.6,20),expand=FALSE) +
  theme(legend.position='none',
        legend.justification=c(1,.5),
        legend.title=element_blank()) + 
  scale_fill_manual(values=c('purple','chartreuse4')) +
  scale_color_manual(values=c('purple','chartreuse4')) +
  annotate(geom="text",x=15,y=16,label="A-type alleles",
           size=8*5/14,fontface='bold',angle=45,vjust=0,
           color='purple') + 
  annotate(geom="text",x=15,y=14,label="C-type alleles",
           size=8*5/14,fontface='bold',angle=45,vjust=1,
           color='chartreuse4') 

# ggsave('Prdm9_alleles_byACtype.png',plot = g,height=5,width=5,units = "in",dpi  = 400)

dfACall <- fread('PrZFA_alleles.details.txt',header=TRUE,skip = 11) %>% 
  mutate(dA = ABDdist,
         dC = CBDdist) %>%
  group_by(dA,dC,ACtype) %>% 
  summarize(alleles = toString(short_ID),
            n=n()) %>%
  select(alleles,dA,dC,ACtype,n) %>%
  mutate(l1    = gsub("^(.{25}.*?),","\\1\n",alleles),
         l2    = gsub("^(.{52}.*?),","\\1\n",l1),
         l3    = gsub("^(.{79}.*?),","\\1\n",l2),
         label = gsub("^(.{111}.*?),","\\1\n",l3))

gAll <- ggplot(dfACall,aes(x=dA,y=dC)) + 
  geom_line(data=data.frame(dA=c(0,23),dC=c(0,23))) +
  geom_hline(lty='dotted',yintercept=0,lwd=.2,color='grey70') + 
  geom_vline(lty='dotted',xintercept=0,lwd=.2,color='grey70') + 
  geom_point(aes(size=n,
                 fill=factor(ACtype)),
             shape=21,
             position=position_jitterdodge()) + 
  xlab(bquote("Distance to "*PRDM9[A]*" binding site")) + 
  ylab(bquote("Distance to "*PRDM9[C]*" binding site")) + 
  coord_cartesian(xlim=c(-0.6,20),ylim=c(-0.6,20),expand=FALSE) +
  theme(legend.position='none',
        legend.justification=c(1,.5),
        legend.title=element_blank()) + 
  scale_fill_manual(values=c('purple','chartreuse4')) +
  scale_color_manual(values=c('purple','chartreuse4')) +
  annotate(geom="text",x=15,y=16,label="A-type alleles",
           size=8*5/14,fontface='bold',angle=45,vjust=0,
           color='purple') + 
  annotate(geom="text",x=15,y=14,label="C-type alleles",
           size=8*5/14,fontface='bold',angle=45,vjust=1,
           color='chartreuse4') 

###############################################################################
gPop <- getPopAlleleSizeDistByType()
lstSperm <- checkSpermVariantDets()

gPopSperm <- ggarrange(gPop,
          lstSperm$gLen,
          lstSperm$gSc,
          ncol=3,nrow=1,
          labels=c('D','E','F'),
          font.label = list(size=8,face='bold'),
          hjust=0,vjust=1)

###############################################################################
gACtypes <- ggarrange(g,gAll,
                      ncol=1,nrow=2,
                      heights=c(1.4,1),
                      labels=c('B','C'),
                      font.label = list(size=8,face='bold'),
                      hjust=0,vjust=1)

##Commented out because we use a manually aligned version
##makeMotifsPlot()

myPNG <- readPNG('motifAlignmentManual.png')
gPic  <- rasterGrob(myPNG, interpolate=TRUE)

gTop <- ggarrange(gPic,gACtypes,
                  ncol=2,nrow=1,widths=c(1,.8),
                  labels=c('A',''),
                  font.label = list(size=8,face='bold'),
                  hjust=0,vjust=1)

gFinal <- ggarrange(gTop,gPopSperm,
                  ncol=1,nrow=2,
                  heights=c(3,.8))

ggsave('Alleva_et_al_ACtypePlot.png',
       plot = gFinal,
       bg='white',
       height=11,width=8)