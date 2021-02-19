source('accessoryFiles/scripts/genericFunctions.R')

library(ggplot2)
library(motifStack)
library(Biostrings)
library('universalmotif')

## Functions -------------------------------------------------------------------
parsePWM <- function(pwmFile,sample=''){
  if (sample == 'A') {pwm <- read.table('selectPRDM9.AA.pwm',skip = 9,nrows = 4)}
  if (sample == 'B') {pwm <- read.table('selectPRDM9.AA.pwm',skip = 14,nrows = 4)}
  if (sample == 'C') {pwm <- read.table('selectPRDM9.AA.pwm',skip = 19,nrows = 4)}
  if (sample == 'L4'){pwm <- read.table('selectPRDM9.AA.pwm',skip = 24,nrows = 4)}
  if (sample == 'N') {pwm <- read.table('selectPRDM9.AA.pwm',skip = 29,nrows = 4)}
  
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

getMotifFromMEMEfile <- function(memeFile,rc=TRUE,n=1,name='motif'){

  meme    <- read_meme(memeFile,readsites = TRUE)
  if (rc){
    motif <- convert_motifs(motif_rc(meme$motifs[[n]]),'motifStack-pfm')
  }else{
    motif <- convert_motifs(meme$motifs[[n]],'motifStack-pfm')
  }
  
  #plot(motif)
  return(list(motif=motif,
              meme=meme,
              rc=rc,
              name=name,
              file=memeFile,
              len=nchar(summarise_motifs(motif)$consensus)))
}

getAAcodes <- function(aaFasta='selectPRDM9.AA.fa'){
  fa <- readAAStringSet(aaFasta)
  
  dfZF <- data.frame(nallele=0,allele=NA,n=rep(0,1000),seq=NA)
  cnt <- 0
  
  for (nSeq in 1:length(fa)){
    len = length(fa[[nSeq]])/28
    
    for (nzf in 1:(length(fa[[nSeq]])/28)){
      cnt <- cnt+1
      zf <- paste0(as.character(Views(fa[[nSeq]],start=((nzf-1)*28)+c(10,13,16),end=((nzf-1)*28)+c(10,13,16))),collapse = '')
      dfZF[cnt,] <- c(nSeq,names(fa)[nSeq],nzf,zf)
    }
  }
  
  dfZF         <- dfZF[!is.na(dfZF$allele),]
  dfZF$n       <- as.numeric(dfZF$n)
  dfZF$nallele <- as.numeric(dfZF$nallele)
  dfZF$from    <- (1 + (dfZF$n-1)*3) +.2
  dfZF$to      <- dfZF$from+2.6
  
  dfZF$xpos    <- dfZF$nallele
  
  return(dfZF)
}

drawMotifsFigure <- function(){
  ## AA codes --------------------------------------------------------------------
  dfZF     <- getAAcodes(aaFasta='selectPRDM9.AA.fa')
  
  dfZF$ypos[dfZF$allele == 'A']  <- 10
  dfZF$ypos[dfZF$allele == 'B']  <- 13
  dfZF$ypos[dfZF$allele == 'N']  <- 7
  dfZF$ypos[dfZF$allele == 'C']  <- 4
  dfZF$ypos[dfZF$allele == 'L4'] <- 1
  
  ## Predictions------------------------------------------------------------------
  lstPFMA  <- getPredMotif(parsePWM("selectPRDM9.AA.pwm",'A'),rc=TRUE,'Prdm9A')
  lstPFMB  <- getPredMotif(parsePWM("selectPRDM9.AA.pwm",'B'),rc=TRUE,'Prdm9B')
  lstPFMN  <- getPredMotif(parsePWM("selectPRDM9.AA.pwm",'N'),rc=TRUE,'Prdm9N')
  lstPFMC  <- getPredMotif(parsePWM("selectPRDM9.AA.pwm",'C'),rc=TRUE,'Prdm9C')
  lstPFML4 <- getPredMotif(parsePWM("selectPRDM9.AA.pwm",'L4'),rc=TRUE,'Prdm9L4')
  
  write_motifs(lstPFMA$motif,  'PRDM9A.predictedBindingSite.Alleva_et_al_2021.universalmotif.yaml',overwrite=TRUE)
  write_motifs(lstPFMB$motif,  'PRDM9B.predictedBindingSite.Alleva_et_al_2021.universalmotif.yaml',overwrite=TRUE)
  write_motifs(lstPFMC$motif,  'PRDM9C.predictedBindingSite.Alleva_et_al_2021.universalmotif.yaml',overwrite=TRUE)
  write_motifs(lstPFML4$motif, 'PRDM9L4.predictedBindingSite.Alleva_et_al_2021.universalmotif.yaml',overwrite=TRUE)
  write_motifs(lstPFMN$motif,  'PRDM9N.predictedBindingSite.Alleva_et_al_2021.universalmotif.yaml',overwrite=TRUE)
  
  # Motifs -----------------------------------------------------------------------
  lstMotA   <- getMotifFromMEMEfile("meme_details_A_Auto.top/meme_out/meme.txt",        n=1,rc=FALSE,name='Prdm9A')
  lstMotB   <- getMotifFromMEMEfile("meme_details_B_biased_Auto.top/meme_out/meme.txt", n=1,rc=TRUE,name='Prdm9B')
  lstMotC   <- getMotifFromMEMEfile("meme_details_C_Auto.top/meme_out/meme.txt",        n=1,rc=FALSE ,name='Prdm9C')
  lstMotL4  <- getMotifFromMEMEfile("meme_details_L4_biased_Auto.top/meme_out/meme.txt",n=1,rc=TRUE ,name='Prdm9L4')
  lstMotN   <- getMotifFromMEMEfile("meme_details_N_biased_Auto.top/meme_out/meme.txt", n=1,rc=TRUE,name='Prdm9N')

  write_motifs(lstMotA$motif,  'PRDM9A.topAutosomalHotspots.Alleva_et_al_2021.universalmotif.yaml',overwrite=TRUE)
  write_motifs(lstMotB$motif,  'PRDM9B.topAutosomalHotspots.Alleva_et_al_2021.universalmotif.yaml',overwrite=TRUE)
  write_motifs(lstMotC$motif,  'PRDM9C.topAutosomalHotspots.Alleva_et_al_2021.universalmotif.yaml',overwrite=TRUE)
  write_motifs(lstMotL4$motif, 'PRDM9L4.topAutosomalHotspots.Alleva_et_al_2021.universalmotif.yaml',overwrite=TRUE)
  write_motifs(lstMotN$motif,  'PRDM9N.topAutosomalHotspots.Alleva_et_al_2021.universalmotif.yaml',overwrite=TRUE)
  
  df <- data.frame(xmin=c(1, 36,
                          1, 19,  
                          1, 20, 
                          1, 20, 
                          1, 15 ), 
                   xmax=c(lstPFML4$len+1, 36+lstMotL4$len, 
                          lstPFMC$len+1 , 19+lstMotC$len,  
                          lstPFMN$len+1 , 20+lstMotN$len, 
                          lstPFMA$len+1 , 20+lstMotA$len,
                          lstPFMB$len+1 , 15+lstMotB$len), 
                   ymin=c(1,2.2, 4,5.2, 7,8.2, 10,11.2, 13,14.2), 
                   ymax=c(1,2.2, 4,5.2, 7,8.2, 10,11.2, 13,14.2)+1,
                   name=rep(c('L4','C','N','A','B'),each=2))
  
  df$mot <- list(lstPFML4$motif, lstMotL4$motif,
                 lstPFMC$motif, lstMotC$motif,
                 lstPFMN$motif, lstMotN$motif,
                 lstPFMA$motif, lstMotA$motif,
                 lstPFMB$motif, lstMotB$motif)
  
  dfPolyCL4 <- data.frame(x=c(19.00, 19.00, 31.00, 31.00, 55.00, 55.00, 43.00, 43.00),                     
                          y=c( 6.30,  3.40,  2.80,  0.40,  0.40,  2.80,  3.40,  6.30),
                          name='CL4',
                          col='royalblue2')
  
  dfPolyAB <- data.frame(x=c( 16.00, 16.00, 20.00, 20.00),                     
                         y=c( 9.40, 15.50, 15.50,  9.40),
                         name='AB',
                         col='forestgreen')
  
  dfPolyABN  <- data.frame(x=c(34.00, 34.00, 40.00, 40.00),                     
                           y=c( 6.40, 12.30, 12.30,  6.40),
                           name='ABN',
                           col='firebrick')
  
  dfPoly <- rbind(dfPolyCL4,
                  dfPolyABN,
                  dfPolyAB)
  
  dfLbl <- data.frame(x    = c(43.00, 43.0, 43.0),
                      xend = c(40.50, 40.5, 40.5),
                      y    = c(12.75, 13.4, 14.6),
                      lbl  = c('PRDM9 ZFs',
                               'Predicted\nbinding site',
                               'Hotspot motif'))
  
  dfArrows <- data.frame(x    = c(43.00, 43.0, 43.0),
                         xend = c(40.50, 40.5, 40.5),
                         y=c(12.75, 13.4, 14.6))
  
  gMotifs <- ggplot() + 
    geom_segment(data=df,aes(x=xmin,xend=xmax,y=ymin,yend=ymin),lwd=.2) + 
    scale_y_continuous(breaks=c(2,5,8,11,14),labels=c('L4','C','N','A','B')) + 
    geom_polygon(data=dfPoly,aes(x=x,y=y,group=name,color=col),
                 alpha=0.1,lty='dotted',lwd=.3, fill='grey50') + 
    geom_motif(data=df,aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, motif=mot)) + 
    geom_text(data=df[seq(2,10,2),],aes(y=ymin,x=1,label=paste0("PRDM9-",name)),size=7*5/14,hjust=0,vjust=0) +
    geom_text(data=dfLbl,aes(x=x,y=y,label=lbl),size=7*5/14,hjust=0,vjust=0.5) +
    geom_segment(data=dfLbl,aes(x=x-0.5,xend=xend,y=y,yend=y),lwd=.2,arrow = arrow(length = unit(0.1,"cm"))) + 
    xlab('Position (bp)') + ylab('') + 
    theme(legend.position='none',
          axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + 
    coord_cartesian(xlim=c(0.8,55),ylim=c(.2,15.6),expand=FALSE)
  
  dfL4 <- data.frame(y   = c( 0.40,  3.40,  0.40),
                     yend= c( 0.40,  3.40,  0.40),
                     x   = c( 7.20,  7.20, 22.20),
                     xend= c(18.80, 18.80, 33.80))
  
  gMotPlus <- gMotifs + 
    geom_rect(data=dfZF,
              aes(xmin=from,xmax=to,
                  ymin=ypos-.5,ymax=ypos-.1,fill=seq),
              color=NA,lwd=.2) + 
    geom_text(data=dfZF,aes(x=(from+to)/2,y=ypos-.3,label=seq),
              hjust=0.5,
              vjust=0.5,
              size=5*5/14) + 
    geom_segment(data=dfL4,aes(x=x,xend=xend,y=y,yend=yend),color='red',lwd=0.4)
  
  return(gMotPlus)
}

binom_P <- function(x, n, p) {
  P <- binom.test(x, n, p, alternative =  "two.sided", )$p.value
  return(-log10(P))
}

drawMAplot <- function(df,col1='AA1',col2='AN',colHighlight='magenta'){
  tot1 = sum(df[[col1]])
  tot2 = sum(df[[col2]])
  
  df$exp <- tot1/(tot1+tot2)

  df$col1 <- df[[col1]]
  df$col2 <- df[[col2]]
  
  dfMA <- df %>% rowwise() %>% mutate(pval=binom_P(col1,(col1+col2),exp))
  
  dfMA$pOK <- dfMA$pval > -log10(0.001/dim(dfMA)[1])
  
  ## Assure there's a similar enrichment in both AA samples
  dfMA$a2ok <- FALSE
  dfMA$a2ok[abs(dfMA$A-dfMA$A2) < 0.5] <- TRUE
  
  lblCount <- paste0('N = ',
                     format(sum(dfMA$pOK),big.mark=','),
                     ' / ',
                     format(dim(dfMA)[1],big.mark=','),
                     ' (',round(sum(dfMA$pOK)/dim(dfMA)[1]*100,0),'%)\n')
  
  gMA <- ggplot(dfMA,aes(x=M,y=A,color=pOK)) + 
    geom_point(size=.1,alpha=.4,fill='black',shape=21,lwd=.2) + 
    scale_x_log10() + 
    scale_color_manual(values=c('grey70',colHighlight)) + 
    geom_hline(yintercept=0,lwd=.2) + 
    theme(legend.position='none') + 
    xlab(paste0('Hotspot strength\nM = ',col1,' + ',col2)) + 
    ylab(paste0('Fold-change\nA = log2(',col1,'/',col2,')')) + 
    annotate(geom='label',
             x=Inf,y=-Inf,
             hjust=1,vjust=0,
             label=lblCount,
             color=colHighlight,
             size=7*5/14,
             label.size=NA) +
    annotation_logticks(sides='b',
                        long=unit(0.2,'cm'),
                        mid=unit(0.1,'cm'),
                        short=unit(0.1,'cm'),
                        size=.2)
  
  return(list(fig=gMA, data=dfMA))
}

toexpr <- function(l) {
  parse(text=l)
}

## Load data and plot ----------------------------------------------------------
theme7point()

df4 <- read.table('hotspotsData.tab',header=TRUE)

df4 <- df4[df4$cs %in% paste0('chr',1:21),]

alleles  <- c('AA1','AA2','AA3','AA4','AB1','AN','AC','CL4')
alleleNm <- c('A/A[1]','A/A[2]','A/A[3]','A/A[4]','A/B','A/N','A/C','C/L4')

mOverlap <- matrix(data=0,nrow=length(alleles),ncol=length(alleles))
colnames(mOverlap) <- alleleNm
rownames(mOverlap) <- alleleNm

mCC      <- mOverlap
plotList <- list()
for (ni in 1:length(alleles)){
  for (nj in 1:length(alleles)){
    
    i <- alleles[ni]
    j <- alleles[nj]
    
    overlap_PC <- sum(df4[[i]]>0 & df4[[j]]> 0) / min(sum(df4[[i]]>0),sum(df4[[j]]>0))
    mOverlap[ni,nj] <- overlap_PC
    
    dfCC <- df4[df4[[i]]>0 & df4[[j]]>0,]
    cc <- cor(log(dfCC[[i]]),log(dfCC[[j]]),method='pearson')
    mCC[ni,nj] <- cc
    
    plotList[[paste0(i,"v",j)]] <- plotHSstrength(dfCC[[i]],dfCC[[j]],i,j,topLeft = TRUE, pointSz = .1, xLimits = c(3,2000), yLimits = c(3,2000)) + 
      theme(legend.position='none')
    
    dfThis <- data.frame(gt1=i,gt2=j,pc=overlap_PC,cc=cc)
    
    if (ni*nj == 1){
      dfComp <- dfThis  
    }else{
      dfComp <- rbind(dfComp,dfThis)
    }
  }
}

## Overlap heatmap -------------------------------------------------------------
mOL       <- apply(mOverlap,2,function(x){mean(x[1:4][x[1:4]<1])})
dfOL      <- as.data.frame(mOL)
dfOL$name <- rownames(dfOL)
dfOL$name <- factor(dfOL$name,levels=alleleNm)

## CC heatmap ------------------------------------------------------------------
gT <- ggCorMat2(mCC,newOrd1 = alleleNm, newOrd2 = alleleNm, keepLeadingZeros = TRUE)
gO <- ggCorMat2(mOverlap,newOrd1 = alleleNm, newOrd2 = alleleNm, keepLeadingZeros = TRUE)

dfHMCC <- gT$data[!is.na(gT$data$value) & gT$data$value < 1,] %>% 
  mutate(x=Var1,y=Var2,val= value, lbl=round(value,2)) %>% 
  select(x,y,val,lbl)

dfHMOL <- gO$data[!is.na(gO$data$value) & gO$data$value < 1,] %>% 
  mutate(x=Var2,y=Var1,val=-value, lbl=paste0(round(value*100,0),'%')) %>% 
  select(x,y,val,lbl)

dfPlot <- rbind(dfHMOL,dfHMCC)

gAllvAll <- ggplot(dfPlot,aes(x=x,y=y,fill=val,label=lbl)) + geom_tile() + geom_text(size=7*5/14) + 
  scale_fill_gradient2(low='#56B4E9',mid='white',high='#D55E00',midpoint=0,na.value = NA) +
  xlab(bquote('Strength correlation at shared hotspots ('*R^2*')')) + ylab('Hotspot overlap (%)') + 
  theme(legend.position='none',
        axis.text.x=element_text(angle=45,hjust=1)) + 
  scale_x_discrete(labels=toexpr) +
  scale_y_discrete(labels=toexpr) 

# AA v AN ----------------------------------------------------------------------
dfAA1vAN <- df4 %>% filter(AA1>0 & AN > 0 & cs %in% paste0('chr',1:21)) %>% 
  mutate(tpmAA1 = toTPM(AA1), 
         tpmAA2 = toTPM(AA2), 
         tpmAN = toTPM(AN), 
         M = (toTPM(AA1) + toTPM(AN))/2, 
         A=log2(toTPM(AA1)/toTPM(AN)), 
         A2=log2(toTPM(AA2)/toTPM(AN)), 
         type='A/N') %>%
  select(cs,from,to,M,A,A2,type,tpmAA1,tpmAA2,tpmAN,AA1,AA2,AN) 

dfAA1vAA2 <- df4 %>% filter(AA1>0 & AA2 > 0 & cs %in% paste0('chr',1:21)) %>% 
  mutate(tpmAA1 = toTPM(AA1), 
         tpmAA2 = toTPM(AA2), 
         tpmAN = toTPM(AN), 
         M = (toTPM(AA1) + toTPM(AA2))/2, 
         A=log2(toTPM(AA1)/toTPM(AA2)), 
         A2=log2(toTPM(AA1)/toTPM(AA2)), 
         type='A/A (2)') %>%
  select(cs,from,to,M,A,A2,type,tpmAA1,tpmAA2,tpmAN,AA1,AA2,AN) 

lMA_AA1_AA2 <- drawMAplot(dfAA1vAA2,'AA1','AA2')
lMA_AA1_AN  <- drawMAplot(dfAA1vAN,'AA1','AN')

gMA2  <- ggarrange(lMA_AA1_AA2$fig + coord_cartesian(ylim=c(-6,6)),
                   lMA_AA1_AN$fig  + coord_cartesian(ylim=c(-6,6)),
                   ncol=1,
                   nrow=2,
                   widths=c(1,1),
                   align='v',
                   labels=c('C',''),
                   font.label = list(size=8,face='bold'))

# Motifs Plot ------------------------------------------------------------------
gMotifs <- drawMotifsFigure()

# Compound figure --------------------------------------------------------------
gX3  <- ggarrange(gAllvAll,
                  gMA2,
                  ncol=1,
                  nrow=2,
                  heights=c(1.1,2),
                  labels=c('B',''),
                  font.label = list(size=8,face='bold'))

ggarrange(gMotifs,gX3,
          ncol=2,
          nrow=1,
          widths=c(3,1.8),
          labels=c('A',''),
          font.label = list(size=8,face='bold'),
          hjust=0,vjust=1)

nScale <- .85
ggsave('Alleva_et_al_Figure4.png',width=8*nScale,height=7*nScale,dpi=400)
ggsave('Alleva_et_al_Figure4.pdf',width=8*nScale,height=7*nScale)
