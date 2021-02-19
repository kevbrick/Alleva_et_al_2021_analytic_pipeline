source('accessoryFiles/scripts/genericFunctions.R')
library(ggplot2)

binom_P <- function(x, n, p) {
  P <- binom.test(x, n, p, alternative =  "two.sided", )$p.value
  return(-log10(P))
}

drawMAplot <- function(dfInit,col1='AA1',col2='AN',name='A/N',colHighlight='magenta'){
  df <- dfInit[dfInit[[col1]] > 0 & dfInit[[col2]] > 0 &dfInit$cs %in% paste0('chr',1:21),]
  df$tpm1 <- toTPM(df[[col1]])
  df$tpm2 <- toTPM(df[[col2]])
  df$M    <- (df$tpm1 + df$tpm2)/2
  df$A    <- log2(df$tpm1/df$tpm2)
  df$type <- name
  
  tot1 = sum(df[[col1]])
  tot2 = sum(df[[col2]])
  
  df$exp <- tot1/(tot1+tot2)
  
  df$col1 <- df[[col1]]
  df$col2 <- df[[col2]]
  
  dfMA <- df %>% rowwise() %>% mutate(pval=binom_P(col1,(col1+col2),exp))
  
  dfMA$pOK <- dfMA$pval > -log10(0.001/dim(dfMA)[1])

  dfMA$a2ok <- FALSE
  
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

lstMAL4 <- drawMAplot(df4,col1='AC' ,col2='CL4',name='ACCL4')
lstMAB  <- drawMAplot(df4,col1='AA1',col2='AB1',name='AAAB')
lstMAN  <- drawMAplot(df4,col1='AA1',col2='AN' ,name='AAAN')

dfL4    <- lstMAL4$data[lstMAL4$data$pOK & lstMAL4$data$A < 0,] %>% mutate(strength=CL4) %>% select(cs,from,to,strength)
dfB     <- lstMAB$data [lstMAB$data$pOK  & lstMAB$data$A < 0,]  %>% mutate(strength=AB1) %>% select(cs,from,to,strength)
dfN     <- lstMAN$data [lstMAN$data$pOK  & lstMAN$data$A < 0,]  %>% mutate(strength=AN)  %>% select(cs,from,to,strength)

write.table(dfL4,'L4_biased_HS.bedgraph',quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(dfB,'B_biased_HS.bedgraph',quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(dfN,'N_biased_HS.bedgraph',quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)