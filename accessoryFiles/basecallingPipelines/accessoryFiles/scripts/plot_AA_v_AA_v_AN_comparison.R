source('genericFunctions.R')
library(ggplot2)

theme7point()
df4 <- read.table('hotspotsData.tab',header=TRUE)

df4 <- df4[df4$cs %in% paste0('chr',1:21),]
alleles  <- c('AA1','AA2','AN','CL4')
alleleNm <- c('A/A(1)','A/A(2)','A/N','C/L4')

mOverlap <- matrix(data=0,nrow=length(alleles),ncol=length(alleles))
colnames(mOverlap) <- alleleNm
rownames(mOverlap) <- alleleNm

mCC      <- mOverlap
plotList <- list()
for (ni in 1:4){
  for (nj in 1:4){
    
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

#gScatAA <- plotList$AA1vAA2 + coord_cartesian(xlim=c(1,10000),ylim=c(1,10000))
#gScatAN <- plotList$AA1vAN + coord_cartesian(xlim=c(1,10000),ylim=c(1,10000))

gCC <- ggCorMat(mCC,
                  noDiagonal = TRUE,
                  keepLeadingZeros = TRUE,
                  decimalPlaces = 2, 
                  flipIt = TRUE, 
                  yOnRight = TRUE,
                  xTilt=0) + 
  theme(axis.text.x=element_text(vjust=1)) +
  annotate(geom='text',x=0.5,y=3,hjust=0,size=8*5/14,
                           label="    Correlation (Pearson's R)") + 
  scale_fill_gradient(low=alpha('darkseagreen1',.3),high='darkolivegreen3',na.value = 'white')

gOL <- ggCorMat(mOverlap,
                  decimalPlaces = 0,
                  asPercentage = TRUE,
                  noDiagonal = TRUE,
                  flipIt = TRUE, 
                  yOnRight = TRUE,
                  xTilt=0) + 
  theme(axis.text.x=element_text(vjust=1)) +
  annotate(geom='text',x=0.5,y=3,hjust=0,size=8*5/14,
           label="    Shared hotspots (%)") + 
  scale_fill_gradient(low=alpha('pink',.3),high='firebrick',na.value = 'white')

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

dfs <- rbind(dfAA1vAN, dfAA1vAA2)
dfs$type <- factor(dfs$type,levels=c('A/A (2)','A/N'))

gDiff1 <- ggplot(dfs,aes(x=abs(A),fill=type,color=type)) + 
  geom_vline(xintercept=c(log2(2)),color='red',lwd=.2,lty='dashed') + 
  geom_density(alpha=.4) + 
  scale_color_manual(values=c('grey30','darkorange')) + 
  scale_fill_manual(values=c('grey50','orange')) + 
  xlab(bquote("Strength change: "*log[2]*(AA1/other))) + 
  ylab('Hotspot density') +
  coord_cartesian(xlim=c(0,5)) + 
  theme(legend.position=c(1,1),
        legend.justification=c(1,1),
        legend.key.size=unit(0.3,'cm'),
        legend.title=element_blank()) 

## Get 
dfHist <- dfs %>% mutate(rA = round(abs(A),2)) %>% 
  group_by(type,rA) %>% 
  count() %>% 
  group_by(type) %>% 
  add_tally(n,name='total') %>% 
  mutate(cN = cumsum(n)) %>% 
  mutate(cumPC = cN/total*100, pc = n/total*100)

gDiffA <- ggplot(dfHist,aes(x=rA,y=pc,color=type)) + 
  geom_vline(xintercept=c(log2(2)),color='red',lwd=.2,lty='dashed') + 
  geom_point(size=.1,alpha=.3) + 
  geom_smooth(lwd=.4,se = FALSE,span=.2) + 
  scale_color_manual(values=c('grey30','darkorange')) + 
  scale_fill_manual(values=c('grey50','orange')) + 
  xlab('') +
  ylab("Hotspots\n(%)") + 
  coord_cartesian(xlim=c(0,5)) + 
  theme(legend.position=c(1,1),
        legend.justification=c(1,1),
        legend.key.size=unit(0.3,'cm'),
        legend.title=element_blank(),
        axis.text.x = element_blank()) 

gDiffC <- ggplot(dfHist,aes(x=rA,y=cumPC,color=type)) + 
  geom_vline(xintercept=c(log2(2)),color='red',lwd=.2,lty='dashed') + 
  geom_line(size=.4) + 
  scale_color_manual(values=c('grey30','darkorange')) + 
  scale_fill_manual(values=c('grey50','orange')) + 
  xlab(bquote("Strength change: "*log[2]*(AA1/other))) + 
  ylab("Hotspots\n(cumulative %)") + 
  coord_cartesian(xlim=c(0,5)) + 
  theme(legend.position="none")

gDiff <- ggarrange(gDiff1 + xlab('') + theme(axis.text.x=element_blank()),
                   gDiffC,ncol=1,nrow=2,heights=c(5,5.5),align='v')

gScatters <- ggarrange(plotList$AA1vAA2 + scale_fill_gradient(low='grey50',high='black'),
                       plotList$AA1vAN  + scale_fill_gradient(low='grey50',high='black'),
                       gDiff,ncol=3,nrow=1)
gHM       <- ggarrange(gCC,gOL,ncol=2,nrow=1)

gX2 <- ggarrange(gScatters,gHM,heights=c(6,4),ncol=1,nrow=2)
ggsave(filename = 'AA_v_AN.png',gX2,width=6,height=4,dpi=400)

#-------------------------------------------------------------------------------------------
aa1tot = sum(dfAA1vAN$AA1)
antot = sum(dfAA1vAN$AN)
pE <- aa1tot/(aa1tot+antot)

dfAA1vAN$exp <- pE

binom_P <- function(x, n, p) {
  P <- binom.test(x, n, p, alternative =  "two.sided", )$p.value
  return(-log10(P))
}

dfMA <- dfAA1vAN %>% rowwise() %>% mutate(pval=binom_P(AA1,(AA1+AN),exp))

dfMA$pOK <- dfMA$pval > -log10(0.001/dim(dfMA)[1])

## Assure there's a similar enrichment in both AA samples
dfMA$a2ok <- FALSE
dfMA$a2ok[abs(dfMA$A-dfMA$A2) < 0.5] <- TRUE
#dfMA$a2ok[dfMA$A*dfMA$A_AA2 < 0]                        <- FALSE

gMA <- ggplot(dfMA,aes(x=M,y=A,color=pOK)) + 
  geom_point(size=.1,alpha=.4,fill='black',shape=21,lwd=.2) + 
  scale_x_log10() + 
  scale_color_manual(values=c('grey70','darkorange')) + 
  geom_hline(yintercept=0,lwd=.2) + 
  theme(legend.position='none') + 
  xlab('Hotspot strength\n(M = AA1 + AN)') + 
  ylab('Fold-change\n(A = log2(AA1/AN))') + 
  annotation_logticks(sides='b',
                      long=unit(0.2,'cm'),
                      mid=unit(0.1,'cm'),
                      short=unit(0.1,'cm'),
                      size=.2)

ggsave(filename = 'MAplot_AA_v_AN.png',gMA,width=3,height=3,dpi=400)

dfMA$set <- 'Others'
dfMA$set[dfMA$A < 0 & dfMA$pOK & dfMA$a2ok] <- 'N-up'
dfMA$set[dfMA$A > 0 & dfMA$pOK & dfMA$a2ok] <- 'N-down'
dfMA$set[(abs(dfMA$A) < log2(1.1)) & (!dfMA$pOK) & dfMA$a2ok] <- 'Stable'

gMASets <- ggplot(dfMA,aes(x=M,y=A,color=set)) + 
  geom_point(data=dfMA[dfMA$set == 'Others',],size=.1,alpha=.4,fill='black',shape=21,lwd=.2) + 
  geom_point(data=dfMA[dfMA$set != 'Others',],size=.1,alpha=.4,fill='black',shape=21,lwd=.2) +
  scale_x_log10() + 
  scale_color_manual('',values=c('green','red','grey80','dodgerblue2')) + 
  geom_hline(yintercept=0,lwd=.2) + 
  theme(legend.position=c(0,1),
        legend.justification=c(0,1),
        legend.background = element_blank(),
        legend.key.size = unit(c(0.2),'cm')) + 
  xlab('Hotspot strength\n(M = AA1 + AN)') + 
  ylab('Fold-change\n(A = log2(AA1/AN))') + 
  annotation_logticks(sides='b',
                      long=unit(0.2,'cm'),
                      mid=unit(0.1,'cm'),
                      short=unit(0.1,'cm'),
                      size=.2)

ggsave(filename = 'MAplot_AA_v_AN_UpDownInN.png',gMASets,width=3,height=3,dpi=400)
fwrite(x = dfMA, file = 'MAdetails_AA_v_AN.tab', sep = "\t")
#-------------------------------------------------------------------------------
gVolcano <- ggplot(dfMA,aes(x=A,y=pval,color=pOK)) + 
  geom_point(size=.1,alpha=.4,fill='black',shape=21,lwd=.2) + 
  scale_color_manual(values=c('grey70','darkorange')) + 
  theme(legend.position='none') + 
  ylab('-log10(P)') + 
  ylab('Fold-change\n(log2(AA1/AN))') + 
  coord_cartesian(ylim=c(0,50),xlim=c(-5,5))


ggsave(filename = 'Volcanoplot_AA_v_AN.png',gVolcano,width=3,height=3,dpi=400)

#-------------------------------------------------------------------------------------------
aa1tot = sum(dfAA1vAA2$AA1)
aa2tot = sum(dfAA1vAA2$AA2)
pE <- aa1tot/(aa1tot+aa2tot)

dfAA1vAA2$exp <- pE

binom_P <- function(x, n, p) {
  P <- binom.test(x, n, p, alternative =  "two.sided", )$p.value
  return(-log10(P))
}

dfMA <- dfAA1vAA2 %>% rowwise() %>% mutate(pval=binom_P(AA1,(AA1+AA2),exp))

dfMA$pOK <- dfMA$pval > -log10(0.001/dim(dfMA)[1])

gMA <- ggplot(dfMA,aes(x=M,y=A,color=pOK)) + 
  geom_point(size=.1,alpha=.4,fill='black',shape=21,lwd=.2) + 
  scale_x_log10() + 
  scale_color_manual(values=c('grey70','darkorange')) + 
  geom_hline(yintercept=0,lwd=.2) + 
  theme(legend.position='none') + 
  xlab('Hotspot strength\n(M = AA1 + AA2)') + 
  ylab('Fold-change\n(A = log2(AA1/AA2))') + 
  annotation_logticks(sides='b',
                      long=unit(0.2,'cm'),
                      mid=unit(0.1,'cm'),
                      short=unit(0.1,'cm'),
                      size=.2)

ggsave(filename = 'MAplot_AA1_v_AA2.png',gMA,width=3,height=3,dpi=400)


gVolcano <- ggplot(dfMA,aes(x=A,y=pval,color=pOK)) + 
  geom_point(size=.1,alpha=.4,fill='black',shape=21,lwd=.2) + 
  scale_color_manual(values=c('grey70','darkorange')) + 
  theme(legend.position='none') + 
  ylab('-log10(P)') + 
  ylab('Fold-change\n(log2(AA1/AN))') + 
  coord_cartesian(ylim=c(0,50),xlim=c(-5,5))

ggsave(filename = 'Volcanoplot_AA1_v_AA2.png',gVolcano,width=3,height=3,dpi=400)