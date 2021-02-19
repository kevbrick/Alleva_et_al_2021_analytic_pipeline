library(ggplot2)
library(tidyverse)
library(reshape2)

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

prdm9_DF$ok <- TRUE
prdm9_DF$ok[prdm9_DF$diploid == '/'] <- FALSE
prdm9_DF$ok[grepl('Unk|noZF',prdm9_DF$diploid)] <- FALSE

prdm9_DF$typeHH  <- 'het'
prdm9_DF$typeHH[prdm9_DF$size_1 == prdm9_DF$size_2 & prdm9_DF$allele_1 == prdm9_DF$allele_2]  <- 'hom'
prdm9_DF$typeHH <- factor(prdm9_DF$type,levels=c('het','hom'))

prdm9_DF$type  <- 'unk'
prdm9_DF$type[prdm9_DF$size_1 != prdm9_DF$size_2]  <- 'easy-het'
prdm9_DF$type[prdm9_DF$size_1 == prdm9_DF$size_2 & prdm9_DF$allele_1 != prdm9_DF$allele_2]  <- 'tough-het'
prdm9_DF$type[prdm9_DF$size_1 == prdm9_DF$size_2 & prdm9_DF$allele_1 == prdm9_DF$allele_2]  <- 'hom'

prdm9_DF$type <- factor(prdm9_DF$type,levels=c('unk','tough-het','easy-het','hom'))

## Start making plots
plotDF <- prdm9_DF %>% filter(ok & pop != 'OTH') %>%
                       select(pop,diploid,allele_1,allele_2,type,typeHH) 

popsDF <- data.frame(pop = c('CHB','FIN','LWK','PEL','PJL','TSI','YRI'),
                     tot = c(120,100,120,70,108,114,120))

###########################################################################
plot1DF <- melt(plotDF[,c('diploid','pop')] %>% group_by(pop) %>% count() %>% 
                  inner_join(popsDF) %>% mutate(yes=n,no=tot-n) %>% select(pop,yes,no))

ggplot(plot1DF,aes(y=pop,x=value,fill=variable)) + 
  geom_bar(stat='identity',width=.75) + 
  geom_text(aes(label=value),
            color='white',
            hjust=2,size=8*5/14) +
  scale_fill_manual('PRDM9 diploid genotype?',
                    values=c('grey40','darkorange')) + 
  xlab('Individuals (#)') + 
  ylab('') + 
  theme(legend.position='top')

ggsave(paste0('PRDM9_diploid_GT_Counts.png'),width=4,height=3)
ggsave(paste0('PRDM9_diploid_GT_Counts.pdf'),width=4,height=3)

###########################################################################

plot2DF <- as.data.frame(plotDF[,c('diploid','pop')] %>% 
                      group_by(pop,diploid) %>% 
                      count() %>% 
                      group_by(pop) %>% 
                      add_count() %>% 
                      group_by(pop,diploid) %>% 
                      summarize(pc=n/nn*100,
                                pcRnd=round(n/nn*100/10)*10))

ggplot(plot2DF,aes(y=pop,x=diploid,fill=pcRnd)) + 
  geom_tile(lwd=.2,color='black') + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5)) + 
  geom_text(aes(label=round(pc,0)),
            size=8*5/14) +
  scale_fill_gradient('Freq(%)',high = 'firebrick1', low='white', na.value = 'white') + 
  xlab('') + ylab('') +
  ggtitle(paste0('Genotype frequency (%); ',dim(plotDF)[1]," genotypes inferred"))

ggsave(paste0('PRDM9_diploid_GT_by_Pop_Freqs.png'),width=14,height=4)
ggsave(paste0('PRDM9_diploid_GT_by_Pop_Freqs.pdf'),width=14,height=4)
###########################################################################
plot3DF <- as.data.frame(plotDF[,c('diploid','pop')] %>% 
                           group_by(pop,diploid) %>% 
                           count() %>% 
                           group_by(pop) %>% 
                           add_count() %>% 
                           group_by(pop,diploid))
                         
ggplot(plot3DF,aes(y=pop,x=diploid,fill=n)) + 
  geom_tile(lwd=.2,color='black') + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5)) + 
  geom_text(aes(label=round(n,0)),
            size=8*5/14) +
  scale_fill_gradient('Freq(#)',high = 'forestgreen', low='white', na.value = 'white') + 
  xlab('') + ylab('') +
  ggtitle(paste0('Genotypes (N); ',dim(plotDF)[1]," genotypes inferred"))

ggsave(paste0('PRDM9_diploid_GT_by_Pop_Counts.png'),width=14,height=4)
ggsave(paste0('PRDM9_diploid_GT_by_Pop_Counts.pdf'),width=14,height=4)
#######################
plot4DF <- melt(plotDF[,c('pop','allele_1','allele_2')],id.vars='pop',value.name='allele') %>%
  group_by(pop,allele) %>% 
  count() %>% 
  group_by(pop) %>% 
  add_count(name='tot') %>% 
  group_by(pop,allele) %>% 
  left_join(plotDF %>% count(pop, name = "totPop")) %>%
  summarize(N=n,
            pc=n/tot*100,
            pcRnd=round(n/tot*100/10)*10,
            popPC=n/totPop*100)

ggplot(plot4DF,aes(y=pop,x=allele,fill=pc)) + 
  geom_tile(lwd=.2,color='black') + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5)) + 
  geom_text(aes(label=round(pc,0)),
            size=8*5/14) +
  scale_fill_gradient('%',high = 'firebrick1', low='white', na.value = 'white') + 
  xlab('') + ylab('') + 
  ggtitle('Allele frequency (%)')

ggsave(paste0('PRDM9_ZFallele_Freq.png'),width=14,height=4)
ggsave(paste0('PRDM9_ZFallele_Freq.pdf'),width=14,height=4)
#######################
ggplot(plot4DF,aes(y=pop,x=allele,fill=N)) + 
  geom_tile(lwd=.2,color='black') + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5)) + 
  geom_text(aes(label=round(N,0)),
            size=8*5/14) +
  scale_fill_gradient('#',high = 'forestgreen', low='white', na.value = 'white') + 
  xlab('') + ylab('') + 
  ggtitle('Allele counts (#))')

ggsave(paste0('PRDM9_ZFallele_Freq_Counts.png'),width=14,height=4)
ggsave(paste0('PRDM9_ZFallele_Freq_Counts.pdf'),width=14,height=4)
###########################
dfA1  <- plotDF[,c('allele_1','pop')]; names(dfA1) <- c('allele','pop')
dfA2  <- plotDF[,c('allele_2','pop')]; names(dfA2) <- c('allele','pop')
dfA3 <- dfA2[dfA2$allele != dfA1$allele,]
plot6DF <- rbind(dfA1,dfA3)  %>% 
  group_by(pop,allele) %>% add_count() %>% inner_join(popsDF) %>%
  summarise(N=n, 
            popPC = n/tot*100)
           
ggplot(plot6DF,aes(y=pop,x=allele,fill=popPC)) + 
  geom_tile(lwd=.2,color='black') + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5)) + 
  geom_text(aes(label=round(popPC,0)),
            size=8*5/14,check_overlap = TRUE) +
  scale_fill_gradient(high = 'firebrick1', low='white', na.value = 'white') + 
  xlab('') + ylab('') + 
  ggtitle('Percentage of population with at least one copy of an allele')

ggsave(paste0('PRDM9_ZFalleles_by_Pop_Freqs.png'),width=14,height=4)
ggsave(paste0('PRDM9_ZFalleles_by_Pop_Freqs.pdf'),width=14,height=4)
###########################################################################

ggplot(plot6DF,aes(y=pop,x=allele,fill=N)) + 
  geom_tile(lwd=.2,color='black') + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5)) + 
  geom_text(aes(label=round(N,0)),
            size=8*5/14,check_overlap = TRUE) +
  scale_fill_gradient('N',high = 'forestgreen', low='white', na.value = 'white') + 
  xlab('') + ylab('') + 
  ggtitle('Number of individuals with at least one copy of an allele')

ggsave(paste0('PRDM9_ZFalleles_by_Pop_Freqs_Count.png'),width=14,height=4)
ggsave(paste0('PRDM9_ZFalleles_by_Pop_Freqs_Count.pdf'),width=14,height=4)

##########################################################################

plot7DF <- plotDF %>%  select(pop,type) %>% group_by(pop,type) %>% count() %>%  
  inner_join(plotDF %>% group_by(pop) %>% count(name = "tot")) %>% 
  summarise(N=tot, pc=n/tot*100)

ggplot(plot7DF,aes(x=pc,y=pop,fill=type)) + geom_bar(stat='identity',width=.8) + 
  scale_fill_manual('PRDM9 genotype',values=c('firebrick3','pink','red','grey40')) + 
  xlab('Individuals (%)') + ylab('') + 
  theme(legend.position='top')

ggsave(paste0('PRDM9_homhetGT_By_Pop.png'),width=4,height=3)
ggsave(paste0('PRDM9_homhetGT_By_Pop.pdf'),width=4,height=3)

##########################################################################

plot8DF <- plotDF %>%  select(pop,typeHH) %>% group_by(pop,typeHH) %>% count() %>%  
  inner_join(plotDF %>% group_by(pop) %>% count(name = "tot")) %>% 
  summarise(N=tot, pc=n/tot*100)

ggplot(plot8DF,aes(x=pc,y=pop,fill=typeHH)) + geom_bar(stat='identity',width=.8) + 
  scale_fill_manual('PRDM9 genotype',values=c('firebrick4','grey40')) + 
  xlab('Individuals (%)') + ylab('') + 
  theme(legend.position='top')

ggsave(paste0('PRDM9_homhetGT_Simple_By_Pop.png'),width=4,height=3)
ggsave(paste0('PRDM9_homhetGT_Simple_By_Pop.pdf'),width=4,height=3)

##########################################################################
#### Missing barcodes
#ggplot(prdm9_DF %>% group_by(bc1,bc2) %>% summarize(N=sum(okpb)) %>% mutate(gt=ifelse(N>0,"YES","NO")), aes(x=bc2,y=bc1)) + geom_tile(aes(fill=factor(gt)),color='white',lwd=.4) + theme(panel.grid =element_line(),axis.text.x=element_text(angle = 90), legend.position='top') + scale_fill_manual('Diploid GT (Pacbio)?',values=c('grey90','firebrick'))+ xlab('') +ylab('') 
#ggsave('missing_individuals.png',height=10, width=2.5,dpi=400)
