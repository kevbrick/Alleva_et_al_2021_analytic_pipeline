library(imager)
library(ggplot2)
library(ggpubr)
library(grid)
library(reshape2)
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

#############################
### fancy_scientific
### KB July 10 2018
### Change number to scientific notation format:
### Most useful for scale_x[y]_log10(labels=fancy_scientific)
### ARGS:
# l   number
## RETURNS: Formatted number as expression
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "e", l)
  # remove +s
  l <- gsub("\\+", "", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "10^", l)
  # return this as an expression
  parse(text=l)
}

theme7point()
quantifyGel <- function(img, cropX, cropY, nLanes, names, ladderLanes = 1, showLadderQuant = TRUE){

  imgRaw <- load.image(img)
  
  ## Cropping distances must be specified a-priori
  ## No autodetection ...
  imgGel <- crop.borders(imgRaw, nx = cropX, ny = cropY, nz = 0)
  
  ## Get the width of the cropped image
  nImageWidth <- dim(imgGel)[2]
  
  nFontNormalizer <- (dim(imgGel)[1])/350 
  
  ## Make Labels
  laneOrder <- vector(length = nLanes)
  laneNums  <- vector(length = nLanes)
  laneNames <- vector(length = nLanes)
  
  nLadders  <- 0
  
  for (n in 1:nLanes){
    if (n %in% ladderLanes){
      nLadders <- nLadders+1
      if (length(ladderLanes) == 1){
        laneOrder[n] <- 'Ladder'
        laneNums[n]  <- 'Ladder'
        laneNames[n] <- 'Ladder'
      }else{
        laneOrder[n] <- paste0('Ladder ',nLadders)
        laneNums[n]  <- paste0('Ladder ',nLadders)
        laneNames[n] <- paste0('Ladder ',nLadders)
      }
    }else{
      if (length(names) > 0){
        laneOrder[n] <- paste0("Lane ",n-nLadders,"\n",names[n])
        laneNums[n]  <- paste0("Lane ",n-nLadders)
        laneNames[n] <- names[n]
      }else{
        laneOrder[n] <- paste0("Lane ",n-nLadders)
        laneNums[n]  <- paste0("Lane ",n-nLadders)
        laneNames[n] <- paste0("Lane ",n-nLadders)
      }
    }
  } 
  
  ## Make DF to with empty intensities 
  dfPlot <- data.frame(lane=rep(laneOrder,each=nImageWidth),
                       laneName=rep(laneNames,each=nImageWidth),
                       pos=1:nImageWidth,
                       intensity=vector(length = nLanes*nImageWidth))
  
  dfPlot$lane <- factor(dfPlot$lane, levels=laneOrder)
  
  ## make Gel Image DF
  dfGels <- data.frame()
  ## Get lane widths
  nLaneWidth <- dim(imgGel)[1]/nLanes
  nW <- round(nLaneWidth)/2
  
  ## Initialize empty list for graphs
  ggList      <- list()
  ggListNoLbl <- list()
  
  ## Loop through lanes to caluclate intensity for each
  for (p in 1:nLanes) {
    
    ## Get x-coordinate position on image for this lane
    xPos <- (p-1)*(nW*2) + nW
  
    ## Get positions for this lane in the dataframe
    from <- (p-1)*nImageWidth+1
    to   <- p*nImageWidth

    ## Get pixel values for this lane
    ## R+G+B in a 10-pixel-wide region around the center
    i <- at(R(imgGel),(xPos-5):(xPos+5)) + 
      at(G(imgGel),(xPos-5):(xPos+5)) + 
      at(B(imgGel),(xPos-5):(xPos+5))
    
    ## Intensity is mean of 11 pixels at each point
    dfPlot$intensity[from:to] <- colMeans(i)
    
    dfGel <- reshape2::melt(i)
    names(dfGel)=c('xpixel','ypixel','intensity')
    dfGel$id <- laneNames[p]
    
    dfGels <- rbind(dfGels,dfGel)
    
    ## Draw the image of the gel segment used for quantitation
    gGel <- ggplot(melt(i),
                   aes(y=Var1,x=Var2,fill=value)) + 
      geom_tile() + 
      scale_fill_gradient(low='black',high='white') + 
      coord_cartesian(xlim=c(1,to-from),expand=FALSE) + 
      theme(legend.position='none',
            axis.line=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank(),
            plot.margin = unit(c(0,0,0,0.1),'cm')) + 
      ylab('') + 
      xlab('')
    
    ## Draw the quantitation graph
    ## NOTE: High intensity == 0 so y-normalization is to assure high 
    ##       intensity is "Up" on the y-axis
    gQuant <- ggplot(dfPlot[from:to,],
                     aes(x=pos,
                         y=-1*intensity + max(intensity))) + 
      geom_line(lwd=.3) + 
      coord_cartesian(xlim=c(1,to-from),expand=FALSE) + 
      xlab('Run time (pixels)') + ylab('Mean intensity (A.U.)') + 
      theme(plot.margin = unit(c(0,0,0,0.1),'cm')) + 
      geom_text(x=Inf,y=Inf,
                label=dfPlot$lane[from],
                size=8*5/14,
                check_overlap = TRUE,
                hjust=1,vjust=1,
                color='red')
    
    ggList[[paste0("Lane",p)]] <- ggarrange(gGel,gQuant,
                                            ncol=1,nrow=2,
                                            heights=c(1,3),
                                            align='v')

    ggListNoLbl[[paste0("Lane",p)]] <- ggarrange(gGel,
                                                 ggplot() + theme_void(), 
                                            gQuant + xlab('') + ylab(''),
                                            ncol=1,nrow=3,
                                            heights=c(10,1,30),
                                            align='v')
    
    ## Add lane number unless this lane is a ladder
    if (grepl(pattern = 'Lane', x = dfPlot$lane[from])){
      imgGel <- imager:::draw_text(im = imgGel, 
                                   x = xPos-4, 
                                   y = 5, 
                                   text = paste0(p-sum(ladderLanes < p)),
                                   color = 'red', 
                                   opacity = 1, 
                                   fsize = 10*nFontNormalizer)  
    }
  }
  
  ## Make merged plots with all subplots
  gAllPlots <- annotate_figure(ggarrange(plotlist = ggListNoLbl),
                              bottom=text_grob(label='Run time (pixels)',
                                               size=8),
                              left=text_grob(label='Mean RGB intensity (A.U.)',
                                             size=8,
                                             rot=90))
  
  ## Convert gel image to geom
  gFullGel <- rasterGrob(imgGel)
  
  ## Full figure: Gel & top - Quant on bottom
  gBigFig <- ggarrange(gFullGel,gAllPlots,
                       heights=c(1,2),
                       ncol=1,nrow=2)
  
  ## Return gel image and list of graphs per lane
  return (list(img  = imgGel, 
               gGel = gFullGel,
               gAll = gAllPlots,
               gBig = gBigFig,
               dets = ggList,
               data = dfPlot,
               gelData=dfGels))
}

### Plot gel vs reads vs ZFs for select samples
plotGelVReadsVZFs <- function(dfIntensity, dfGel, ids){
  dfZF       <- fread('allZFAsizes.txt',header=TRUE)
  dfReadLens <- fread('allReadLengths.txt',header=TRUE)
  
  dfReadLens$kb <- as.numeric(dfReadLens$length)/1000
  
  ggList <- list()
  
  for (thisID in ids){
    if (!grepl('Ladder',thisID)){
      noMarg <- theme(plot.margin=unit(c(0,0,0.5,0.5),'cm'))
      
      gI <- ggplot(dfIntensity[dfIntensity$laneName == thisID,],
                   aes(x=pos,
                       y=-1*intensity + max(intensity))) + 
        geom_line() +
        xlab('') + ylab('Intensity (A.U.)') + 
        geom_text(x=Inf,y=Inf,
                  hjust=1,vjust=1,
                  aes(label=lane),
                  size=8*5/14,
                  color='red',
                  check_overlap = TRUE) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              plot.margin=unit(c(0,0,0,0.5),'cm')) + 
        coord_cartesian(expand=FALSE,clip='off')
      
      gGel <- ggplot(dfGel[dfGel$id == thisID,],
                     aes(x=ypixel,y=xpixel,fill=intensity)) + 
        geom_tile() +
        scale_fill_gradient(low='black',high='grey99') +
        xlab('Run time (pixels)') + ylab('') +
        theme(legend.position='none',
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank(),
              plot.margin=unit(c(0,0,0,0.5),'cm')) + 
        coord_cartesian(expand=FALSE,clip='off')
      
      if (thisID == 'PrZFA A/N ctrl'){
        thisID <- 'AN-472'
      }
      
      if (thisID == 'PrZFA A/A ctrl'){
        thisID <- 'AA-428'
      }
      
      if (thisID == 'LCL gDNA'){
        thisID <- 'LCL'
      }
      
      gR <- ggplot(dfReadLens[dfReadLens$id == thisID,],
                   aes(x=length/1000)) + 
        geom_histogram(binwidth = .042) + 
        xlab('Read length (Kbp)') + ylab('Count') + 
        geom_text(x=-Inf,y=Inf,
                  hjust=0,vjust=1,
                  label=paste0("  # reads = ",dim(dfReadLens[dfReadLens$id == thisID,])[1]),
                  size=7*5/14,
                  color='forestgreen',
                  check_overlap = TRUE) + 
        theme(plot.margin=unit(c(0,0,0,0.5),'cm')) + 
        coord_cartesian(xlim=c(0,2.8),expand=FALSE,clip='off')
        
      
      gZ <- ggplot(dfZF[dfZF$id == thisID,],
                   aes(x=size)) + 
        geom_bar() +
        xlab('ZF array length (#ZFs)') + ylab('Count') +
        scale_x_continuous(breaks = seq(0,22,1),
                           labels = c(0,rep('',4),
                                      5,rep('',4),
                                      10,rep('',4),
                                      15,rep('',4),
                                      20,rep('',2))) +
        geom_text(x=-Inf,y=Inf,
                  hjust=0,vjust=1,
                  label=paste0("  # ZFAs = ",dim(dfZF[dfZF$id == thisID,])[1]),
                  size=7*5/14,
                  color='forestgreen',
                  check_overlap = TRUE) + theme(plot.margin=unit(c(0,0,0.8,0.5),'cm')) + 
        coord_cartesian(xlim=c(0,22),expand=FALSE,clip='off')

      
      ggList[[thisID]] <- ggarrange(gI,gGel,gR,gZ,
                                ncol=1,nrow=4,
                                heights=c(1,.6,1,1.4),
                                align='v')
    }
  }
  
  ## Make merged plots with all subplots
  gAllPlots <- ggarrange(plotlist = ggList)
  
  nHalf <- floor(length(ggList)/2)
  nTot <- length(ggList)
  
  gPlotsA <- ggarrange(plotlist = ggList[1:nHalf],ncol=4,nrow = ceiling(nHalf/4))
  gPlotsB <- ggarrange(plotlist = ggList[(nHalf+1):nTot],ncol=4, nrow = ceiling((nTot-nHalf)/4))
  
  return(list(plotsList = ggList,
              gAll = gAllPlots,
              gSetA = gPlotsA,
              gSetB = gPlotsB))
}

## IMG 1 ####################################################################
gQuantImg1 <- quantifyGel(img = "SuppGelFigure_justGel.jpg",
                          cropX = 6, 
                          cropY = 40, 
                          nLanes = 17,
                          names = c('Ladder',
                                    'NA18877',
                                    'NA18881',
                                    'NA18908',
                                    'NA18909',
                                    'NA18910',
                                    'NA18912',
                                    'NA18913',
                                    'NA18915',
                                    'NA18923',
                                    'NA18933',
                                    'NA18934',
                                    'NA19092',
                                    'NA19093',
                                    'NA19095',
                                    'NA19096',
                                    'NA19099'),
                          ladderLanes = 1)

### old order
# 'NA18871',
# 'NA18923',
# 'NA18873',
# 'NA18924', 
# 'NA18877', 
# 'NA18933',
# 'NA18881', 
# 'NA18934', 
# 'NA18909', 
# 'NA19092',
# 'NA18915', 
# 'NA19095',
# 'NA18916',
# 'NA19096', 
# 'NA18917', 
# 'NA19099'),

gFigI1 <- ggarrange(gQuantImg1$gGel,
                    ggplot() + theme_void(),
                    gQuantImg1$gAll,
          labels=c('A',' ','B'),
          ncol=1,nrow=3,
          heights=c(20,1,40),
          font.label = list(size=8,face='bold'))

ggsave('Alleva_et_al_PrZFA_Img1Quantified_JustGelData.png',gFigI1,height=8,width=8,dpi=500)
ggsave('Alleva_et_al_PrZFA_Img1Quantified_JustGelData.pdf',gFigI1,height=8,width=8)

## IMG 2 ####################################################################
gQuantImg2 <- quantifyGel(img = "SuppGelFigure_justGel2.jpg",
                          cropX = 2,
                          cropY = 2,
                          nLanes = 17,
                          names = c('Ladder',
                                    'NA20818',
                                    'NA20831',
                                    'NA20819', 
                                    'NA20832',
                                    'NA20820', 
                                    'PrZFA A/A ctrl',
                                    'NA20821', 
                                    'PrZFA A/N ctrl',
                                    'NA20822', 
                                    'NA08873', 
                                    'NA20826', 
                                    'NA17221', 
                                    'NA20827', 
                                    'NA17236', 
                                    'NA20828', 
                                    'LCL gDNA'),
                          ladderLanes = 1)

gFigI2 <-  ggarrange(gQuantImg2$gGel,
                     ggplot() + theme_void(),
                     gQuantImg2$gAll,
                     labels=c('C',' ','D'),
                     ncol=1,nrow=3,
                     heights=c(20,1,40),
                     font.label = list(size=8,face='bold'))

ggsave('Alleva_et_al_PrZFA_Img2Quantified_JustGelData.png',gFigI2,height=8,width=8,dpi=500)
ggsave('Alleva_et_al_PrZFA_Img2Quantified_JustGelData.pdf',gFigI2,height=8,width=8)

## ADD OTHER DATA
gCompleteGels <- ggarrange(gQuantImg1$gGel,gQuantImg2$gGel,ncol=1,nrow=2,labels=c('A','B'),font.label=list(size=8,face='bold'))

idList <- unique(gQuantImg1$data$laneName)
lstQuants1 <- plotGelVReadsVZFs(gQuantImg1$data,gQuantImg1$gelData,idList)

idList <- unique(gQuantImg2$data$laneName)
lstQuants2 <- plotGelVReadsVZFs(gQuantImg2$data,gQuantImg2$gelData,idList)

## Make PNGs ##############################################################
nHeight <- 9.5
ggsave('Alleva_et_al_PrZFA_Gel1.png',
       ggarrange(gQuantImg1$gGel,
                 ggplot() + theme_void(),
                 lstQuants1$gSetA,
                 labels=c("A",'','B'),
                 ncol=1,nrow=3,
                 heights=c(1,0.2,4),
                 font.label = list(size=8,face='bold')),
       height=nHeight,width=7,dpi=500)

ggsave('Alleva_et_al_PrZFA_Gel1_Part2.png',
       lstQuants1$gSetB,
       height=nHeight*4/5.2,width=7,dpi=500)

ggsave('Alleva_et_al_PrZFA_Gel2.png',
       ggarrange(gQuantImg2$gGel,
                 ggplot() + theme_void(),
                 lstQuants2$gSetA,
                 labels=c("C",'','D'),
                 ncol=1,nrow=3,
                 heights=c(1,0.2,4),
                 font.label = list(size=8,face='bold')),
       height=nHeight,width=7,dpi=500)

ggsave('Alleva_et_al_PrZFA_Gel2_Part2.png',
       lstQuants2$gSetB,
       height=nHeight*4/5.2,width=7,dpi=500)

## PDFs TOO #####################################################
ggsave('Alleva_et_al_PrZFA_Gel1.pdf',
       ggarrange(gQuantImg1$gGel,
                 ggplot() + theme_void(),
                 lstQuants1$gSetA,
                 labels=c("A",'','B'),
                 ncol=1,nrow=3,
                 heights=c(1,0.2,4),
                 font.label = list(size=8,face='bold')),
       height=nHeight,width=7)

ggsave('Alleva_et_al_PrZFA_Gel1_Part2.pdf',
       lstQuants1$gSetB,
       height=nHeight*4/5.2,width=7)

ggsave('Alleva_et_al_PrZFA_Gel2.pdf',
       ggarrange(gQuantImg2$gGel,
                 ggplot() + theme_void(),
                 lstQuants2$gSetA,
                 labels=c("C",'','D'),
                 ncol=1,nrow=3,
                 heights=c(1,0.2,4),
                 font.label = list(size=8,face='bold')),
       height=nHeight,width=7)

ggsave('Alleva_et_al_PrZFA_Gel2_Part2.pdf',
       lstQuants2$gSetB,
       height=nHeight*4/5.2,width=7)