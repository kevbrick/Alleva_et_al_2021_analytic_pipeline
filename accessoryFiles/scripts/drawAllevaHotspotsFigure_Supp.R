  source('accessoryFiles/scripts/genericFunctions.R')
  
  library(ggplot2)
  library(treemapify)
  
  ################################################################################
  theme7point()
  
  hsComp <- function(dfX,col1,col2){
    dfX <- dfScatters
    dfX <- dfX[,c(col1,col2)]
    names(dfX) <- c('C1','C2')
    
    dfX <- dfX %>% 
      dplyr:::filter(!is.na(C1) | !is.na(C2)) %>%
      dplyr:::filter(C1 > 0 | C2 > 0)
    
    dfPlot <- dfX %>% 
      rowwise %>% 
      mutate(comp     = ifelse(is.na(C1),col2,ifelse(is.na(C2),col1,'Both')),
             filltype = ifelse(is.na(C1),'1',ifelse(is.na(C2),'2','Both'))) %>%
      group_by(comp,filltype) %>%
      dplyr:::count() %>%
      mutate(n=n/1000,lbl=paste0(comp," [",round(n),"]"))
    
    g <- ggplot(dfPlot,aes(x=comp, y=n)) + 
      geom_bar(aes(fill=filltype),width=0.8,color='black',lwd=.2,stat='identity') + 
      xlab('') + ylab('') +
      coord_cartesian(xlim=c(0.55,3.45),ylim=c(0,max(dfPlot$n)*1.1),expand=FALSE) + 
      scale_y_continuous(position='right') +
      theme(legend.position='none',
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            plot.margin=unit(c(0,0,0,0),'cm'))
    
    g <- ggplot(dfPlot,aes(y="", fill=comp, x=n)) + 
      geom_bar(aes(fill=filltype),width=1,color='black',lwd=.2,stat='identity') +
      coord_polar()+
      xlab('') + ylab('') +
      theme_void() + 
      theme(legend.position='none',
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            plot.margin=unit(c(0,0,0,0),'cm')) + 
      scale_fill_manual(values=c('violet','orange','grey30'))
    
    g <- ggplot(dfPlot,aes(area=n, label=comp, fill=filltype)) + 
      geom_treemap(color='black',lwd=.4) + 
      geom_treemap_text(colour = "white", 
                        place = "centre",
                        grow = TRUE,) + 
      theme(legend.position='none') + 
      scale_fill_manual(values=c('violet','orange','grey30'))
    
    #ylab(bquote('Hotspots (x'*10^3*')')) + 
    #geom_text(aes(label=lbl),size=7*5/14,vjust=-0.5,color='black') + 
    return(g)
  }
  
  hsScatter <- function(dfX,col1,col2){
    
    dfX <- dfX[,c(col1,col2)]
    names(dfX) <- c('C1','C2')
    cc <- cor(dfX, use = "pairwise.complete.obs")
    cc <- cc[1,2]
    
    dfX <- dfX %>% 
      dplyr:::filter(!is.na(C1) | !is.na(C2)) %>%
      dplyr:::filter(C1 > 0 | C2 > 0)
    
    dfPlot <- dfX %>% 
      rowwise %>% 
      mutate(comp     = ifelse(is.na(C1),col2,ifelse(is.na(C2),col1,'Both')),
             filltype = ifelse(is.na(C1),'1',ifelse(is.na(C2),'2','Both'))) 
    
    g <- ggplot(dfPlot,aes(x=C1, y=C2)) + 
      geom_point(color='black',alpha=.1,lwd=.2) + 
      xlab(col1) + 
      ylab(col2) + 
      coord_cartesian(expand=FALSE) + 
      theme(legend.position='none',
            plot.margin=unit(c(0,0,0,0),'cm')) + 
      annotate(geom='text',x=Inf,y=-Inf,hjust=1,vjust=0,
               color='red',
               fontface='bold',
               size=6*5/14,
               label=paste0("R=",round(cc,2))) + 
      annotate(geom='text',x=-Inf,y=Inf,hjust=0,vjust=1,
               color='forestgreen',
               fontface='bold',
               size=6*5/14,
               label=paste0(" N=",format(sum(dfPlot$comp == 'Both'),big.mark=",")))
    
    return(g)
  }
  
  logTPM <- function(x){
    xl <- log10(toTPM(x))
    xl[is.infinite(xl)] <- NA
    return(xl)
  }
  
  getHotspotSrc <- function(A,B,C,N,L4){
    if (sum(A,B,C,N,L4)>1){return("X")}
    if (A){return("A")}
    if (B){return("B")}
    if (C){return("C")}
    if (N){return("N")}
    if (L4){return("L4")}
    return("X")
  }
  
  drawHSCounts_and_Strengths <- function(dfInit){
    
    dfHeat <- dfInit %>%
      mutate(is_A_HS=ifelse(AA1+AA2+AA3+AA4 > 0,TRUE,FALSE),
             is_C_HS=ifelse(is_A_HS == 0 & AC > 0,TRUE,FALSE),
             is_B_HS=ifelse(is_A_HS == 0 & AB > 0,TRUE,FALSE),
             is_N_HS=ifelse(is_A_HS == 0 & AN > 0,TRUE,FALSE),
             is_L4_HS=ifelse(is_C_HS == 0 & CL4 > 0,TRUE,FALSE))%>%
      rowwise %>%
      mutate(hsSrc = getHotspotSrc(is_A_HS,
                                   is_B_HS,
                                   is_C_HS,
                                   is_N_HS,
                                   is_L4_HS)) %>%
      select(cs,from,to,
             AA1,AA2,AA3,AA4,AB,AN,AC,CL4,
             is_A_HS,
             is_C_HS,
             is_B_HS,
             is_N_HS,
             is_L4_HS,hsSrc)
    
    dfMeltHeat <- reshape2:::melt.data.frame(data = dfHeat, 
                                             id.vars = c('cs','from','to',"hsSrc",
                                                         paste0("is_",c("A","B","N","C","L4"),"_HS")),
                                             variable.name = 'sample',
                                             value.name = 'strength') %>%
      mutate(sample = factor(sample, 
                             levels=c('AA1','AA2','AA3','AA4',
                                      'AA4a','AA4b','AB',
                                      'AN','AC','CL4')),
             stype = ifelse(grepl("A",sample),
                            ifelse(grepl("AA",sample),
                                   "AA",
                                   "A-het"),
                            "Oth")) %>%
      dplyr:::filter(sample %in% c('AA1','AA2','AA3','AA4','AB','AN','AC','CL4'),
                     strength > 0)
    
    ###########################################################
    dfHSCounts <- dfMeltHeat %>% 
      group_by(sample) %>% 
      add_count(name = "tot") %>%
      group_by(sample,hsSrc,tot) %>% 
      count(name = "n") %>%
      mutate(pc=n/tot*100)
    
    gE <- ggplot(dfHSCounts,
                 aes(y=sample,x=pc,fill=hsSrc)) + 
      geom_bar(stat='identity',color='black',lwd=.3) + 
      xlab('Hotspots (%)') + ylab('') + 
      scale_fill_manual("Hotspots defined by: ",
                        values=c("X"="grey80",
                                 "A"="salmon",
                                 "B"="darkorange2",
                                 "N"="brown2",
                                 "C"="slateblue2",
                                 "L4"="orchid"))
    
    ###########################################################
    dfPerHS <- dfMeltHeat %>% 
      group_by(sample) %>% 
      add_tally(strength,name = "tot") %>%
      group_by(sample,hsSrc,tot) %>% 
      tally(strength,name = "n") %>%
      mutate(pc=n/tot*100)
    
    gF <- ggplot(dfPerHS,aes(y=sample,x=pc,fill=hsSrc)) + 
      geom_bar(stat='identity',color='black',lwd=.3) + 
      xlab('Total strength (%)') + ylab('') + 
      scale_fill_manual("Hotspots defined by: ",
                        values=c("X"="grey80",
                                 "A"="salmon",
                                 "B"="darkorange2",
                                 "N"="brown2",
                                 "C"="slateblue2",
                                 "L4"="orchid"))
    
    dfPerHSBP <- dfMeltHeat %>% 
      group_by(sample) %>% 
      add_tally(strength,name = "tot") %>%
      dplyr:::filter(hsSrc != "X")
    
    gG <- ggplot(dfPerHSBP,aes(y=sample,fill=hsSrc,x=strength/tot*100)) + 
      geom_boxplot(outlier.alpha=0,notch=TRUE,lwd=.3,
                   position = position_dodge(preserve = "single")) + 
      scale_fill_manual("Hotspots defined by: ",
                        values=c("A"="salmon",
                                 "B"="darkorange2",
                                 "N"="brown2",
                                 "C"="slateblue2",
                                 "L4"="orchid")) + 
      scale_x_log10(labels=fancy_scientific) + 
      coord_cartesian(xlim=c(1e-4,1e-1),ylim=c(.4,8.6),expand=FALSE) + 
      xlab('Strength (%)') + ylab('')
    
    noLeg <- theme(legend.position='none')
    noM <- theme(plot.margin=unit(c(0,0,0,0),'cm'))
    
    gLeg <- get_legend(gE + 
                         theme(legend.position='top',
                               legend.direction='horizontal',
                               legend.key.size=unit(0.3,'cm')) + 
                         guides(fill=guide_legend(nrow=1,byrow=TRUE)))
    
    gTop <- ggarrange(ggplot() + theme_void() + noM,
                      ggplot() + theme_void() + noM,
                      ggplot() + theme_void() + noM,
                      labels=c('F','G','H'),
                      ncol=3,nrow=1,
                      font.label = list(size=8,fontface='bold'))
    
    gData <- ggarrange(gE+noM+noLeg,
                       gF+noM+noLeg,
                       gG+noLeg+noM,
                       nrow=1)
    
    
    return(list(gE=gE,
                gF=gF,
                gG=gG,
                gAll=ggarrange(gTop,gLeg,gData,ncol=1,heights=c(1,1,10))))
  }
  
  ################################################################################
  dfHotspots <- read.table('hotspotsData.tab',header=TRUE)
  
  dfHotspots <- dfHotspots[dfHotspots$cs %in% paste0('chr',1:21),]  %>%
    rename(AB=AB1) %>%
    mutate(anyHS = AA1+AA2+AA3+AA4+AB+AN+AC+CL4)
  
  ################################################################################
  ## Make Panel A
  dfScatters <- dfHotspots %>%
    mutate(AA1=logTPM(AA1),
           AA2=logTPM(AA2),
           AA3=logTPM(AA3),
           AA4=logTPM(AA4),
           AA4a=logTPM(AA4a),
           AA4b=logTPM(AA4b),
           AB=logTPM(AB),
           AN=logTPM(AN),
           AC=logTPM(AC),
           CL4=logTPM(CL4)) %>%
    select(AA1,AA2,AA3,AA4,AB,AN,AC,CL4)
  
  allIDs <- c('AA1','AA2','AA3','AA4','AB','AN','AC','CL4')
  plotList <- list()
  
  numIDs <- length(allIDs)
  for (i in allIDs){
    for (j in allIDs){
      iPos <- which(allIDs==i)
      jPos <- which(allIDs==j)
      
      barPos     <- iPos + (numIDs*(jPos - 1))
      scatterPos <- jPos + (numIDs*(iPos - 1))
        
      if (i == j){
        plotList[[scatterPos]] <- ggplot() + theme_void() + 
          annotate(geom='text',x=.5,y=.5,label=i,size=8*5/14,fontface='bold')
      }else{
        plotList[[barPos]]     <- hsComp(dfScatters,i,j)
        plotList[[scatterPos]] <- hsScatter(dfScatters,j,i)
      }
  
      if ((barPos)%%numIDs > 0){
        plotList[[barPos]] <- plotList[[barPos]] + ylab('')
      }
      
      if ((scatterPos-1)%%numIDs > 0){
        plotList[[scatterPos]] <- plotList[[scatterPos]] + ylab('')
      }
      
      if (iPos < numIDs){
        plotList[[scatterPos]] <- plotList[[scatterPos]] + xlab('')
      }
    }
  }
  
  gAllvAll  <- ggarrange(plotlist = plotList)
  
  yRight  <- textGrob(expression(paste("Hotspots (thousands)")), rot=90, gp = gpar(fontsize = 8))
  yLeft   <- textGrob(expression(paste("Hotspot Strength (log)")), rot=90,gp = gpar(fontsize = 8))
  xBottom <- textGrob(expression(paste("Hotspot Strength (log)")), gp = gpar(fontsize = 8))
  
  gSave <- grid.arrange(gAllvAll,
                      right=yRight,
                      bottom=xBottom,
                      left=yLeft)
  
  ggsave('tempPanelA.png',gSave,bg='white',height=6,width=6)
  
  # gP <- ggpairs(dfScatters , switch="both",
  #               lower = list(continuous = wrap("points", alpha = 0.1, size=0.05)), 
  #               diag  = list(continuous = wrap("blank")), 
  #               upper = list(continuous = wrap("cor",stars=FALSE,method='pearson',digits=2,title='R',size=7*5/14,color='red'))) + 
  #   theme(strip.placement = "inside",
  #         strip.background=element_rect(fill='grey90'),
  #         strip.text=element_text(size=7,face='bold'),
  #         panel.border=element_rect(fill=NA)) + 
  #     scale_x_continuous(breaks=0:3,labels=c('',1,'',3)) +
  #     scale_y_continuous(breaks=0:3,labels=c('',1,'',3)) +
  #     xlab(bquote('Hotspot strength ('*log[10]*')')) + 
  #     ylab(bquote('Hotspot strength ('*log[10]*')')) 
  
  ################################################################################
  ## Make other panels
  dfMeltHS <- reshape2:::melt.data.frame(data = dfHotspots, 
                                       id.vars = c('cs','from','to','anyHS'),
                                       variable.name = 'sample',
                                       value.name = 'strength') %>%
  mutate(sample = factor(sample, 
                         levels=c('AA1','AA2','AA3','AA4',
                                  'AA4a','AA4b','AB',
                                  'AN','AC','CL4')),
         stype = ifelse(grepl("A",sample),
                        ifelse(grepl("AA",sample),
                               "AA",
                               "A-het"),
                        "Oth")) %>%
  dplyr:::filter(sample %in% c('AA1','AA2','AA3','AA4','AB','AN','AC','CL4'))
  
  
  #### Hotspot Counts
  dfCounts <- dfMeltHS %>% 
    group_by(sample,stype) %>% 
    filter(strength>0) %>% 
    count() 
  
  gHSCount <- ggplot(dfCounts,aes(y=sample,x=n/1000,fill=stype)) +
    geom_bar(color='black',lwd=.2,stat='identity') +  
    xlab(bquote('Hotspots (x'*10^3*')')) + 
    ylab('') + theme(legend.position='none') + 
    coord_cartesian(expand=FALSE,clip='off')
  
  ## Unique hotspots per sample
  dfUnique <- dfMeltHS %>% 
    mutate(unique=ifelse(strength==anyHS,1,0)) %>%
    group_by(sample,stype) %>% 
    summarize(tot=sum(strength>0), 
              unique=sum(unique>0)) %>%
    mutate(uniquePC=unique/tot*100)
  
  gUnique <- ggplot(dfUnique,
                    aes(x=uniquePC,
                        y=sample, 
                        fill=stype)) + 
    geom_bar(color='black',lwd=.2,stat='identity') + 
    xlab('Unique hotspots\n(% of total)') + 
    ylab('') + theme(legend.position='none') +
    coord_cartesian(expand=FALSE,clip='off')
  
  ## Unique hotspots per sample (downsampled to smallest)
  dfTopHS <- dfMeltHS %>%
    group_by(sample,stype) %>%
    slice_max(strength, 
              n = min(dfCounts$n), 
              with_ties = FALSE) %>%
    mutate(unique=ifelse(strength==anyHS,1,0)) %>%
    group_by(sample,stype) %>% 
    summarize(tot=sum(strength>0), 
              unique=sum(unique>0)) %>%
    mutate(uniquePC=unique/tot*100)
  
  gUniqueTop <- ggplot(dfTopHS,
                       aes(x=uniquePC,
                           y=sample, 
                           fill=stype)) + 
    geom_bar(color='black',lwd=.2,stat='identity') + 
    xlab(paste0('Unique hotspots\n(% of top ',min(dfCounts$n),')')) + 
    ylab('') + theme(legend.position='none') +
    coord_cartesian(expand=FALSE,clip='off')
  
  #Get all CCs
  dfCCs <- reshape2:::melt.data.frame(as.data.frame(cor(dfScatters, 
                                                        method='pearson', 
                                                        use = "pairwise.complete.obs")),
                                      variable.name = "sample",
                                      value.name = "CC") %>% 
    dplyr:::filter(CC<1) %>%
    mutate(stype = ifelse(grepl("A",sample),ifelse(grepl("AA",sample),"AA","A-het"),"Oth"))
  
  gCorBoxplot <- ggplot(dfCCs,aes(y=sample,x=CC)) + 
    geom_boxplot(aes(fill=stype),
                 color='grey50',
                 lwd=.3,
                 outlier.alpha=0) + 
    geom_point(size=.3,
               position=position_jitter(height=0,width=.2)) + 
    ylab('') + xlab('Strength correlation with\nother samples (Pearson R)') + 
    theme(legend.position='none',
          axis.text.y=element_text(angle=0,hjust=1)) 
  
  noX    <- theme(axis.text.x = element_blank())
  noY    <- theme(axis.text.y = element_blank())
  noMarg <- theme(plot.margin = margin(c(0,0,0,0),'cm'))
  
  fillColz <- scale_fill_manual(values=c('dodgerblue3','grey30','dodgerblue1'))
  gStats <- ggarrange(gHSCount + noMarg + fillColz,
                      gCorBoxplot + noY + noMarg + ylab('') + fillColz,
                      gUnique+noY+noMarg + fillColz,
                      gUniqueTop+noY+noMarg + fillColz,
                      nrow=1,ncol=4,
                      align='h',
                      labels=c('B','C','D','E'),
                      font.label = list(size=8,fontface='bold',
                                        hjust=0,vjust=1))
  
  gFGH <- drawHSCounts_and_Strengths(dfHotspots)
  
  gOut <- ggarrange(gSave,gStats,gFGH$gAll,
                    ncol=1,nrow=3,
                    heights=c(1,.25,.35),
                    labels=c('A',''),
                    font.label = list(size=8,fontface='bold',hjust=0,vjust=1))
  
  ggsave('Alleva_et_al_Supplemental_Hotspot_Analyses.png',gOut,
         bg='white',height=10,width=6.5)
