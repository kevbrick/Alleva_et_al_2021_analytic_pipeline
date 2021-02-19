#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("--distmatrix"), 
              type="character", 
              default=NULL,  
              help="Dissimilarity matrix csv",        
              metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

c_counts <- data.frame(numclusters=999,clust=999,Freq=999);

#################################################################################################
getOptimalClusterNumber <- function(hc,dm,max_clusters=5){
  ## Vector to store cluster mean distances
  means <- rep(999,max_clusters)

  ## check each number of clusters and find optimum
  ## Optimum = minimum within cluster mean distance 
  ##           (only consider two largest clusters; 1&2)
  for (num_clusters in 1:max_clusters){
    cdata                <- cutree(hc,num_clusters)
    
    ## per cluster mean distances
    means[num_clusters]  <- mean(dm[which(cdata == 1),which(cdata == 1)])
    
    if (num_clusters > 1){
      means[num_clusters] <- means[num_clusters] + mean(dm[which(cdata == 2),which(cdata == 2)])
    }
  }
  
  optimal_clusters <- min(which(means==min(means)))
  cdata            <- cutree(hc,optimal_clusters)
  
  ## mark unusable sequences as being in cluster 9
  ## 1. sequences in cluster with one member sequence
  ## 2. sequences in cluster 3+ (can only have two haplotypes)
  for (i in 1:optimal_clusters){if (sum(cdata == i) == 1){cdata[cdata==i] <- 9}}
  cdata[cdata > 2] <- 9
  
  return(cdata)
}

##################################################################
full_dist_matrix <- as.matrix(read.csv(opt$distmatrix,header=FALSE))
pngOut           <- paste0(opt$distmatrix,'.png')
clustOut         <- paste0(opt$distmatrix,'.clusters.csv')

dist_matrix      <- full_dist_matrix
dist_matrix[upper.tri(dist_matrix,diag=TRUE)] <- NA
diag(dist_matrix) <- 0

dist_matrix <- as.dist(dist_matrix, diag = TRUE)

h <-hclust(dist_matrix)
png(filename = pngOut)
print(plot(h))
dev.off()

#clust <- cutree(h,k=2)
clust <- getOptimalClusterNumber(h,full_dist_matrix,5)

write.table(clust,
            file=clustOut, 
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = ",")

