
require(cluster)
require(pamr)
require(gplots)
require(ROCR)
require(kmed)

check_max_k=function(distanceMatrices1,k)
{
  distdata=distanceMatrices1
  ncluster=k
  
  sorted_object <- order(unique(colSums(distdata/sum(distdata))))
  medoid_init <- sorted_object[1:ncluster]
  
  dist_0 <- distdata[, medoid_init]
  member_0 <- apply(dist_0, 1, which.min)
  
  check_k=max(member_0)
  check_k
}

#evaluate class number
selectK=function(expr,cluster_method,kmax=round(sqrt(dim(expr)[1])),distanceMatrices,fileName,distance_method=NA)
{
  
  if(cluster_method=="kmedoids")
  {
    for(k in 2:kmax)
    {
      check_k=check_max_k(distanceMatrices,k)
      if(check_k<k)
      {
        kmax=k-1
        print("break")
        break
        
      }
    }
  }

  print(paste("kmax",kmax))
  
  #if(kmax <2)
  #{
   # return(NA)
  #}
  
  
  criteria <- CHCriterion( data =expr,cluster_method, kmax = kmax,distanceMatrices)
  # result 
  #head(criteria$data)
  
  if(!is.na(criteria))
  {
    pdf(fileName)
    print(criteria$plot)
    dev.off()
  }
  
  
}


selectK_gap=function(expr,cluster_method,kmax=round(sqrt(dim(expr)[1])),distanceMatrices,fileName,distance_method)
{
  
  if(cluster_method=="kmedoids")
  {
    for(k in 2:kmax)
    {
      check_k=check_max_k(distanceMatrices,k)
      if(check_k<k)
      {
        kmax=k-1
        print("break")
        break
        
      }
    }
  }
  
  print(paste("kmax",kmax))
  
  
  if(kmax>1)
  {

  gaps <- compGapStats(expr, K.max=kmax,nboot=3,distance_method=distance_method,cluster_method=cluster_method)
  gapsmat <- gaps$gapsmat # gap value
  gapsSE <- gaps$gapsSE # standard error of gap
  
  # use the best K in top 300
  bestK <-  which(max(gapsmat[, 3]) == gapsmat[, 3])
  print(paste("bestk",bestK))
  pdf(file=fileName, width=14, height=5)
  figGAP(gapsmat, gapsSE)
  graphics.off()
  
  }
  
  bestK
}

doCluster=function(distanceMatrices,cluster_method,classNum,fileName)
{
  if(!is.na(fileName))
  {
    pdf(fileName)
  }
  
  distObject=as.dist(distanceMatrices)
  #hcluster
  if(cluster_method=="average")
  {
    hc=hclust(distObject,method = "average")
    plot(hc)
    memb=cutree(hc,k=classNum)
  }
  if(cluster_method=="complete")
  {
    hc=hclust(distObject,method = "complete")
    plot(hc)
    memb=cutree(hc,k=classNum)
  }
  if(cluster_method=="single")
  {
    hc=hclust(distObject,method = "single")
    plot(hc)
    memb=cutree(hc,k=classNum)
  }
  if(cluster_method=="kmedoids")
  {
    res=fastkmed(distanceMatrices, ncluster=classNum)
    memb=res$cluster
    
    
  }
  if(cluster_method=="pam")
  {
    memb=pam(distanceMatrices,k=classNum,diss=T, cluster.only=TRUE)

  }
  if(cluster_method=="hcluster_dynamicK")
  {
    hc=hclust(distObject,method = "average")
    k=classNum
    
    runagin=T
    while(runagin)
    {
      
      memb=cutree(hc,k=k)
      
      runagin=F
      noise_class_num=0
      for(classi in unique(memb))
      {
        A=which(memb==classi,arr.ind = T)
        if(length(A)<4)
        {
          noise_class_num=noise_class_num+1
        }
      }
      if(noise_class_num>(k-classNum))
      {
        runagin=T
        k=k+1
      }
      if(k>sqrt(dim(distanceMatrices)[1]))
      {
        runagin=F
      }
      
    }
    memb=cutree(hc,k=k)
    print(paste("assigned class Num",classNum,"dynamic K", k))
    
  }
  
  if(!is.na(fileName))
  {
    dev.off()
  }
  
  #return hcluster memb
  memb
}

convertID=function(geneName,fromID,toID)
{
  index=1
  geneSymbol=c()
  for(i in c(1:length(geneName)))
  {
    A=which(geneName[i]==fromID)
    if(length(A)>0)
    {
      geneSymbol[index]=toID[A]
      index=index+1
    }
    
  }
  geneSymbol
}







# [Distance] : calculates the sum squared distance of a given cluster of points,
#              note that "sum squared distance" is used here for measuring variance 
Distance <- function(cluster)
{
  # the center of the cluster, mean of all the points
  center <- colMeans(cluster)
  
  # calculate the summed squared error between every point and 
  # the center of that cluster 
  distance <- apply( cluster, 1, function(row)
  {
    sum( ( row - center )^2 )
  }) %>% sum()
  
  return(distance)
}

# calculate the within sum squared error manually for hierarchical clustering 
# [WSS] : pass in the dataset, and the resulting groups(cluster)
WSS <- function( data, groups )
{
  k <- max(groups)
  
  # loop through each groups (clusters) and obtain its 
  # within sum squared error 
  total <- lapply( 1:k, function(k)
  {
    # extract the data point within the cluster
    cluster <- subset( data, groups == k )
    
    distance <- Distance(cluster)
    return(distance)
  }) %>% unlist()
  
  return( sum(total) )
}

# testing 
# sum_squared_error <- WSS( data = mtcars_scaled, groups =  groups )

# this value will will decrease as the number of clusters increases, 
# because each cluster will be smaller and tighter.
# And the rate of the decrease will slow down after the optimal cluster number


## method 2 : Calinski-Harabasz index, ratio of the between cluster variance
#			  to the total within cluster variance
# http://www.mathworks.com/help/stats/clustering.evaluation.calinskiharabaszevaluation-class.html 

# TSS (total sum of square) : the squared distance of all the data points from 
# the dataset's centroid 

# BSS (between sum of square) = TSS - WSS, measures how far apart are the clusters
# from each other 
# !! a good clustering has a small WSS and a high BSS

# CHIndex = B / W, the ratio should be maximized at the optimal k
# B = BSS(k) / (k-1) ; k = # of cluster
# W = WSS(k) / (n-k) ; n = # of data points

# [CHCriterion] : calculates both Calinski-Harabasz index and within sum squared error
# @kmax          = maximum cluster number, caculates the CH index from 2 cluster to kmax
# @clustermethod = "hclust"

CHCriterion <- function( data,method, kmax, distanceMatrices, ...  )
{
  
  
  # total sum squared error (independent with the number of cluster k)
  tss <- Distance( cluster = data )
  
  # initialize a numeric vector storing the score
  wss <- numeric(kmax)
  
  # k starts from 2, cluster 1 is meaningless
  if(method=="kmedoids")
  {
    if(kmax<2)
    {
      return(NA)
    }
    for( k in 2:kmax )
    {
      print(paste("k",k))
      #write.csv(distanceMatrices,"test_distMat.csv",quote=F)
      res=fastkmed(distanceMatrices,k)
      print(paste("k",k,"finished"))
      groups=res$cluster
      wss[k] <- WSS( data = data, groups =  groups )
    }
  }else if(method=="pam")
  {
    
    for( k in 2:kmax )
    {
      memb=pam(distanceMatrices,k=k,diss=T, cluster.only=TRUE)
      groups <- memb
      wss[k] <- WSS( data = data, groups =  groups )
    }
  }else{
    d=as.dist(distanceMatrices)
    hc <- hclust(d, method=method )
    for( k in 2:kmax )
    {
      
      groups <- cutree( hc, k )
      wss[k] <- WSS( data = data, groups =  groups )
    }
  }

  
  
  # between sum of square
  bss <- tss - wss[-1]
  
  # cluster count start from 2! 
  numerator <- bss / ( 1:(kmax-1) )
  denominator <- wss[-1] / ( nrow(data) - 2:kmax )
  
  criteria <- data.frame( k = 2:kmax,
                          CHIndex = numerator / denominator,
                          wss = wss[-1] )
  
  # convert to long format for plotting 
  criteria_long <- gather( criteria, "index", "value", -1 )
  
  
  plot <- ggplot( criteria_long, aes( k, value, color = index ) ) + 
    geom_line() + geom_point( aes( shape = index ), size = 3 ) +
    facet_wrap( ~ index, scale = "free_y" ) + 
    guides( color = FALSE, shape = FALSE )

  
  return( list( data = criteria, 
                plot = plot ) )
}




ClusterMethod_fromDist=function(distanceMatrices,classNum, noise.cut=noise.cut,cluster_method)
{
  
  distObject=as.dist(distanceMatrices)
  #hcluster
  if(cluster_method=="average")
  {
    hc=hclust(distObject,method = "average")
    #plot(hc)
    partition=cutree(hc,k=classNum)
  }
  if(cluster_method=="complete")
  {
    hc=hclust(distObject,method = "complete")
    #plot(hc)
    partition=cutree(hc,k=classNum)
  }
  if(cluster_method=="single")
  {
    hc=hclust(distObject,method = "single")
    #plot(hc)
    partition=cutree(hc,k=classNum)
  }
  if(cluster_method=="kmedoids")
  {
    hc=NA
    res=fastkmed(distanceMatrices, ncluster=classNum)
    partition=res$cluster
    
    
  }
  if(cluster_method=="pam")
  {
    hc=NA
    res=pam(distanceMatrices,k=classNum,diss=T)
    partition=res$clustering

    
  }
  if(cluster_method=="hcluster_dynamicK")
  {
    hc=hclust(distObject,method = "average")
    k=classNum
    
    runagin=T
    while(runagin)
    {
     
      #print(paste("k",k)) 
      memb=cutree(hc,k=k)
      #print(paste("memb",memb)) 
      runagin=F
      noise_class_num=0
      for(classi in unique(memb))
      {
        A=which(memb==classi,arr.ind = T)
        if(length(A)<4)
        {
          noise_class_num=noise_class_num+1
        }
      }
      if(noise_class_num>(k-classNum))
      {
        runagin=T
        k=k+1
      }
      if(k>sqrt(dim(distanceMatrices)[1]))
      {
        runagin=F
      }
      
    }
    partition=cutree(hc,k=k)
    
    
  }
  
  ###################modify begins here
  # equivalent to k
  cluster_num <- max(partition) 
  
  # calculate each cluster's size 
  cluster_size <- numeric(cluster_num)
  for( i in 1:cluster_num )
    cluster_size[i] <- sum( partition == i )
  
  # if there're cluster size smaller than the specified noise.cut, do :
  not_noise_num <- sum( cluster_size > noise.cut )
  
  if( cluster_num > not_noise_num )
  {
    # extract the cluster whose size is larger than noise.cut
    cluster_new <- (1:cluster_num)[ cluster_size > noise.cut ]
    
    # all the data points whose original cluster is smaller than the noise.cut
    # will be assigned to the same new cluster
    cluster_num <- not_noise_num + 1
    
    # new clustering number, assign the noise cluster's number first
    # then adjust the original cluster's number
    new <- rep( cluster_num, nrow(data) )
    
    for( i in 1:not_noise_num )
      new[ ( partition == cluster_new[i] ) ] <- i
    
    partition <- new
  }
  

  # boolean vector indicating which data point belongs to which cluster
  cluster_list <- lapply( 1:cluster_num, function(x)
  {
    return( partition == x )
  })

  cluster_result <- list( result      = hc,	                        
                          partition   = partition,
                          clusternum  = cluster_num,
                          clusterlist = cluster_list )
  return(cluster_result)
}

ClusterBootstrap_fromDist <- function( distanceMatrices, k, noise.cut = 0, bootstrap = 100, 
                                       dissolve = .5, cluster_method)
{
  # step 1
  cluster_result <- ClusterMethod_fromDist(distanceMatrices, k, noise.cut, 
                                           cluster_method)
  
  
  
  cluster_num  <- cluster_result$clusternum
  boot_jaccard <- matrix( 0, nrow = bootstrap, ncol = cluster_num )
  
  # pass in two vectors containing TRUE and FALSE
  # ( do not use built in intersect or union ! )
  jaccardSimilarity <- function( x, y )
  {
    jaccard <- sum( x & y ) / ( sum(x) + sum(y) - sum( x & y ) )
    return(jaccard)
  }
  
  n <- nrow(distanceMatrices)
  for( i in 1:bootstrap )
  {
    # step 2, cluster the new sampled data 
    sampling  <- sample( n, round(n*4/5), replace = F )
    #sampling  <- sample( n, n, replace = T )
    boot_distanceMatrices <- distanceMatrices[ sampling,sampling ]
    
    boot_result <- ClusterMethod_fromDist(boot_distanceMatrices, k, noise.cut, 
                                          cluster_method)
    boot_num <- boot_result$clusternum
    
    # step 3
    for( j in 1:cluster_num )
    {
      # compare the original cluster with every other bootstrapped cluster
      similarity <- lapply( 1:boot_num, function(k)
      {
        jaccard <- jaccardSimilarity( x = cluster_result$clusterlist[[j]][sampling],
                                      y = boot_result$clusterlist[[k]] )
      }) 
      similarity=unlist(similarity)
      
      # return the largest jaccard similarity
      boot_jaccard[ i, j ] <- max(similarity)
    }	
  }
  
  # cluster's stability, mean of all the boostrapped jaccard similarity 
  boot_mean <- colMeans(boot_jaccard)
  
  # how many times are each cluster's jaccard similarity below the 
  # specified "dissolved" value  
  boot_dissolved <- apply( boot_jaccard, 2, function(x)
  {
    sum( x < dissolve, na.rm = TRUE )
    
  })
  
  boot_result <- list( result        = cluster_result$result,
                       bootmean      = boot_mean,
                       partition     = cluster_result$partition,
                       clusternum    = cluster_num,                     
                       bootdissolved = boot_dissolved,
                       clusterlist   = cluster_result$clusterlist)
  return(boot_result)
}


ClusterBootstrap_fromDist_sampleTogetherMat <- function( distanceMatrices, k, noise.cut = 0, bootstrap = 100, 
                                                         dissolve = .5, cluster_method)
{
  # step 1
  cluster_result <- ClusterMethod_fromDist(distanceMatrices, k, noise.cut, 
                                           cluster_method)
  cluster_num  <- cluster_result$clusternum
  sample_num <- length(cluster_result$partition)
  boot_jaccard <- matrix( 0, nrow = bootstrap, ncol = cluster_num )
  sampleTogetherMat <- matrix(0, nrow=sample_num, sample_num)
  sampleTogether_clusteredTimes_Mat <- matrix(0, nrow=sample_num, sample_num)
  ref_sampleTogetherMat <- matrix(0, nrow=sample_num, sample_num)
  
  # pass in two vectors containing TRUE and FALSE
  # ( do not use built in intersect or union ! )
  jaccardSimilarity <- function( x, y )
  {
    jaccard <- sum( x & y ) / ( sum(x) + sum(y) - sum( x & y ) )
    return(jaccard)
  }
  
  n <- nrow(distanceMatrices)
  for( i in 1:bootstrap )
  {
    # step 2, cluster the new sampled data 
    #sampling  <- sample( n, round(n*4/5), replace = F )
    sampling  <- sample( n, n, replace = T )
    boot_distanceMatrices <- distanceMatrices[ sampling,sampling ]
    
    boot_result <- ClusterMethod_fromDist(boot_distanceMatrices, k, noise.cut, 
                                          cluster_method)
    boot_num <- boot_result$clusternum
    
    # step 3
    for( j in 1:cluster_num )
    {
      # compare the original cluster with every other bootstrapped cluster
      similarity <- lapply( 1:boot_num, function(k)
      {
        jaccard <- jaccardSimilarity( x = cluster_result$clusterlist[[j]][sampling],
                                      y = boot_result$clusterlist[[k]] )
      }) 
      similarity=unlist(similarity)
      
      # return the largest jaccard similarity
      boot_jaccard[ i, j ] <- max(similarity)
      
      
    }	
    
    for( q in 1:boot_num)
    {
      in_same_cluster_array=boot_result$clusterlist[[q]]
      A=which(in_same_cluster_array==1, arr.ind = T)
      
      sampleTogetherMat[sampling[A],sampling[A]]=sampleTogetherMat[sampling[A],sampling[A]]+1
      
    }
    sampleTogether_clusteredTimes_Mat[sampling,sampling]=sampleTogether_clusteredTimes_Mat[sampling,sampling]+1
  }
  
  for(q in 1:cluster_num)
  {
    in_same_cluster_array=cluster_result$clusterlist[[q]]
    A=which(in_same_cluster_array==1, arr.ind = T)
    ref_sampleTogetherMat[A,A]=ref_sampleTogetherMat[A,A]+1
  }
  
  # cluster's stability, mean of all the boostrapped jaccard similarity 
  boot_mean <- colMeans(boot_jaccard)
  
  # how many times are each cluster's jaccard similarity below the 
  # specified "dissolved" value  
  boot_dissolved <- apply( boot_jaccard, 2, function(x)
  {
    sum( x < dissolve, na.rm = TRUE )
  })
  
  boot_result <- list( result        = cluster_result$result,
                       bootmean      = boot_mean,
                       partition     = cluster_result$partition,
                       clusternum    = cluster_num,                     
                       bootdissolved = boot_dissolved,
                       clusterlist   = cluster_result$clusterlist,
                       sampleTogetherMat = sampleTogetherMat,
                       sampleTogether_clusteredTimes_Mat = sampleTogether_clusteredTimes_Mat,
                       ref_sampleTogetherMat = ref_sampleTogetherMat)
  return(boot_result)
}




compGapStats <- function(ge.CRC, ntops=c(1, 2, 3)*100, K.max=6, nboot=100,distance_method,cluster_method) {
  
  MAD <- apply(ge.CRC, 1, mad)
  names(MAD)=c(1:length(MAD))
  ords <- names(sort(MAD, decreasing=TRUE))
  
  if(length(ords) < 300) {
    ntops <- c(2, 3, length(ords));
  }
  
  fun <- function(x, k) {
    
    distanceMatrices=produceDistanceMatrices(ge.CRC,distance_method)
    clus=doCluster(distanceMatrices,cluster_method,k,NA)
    return(list(cluster=clus))
  }
  gaps <- list()
  for(i in 1:length(ntops)) {
    ntop <- ntops[i]
    sdat <- ge.CRC[ords[1:ntop], ]
    sdat <- sweep(sdat,1, apply(sdat,1,median))
    gaps[[i]] <- clusGap(t(sdat), FUNcluster=fun, K.max=K.max, B=nboot)
  }
  gapsmat <- matrix(0, K.max, length(ntops))
  gapsSE <- matrix(0, K.max, length(ntops))
  for(i in 1:length(ntops)) {
    gapsmat[, i] <- gaps[[i]]$Tab[, 3]
    gapsSE[, i] <- gaps[[i]]$Tab[, 4]
  }
  
  colnames(gapsmat) <- colnames(gapsSE) <- ntops
  rownames(gapsmat) <- rownames(gapsSE) <- 1:K.max
  return(list(gapsmat=gapsmat, gapsSE=gapsSE))
}


figGAP <- function(gapsmat, gapsSE) {
  ##	gapsmat <- gapsSE <- NULL
  ##	data("gaps", package="DeSousa2013", envir = environment())
  n=min(dim(gapsmat)[1],dim(gapsSE)[1])
  par(mfrow=c(2, 3))
  for(i in 1:ncol(gapsmat)) {
    plotCI(2:n, gapsmat[2:n, i], uiw = gapsSE[2:n, i],
           liw = gapsSE[2:n, i], xlab="No. of clusters", ylab="GAP")
    lines(2:n, gapsmat[2:n, i])
    title(main=paste("top ", colnames(gapsmat)[i], sep=""))
  }
}

bootstrap_scalePCA=function (x, FUNcluster, K.max, B = 100, d.power = 1) 
{
  stopifnot(is.function(FUNcluster), length(dim(x)) == 2, 
            K.max >= 2, (n <- nrow(x)) >= 1, ncol(x) >= 1)
  
  n <- nrow(x)
  
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  
  spaceH0 = c("scaledPCA","original")
  spaceH0 <- match.arg(spaceH0)
  xs <- scale(x, center = TRUE, scale = FALSE)
  m.x <- rep(attr(xs, "scaled:center"), each = n)
  
  switch(spaceH0, scaledPCA = {
    V.sx <- svd(xs, nu = 0)$v
    xs <- xs %*% V.sx
  }, original = {})
  
  rng.x1 <- apply(xs, 2L, range)
  
  
  for (b in 1:B) {
    z1 <- apply(rng.x1, 2, function(M, nn) runif(nn, min = M[1], max = M[2]), nn = n)
    
    z <- switch(spaceH0, scaledPCA = tcrossprod(z1, V.sx), original = z1) + m.x
    
    for (k in 1:K.max) {
      FUNcluster(z, k)$cluster
      #evaluate boot
    }
    
  }
}
