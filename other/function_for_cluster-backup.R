require("cluster")

#evaluate class number
selectK=function(expr,method,kmax=round(sqrt(dim(expr)[1])),distanceMatrices,fileName)
{
  
  #pdf(fileName)
  criteria <- CHCriterion( data =expr,method, kmax = kmax,distanceMatrices)
  # result 
  head(criteria$data)
  
  criteria$plot
  #dev.off()
}


doCluster=function(distanceMatrices,cluster_method,classNum,fileName)
{
  pdf(fileName)
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
  #if(cluster_method=="kmeans")
  #{
  #  (cl <- kmeans(x, 5, nstart = 25))
  #  memb=cl$cluster
  #}
  if(cluster_method=="kmedoids")
  {
    res=fastkmed(distanceMatrices, ncluster=classNum)
    memb=res$cluster
  }
  dev.off()
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
    for( k in 2:kmax )
    {
      res=fastkmed(distanceMatrices,k)
      groups=res$cluster
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
    res=hclust(distObject,method = "average")
    #plot(hc)
    partition=cutree(res,k=classNum)
  }
  if(cluster_method=="complete")
  {
    res=hclust(distObject,method = "complete")
    #plot(hc)
    partition=cutree(res,k=classNum)
  }
  if(cluster_method=="single")
  {
    res=hclust(distObject,method = "single")
    #plot(hc)
    partition=cutree(res,k=classNum)
  }

  if(cluster_method=="kmedoids")
  {
    res=fastkmed(distObject, ncluster=classNum)
    partition=res$cluster
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
  
  cluster_result <- list( result      = res,	                        
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
    sampling  <- sample( n, round(n*4/5), replace = F )
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

