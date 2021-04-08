#source and library
library("kmed")
library("clues")
library("ggplot2")
library("abind")
source("function_for_cluster.R")
source("function_for_distance.R")
source("function_for_test_difference.R")

#input and parameter
input_dir="./CompCancer/"
list_file=paste(input_dir,"list.txt",sep="")

#load data and reference cluster
fileList=as.matrix(read.table(list_file,header = F))


do_sampleCluster_oneTest_bootstrapTest_forV=function(fileName,distance_method,cluster_method,output_dir,ith,bootstrap=10,noise.cut=0)
{
  #load data and reference cluster
  data=load_data(fileName)
  expr=data$expr
  geneID=data$geneID
  sampleID=data$sampleID
  ref_cluster_array=data$ref_cluster_array
  ref_cluster=data$ref_cluster
  
  unique_ref_cluster=unique(as.character(ref_cluster))
  k=length(unique_ref_cluster)
  #print("estimate k")
  #print(k)
  
  distanceMatrices=produceDistanceMatrices(expr,distance_method)
  ###############################################################
  cluster_result <- ClusterMethod_fromDist(distanceMatrices, k, noise.cut, 
                                           cluster_method)
  
  
  
  cluster_num  <- cluster_result$clusternum
  boot_jaccard <- matrix( 0, nrow = bootstrap, ncol = cluster_num )
  
  # pass in two vectors containing TRUE and FALSE
  # ( do not use built in intersect or union ! )

  n <- nrow(distanceMatrices)
  boot_partitions=matrix(0,nrow = bootstrap, ncol=n)
  boot_samplings=matrix(0,nrow=bootstrap,ncol=n)
  #boot_hc_merges=matrix(NA,nrow=bootstrap,ncol=1)
  boot_hc_merges = c()
  for( i in 1:bootstrap )
  {
    # step 2, cluster the new sampled data 
    #sampling  <- sample( n, round(n*4/5), replace = F )
    sampling  <- sample( n, n, replace = T )
    boot_distanceMatrices <- distanceMatrices[ sampling,sampling ]
    
    boot_result <- ClusterMethod_fromDist(boot_distanceMatrices, k, noise.cut, 
                                          cluster_method)
    
    # step 3
    boot_partitions[i,]=boot_result$partition
    boot_samplings[i,]=sampling
    #print(summary(boot_result$result))
    #boot_hcs=cbind(boot_hcs,list(boot_result$result))
    boot_hc_merges = abind(boot_hc_merges,boot_result$result$merge,along=3)
  }
  

  
  t_indices=adjustedRand(cluster_result$partition,ref_cluster_array)
  ARI=t_indices[1]


  
  
  boot_result <- list( result        = cluster_result$result,
                       partition     = cluster_result$partition,
                       clusternum    = cluster_num,
                       clusterlist   = cluster_result$clusterlist,
                       boot_partitions = boot_partitions,
                       boot_hc_merges = boot_hc_merges,
                       boot_samplings = boot_samplings,
                       ARI = ARI,
                       data = data,
                       distanceMatrices = distanceMatrices)
  return(boot_result)
  
  
  ###############################################################
}

compare_two_bootstrap=function(boot_result1, boot_result2)
{
  
  cluster_num1=boot_result1$clusternum
  cluster_num2=boot_result2$clusternum
  diff_array=rep(0,cluster_num1)
  
  mat_jaccard <- matrix( 0, nrow = cluster_num1, ncol = cluster_num2 )
  
  # pass in two vectors containing TRUE and FALSE
  # ( do not use built in intersect or union ! )
  jaccardSimilarity <- function( x, y )
  {
    jaccard <- sum( x & y ) / ( sum(x) + sum(y) - sum( x & y ) )
    return(jaccard)
  }
  
  
  for( i in 1:cluster_num1 )
  {
    for(j in 1:cluster_num2)
    {
      jaccard <- jaccardSimilarity( x = boot_result1$clusterlist[[i]],
                                    y = boot_result2$clusterlist[[j]] )
      mat_jaccard[i,j]=jaccard
    }
  }
  
  for(i in 1:cluster_num1)
  {
    array_jaccard=mat_jaccard[i,]
    A=which(array_jaccard==max(array_jaccard),arr.ind = T)
    
    diff=boot_result1$bootmean[i]-boot_result2$bootmean[A]
    
    diff_array[i]=diff
  }
  
  diff_array
}






jaccardSimilarity <- function( x, y )
{
  jaccard <- sum( x & y ) / ( sum(x) + sum(y) - sum( x & y ) )
  return(jaccard)
}


load_data=function(fileName)
{
  #load data and reference cluster
  data=read.table(fileName,header = T,sep="\t")
  
  
  geneID=as.matrix(data[2:dim(data)[1],1])
  ref_cluster=c(as.matrix(data[1,c(2:dim(data)[2])]))
  ref_cluster_array=rep(0,length(ref_cluster))
  unique_ref_cluster=unique(as.character(ref_cluster))
  for(j in c(1:length(unique_ref_cluster)))
  {
    A=which(ref_cluster==unique_ref_cluster[j],arr.ind = T)
    ref_cluster_array[A]=j
  }
  
  expr=as.matrix(data[c(2:dim(data)[1]),c(2:dim(data)[2])])
  class(expr)="numeric"
  sampleID=as.matrix(colnames(data)[2:dim(data)[2]])
  B=which(is.na(expr),arr.ind = T)
  if(dim(B)[1]>0)
  {
    selectSample=setdiff(c(1:length(sampleID)),B[,2])
    selectGene=setdiff(c(1:length(geneID)),B[,1])
    ref_cluster_array=ref_cluster_array[selectSample]
    geneID=geneID[selectGene]
    sampleID=sampleID[selectSample]
    expr=expr[selectGene,selectSample]
  }
  expr=t(expr)
  expr=scale(expr)

  list(expr=expr,
       sampleID=sampleID,
       geneID=geneID,
       ref_cluster_array=ref_cluster_array,
       ref_cluster=ref_cluster
       )
}


visualize_with_partition=function(abundance, partition,fileName,ref_cluster_array,sampling=NA)
{
  #do pca
  #abundance=expr
  #partition=boot_result$partition
  ir.pca <- prcomp(abundance,center = TRUE, scale. = TRUE) 
  x=ir.pca$x[,1]
  y=ir.pca$x[,2]
  pdf(fileName)
  
  if(is.na(sampling))
  {
    #ggplot(mtcars, aes(wt, mpg)) + geom_point(colour = "red", size = 3)
    df=data.frame(x,y,partition,ref_cluster_array)
    
    
  }else{
    tmp_partition=rep((max(ref_cluster_array)+1),length(ref_cluster_array))
    tmp_partition[sampling]=partition
    df=data.frame(x,y,partition=tmp_partition,ref_cluster_array=ref_cluster_array)
    #df=data.frame(x=x[sampling],y=y[sampling],partition=partition,ref_cluster_array=ref_cluster_array[sampling])
    #plot(x[sampling],y[sampling],col=partition,pch=ref_cluster_array[sampling])
  }
  p1=ggplot(df)+geom_point(aes(x,y,color=factor(partition),shape=factor(ref_cluster_array)))
  print(p1)
  dev.off()
  #color according to partition
}


visualize_with_partition_two=function(abundance, partition1, partition2)
{
  #do pca
  
  #color according to partition
  #one with color, one with shape
  
  
  
}





visualize_compare=function(data,boot_result1,boot_result2,prefix)
{

  expr=data$expr
  geneID=data$geneID
  sampleID=data$sampleID
  ref_cluster_array=data$ref_cluster_array
  
  print("k")
  print(length(unique(ref_cluster_array)))

  visualize_with_partition(expr,boot_result1$partition,paste(prefix,"cluster1.pdf",sep="_"),ref_cluster_array)
  for(i in c(1:dim(boot_result1$boot_partitions)[1]))
  {
    visualize_with_partition(expr,boot_result1$boot_partitions[i,],paste(prefix,i,"boot1.pdf",sep="_"),ref_cluster_array,boot_result1$boot_samplings[i,])
  }
  
  visualize_with_partition(expr,boot_result2$partition,paste(prefix,"cluster2.pdf",sep="_"),ref_cluster_array)
  for(i in c(1:dim(boot_result2$boot_partitions)[1]))
  {
    visualize_with_partition(expr,boot_result2$boot_partitions[i,],paste(prefix,i,"boot2.pdf",sep="_"),ref_cluster_array,boot_result2$boot_samplings[i,])
  }
  
  visualize_with_partition_two(boot_result1$partition,boot_result2$partition)
  
}

estimate_by_sampleTogetherAndNotTogether=function(expr,boot_result)
{
  boot_result$partition
  sample_n=dim(boot_result$boot_partitions)[2]
  #compute sample sampled together by pair 
  sampled_together_mat=compute_sample_sampled_together_mat(boot_result$boot_samplings)
  #compute sample clustered together by pair, boot#matrix (n*n), n refers to number of samples
  clustered_together_mat=compute_sample_clustered_together_mat(boot_result$boot_partitions,boot_result$boot_samplings)
  #compute sample clustered not together by pair, boot
  clustered_not_together_mat=compute_sample_not_clustered_together_mat(boot_result$boot_partitions,boot_result$boot_samplings)
  

  together_ratio_mat=clustered_together_mat/sampled_together_mat
  A=which(is.nan(together_ratio_mat),arr.ind = T)
  together_ratio_mat[A]=0
  not_together_ratio_mat=clustered_not_together_mat/sampled_together_mat
  A=which(is.nan(not_together_ratio_mat),arr.ind = T)
  not_together_ratio_mat[A]=0
  
  decision_together_mat=matrix(0,sample_n,sample_n)
  decision_notTogether_mat=matrix(0,sample_n,sample_n)
  
  A=which(together_ratio_mat>not_together_ratio_mat,arr.ind = T)
  decision_together_mat[A]=together_ratio_mat[A]
  A=which(not_together_ratio_mat>together_ratio_mat,arr.ind = T)
  decision_notTogether_mat[A]=not_together_ratio_mat[A]
  
  list(sampled_together_mat=sampled_together_mat,
       clustered_together_mat=clustered_together_mat,
       clustered_not_together_mat=clustered_not_together_mat,
       together_ratio_mat=together_ratio_mat,
       not_together_ratio_mat=not_together_ratio_mat,
       decision_together_mat=decision_together_mat,
       decision_notTogether_mat=decision_notTogether_mat)
  
}





compute_sample_not_clustered_together_mat=function(partition_mat,sampling_mat)
{
  boot_n=dim(partition_mat)[1]
  sample_n=dim(partition_mat)[2]
  result_mat=matrix(0,sample_n,sample_n)
  for(i in c(1:boot_n))
  {
    
    vector_one=rep(0,sample_n)
    sampling=sampling_mat[i,]
    vector_one[sampling]=partition_mat[i,]
    
    unique_term=unique(partition_mat[i,])
    
    
    for(unique_term_i in unique_term)
    {
      A=which(vector_one==unique_term_i,arr.ind = T)
      B=setdiff(unique(sampling),A)
      result_mat[A,B]=result_mat[A,B]+1
      #result_mat[B,A]=result_mat[B,A]+1
    }
  }
  result_mat
}

compute_sample_clustered_together_mat=function(partition_mat,sampling_mat)
{
  boot_n=dim(partition_mat)[1]
  sample_n=dim(partition_mat)[2]
  result_mat=matrix(0,sample_n,sample_n)
  for(i in c(1:boot_n))
  {
    vector_one=rep(0,sample_n)
    sampling=sampling_mat[i,]
    vector_one[sampling]=partition_mat[i,]
    
    unique_term=unique(partition_mat[i,])
    for(unique_term_i in unique_term)
    {
      A=which(vector_one==unique_term_i,arr.ind = T)
      result_mat[A,A]=result_mat[A,A]+1
    }
  }
  result_mat
}

compute_sample_sampled_together_mat=function(sampling_mat)
{
  boot_n=dim(sampling_mat)[1]
  sample_n=dim(sampling_mat)[2]
  result_mat=matrix(0,sample_n,sample_n)
  for(i in c(1:boot_n))
  {
    vector_one=sampling_mat[i,]
    unique_term=unique(vector_one)
    result_mat[unique_term,unique_term]=result_mat[unique_term,unique_term]+1
  }
  result_mat
}



###################################################
estimate_by_compareWithCompleteAnswer=function(expr, boot_result,cluster_method="average")
{
  ######for test, comment it after testing!!!!!!
  boot_result=do_sampleCluster_oneTest_bootstrapTest_forV(fileName,"root absolute","average",output_dir,i,bootstrap = 100)
  ######
  #boot_result$partition
  #boot_result$boot_samplings
  sample_n=dim(boot_result$boot_partitions)[2]
  boot_n=dim(boot_result$boot_samplings)[1]
  sampled_together_mat=compute_sample_sampled_together_mat(boot_result$boot_samplings)
  
  
  mat_samplePairConsistWithComplete_boots=matrix(0,sample_n,sample_n)
  mat_samplePairNotConsistWithComplete_boots=matrix(0,sample_n,sample_n)
  
  #score scheme, get the difference?, binary scheme, count the same?
  for(booti in c(1:boot_n))
  {
    sampling=boot_result$boot_samplings[booti,]
    hc_merge=boot_result$boot_hc_merges[,,booti]
    mat_node_in_cluster=get_internal_node_clusters(hc_merge)
    
    #for each boot, compute the number of internal node
    num_of_internal_node=dim(mat_node_in_cluster)[1]-1
    size_of_internal_cluster=apply(mat_node_in_cluster,1,sum)
    
    #for each boot, consider the internal node that a sample pair clustered together.
    for(sample_i in c(1:length(sampling)))
    {
      sample_i_classes=mat_node_in_cluster[,sample_i]
      for(sample_j in c(1:length(sampling)))
      {
        if(sample_i != sample_j)
        {
          sample_j_classes=mat_node_in_cluster[,sample_j]
          tmp=sample_i_classes+sample_j_classes
          A=which(tmp==2,arr.ind = T)
          length(A)
          #consider the size of cluster
          #score for the sample pair, consider all the internal node
          #for this pair in each internal cluster, score refers to 2/(size of cluster)
          #the score for this pair refers to the max 2/(size of cluster) it can get
          score_for_this_pair=0
          for(a in A)
          {
            if(size_of_internal_cluster[a]<sample_n)
            {
              score_for_this_pair=max(score_for_this_pair,(2/size_of_internal_cluster[a]))
            }
          }
          
          
          #score scheme, get the difference?, binary scheme, count the same?
          
          
          
        }
      }
    }
    
    
  }
}

estimate_by_sampleTogetherAndNotTogether_internalNode=function(expr,boot_result,cluster_method=
                                                                 "average")
{
  #boot_result=do_sampleCluster_oneTest_bootstrapTest_forV(fileName,"root absolute","average",output_dir,i,bootstrap = 100)
  boot_result$partition
  boot_result$boot_samplings
  sample_n=dim(boot_result$boot_partitions)[2]
  boot_n=dim(boot_result$boot_samplings)[1]
  sampled_together_mat=compute_sample_sampled_together_mat(boot_result$boot_samplings)
  
  nboot_mat_sampleTogether_oneboot=c()
  mat_sampleTogether_boots=matrix(0,sample_n,sample_n)
  mat_sampleNotTogether_boots=matrix(0,sample_n,sample_n)
  for(booti in c(1:boot_n))
  {
    mat_sampleTogether_oneboot=matrix(0,sample_n,sample_n)
    mat_sampleTogether_oneboot_without_order=matrix(0,sample_n,sample_n)
    mat_sampleNotTogether_oneboot=matrix(0,sample_n,sample_n)
    mat_sampleNotTogether_oneboot_without_order=matrix(0,sample_n,sample_n)
    sampling=boot_result$boot_samplings[booti,]
    
    hc_merge=boot_result$boot_hc_merges[,,booti]
    
    mat_node_in_cluster=get_internal_node_clusters(hc_merge)
    #compute sample sampled together by pair 
    
    
    #for each boot, compute the number of internal node
    num_of_internal_node=dim(mat_node_in_cluster)[1]-1
    size_of_internal_cluster=apply(mat_node_in_cluster,1,sum)
    
    #for each boot, consider the internal node that a sample pair clustered together.
    #
    for(sample_i in c(1:length(sampling)))
    {
      sample_i_classes=mat_node_in_cluster[,sample_i]
      for(sample_j in c(1:length(sampling)))
      {
        if(sample_i != sample_j)
        {
          sample_j_classes=mat_node_in_cluster[,sample_j]
          tmp=sample_i_classes+sample_j_classes
          A=which(tmp==2,arr.ind = T)
          length(A)#number of internal nodes that sampling[sample_i] and sampling[samplej] were clustered together 
          
          #consider the size of cluster
          #score for the sample pair, consider all the internal node
          #for this pair in each internal cluster, score refers to 2/(size of cluster)
          #the score for this pair refers to the max 2/(size of cluster) it can get
          score_for_this_pair=0
          for(a in A)
          {
            if(size_of_internal_cluster[a]<sample_n)
            {
              score_for_this_pair=max(score_for_this_pair,(2/size_of_internal_cluster[a]))
            }
          }
          
          score_for_this_pair_planB=0
          score_for_this_pair_planB=length(A)
          
          B=which(tmp!=2,arr.ind=T)
          #print(length(B))#number of internal nodes that sampling[sample_i] and sampling[samplej] were clustered together 
          
          
          #mat_sampleTogether_oneboot_without_order[sample_i,sample_j]=length(A)/num_of_internal_node
          mat_sampleTogether_oneboot_without_order[sample_i,sample_j]=score_for_this_pair
          mat_sampleNotTogether_oneboot_without_order[sample_i,sample_j]=length(B)/num_of_internal_node
          
        }
      }
    }
    
    mat_sampleTogether_oneboot[sampling,sampling]=mat_sampleTogether_oneboot_without_order
    mat_sampleNotTogether_oneboot[sampling,sampling]=mat_sampleNotTogether_oneboot_without_order
    
    ##for calculate mean value for all boot 
    mat_sampleTogether_boots=mat_sampleTogether_boots+mat_sampleTogether_oneboot
    mat_sampleNotTogether_boots=mat_sampleNotTogether_boots+mat_sampleNotTogether_oneboot
    nboot_mat_sampleTogether_oneboot=abind(nboot_mat_sampleTogether_oneboot,mat_sampleTogether_oneboot,along=3)

  }
  ##calculate mean value for all boot 
  
  mat_sampleTogether_avg=mat_sampleTogether_boots/sampled_together_mat
  mat_sampleNotTogether_avg=mat_sampleNotTogether_boots/sampled_together_mat
  #then measure the variance of mat_sampleTogether_avg in multiple iteration
  std_mat=matrix(NA,sample_n,sample_n)
  for(i in c(1:sample_n))
  {
    for(j in c(1:sample_n))
    {
      if(i != j)
      {
        score_vector=nboot_mat_sampleTogether_oneboot[i,j,]
        together_vector=compute_sample_sampled_together_vector(boot_result$boot_samplings,i,j)
        mat_sampleTogether_avg[i,j]
        
        variance_i_j=sum((score_vector[together_vector]-mat_sampleTogether_avg[i,j])^2)
        std_mat[i,j]=sqrt(variance_i_j/(length(which(together_vector,arr.ind = T))-1))
      }
      
    }
  }
  mean_of_variance_of_sampleTogetherScore_in_all_boot=mean(std_mat,na.rm=T)
  
  
  
  list(std_mat=std_mat,
       mat_sampleTogether_avg=mat_sampleTogether_avg,
       mat_sampleNotTogether_avg=mat_sampleNotTogether_avg,
       #nboot_mat_sampleTogether_oneboot=nboot_mat_sampleTogether_oneboot
       mean_of_variance_of_sampleTogetherScore_in_all_boot=mean_of_variance_of_sampleTogetherScore_in_all_boot
       )
}



dissolve_test_hcluster_internal_node=function(expr,boot_result,cluster_method=
                                                "average")
{
  
  
  completed_mat_node_in_cluster=get_internal_node_clusters(boot_result$result$merge)
  completed_num_of_internal_node=dim(completed_mat_node_in_cluster)[1]-1
  
  sample_n=dim(boot_result$boot_partitions)[2]
  boot_n=dim(boot_result$boot_samplings)[1]
  
  max_jaccard_for_each_internal_for_complete=matrix(NA,completed_num_of_internal_node,boot_n)
  
  
  nboot_mat_sampleTogether_oneboot=c()
  mat_sampleTogether_boots=matrix(0,sample_n,sample_n)
  mat_sampleNotTogether_boots=matrix(0,sample_n,sample_n)
  for(booti in c(1:boot_n))
  {
    sampling=boot_result$boot_samplings[booti,]
    
    hc_merge=boot_result$boot_hc_merges[,,booti]
    
    mat_node_in_cluster=get_internal_node_clusters(hc_merge)
    #compute sample sampled together by pair 
    
    
    #for each boot, compute the number of internal node
    num_of_internal_node=dim(mat_node_in_cluster)[1]-1
    size_of_internal_cluster=apply(mat_node_in_cluster,1,sum)
    
    for(i in c(1:completed_num_of_internal_node))
    {
      max_jaccard_for_this_internal_node=0
      for(j in c(1:num_of_internal_node))
      {
        jaccard_for_completei_bootj_node=jaccardSimilarity(mat_node_in_cluster[j,],completed_mat_node_in_cluster[i,][sampling]) 
        max_jaccard_for_this_internal_node=max(max_jaccard_for_this_internal_node,jaccard_for_completei_bootj_node)
      }
      max_jaccard_for_each_internal_for_complete[i,booti]=max_jaccard_for_this_internal_node
    }
  }
  max_jaccard_for_each_internal_for_complete
}


compute_sample_sampled_together_vector=function(sampling_mat,samplei,samplej)
{
  boot_n=dim(sampling_mat)[1]

  sample_n=dim(sampling_mat)[2]
  result_together_vector=rep(F,boot_n)
  for(i in c(1:boot_n))
  {
    vector_one=sampling_mat[i,]
    A=which(vector_one==samplei,arr.ind = T)
    B=which(vector_one==samplej,arr.ind = T)

    if(length(A)>0)
    {
      if(length(B)>0)
      {
        result_together_vector[i]=T
      }
    }
    
  }
  result_together_vector
}

get_internal_node_clusters=function(hc_merge)
{
  
  mat_node_in_cluster=matrix(0,dim(hc_merge)[1],(-min(hc_merge)))
  for(i in c(1:dim(hc_merge)[1]))
  {
    internal_cluster_index=i
    for(j in c(1:dim(hc_merge)[2]))
    {
      
      if(hc_merge[i,j]>0)
      {
        internal_node=hc_merge[i,j]
        internal_node_childs=mat_node_in_cluster[internal_node,]
        A=which(internal_node_childs==1,arr.ind = T)
        mat_node_in_cluster[internal_cluster_index,A]=1
      }else{
        leaf=-hc_merge[i,j]
        mat_node_in_cluster[internal_cluster_index,leaf]=1
      }
    }
    
  }
  mat_node_in_cluster
}





####################################################

convert_file_for_plot=function(output_dir,fileList)
{
  mat_for1_pvalues=rep(0,length(fileList))
  mat_for2_method1=c()
  mat_for2_method2=c()
  mat_for3_win=rep(0,length(fileList))
  mat_for3_lose=rep(0,length(fileList))
  for(i in c(1:length(fileList)))
  {
    filei=paste(output_dir,"file",i,"bootresult.csv",sep=" ")
    data=read.csv(filei,header = T)
    data=data[,2:dim(data)[2]]
    
    #for figure2(boxplot), convert to two matrix, dataset*boot_n

    mat_for2_method1=rbind(mat_for2_method1,data[,13])
    mat_for2_method2=rbind(mat_for2_method2,data[,14])
    
    #for figure3, win and lose times
    #calculate win times and lose times 
    #convert to two vectors, vector length 35 (value record win times for each datafile)
    A=which(data[,13]>=data[,14],arr.ind = T)
    B=which(data[,13]<data[,14],arr.ind = T)
    mat_for3_win[i]=length(A)
    mat_for3_lose[i]=length(B)
    
    #for figure1 (pvalue lines over threshold abline)
    wilcox_test_r1=wilcox.test(data[,13],data[,14])
    wilcox_test_r1$p.value
    mat_for1_pvalues[i]=wilcox_test_r1$p.value
    
  }
  wilcox_test_winRate = wilcox.test(mat_for3_win,mat_for3_lose)
  print(wilcox_test_winRate)
  print(wilcox_test_winRate$p.value)
  
  t_test_winRate = t.test(mat_for3_win,mat_for3_lose)
  print(t_test_winRate)
  
  list(
    mat_for1_pvalues=mat_for1_pvalues,
    mat_for2_method1=mat_for2_method1,
    mat_for2_method2=mat_for2_method2,
    mat_for3_win=mat_for3_win,
    mat_for3_lose=mat_for3_lose
  )
}


plot_figure3_win_vs_lose=function(root_wins,root_lose,fileName)
{
  times=c(root_wins,root_lose)
  dataset=c(1:length(root_wins),1:length(root_lose))
  type=c(rep("win",length(root_wins)),rep("lose",length(root_lose)))
  colors=c(rep('deepskyblue2',length(root_wins)),rep("firebrick3",length(root_lose)))
  df=data.frame(times,dataset,type,colors)
  
  avg_count=mean(root_wins)
  print(avg_count)
  
  pdf(fileName,width = 12,height = 12)
  p1=
    ggplot(df,aes(x=dataset,y=times,fill=type))+
    theme_light(base_size = 25)+
    scale_colour_manual(values = c("deepskyblue2", "firebrick3"))+
    geom_bar(stat="identity")+coord_flip()+xlab("dataset count")+ylab("bootstrap iteration")#+geom_segment(aes(y=avg_count,x=0,yend=avg_count,xend=length(root_wins)),linetype=3)
  print(p1)
  dev.off()
  
}

plot_figure4_win_vs_lose_vs_equal=function(root_wins,root_lose,root_equal,fileName)
{
  times=c(root_wins,root_equal,root_lose)
  dataset=c(1:length(root_wins),1:length(root_equal),1:length(root_lose))
  type=c(rep("3win",length(root_wins)),rep("2equal",length(root_equal)),rep("1lose",length(root_lose)))
  colors=c(rep('deepskyblue2',length(root_wins)),rep("grey",length(root_equal)),rep("firebrick3",length(root_lose)))
  df=data.frame(times,dataset,type,colors)
  
  avg_count=mean(root_wins)
  print(avg_count)
  
  pdf(fileName,width = 12,height = 12)
  p1=
    ggplot(df,aes(x=dataset,y=times,fill=type))+
    theme_light(base_size = 25)+
    scale_colour_manual(values = c("deepskyblue2", "grey","firebrick3"))+
    geom_bar(stat="identity")+coord_flip()+xlab("dataset count")+ylab("bootstrap iteration")#+geom_segment(aes(y=avg_count,x=0,yend=avg_count,xend=length(root_wins)),linetype=3)
  print(p1)
  dev.off()
  
}

plot_figure5_win_vs_lose_vs_equal=function(root_wins,root_lose,root_equal,fileName)
{
  times=c(root_wins,root_equal,root_lose)
  dataset=c(1:length(root_wins),1:length(root_equal),1:length(root_lose))
  type=c(rep("2win",length(root_wins)),rep("3equal",length(root_equal)),rep("1lose",length(root_lose)))
  colors=c(rep('deepskyblue2',length(root_wins)),rep("grey",length(root_equal)),rep("firebrick3",length(root_lose)))
  df=data.frame(times,dataset,type,colors)
  
  wilcox_test_dissolve=wilcox.test(root_wins,root_lose)
  print(wilcox_test_dissolve)
  print(wilcox_test_dissolve$p.value)
  
  t_test_dissolve=t.test(root_wins,root_lose)
  print(t_test_dissolve)
  print(t_test_dissolve$p.value)
  
  avg_root_win=mean(root_wins)
  print(paste("root win avg",avg_root_win))
  avg_root_lose=mean(root_lose)
  print(paste("root lose avg",avg_root_lose))
  
  pdf(fileName,width = 12,height = 12)
  p1=
    ggplot(df,aes(x=dataset,y=times,fill=type))+
    theme_light(base_size = 25)+
    scale_colour_manual(values = c("deepskyblue2", "grey","firebrick3"))+
    geom_bar(stat="identity", position=position_dodge())+xlab("iteration")+ylab("number of dataset")+
    geom_segment(aes(y=avg_root_win,x=0,yend=avg_root_win,xend=length(root_wins)),linetype=2,color=3)+
    geom_segment(aes(y=avg_root_lose,x=0,yend=avg_root_lose,xend=length(root_wins)),linetype=3,color=2)
  print(p1)
  dev.off()
  
}


plot_figure2_boxplot=function(mean_select_ref_ratio_ab,mean_select_ref_ratio_root,fileName,ylim1=0.2,ylim2=0.95,reorder_flag=T)
{
  #dataset * 100
  datasets_n=dim(mean_select_ref_ratio_ab)[1]
  boot_n=dim(mean_select_ref_ratio_ab)[2]
  
  selected_dataset=c()
  for(i in c(1:datasets_n))
  {
    wtr=wilcox.test(as.numeric(mean_select_ref_ratio_ab[i,]),as.numeric(mean_select_ref_ratio_root[i,]))
    print(wtr$p.value)
    #if(wtr$p.value<0.05)
    
    #{
      selected_dataset=c(selected_dataset,i)
    
    
      #}
  }
  print(selected_dataset)
  datasets_n=length(selected_dataset)
  mean_select_ref_ratio_ab=mean_select_ref_ratio_ab[selected_dataset,]
  mean_select_ref_ratio_root=mean_select_ref_ratio_root[selected_dataset,]
  
  
  
  x=c()
  y=c()
  type=c()
  means=apply(mean_select_ref_ratio_ab,1,mean)
  tmp=sort(means,index.return=T,decreasing = T)
  
  if(reorder_flag)
  {
    reorder_dataset=tmp$ix
  }else{
    reorder_dataset=c(1:datasets_n)
  }
  
  
  for(j in c(1:datasets_n))
  {
    
  }
  index=1
  for(j in c(1:datasets_n))
  {
    i=reorder_dataset[j]
    list1=as.numeric(mean_select_ref_ratio_root[i,])
    list2=as.numeric(mean_select_ref_ratio_ab[i,])
    
    
    
    #if((max(list1)-min(list1)>0.01) & (max(list2)-min(list2)>0.01))
    #{
      
      x=c(x,as.numeric(rep(index,boot_n)))
      y=c(y,as.numeric(mean_select_ref_ratio_root[i,]))
      type=c(type,rep("root",boot_n))
      
      
      x=c(x,as.numeric(rep(index,boot_n)))
      y=c(y,as.numeric(mean_select_ref_ratio_ab[i,]))
      type=c(type,rep("ab",boot_n))
      
      
      index=index+1
    #}
    
  }
  x=as.character.Date(x)
  pdf(fileName,width = 17,height = 15)
  df=data.frame(x,y,type)
  p1=ggplot(df, aes(x=x, y=y, fill=type)) +
    geom_boxplot(notch=F,outlier.size = 0.5)+ylim(ylim1,ylim2)+theme_light(base_size = 27)+
    xlab("dataset")+ylab("v_a * v_b")+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())
  print(p1)
  dev.off()
}

plot_figure1_pvalues=function(pvalues,fileName)
{

  log_ratio=log(pvalues)
  dataset=c(1:length(pvalues))
  A=which(pvalues<0.05,arr.ind = T)
  type=rep('not significant',length(pvalues))
  type[A]='significant'
  
  df=data.frame(log_ratio,dataset,type)
  
  pdf(fileName,width = 12,height = 12)

  
  p2=ggplot(df,aes(x=dataset,y=log_ratio))+
    theme_classic(base_size = 25)+
    ylab("log (Pvalue)")+xlab("dataset")+
    geom_point(aes(x=dataset,y=log_ratio),shape=1)+
    geom_segment(aes(x=dataset,y=log_ratio,xend=dataset,yend=0,color=type))+
    geom_segment(aes(x=0,y=log(0.05),xend=dataset,yend=log(0.05)))
  
  print(p2)
  dev.off()
  
  
}


##########################end function define####################################


###########################begin main############################################


  #for files
  #do cluster
  #do cluster result visualize
work=function(output_dir,cluster_method,distance_method_root,distance_method_ab,boot_n=20)
{
  input_dir="./CompCancer/"
  #output_dir="./visualize_test_hcluster/"
  list_file=paste(input_dir,"list.txt",sep="")
  fileList=as.matrix(read.table(list_file,header = F))
  
  
  
  
  for(i in c(1:length(fileList)))
  {
    boot_for_one_file=c()
    
    for(booti in c(1:boot_n))
    {
      
      #loop
      #print(i)
      filei=fileList[i]
      fileName=paste(input_dir,filei,sep="")
      
      boot_result_root=do_sampleCluster_oneTest_bootstrapTest_forV(fileName,distance_method_root,cluster_method,output_dir,i,bootstrap = 100)
      boot_result_ab=do_sampleCluster_oneTest_bootstrapTest_forV(fileName,distance_method_ab,cluster_method,output_dir,i,bootstrap = 100)
      
      #estimate_by_sampleTogetherAndNotTogether_internalNode(boot_result_ab$data,boot_result_ab)
      #estimate_by_sampleTogetherAndNotTogether_internalNode(boot_result_root$data,boot_result_root)
      
      
      together_result1=estimate_by_sampleTogetherAndNotTogether(boot_result_ab$data,boot_result_ab)
      together_result2=estimate_by_sampleTogetherAndNotTogether(boot_result_root$data,boot_result_root)
      
      #visualize_compare(boot_result_ab$data,boot_result_root,boot_result_ab,paste(output_dir,"file",i,sep="_"))
      
    
      #write.csv(together_result1$sampled_together_mat,file="together_result_root_sampled_together_mat",quote=F)
      #write.csv(together_result1$clustered_together_mat,file="together_result_root_clustered_together_mat",quote=F)
      #write.csv(together_result1$clustered_not_together_mat,file="together_result_root_clustered_not_together_mat",quote=F)
      #write.csv(together_result1$together_ratio_mat,file="together_result_root_together_ratio_mat",quote=F)
      #write.csv(together_result1$not_together_ratio_mat,file="together_result_root_not_together_ratio_mat",quote=F)
      
      
      #the distribution should be enrichment on the two side
      
      #for hist to see the trend, the enrich more close to 1 or 0 the better
      hist(together_result1$together_ratio_mat)
      hist(together_result1$not_together_ratio_mat)
      
      hist(together_result2$together_ratio_mat)
      hist(together_result2$not_together_ratio_mat)
      
      #the number of stable pair (together and not together, use 0.7 as threshold)
      
      A=which(together_result1$together_ratio_mat>0.7)
      B=which(together_result1$not_together_ratio_mat>0.7)
      
      C=which(together_result2$together_ratio_mat>0.7)
      D=which(together_result2$not_together_ratio_mat>0.7)
      
      length(A)
      length(B)
      length(C)
      length(D)
      
      
      #for pmax mat
      
      mat1=pmax(together_result1$together_ratio_mat,together_result1$not_together_ratio_mat)
      mat2=pmax(together_result2$together_ratio_mat,together_result2$not_together_ratio_mat)
      
      hist(mat1)
      hist(mat2)
      
      E=which(mat1>0.7,arr.ind = T)
      H=which(mat2>0.7,arr.ind = T)
      length(E)
      length(H)
      
      non_zero_index=which(mat1>0,arr.ind = T)
      non_zero_values_mat1=mat1[non_zero_index]
      median_mat1=median(non_zero_values_mat1)

      non_zero_index=which(mat2>0,arr.ind = T)
      non_zero_values_mat2=mat2[non_zero_index]
      median_mat2=median(non_zero_values_mat2)

      
      
      #for decision mat, divide consider pmax mat
      non_zero_index=which(together_result1$decision_together_mat>0,arr.ind = T)
      non_zero_values_mat1=together_result1$decision_together_mat[non_zero_index]
      v_a=median(non_zero_values_mat1)

      
      non_zero_index=which(together_result1$decision_notTogether_mat>0,arr.ind = T)
      non_zero_values_mat1=together_result1$decision_notTogether_mat[non_zero_index]
      v_b=median(non_zero_values_mat1)

      
      non_zero_index=which(together_result2$decision_together_mat>0,arr.ind = T)
      non_zero_values_mat1=together_result2$decision_together_mat[non_zero_index]
      v_c=median(non_zero_values_mat1)

      
      non_zero_index=which(together_result2$decision_notTogether_mat>0,arr.ind = T)
      non_zero_values_mat1=together_result2$decision_notTogether_mat[non_zero_index]
      v_d=median(non_zero_values_mat1)

      
      one_boot=c(length(A)/2,
                           length(B)/2,
                           length(C)/2,
                           length(D)/2,
                           length(E)/2,
                           length(H)/2,
                           median_mat1,
                           median_mat2,
                           v_a,v_b,v_c,v_d,
                            v_a*v_b,
                            v_c*v_d)
      
      boot_for_one_file=rbind(boot_for_one_file,one_boot)
      
    }
    colnames(boot_for_one_file)=c("num_boot1_together_than0.7_pairs",
                                  "num_boot1_notTogether_than0.7_pairs",
                                  "num_boot2_together_than0.7_pairs",
                                  "num_boot2_notTogether_than0.7_pairs",
                                  "num_pmax_mat1_than0.7",
                                  "num_pmax_mat1_than0.7",
                                  "median_pmax_mat1",
                                  "median_pmax_mat2",
                                  "median_non_zero_decision_together_mat1",
                                  "median_non_zero_decision_notTogether_mat1",
                                  "median_non_zero_decision_together_mat2",
                                  "median_non_zero_decision_notTogether_mat2",
                                  "median_non_zero_decision_together_mat1*itnotTogether",
                                  "median_non_zero_decision_together_mat2*itnotTogether"
                                  )
    write.csv(boot_for_one_file,file=paste(output_dir,"file",i,"bootresult.csv",sep=""),quote=F)
    
    A=which(boot_for_one_file[,13]>boot_for_one_file[,14],arr.ind = T)
    print(length(A))
    wilcox_test_r1=wilcox.test(boot_for_one_file[,13],boot_for_one_file[,14])

    print(wilcox_test_r1$p.value)
  }
  
}

tmp=function()
{
work(output_dir="./visualize_test_pam/",cluster_method = "pam",distance_method_root="root absolute",distance_method_ab="absolute",boot_n=20)
convert_data=convert_file_for_plot(output_dir="./visualize_test_pam/",fileList)
plot_figure2_boxplot(convert_data$mat_for2_method2,convert_data$mat_for2_method1,"new_boot_figure2_box_pam_pearson_noreorder.pdf",ylim1=0.2,ylim2=0.95,reorder_flag=F)
plot_figure2_boxplot(convert_data$mat_for2_method2,convert_data$mat_for2_method1,"new_boot_figure2_box_pam_pearson_reorder.pdf",ylim1=0.2,ylim2=0.95,reorder_flag=T)
plot_figure3_win_vs_lose(convert_data$mat_for3_win,convert_data$mat_for3_lose,"new_boot_figure3_pearson_pam.pdf")
plot_figure1_pvalues(convert_data$mat_for1_pvalues,"new_boot_figure1_pvalues_pearson_pam.pdf")


work(output_dir="./visualize_test_pam_spearman/",cluster_method = "pam",distance_method_root="root spearman",distance_method_ab="spearman",boot_n=20)
convert_data=convert_file_for_plot(output_dir="./visualize_test_pam_spearman/",fileList)
plot_figure2_boxplot(convert_data$mat_for2_method2,convert_data$mat_for2_method1,"new_boot_figure2_box_pam_speaman_noreorder.pdf",ylim1=0.2,ylim2=0.95,reorder_flag=F)
plot_figure2_boxplot(convert_data$mat_for2_method2,convert_data$mat_for2_method1,"new_boot_figure2_box_pam_spearman_reorder.pdf",ylim1=0.2,ylim2=0.95,reorder_flag=T)
plot_figure3_win_vs_lose(convert_data$mat_for3_win,convert_data$mat_for3_lose,"new_boot_figure3_spearman_pam.pdf")
plot_figure1_pvalues(convert_data$mat_for1_pvalues,"new_boot_figure1_pvalues_spearman_pam.pdf")

work(output_dir="./visualize_test_pam_uncentered/",cluster_method = "pam",distance_method_root="root uncentered",distance_method_ab="uncentered",boot_n=20)
convert_data=convert_file_for_plot(output_dir="./visualize_test_pam_uncentered/",fileList)
plot_figure2_boxplot(convert_data$mat_for2_method2,convert_data$mat_for2_method1,"new_boot_figure2_box_pam_uncentered_noreorder.pdf",ylim1=0.2,ylim2=0.95,reorder_flag=F)
plot_figure2_boxplot(convert_data$mat_for2_method2,convert_data$mat_for2_method1,"new_boot_figure2_box_pam_uncentered_reorder.pdf",ylim1=0.2,ylim2=0.95,reorder_flag=T)
plot_figure3_win_vs_lose(convert_data$mat_for3_win,convert_data$mat_for3_lose,"new_boot_figure3_uncentered_pam.pdf")
plot_figure1_pvalues(convert_data$mat_for1_pvalues,"new_boot_figure1_pvalues_uncentered_pam.pdf")


work(output_dir="./visualize_test_hcluster_dynamicK_pearson/",cluster_method = "hcluster_dynamicK",distance_method_root="root absolute",distance_method_ab="absolute",boot_n=20)
convert_data=convert_file_for_plot(output_dir="./visualize_test_hcluster_dynamicK_pearson/",fileList)
plot_figure2_boxplot(convert_data$mat_for2_method2,convert_data$mat_for2_method1,"new_boot_figure2_box_hcluster_dynamicK_pearson_noreorder.pdf",ylim1=0.2,ylim2=0.95,reorder_flag=F)
plot_figure2_boxplot(convert_data$mat_for2_method2,convert_data$mat_for2_method1,"new_boot_figure2_box_hcluster_dynamicK_pearson_reorder.pdf",ylim1=0.2,ylim2=0.95,reorder_flag=T)
plot_figure3_win_vs_lose(convert_data$mat_for3_win,convert_data$mat_for3_lose,"new_boot_figure3_pearson_hcluster_dynamicK.pdf")
plot_figure1_pvalues(convert_data$mat_for1_pvalues,"new_boot_figure1_pvalues_pearson_hcluster_dynamicK.pdf")
}

work_hcluster_samplePair_internal=function(output_dir,cluster_method,distance_method_root,distance_method_ab,boot_n=20)
{
  input_dir="./CompCancer/"
  #output_dir="./visualize_test_hcluster/"
  list_file=paste(input_dir,"list.txt",sep="")
  fileList=as.matrix(read.table(list_file,header = F))
  
  mean_of_std_mat_ab=c()
  mean_of_std_mat_root=c()
  
  
  for(i in c(1:length(fileList)))
  {
    boot_for_one_file=c()
    mean_of_std_vector_ab=c()
    mean_of_std_vecotr_root=c()
    
    for(booti in c(1:boot_n))
    {
      
      #loop
      #print(i)
      filei=fileList[i]
      fileName=paste(input_dir,filei,sep="")
      
      boot_result_root=do_sampleCluster_oneTest_bootstrapTest_forV(fileName,distance_method_root,cluster_method,output_dir,i,bootstrap = 100)
      boot_result_ab=do_sampleCluster_oneTest_bootstrapTest_forV(fileName,distance_method_ab,cluster_method,output_dir,i,bootstrap = 100)
      
      ab_result=estimate_by_sampleTogetherAndNotTogether_internalNode(boot_result_ab$data,boot_result_ab)
      root_result=estimate_by_sampleTogetherAndNotTogether_internalNode(boot_result_root$data,boot_result_root)
      
      
      print(" ab mean of std")
      print(mean(ab_result$std_mat,na.rm=T))
      mean_of_std_vector_ab=c(mean_of_std_vector_ab,mean(ab_result$std_mat,na.rm=T))
      print("root mean of std")
      print(mean(root_result$std_mat,na.rm=T))
      mean_of_std_vecotr_root=c(mean_of_std_vecotr_root,mean(root_result$std_mat,na.rm=T))
      
      #hist(ab_result$mat_sampleTogether_avg)
      #hist(root_result$mat_sampleTogether_avg)
      
      A=which(ab_result$mat_sampleTogether_avg>0.5,arr.ind = T)
      print("ab mean avg")
      print(mean(ab_result$mat_sampleTogether_avg[A],na.rm=T))
      print("ab selected mean")
      print(mean(ab_result$std_mat[A],na.rm=T))
      #hist(ab_result$std_mat[A])
      
      B=which(root_result$mat_sampleTogether_avg>0.5,arr.ind = T)
      print("root mean avg")
      print(mean(root_result$mat_sampleTogether_avg[B],na.rm=T))
      print("root selected mean")
      print(mean(root_result$std_mat[B],na.rm=T))
      #hist(root_result$std_mat[B])
      
      
      tmp=function(){
        for(ii in c(1:sample_n))
        {
          for(jj in c(1:sample_n))
          {
            if(ii != jj)
            {
              score_vector_root=root_result$nboot_mat_sampleTogether_oneboot[ii,jj,]
              score_vector_ab=ab_result$nboot_mat_sampleTogether_oneboot[ii,jj,]
              
              together_vector_root=compute_sample_sampled_together_vector(boot_result_root$boot_samplings,ii,jj)
              together_vector_ab=compute_sample_sampled_together_vector(boot_result_ab$boot_samplings,ii,jj)
              
              score_vector_root[together_vector_root]
              score_vector_ab[together_vector_ab]
              
              hist(score_vector_root[together_vector_root])
              hist(score_vector_ab[together_vector_ab])
              
              root_result$mat_sampleTogether_avg[ii,jj]
              ab_result$mat_sampleTogether_avg[ii,jj]
              
              #variance_i_j=sum((score_vector[together_vector]-mat_sampleTogether_avg[i,j])^2)
              #std_mat[i,j]=sqrt(variance_i_j/(length(which(together_vector,arr.ind = T))-1))
            }
            
          }
        }
        
      }
      
    }
    test=mean_of_std_vector_ab-mean_of_std_vecotr_root
    A=which(test<0,arr.ind = T)
  }
  mean_of_std_mat_ab=rbind(mean_of_std_mat_ab,mean_of_std_vector_ab)
  mean_of_std_mat_root=rbind(mean_of_std_mat_root,mean_of_std_vector_root)
  write.csv(mean_of_std_mat_ab,file=paste(output_dir,"mean_of_std_mat_ab.csv",sep=""),quote=F)
  write.csv(mean_of_std_mat_root,file=paste(output_dir,"mean_of_std_mat_root.csv",sep=""),quote=F)
  
  list(mean_of_std_mat_ab=mean_of_std_mat_ab,
       mean_of_std_mat_root=mean_of_std_mat_root)
  
}
  
tmp=function(){
r1=work_hcluster_samplePair_internal(output_dir="./visualize_test_hcluster_pearson/",cluster_method = "average",distance_method_root="root absolute",distance_method_ab="absolute",boot_n=20)
}

work_hcluster_dissolve_internal=function(output_dir,cluster_method,distance_method_root,distance_method_ab,boot_n=20)
{
  output_dir="./visualize_test_hcluster/"
  cluster_method="average"
  distance_method_root="root absolute"
  distance_method_ab="absolute"
  boot_n=20
  
  input_dir="./CompCancer/"
  list_file=paste(input_dir,"list.txt",sep="")
  fileList=as.matrix(read.table(list_file,header = F))
  
  win_count_mat=matrix(0,length(fileList),boot_n)
  lose_count_mat=matrix(0,length(fileList),boot_n)
  
  for(i in c(1:length(fileList)))
  {
    for(booti in c(1:boot_n))
    {
      filei=fileList[i]
      fileName=paste(input_dir,filei,sep="")
      
      boot_result_root=do_sampleCluster_oneTest_bootstrapTest_forV(fileName,distance_method_root,cluster_method,output_dir,i,bootstrap = 20)
      boot_result_ab=do_sampleCluster_oneTest_bootstrapTest_forV(fileName,distance_method_ab,cluster_method,output_dir,i,bootstrap = 20)
    
      jaccard_mat_root=dissolve_test_hcluster_internal_node(boot_result_root$data,boot_result_root,cluster_method="average")
      jaccard_mat_ab=dissolve_test_hcluster_internal_node(boot_result_ab$data,boot_result_ab,cluster_method="average")
      
      
      mean_root=apply(jaccard_mat_root,1,mean)
      mean_ab=apply(jaccard_mat_ab,1,mean)
      
      var_root=apply(jaccard_mat_root,1,var)
      var_ab=apply(jaccard_mat_ab,1,var)
      
      A=which(mean_root>mean_ab,arr.ind=T)
      B=which(mean_root<mean_ab,arr.ind=T)
      (mean_root-mean_ab)[A]
      (mean_root-mean_ab)[B]
      
      
      A=which(jaccard_mat_root<0.5,arr.ind = T)
      mat_dissolveTime_root=matrix(0,dim(jaccard_mat_root)[1],dim(jaccard_mat_root)[2])
      mat_dissolveTime_root[A]=1
      dissolveTime_root=apply(mat_dissolveTime_root,1,sum)
      #print("dissolveTime root")
      #print(dissolveTime_root)
      
      B=which(jaccard_mat_ab<0.5,arr.ind = T)
      mat_dissolveTime_ab=matrix(0,dim(jaccard_mat_ab)[1],dim(jaccard_mat_ab)[2])
      mat_dissolveTime_ab[B]=1
      dissolveTime_ab=apply(mat_dissolveTime_ab,1,sum)
      #print("dissolveTime ab")
      #print(dissolveTime_ab)
      
      win_count=which((dissolveTime_root-dissolveTime_ab)<0,arr.ind=T)
      print("win_count")
      print(length(win_count))
      lose_count=which((dissolveTime_root-dissolveTime_ab)>0,arr.ind=T)
      print("lose count")
      print(length(lose_count))
      
      win_count_mat[i,booti]=length(win_count)
      lose_count_mat[i,booti]=length(lose_count)
    }
  }
  
  write.csv(win_count_mat,file=paste(output_dir,"dissolve_internal_node_win_count_mat.csv"),quote=F)
  write.csv(lose_count_mat,file=paste(output_dir,"dissolve_internal_node_lose_count_mat.csv"),quote=F)
  
  #win_count_mat = read.csv(paste(output_dir,"dissolve_internal_node_win_count_mat.csv"),row.names = 1)
  #lose_count_mat = read.csv(paste(output_dir,"dissolve_internal_node_lose_count_mat.csv"),row.names = 1)
  
  test=win_count_mat-lose_count_mat
  
  A=which(test>0,arr.ind=T)
  B=which(test<0,arr.ind=T)
  
  win_files=0
  lose_files=0
  root_win=c()
  root_lose=c()
  root_equal=c()
  for(i in c(1:dim(win_count_mat)[1]))
  {
    test=win_count_mat[i,]-lose_count_mat[i,]
    A=which(test>0,arr.ind=T)
    B=which(test<0,arr.ind=T)
    D=which(test==0,arr.ind = T)
    root_win=c(root_win,length(A))
    root_lose=c(root_lose,length(B))
    root_equal=c(root_equal,length(D))
    print(paste("win times",length(A),"lose times",length(B)))
    if(length(A)>length(B))
    {
      win_files=win_files+1
    }else if(length(A)<length(B))
    {
      lose_files=lose_files+1
    }
    
  }
  plot_figure4_win_vs_lose_vs_equal(root_win,root_lose,root_equal,"dissolved_internal_node_by_file.pdf")
  plot_figure5_win_vs_lose_vs_equal(root_win,root_lose,root_equal,"dissolved_internal_node_by_file_figure5.pdf")
  

  win_iteration=0
  lose_iteration=0
  root_win=c()
  root_lose=c()
  root_equal=c()
  for(i in c(1:dim(win_count_mat)[2]))
  {
    test=win_count_mat[,i]-lose_count_mat[,i]
    A=which(test>0,arr.ind=T)
    B=which(test<0,arr.ind=T)
    D=which(test==0,arr.ind = T)
    root_win=c(root_win,length(A))
    root_lose=c(root_lose,length(B))
    root_equal=c(root_equal,length(D))
    print(paste("win times",length(A),"lose times",length(B)))
    if(length(A)>length(B))
    {
      win_iteration=win_iteration+1
    }else if(length(A)<length(B))
    {
      lose_iteration=lose_iteration+1
    }
    
  }
  
  plot_figure4_win_vs_lose_vs_equal(root_win,root_lose,root_equal,"dissolved_internal_node_by_iteration.pdf")
  plot_figure5_win_vs_lose_vs_equal(root_win,root_lose,root_equal,"dissolved_internal_node_by_iteration_figure5.pdf")
  
}
