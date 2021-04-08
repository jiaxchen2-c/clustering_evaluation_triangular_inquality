#source and library
library("kmed")
library("clues")
library(ggplot2)
source("function_for_cluster.R")
source("function_for_distance.R")

#input and parameter
input_dir="./CompCancer/"
list_file=paste(input_dir,"list.txt",sep="")
output_dir="./result_kmedoids_bootstrap/"
distance_method="root absolute"
cluster_method="kmedoids"

#load data and reference cluster
fileList=as.matrix(read.table(list_file,header = F))


do_sampleCluster_oneTest_bootstrapTest=function(fileName,distance_method,cluster_method,output_dir,ith)
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
  
  #compute distance
  prefix=paste(ith,distance_method,cluster_method,sep="_")
  distanceMatrices1=produceDistanceMatrices(expr,distance_method)
  
  write.table(distanceMatrices1,file=paste(output_dir,prefix,"distanceMatrix.csv",sep=""),quote=F,row.names = F,col.names = F,sep=",")
  
  #sample cluster
  #selectK(expr,cluster_method,kmax=round(sqrt(dim(expr)[1])),distanceMatrices1,paste(output_dir,prefix,"CHIndex.pdf",sep=""))
  
  k=length(unique_ref_cluster)
  
  boot_result=ClusterBootstrap_fromDist(distanceMatrices1, k, noise.cut = 0, bootstrap = 100, 
                                        dissolve = .5, cluster_method)
  #boot_result$partition
  
  write.table(cbind(as.character(sampleID),boot_result$partition),file = paste(output_dir,prefix,".csv",sep=""),quote=F,sep=",",row.names = F,col.names = F)
  
  
  #evaluation
  #ARI
  #agreement_indices=adjustedRand(memb1,ref_cluster_array)
  #ARI=agreement_indices[1]
  #print(prefix)
  #print(paste("bootmean",boot_result$bootmean,"bootdissolved", boot_result$bootdissolved))
  
  
  #ISA
  boot_result
  
}

get_k_for_all=function(input_dir,fileList)
{
  k_for_all=rep(0,length(fileList))
  for(i in c(1:length(fileList)))
  {
    filei=fileList[i]
    fileName=paste(input_dir,filei,sep="")
    
    data=read.table(fileName,header = T,sep="\t")
    
    
    ref_cluster=c(as.matrix(data[1,c(2:dim(data)[2])]))
    unique_ref_cluster=unique(as.character(ref_cluster))

    
    k=length(unique_ref_cluster)
    
    k_for_all[i]=k
  }
  k_for_all
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
      
      x = boot_result1$clusterlist[[i]]
      y = boot_result2$clusterlist[[j]]
      
      jaccard <- jaccardSimilarity(x,y)
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

plot_figure_dissolveTimes_root_vs_ab=function(root_dissolve_times,ab_dissolve_time,fileName)
{
  times=c(root_dissolve_times,ab_dissolve_time)
  dataset=c(1:length(root_dissolve_times),1:length(ab_dissolve_time))
  type=c(rep("root",length(root_dissolve_times)),rep("ab",length(ab_dissolve_time)))
  colors=c(rep('deepskyblue2',length(root_dissolve_times)),rep("firebrick3",length(ab_dissolve_time)))
  df=data.frame(times,dataset,type,colors)
  
  wilcon_test_dissolved_times = wilcox.test(root_dissolve_times,ab_dissolve_time)
  print(wilcon_test_dissolved_times)
  
  t_test_dissolved_times = t.test(root_dissolve_times,ab_dissolve_time)
  print(t_test_dissolved_times)
  
  
  pdf(fileName,width = 12,height = 12)
  p1=
    ggplot(df,aes(x=dataset,y=times,fill=type))+
    theme_light(base_size = 25)+
    scale_colour_manual(values = c("deepskyblue2", "grey","firebrick3"))+
    geom_bar(stat="identity",position=position_dodge())+xlab("iteration")+ylab("dissolve times")
  print(p1)
  dev.off()
  
}

##############################end iteration#####################################


##############################start iteration###################################
##############################test by bootstrap###################################

work_by_file=function(root_distance_method,ab_distance_method,plot_fileName)
{
  root_win_count=0
  root_lose_count=0
  root_equal_count=0
  
  dissolvedTimes_across20IterPerFile_allFile_root=rep(0,length(fileList))
  dissolvedTimes_across20IterPerFile_allFile_ab=rep(0,length(fileList))
  
  
  for(i in c(1:length(fileList)))
  {
    print(i)
    bootmean_classes_root=c()
    bootmean_classes_ab=c()
    bootdissolved_classes_root=c()
    bootdissolved_classes_ab=c()
    
    for(boot_i in c(1:20))
    {
      #loop
      #print(i)
      filei=fileList[i]
      fileName=paste(input_dir,filei,sep="")
      boot1=do_sampleCluster_oneTest_bootstrapTest(fileName,root_distance_method,"pam",output_dir,i)
      boot2=do_sampleCluster_oneTest_bootstrapTest(fileName,ab_distance_method,"pam",output_dir,i)
      bootmean_classes_root=c(bootmean_classes_root,boot1$bootmean)
      bootmean_classes_ab=c(bootmean_classes_ab,boot2$bootmean)
      
      #bootdissolved: times the class is dissolved
      bootdissolved_classes_root=c(bootdissolved_classes_root,boot1$bootdissolved)
      bootdissolved_classes_ab=c(bootdissolved_classes_ab,boot2$bootdissolved)
  
    }
    #A=which(bootmean_classes_ab>0.5,arr.ind = T)
    #B=which(bootmean_classes_root>0.5,arr.ind = T)
    
    C=which(bootdissolved_classes_ab>40,arr.ind = T)
    D=which(bootdissolved_classes_root>40,arr.ind = T)
    
    #summary classes which dissolved many times
    length(C)
    length(D)
    
    dissolvedTimes_across20IterPerFile_allFile_root[i]=length(D)
    dissolvedTimes_across20IterPerFile_allFile_ab[i]=length(C)
    
    if(length(D)<length(C))
    {
      root_win_count=root_win_count+1
    }else if(length(D)>length(C))
    {
      root_lose_count=root_lose_count+1
    }else if(length(D)==length(C))
    {
      root_equal_count=root_equal_count+1
    }
    #length D smaller the better
  }
  
  plot_figure_dissolveTimes_root_vs_ab(dissolvedTimes_across20IterPerFile_allFile_root,dissolvedTimes_across20IterPerFile_allFile_ab,plot_fileName)
  
}


work_by_iteration=function(root_distance_method,ab_distance_method,print_fileName)
{
  root_win_count=0
  root_lose_count=0
  
  dissolvedTimes_across20IterPerFile_allFile_root=rep(0,length(20))
  dissolvedTimes_across20IterPerFile_allFile_ab=rep(0,length(20))
  
  for(boot_i in c(1:20))
  {
    print(boot_i)
    bootmean_classes_root=c()
    bootmean_classes_ab=c()
    bootdissolved_classes_root=c()
    bootdissolved_classes_ab=c()
    
    for(i in c(1:length(fileList)))
    {
      #loop
      #print(i)
      filei=fileList[i]
      fileName=paste(input_dir,filei,sep="")
      boot1=do_sampleCluster_oneTest_bootstrapTest(fileName,root_distance_method,"pam",output_dir,i)
      boot2=do_sampleCluster_oneTest_bootstrapTest(fileName,ab_distance_method,"pam",output_dir,i)
      bootmean_classes_root=c(bootmean_classes_root,boot1$bootmean)
      bootmean_classes_ab=c(bootmean_classes_ab,boot2$bootmean)
      
      #bootdissolved: times the class is dissolved
      bootdissolved_classes_root=c(bootdissolved_classes_root,boot1$bootdissolved)
      bootdissolved_classes_ab=c(bootdissolved_classes_ab,boot2$bootdissolved)

    }
    #A=which(bootmean_classes_ab>0.5,arr.ind = T)
    #B=which(bootmean_classes_root>0.5,arr.ind = T)
    
    C=which(bootdissolved_classes_ab>40,arr.ind = T)
    D=which(bootdissolved_classes_root>40,arr.ind = T)
    
    #summary classes which dissolved many times
    length(C)
    length(D)
    
    dissolvedTimes_across20IterPerFile_allFile_root[boot_i]=length(D)
    dissolvedTimes_across20IterPerFile_allFile_ab[boot_i]=length(C)
    
    if(length(D)<length(C))
    {
      root_win_count=root_win_count+1
    }else if(length(D)>length(C))
    {
      root_lose_count=root_lose_count+1
    }
    #length D smaller the better
    
  }
  
  plot_figure_dissolveTimes_root_vs_ab(dissolvedTimes_across20IterPerFile_allFile_root,dissolvedTimes_across20IterPerFile_allFile_ab,print_fileName)
  
}

work_by_iteration("root absolute","absolute","dissolve_times_pam_pearson_root_vs_ab_byIteration.pdf")
#work_by_file("root absolute","absolute","dissolve_times_pam_pearson_root_vs_ab_byFile.pdf")

#work_by_iteration("root spearman","spearman","dissolve_times_pam_spearman_root_vs_ab_byIteration.pdf")
#work_by_iteration("root uncentered","uncentered","dissolve_times_pam_uncentered_root_vs_ab_byIteration.pdf")
#work_by_file("root uncentered","uncentered","dissolve_times_pam_uncentered_root_vs_ab_byFile.pdf")

#work_by_file("root spearman","spearman","dissolve_times_pam_spearman_root_vs_ab_byFile.pdf")

