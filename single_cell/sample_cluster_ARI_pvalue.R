library("clues")
library("kmed")
source("function_for_cluster.R")
source("function_for_distance.R")


#######the input format#########
#########file format 1###############
#expression file .txt (sep='\t'), with header (sample ID), rownames=1 (gene name), 
#ref_cluster .txt(sep='\t'), without headeer, rownames=1 (sampleID), keep the sampleID order consistent with the expr file

#########file format 1###############
#expression file .txt(sep='\t'), with header (sample ID), rownames=1 (gene name), row 1 (except header) is the ref_cluster




do_sampleCluster_oneTest=function(fileName,distance_method,cluster_method,output_dir,ith,multiple_k=FALSE,flag_normalize=F,partition_file=NA,t_expr=F)
{
  print(fileName)
  #load data and reference cluster
  
  if(!is.na(partition_file ))
  {
    
    data=read.table(fileName,header = T,sep="\t",row.names = 1)
    if(t_expr)
    {
      data=t(data)
    }
    geneID=rownames(data)
    #print(paste("geneID",geneID))
    sampleID=colnames(data)
    #print(paste("sampleID",sampleID))
    expr=as.matrix(data)
    class(expr)="numeric"
    print(paste("dim(expr)",dim(expr)))
    
    if(partition_file=='random')
    {
      ref_cluster=c(rep('1',1000), rep('2',(length(sampleID)-1000)))
      ref_cluster_array=rep(0,length(ref_cluster))
      unique_ref_cluster=unique(as.character(ref_cluster))
    }else{
      ref_cluster=read.table(partition_file,header = F,row.names = 1)
      ref_cluster=ref_cluster[,1]
      #print(paste("ref_cluster",ref_cluster))
      print(paste("len(sampleID)",length(sampleID),"len(ref_cluster)",length(ref_cluster)))
      ref_cluster_array=rep(0,length(ref_cluster))
      unique_ref_cluster=unique(as.character(ref_cluster))
    }
    
    
  }else{
    data=read.table(fileName,header = T,sep="\t")
    geneID=as.matrix(data[2:dim(data)[1],1])
    ref_cluster=c(as.matrix(data[1,c(2:dim(data)[2])]))
    ref_cluster_array=rep(0,length(ref_cluster))
    unique_ref_cluster=unique(as.character(ref_cluster))
    
    expr=as.matrix(data[c(2:dim(data)[1]),c(2:dim(data)[2])])
    class(expr)="numeric"
    sampleID=as.matrix(colnames(data)[2:dim(data)[2]])
  }
  
  for(j in c(1:length(unique_ref_cluster)))
  {
    A=which(ref_cluster==unique_ref_cluster[j],arr.ind = T)
    ref_cluster_array[A]=j
  }
  
  
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
  
  if(flag_normalize)
  {
    expr=scale(expr)
  }
  
  #compute distance
  prefix=paste(ith,distance_method,cluster_method,sep="_")
  distanceMatrices1=produceDistanceMatrices(expr,distance_method)
  summary(distanceMatrices1)
  write.table(distanceMatrices1,file=paste(output_dir,prefix,"distanceMatrix.csv",sep=""),quote=F,row.names = F,col.names = F,sep=",")
  
  #sample cluster
  #selectK(expr,cluster_method,kmax=round(sqrt(dim(expr)[1])),distanceMatrices1,paste(output_dir,prefix,"CHIndex.pdf",sep=""))
  
  if(multiple_k==F)
  {
    k=length(unique_ref_cluster)
    
    #ClusterBootstrap_fromDist(distanceMatrices, k, noise.cut = 0, bootstrap = 100, 
    #                         dissolve = .5, cluster_method, ... )
    
    memb1=doCluster(distanceMatrices1,cluster_method,k,paste(output_dir,prefix,"dendrogram.pdf",sep=""))
    write.table(cbind(as.character(sampleID),memb1),file = paste(output_dir,prefix,".csv",sep=""),quote=F,sep=",",row.names = F,col.names = F)
    
    
    #evaluation
    #ARI
    agreement_indices=adjustedRand(memb1,ref_cluster_array)
    ARI=agreement_indices[1]
    #print(prefix)
    #print(agreement_indices)
  }else{
    
    max_k = round(sqrt(length(ref_cluster_array)/2))
    #max_k=length(ref_cluster_array)-1
    ARI=c()
    for(k in seq(1,max_k))
    {
      memb1=doCluster(distanceMatrices1,cluster_method,k,paste(output_dir,prefix,"dendrogram.pdf",sep=""))
      write.table(cbind(as.character(sampleID),memb1),file = paste(output_dir,prefix,".csv",sep=""),quote=F,sep=",",row.names = F,col.names = F)
      
      
      #evaluation
      #ARI
      agreement_indices=adjustedRand(memb1,ref_cluster_array)
      ARI=c(ARI,agreement_indices[1])
      #print(prefix)
      #print(agreement_indices)
    }
    ARI=max(ARI)
  }
  print(paste("ARI",ARI))
  ARI
  #ISA
}

##############################test by ARI###################################

plot_diff=function(tmp,file="tmp6.18.pdf")
{
  library(ggplot2)
  x=c(1:length(tmp))
  y=tmp
  cate=sign(y)
  df=data.frame(x,y,cate)
  p1=ggplot(df)+geom_col(aes(x=x,y=y,fill=cate))+xlab("")+ylab("ARI(dr)-ARI(ds)")
  pdf(file)
  plot(p1)
  dev.off()
}

test_by_ARI_for_fileList=function(fileList,distance_method1="root spearman",distance_method2="square spearman",cluster_method1="hcluster_dynamicK",cluster_method2="pam",
                                  flag_normalize=F,partition_file_list=NA,t_expr=F,run_cluster_method2=T)
{
  ARI_dr_hclust=c()
  ARI_ds_hclust=c()
  ARI_dr_pam=c()
  ARI_ds_pam=c()
  ARI_mat=matrix(NA,length(fileList),4)
  colnames(ARI_mat)=c(paste(distance_method1,cluster_method1,sep="_"),
                      paste(distance_method2,cluster_method1,sep="_"),
                      paste(distance_method1,cluster_method2,sep="_"),
                      paste(distance_method2,cluster_method2,sep="_"))
  
  for(i in c(1:length(fileList)))
  {
    #loop
    filei=fileList[i]
    
    fileName=paste(input_dir,filei,sep="")
    partition_file=NA
    if(!is.na(partition_file_list))
    {
      partition_filei=partition_file_list[i]
      partition_file=paste(input_dir,partition_filei,sep="")
    }else{
      partition_file="random"
    }
    
    
    ARI=do_sampleCluster_oneTest(fileName,distance_method1,cluster_method1,output_dir,i,F,flag_normalize=flag_normalize,partition_file=partition_file,t_expr=t_expr)
    ARI_mat[i,1]=ARI[1]
    ARI=do_sampleCluster_oneTest(fileName,distance_method2,cluster_method1,output_dir,i,F,flag_normalize=flag_normalize,partition_file=partition_file,t_expr=t_expr)
    ARI_mat[i,2]=ARI[1]
    if(run_cluster_method2)
    {
      ARI=do_sampleCluster_oneTest(fileName,distance_method1,cluster_method2,output_dir,i,F,flag_normalize=flag_normalize,partition_file=partition_file,t_expr=t_expr)
      ARI_mat[i,3]=ARI[1]
      ARI=do_sampleCluster_oneTest(fileName,distance_method2,cluster_method2,output_dir,i,F,flag_normalize=flag_normalize,partition_file=partition_file,t_expr=t_expr)
      ARI_mat[i,4]=ARI[1]
    }else{
      ARI_mat[i,3]=NA
      ARI_mat[i,4]=NA
    }
  }
  write.csv(ARI_mat,file=paste(output_dir,paste("multik_ARI_mat",distance_method1,"_vs_",distance_method2,"_normalized.csv",sep=""),sep=""),quote=F)
  
  
  if(length(fileList)>2)
  {
    A=which((ARI_mat[,1]-ARI_mat[,2])>0)
    B=which((ARI_mat[,1]-ARI_mat[,2])<0)
    C=which((ARI_mat[,3]-ARI_mat[,4])>0)
    D=which((ARI_mat[,3]-ARI_mat[,4])<0)
    E=which((ARI_mat[,1]-ARI_mat[,2])!=0)
    G=which((ARI_mat[,3]-ARI_mat[,4])!=0)
    ARI_dr=ARI_mat[E,1]
    ARI_ds=ARI_mat[E,2]
    ARI_dr_hclust=c(ARI_dr_hclust,ARI_mat[E,1])
    ARI_ds_hclust=c(ARI_ds_hclust,ARI_mat[E,2])
    ARI_dr_pam=c(ARI_dr_pam,ARI_mat[G,3])
    ARI_ds_pam=c(ARI_ds_pam,ARI_mat[G,4])
    
    t.test(ARI_dr,ARI_ds)
    t.test(ARI_dr_hclust,ARI_ds_hclust)
    wilcox.test(ARI_dr_hclust,ARI_ds_hclust)#significant!!!! pavlue<0.1
    
    tmp=ARI_dr_hclust-ARI_ds_hclust
    tmp2=ARI_dr_pam-ARI_ds_pam
    
    plot_diff(tmp,"diff_ARI_dr-ds_hclust.pdf")
    plot_diff(tmp2,"diff_ARI_dr-ds_pam.pdf")
  
  }
  
}


####################end function definition##########################



#input and parameter
input_dir="./2data_human_neuron3kcell/"
list_file=paste(input_dir,"list.txt",sep="")
partition_list_file=paste(input_dir,"partition_list.txt",sep="")
output_dir="./otherDistance_ARI_2data_human_neuron_dynamichclust/"

#load data and reference cluster
fileList=as.matrix(read.table(list_file,header = F))
partition_file_list=as.matrix(read.table(partition_list_file,header = F))

#test_by_ARI_for_fileList(fileList,"root absolute","root square","hcluster_dynamicK","pam",partition_file=partition_file_list,t_expr=F,run_cluster_method2 = F)
#test_by_ARI_for_fileList(fileList,"root absolute","absolute","hcluster_dynamicK","pam",partition_file=partition_file_list,t_expr=F,run_cluster_method2 = F)

#test_by_ARI_for_fileList(fileList,"root spearman","square spearman","hcluster_dynamicK","pam",partition_file=partition_file_list,t_expr=F,run_cluster_method2 = F)
#test_by_ARI_for_fileList(fileList,"root spearman","spearman","hcluster_dynamicK","pam",partition_file=partition_file_list,t_expr=F,run_cluster_method2 = F)

#test_by_ARI_for_fileList(fileList,"root uncentered","square uncentered","hcluster_dynamicK","pam",partition_file=partition_file_list,t_expr=F,run_cluster_method2 = F)
test_by_ARI_for_fileList(fileList,"root uncentered","uncentered","hcluster_dynamicK","pam",partition_file=partition_file_list,t_expr=F,run_cluster_method2 = F)


