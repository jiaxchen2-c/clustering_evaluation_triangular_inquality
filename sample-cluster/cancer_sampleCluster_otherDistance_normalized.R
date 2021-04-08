#Before each clustering experiment, the expression values of each gene was standardized
#(mean 0, variance 1) across all samples being considered in the clustering. 

#source and library
library("clues")
library("kmed")
source("function_for_cluster.R")
source("function_for_distance.R")


#input and parameter
input_dir="./CompCancer/"
list_file=paste(input_dir,"list.txt",sep="")
output_dir="./otherDistance_ARI/"
distance_method="root absolute"
cluster_method="single"

#load data and reference cluster
fileList=as.matrix(read.table(list_file,header = F))


#######
do_sampleCluster_oneTest_withNormalize=function(fileName,distance_method,cluster_method,output_dir,ith)
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
  
  #normalize2
  #expr=scale(expr)
  expr=t(expr)
  
  #do normalize here
  #normalize1
  expr=scale(expr)
  
  #compute distance
  prefix=paste(ith,distance_method,cluster_method,sep="_")
  distanceMatrices1=produceDistanceMatrices(expr,distance_method)
  
  write.table(distanceMatrices1,file=paste(output_dir,prefix,"distanceMatrix.csv",sep=""),quote=F,row.names = F,col.names = F,sep=",")
  
  #sample cluster
  #selectK(expr,cluster_method,kmax=round(sqrt(dim(expr)[1])),distanceMatrices1,paste(output_dir,prefix,"CHIndex.pdf",sep=""))
  
  k=length(unique_ref_cluster)
  
  #ClusterBootstrap_fromDist(distanceMatrices, k, noise.cut = 0, bootstrap = 100, 
  #                         dissolve = .5, cluster_method, ... )
  
  memb1=doCluster(distanceMatrices1,cluster_method,k,paste(output_dir,prefix,"dendrogram.pdf",sep=""))
  write.table(cbind(as.character(sampleID),memb1),file = paste(output_dir,prefix,".csv",sep=""),quote=F,sep=",",row.names = F,col.names = F)
  
  
  #evaluation
  #ARI
  agreement_indices=adjustedRand(memb1,ref_cluster_array)
  ARI=agreement_indices[1]
  print(prefix)
  print(agreement_indices)
  
  ARI
  #ISA
}



##############################start iteration###################################
##############################test by ARI###################################
ARI_mat=matrix(NA,length(fileList),4)
colnames(ARI_mat)=c("root_average","ab_average","root_pam","ab_pam")
for(i in c(1:length(fileList)))
{
  #loop
  filei=fileList[i]
  fileName=paste(input_dir,filei,sep="")
  
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"root absolute","average",output_dir,i)
  ARI_mat[i,1]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"absolute","average",output_dir,i)
  ARI_mat[i,2]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"root absolute","pam",output_dir,i)
  ARI_mat[i,3]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"absolute","pam",output_dir,i)
  ARI_mat[i,4]=ARI[1]
}
write.csv(ARI_mat,file=paste(output_dir,"ARI_mat_pearson_hclusterAndPam_normalized.csv",sep=""),quote=F)
##############################################################################
##############################test by ARI###################################
ARI_mat=matrix(NA,length(fileList),4)
colnames(ARI_mat)=c("root_spearman_average","ab_spearman_average","root_spearman_pam","ab_spearman_pam")
for(i in c(1:length(fileList)))
{
  #loop
  filei=fileList[i]
  fileName=paste(input_dir,filei,sep="")
  
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"root spearman","average",output_dir,i)
  ARI_mat[i,1]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"spearman","average",output_dir,i)
  ARI_mat[i,2]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"root spearman","pam",output_dir,i)
  ARI_mat[i,3]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"spearman","pam",output_dir,i)
  ARI_mat[i,4]=ARI[1]
}
write.csv(ARI_mat,file=paste(output_dir,"ARI_mat_spearman_hclusterAndPam_normalized.csv",sep=""),quote=F)
##############################################################################
##############################test by ARI###################################
ARI_mat=matrix(NA,length(fileList),4)
colnames(ARI_mat)=c("root_uncentered_average","ab_uncentered_average","root_uncentered_pam","ab_uncentered_pam")
for(i in c(1:length(fileList)))
{
  #loop
  filei=fileList[i]
  fileName=paste(input_dir,filei,sep="")
  
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"root uncentered","average",output_dir,i)
  ARI_mat[i,1]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"uncentered","average",output_dir,i)
  ARI_mat[i,2]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"root uncentered","pam",output_dir,i)
  ARI_mat[i,3]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"uncentered","pam",output_dir,i)
  ARI_mat[i,4]=ARI[1]
}
write.csv(ARI_mat,file=paste(output_dir,"ARI_mat_uncentered_hclusterAndPam_normalized.csv",sep=""),quote=F)
##############################################################################

##############################start iteration###################################
##############################test by ARI###################################
ARI_mat=matrix(NA,length(fileList),4)
colnames(ARI_mat)=c("root_average","square_average","root_pam","square_pam")
for(i in c(1:length(fileList)))
{
  #loop
  filei=fileList[i]
  fileName=paste(input_dir,filei,sep="")
  
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"root absolute","average",output_dir,i)
  ARI_mat[i,1]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"root square","average",output_dir,i)
  ARI_mat[i,2]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"root absolute","pam",output_dir,i)
  ARI_mat[i,3]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"root square","pam",output_dir,i)
  ARI_mat[i,4]=ARI[1]
}
write.csv(ARI_mat,file=paste(output_dir,"ARI_mat_pearson_compareSquare_hclusterAndPam_normalized.csv",sep=""),quote=F)

##############################################################################################
######################following add in 2019.4.24##############################################
########################################################################################

##############################start iteration###################################
##############################test by ARI###################################
ARI_mat=matrix(NA,length(fileList),4)
colnames(ARI_mat)=c("root_average","square_average","root_pam","square_pam")
for(i in c(1:length(fileList)))
{
  #loop
  filei=fileList[i]
  fileName=paste(input_dir,filei,sep="")
  
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"root absolute","average",output_dir,i)
  ARI_mat[i,1]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"root square","average",output_dir,i)
  ARI_mat[i,2]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"root absolute","pam",output_dir,i)
  ARI_mat[i,3]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"root square","pam",output_dir,i)
  ARI_mat[i,4]=ARI[1]
}
write.csv(ARI_mat,file=paste(output_dir,"ARI_mat_pearson_compareSquare_hclusterAndPam_normalized_checkAgain.csv",sep=""),quote=F)
##############################start iteration###################################
##############################test by ARI###################################
ARI_mat=matrix(NA,length(fileList),4)
colnames(ARI_mat)=c("root_spearman_average","square_spearman_average","root_spearman_pam","square_spearman_pam")
for(i in c(1:length(fileList)))
{
  #loop
  filei=fileList[i]
  fileName=paste(input_dir,filei,sep="")
  
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"root spearman","average",output_dir,i)
  ARI_mat[i,1]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"square spearman","average",output_dir,i)
  ARI_mat[i,2]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"root spearman","pam",output_dir,i)
  ARI_mat[i,3]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"square spearman","pam",output_dir,i)
  ARI_mat[i,4]=ARI[1]
}
write.csv(ARI_mat,file=paste(output_dir,"ARI_mat_spearman_compareSquare_hclusterAndPam_normalized_2019.4.24.csv",sep=""),quote=F)

##############################start iteration###################################
##############################test by ARI###################################
ARI_mat=matrix(NA,length(fileList),4)
colnames(ARI_mat)=c("root_uncentered_average","square_uncentered_average","root_uncentered_pam","square_uncentered_pam")
for(i in c(1:length(fileList)))
{
  #loop
  filei=fileList[i]
  fileName=paste(input_dir,filei,sep="")
  
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"root uncentered","average",output_dir,i)
  ARI_mat[i,1]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"square uncentered","average",output_dir,i)
  ARI_mat[i,2]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"root uncentered","pam",output_dir,i)
  ARI_mat[i,3]=ARI[1]
  ARI=do_sampleCluster_oneTest_withNormalize(fileName,"square uncentered","pam",output_dir,i)
  ARI_mat[i,4]=ARI[1]
}
write.csv(ARI_mat,file=paste(output_dir,"ARI_mat_uncentered_compareSquare_hclusterAndPam_normalized_2019.4.24.csv",sep=""),quote=F)


