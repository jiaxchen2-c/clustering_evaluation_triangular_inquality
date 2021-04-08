#source and library###########################
library(ggplot2)
library(ggrepel)
library(reshape2)
source("support_function.R")

work1=function(){
  #load data from file #########################
  fileList=rep("file",8)
  label_distance_pair=rep(NA,8)
  label_cluster=rep(NA,8)
  
  fileList[1]="/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-pearson-pam-normalized/absolute_vs_root.csv"
  fileList[2]="/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-spearman-pam-normalized/absolute_vs_root.csv"
  fileList[3]="/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-uncentered-pam-normalized/absolute_vs_root.csv"
  fileList[4]="/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-pam-normalized/square_vs_root.csv"
  
  fileList[5]="/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-pearson-hcluster-normalized/absolute_vs_root.csv"
  fileList[6]="/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-spearman-hcluster-normalized/absolute_vs_root.csv"
  fileList[7]="/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-uncentered-hcluster-normalized/absolute_vs_root.csv"
  fileList[8]="/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-hcluster-normalized/square_vs_root.csv"
  
  
  label_cluster[1:4]="PAM"
  label_cluster[5:8]="Hierarchical cluster"
  label_distance_pair[1]="Pearson correlation root vs absolute"
  label_distance_pair[2]="Spearman correlation root vs absolute"
  label_distance_pair[3]="Uncentered Pearson correlation root vs absolute"
  label_distance_pair[4]="Pearson correlation root vs square"
  label_distance_pair[5]="Pearson correlation root vs absolute"
  label_distance_pair[6]="Spearman correlation root vs absolute"
  label_distance_pair[7]="Uncentered Pearson correlation root vs absolute"
  label_distance_pair[8]="Pearson correlation root vs square"
  
  
  #work############################################
  #convert from file list to mat 
  expr_fileList=get_gene_fileName()
  mat=load_and_form_mat(fileList)
  colnames(mat)=expr_fileList
  
  pdf("figure_GOenrichmentEvaluation_summary_heatmap2.pdf")
  plot_heatmap(mat,label_distance_pair,expr_fileList)
  dev.off()
  
  select_row=c(1,2,3,5,6,7)
  mat=mat[select_row,]
  label_distance_pair=label_distance_pair[select_row]
  pdf("figure_GOenrichmentEvaluation_root_vs_absolute_only_summary_heatmap2.pdf")
  plot_heatmap(mat,label_distance_pair,expr_fileList)
  dev.off()
}

work_5.6=function()
{
  #load data from file #########################
  fileList=rep("file",8)
  label_distance_pair=rep(NA,8)
  label_cluster=rep(NA,8)
  
  fileList[1]="/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-pam-normalized/square_vs_root.csv"
  fileList[2]="/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-spearman-pam-normalized-2019.5.6/square_vs_root.csv"
  fileList[3]="/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-uncentered-pam-normalized-2019.5.8/square_vs_root.csv"
  fileList[4]="/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-pam-normalized/square_vs_root.csv"
  
  fileList[5]="/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-hcluster-normalized/square_vs_root.csv"
  fileList[6]="/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-spearman-hcluster-normalized-2019.5.6/square_vs_root.csv"
  fileList[7]="/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-uncentered-hcluster-normalized-2019.5.6/square_vs_root.csv"
  fileList[8]="/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-hcluster-normalized/square_vs_root.csv"
  
  
  label_cluster[1:4]="PAM"
  label_cluster[5:8]="Hierarchical cluster"
  label_distance_pair[1]="Pearson correlation root vs square"
  label_distance_pair[2]="Spearman correlation root vs square"
  label_distance_pair[3]="Uncentered Pearson correlation root vs square"
  label_distance_pair[4]="Pearson correlation root vs square"
  label_distance_pair[5]="Pearson correlation root vs square"
  label_distance_pair[6]="Spearman correlation root vs square"
  label_distance_pair[7]="Uncentered Pearson correlation root vs square"
  label_distance_pair[8]="Pearson correlation root vs square"
  
  
  #work############################################
  #convert from file list to mat 
  expr_fileList=get_gene_fileName()
  mat=load_and_form_mat(fileList)
  colnames(mat)=expr_fileList
  

  
  select_row=c(1,2,5,6,7)
  mat=mat[select_row,]
  label_distance_pair=label_distance_pair[select_row]
  pdf("figure_GOenrichmentEvaluation_root_vs_absolute_only_summary_heatmap2_root_vs_square_normalized.pdf")
  plot_heatmap(mat,label_distance_pair,expr_fileList)
  dev.off()
}

work_5.6()

######################################define functions##################################
########################################################################################
get_gene_fileName=function()
{
  input_dir="benchmark_datasets_normalized/"
  list_file=paste(input_dir,"list.txt",sep="")
  fileList=as.matrix(read.table(list_file,header = F))
  labels=convertToID(fileList,".csv")
  labels
}


load_enrichmentEvaluationScore_from_file=function(fileName,row_index=1)
{
  data=read.table(fileName,header = T,row.names = 1,sep=",")
  data[row_index,]
}

load_and_form_mat=function(fileList)
{
  mat=c()
  for(i in c(1:length(fileList)))
  {
    score_one_experiment=load_enrichmentEvaluationScore_from_file(fileList[i],1)
    mat=rbind(mat,score_one_experiment)
  }

  A=which(is.na(mat),arr.ind = T)
  #B=which(is.nan(mat),arr.ind = T)
  mat[A]=0
  #mat[B]=0
  mat
}


plot_heatmap=function(mat,row_label,col_label)
{
  colnames(mat)=col_label
  test=mat
  test=as.data.frame(test)
  test$id=c(1:length(rownames(mat)))
  mat_melt <- melt(test, id.var="id")
  
  A=which(mat_melt[,3]!=0, arr.ind = T)
  #mat_melt=mat_melt[A,]
  
  p1=ggplot(mat_melt,aes(x=as.character(variable),y=id))+
    geom_tile(aes(fill=mat_melt$value),colour="black")+theme_light()+
    #scale_fill_gradient("ratio of together",low="steelblue",medium="white",high="firebrick")
    scale_fill_gradientn(colours = c("steelblue","white", "firebrick"),
                         values = scales::rescale(c(-1, 0, 1)),lim=c(-2,2))+
    #scale_fill_continuous(colours = c("steelblue","white", "firebrick"),limits=c(-2, 2), breaks=seq(-1,1,by=0.25))
    theme(axis.text.x = element_text(angle=75, hjust=1)
          ,axis.title.x=element_blank(),axis.title.y = element_blank())#+theme_void()
  print(p1)
  
}
