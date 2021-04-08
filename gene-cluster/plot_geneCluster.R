library(ggplot2)
library(ggrepel)

###################################################################################
##################################define functions#################################
convertToID=function(toID,sep)
{
  resultID=c()
  for(i in c(1:dim(toID)[1]))
  {
    resultID[i]=unlist(strsplit(toID[i],sep))[1]
  }
  resultID
}
##################################end define functions#############################
###################################################################################



input_dir="benchmark_datasets/"
list_file=paste(input_dir,"list.txt",sep="")
output_dir="./"
fileList=as.matrix(read.table(list_file,header = F))
labels=convertToID(fileList,".csv")

############################################################################
data=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-hcluster-result/result_mat_geneCluster_absolute_vs_root.csv",header = T, row.names = 1)
data_kmedoids=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-hcluster-result-kmedoids/absolute_vs_root.csv",header = T,row.names = 1)

plot_figure1(data,labels,paste(output_dir,"geneCluster_figure1.pdf",sep=""))
plot_figure2(data,labels,paste(output_dir,"geneCluster_figure2.pdf",sep=""))

plot_figure1(data_kmedoids,labels,paste(output_dir,"geneCluster_kmedoids_figure1.pdf",sep=""))
plot_figure2(data_kmedoids,labels,paste(output_dir,"geneCluster_kmedoids_figure2.pdf",sep=""))

############################################################################
data_checkedCHI=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-hcluster-result-checkedCHI-Nov/absolute_vs_root.csv",header = T, row.names = 1)
data_kmedoids_checkedCHI=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-hcluster-result-kmedoids-checkedCHI/absolute_vs_root.csv",header = T,row.names = 1)

plot_figure1(data_checkedCHI,labels,paste(output_dir,"geneCluster_figure1_checkedCHI.pdf",sep=""))
plot_figure2(data_checkedCHI,labels,paste(output_dir,"geneCluster_figure2_checkedCHI.pdf",sep=""))

plot_figure1(data_kmedoids_checkedCHI,labels,paste(output_dir,"geneCluster_kmedoids_figure1_checkedCHI.pdf",sep=""))
plot_figure2(data_kmedoids_checkedCHI,labels,paste(output_dir,"geneCluster_kmedoids_figure2_checkedCHI.pdf",sep=""))

############################################################################
data_spearman_hcluster=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-hcluster-result-spearman-hcluster/absolute_vs_root.csv",header = T, row.names = 1)
data_uncentered_hcluster=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-hcluster-result-uncentered-hcluster/absolute_vs_root.csv",header = T,row.names = 1)

plot_figure1(data_spearman_hcluster[,2:17],labels[2:17],paste(output_dir,"geneCluster_figure1_spearman_hcluster_600.pdf",sep=""),600,600)

plot_figure1(data_spearman_hcluster,labels,paste(output_dir,"geneCluster_figure1_spearman_hcluster_1000.pdf",sep=""),1000,1000)
plot_figure2(data_spearman_hcluster,labels,paste(output_dir,"geneCluster_figure2_spearman_hcluster.pdf",sep=""),-1.7,1.7)

plot_figure1(data_uncentered_hcluster,labels,paste(output_dir,"geneCluster_kmedoids_figure1_uncentered_hcluster.pdf",sep=""))
plot_figure2(data_uncentered_hcluster,labels,paste(output_dir,"geneCluster_kmedoids_figure2_uncentered_hcluster.pdf",sep=""),-3,3)

###############################################################################################
data_compareSquare_hcluster=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-hcluster-result-compareSquare-average/square_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_compareSquare_hcluster,labels,paste(output_dir,"geneCluster_kmedoids_figure1_compareSquare_hcluster.pdf",sep=""))
plot_figure2(data_compareSquare_hcluster,labels,paste(output_dir,"geneCluster_kmedoids_figure2_compareSquare_hcluster.pdf",sep=""),-3,3)


###############################################################################################
data_pearson_kmedoids1121=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-hcluster-result-pearson-kmedoids/absolute_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_pearson_kmedoids1121,labels,paste(output_dir,"geneCluster_figure1_pearson_kmedoids1121.pdf",sep=""))
plot_figure2(data_pearson_kmedoids1121,labels,paste(output_dir,"geneCluster_figure2_pearson_kmedoids1121.pdf",sep=""))

data_uncentered_kmedoids1121=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-hcluster-result-uncentered-kmedoids/absolute_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_uncentered_kmedoids1121,labels,paste(output_dir,"geneCluster_figure1_uncentered_kmedoids1121.pdf",sep=""))
plot_figure2(data_uncentered_kmedoids1121,labels,paste(output_dir,"geneCluster_figure2_uncentered_kmedoids1121.pdf",sep=""))

###############################################################################################
data_pearson_kmedoids1123=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-hcluster-result-pearson-kmedoids-gap/absolute_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_pearson_kmedoids1123,labels,paste("geneCluster_figure1_pearson_kmedoids_gap1123.pdf",sep=""))
plot_figure2(data_pearson_kmedoids1123,labels,paste("geneCluster_figure2_pearson_kmedoids_gap1123.pdf",sep=""))


###############################################################################################

data_pearson_kmedoids1123=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-hcluster-result-pearson-kmedoids-gap-samewithab/absolute_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_pearson_kmedoids1123,labels,paste("geneCluster_figure1_pearson_kmedoids_gap_samewithab_1125.pdf",sep=""))
plot_figure2(data_pearson_kmedoids1123,labels,paste("geneCluster_figure2_pearson_kmedoids_gap_samewithab_1125.pdf",sep=""))

data_pearson_kmedoids1123=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-hcluster-result-pearson-kmedoids-gap-samewithroot/absolute_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_pearson_kmedoids1123,labels,paste("geneCluster_figure1_pearson_kmedoids_gap_samewithroot1125.pdf",sep=""))
plot_figure2(data_pearson_kmedoids1123,labels,paste("geneCluster_figure2_pearson_kmedoids_gap_samewithroot1125.pdf",sep=""))

###################################################################################################
data_pearson_kmedoids1123=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-pearson-pam/absolute_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_pearson_kmedoids1123,labels,paste("geneCluster_figure1_pearson_pam.pdf",sep=""))
plot_figure2(data_pearson_kmedoids1123,labels,paste("geneCluster_figure2_pearson_pam.pdf",sep=""))

data_pearson_kmedoids1123=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-spearman-pam/absolute_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_pearson_kmedoids1123,labels,paste("geneCluster_figure1_spearman_pam.pdf",sep=""))
plot_figure2(data_pearson_kmedoids1123,labels,paste("geneCluster_figure2_spearman_pam.pdf",sep=""))

data_pearson_kmedoids1123=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-uncentered-pam/absolute_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_pearson_kmedoids1123,labels,paste("geneCluster_figure1_uncentered_pam.pdf",sep=""))
plot_figure2(data_pearson_kmedoids1123,labels,paste("geneCluster_figure2_uncentered_pam.pdf",sep=""),-2,2)

###################################################################################################
input_dir="benchmark_datasets_normalized/"
list_file=paste(input_dir,"list.txt",sep="")
output_dir="./"
fileList=as.matrix(read.table(list_file,header = F))
labels=convertToID(fileList,".csv")

data_pearson_kmedoids1123=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-pearson-pam-normalized/absolute_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_pearson_kmedoids1123,labels,paste("geneCluster_figure1_pearson_pam_normalized.pdf",sep=""),850,850)
plot_figure2(data_pearson_kmedoids1123,labels,paste("geneCluster_figure2_pearson_pam_normalized.pdf",sep=""),-1.5,1.5)

data_pearson_kmedoids1123=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-spearman-pam-normalized/absolute_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_pearson_kmedoids1123,labels,paste("geneCluster_figure1_spearman_pam_normalized.pdf",sep=""),1000,1000)
plot_figure2(data_pearson_kmedoids1123,labels,paste("geneCluster_figure2_spearman_pam_normalized.pdf",sep=""))

data_pearson_kmedoids1123=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-uncentered-pam-normalized/absolute_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_pearson_kmedoids1123,labels,paste("geneCluster_figure1_uncentered_pam_normalized.pdf",sep=""))
plot_figure2(data_pearson_kmedoids1123,labels,paste("geneCluster_figure2_uncentered_pam_normalized.pdf",sep=""))


data_pearson_kmedoids1123=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-pearson-hcluster-normalized/absolute_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_pearson_kmedoids1123,labels,paste("geneCluster_figure1_pearson_hcluster_normalized.pdf",sep=""),850,850)
plot_figure2(data_pearson_kmedoids1123,labels,paste("geneCluster_figure2_pearson_hcluster_normalized.pdf",sep=""))

data_pearson_kmedoids1123=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-spearman-hcluster-normalized/absolute_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_pearson_kmedoids1123,labels,paste("geneCluster_figure1_spearman_hcluster_normalized.pdf",sep=""),750,750)
plot_figure2(data_pearson_kmedoids1123,labels,paste("geneCluster_figure2_spearman_hcluster_normalized.pdf",sep=""))

data_pearson_kmedoids1123=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-uncentered-hcluster-normalized/absolute_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_pearson_kmedoids1123,labels,paste("geneCluster_figure1_uncentered_hcluster_normalized.pdf",sep=""))
plot_figure2(data_pearson_kmedoids1123,labels,paste("geneCluster_figure2_uncentered_hcluster_normalized.pdf",sep=""),-1.8,1.8)

data_pearson_kmedoids1123=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-hcluster-normalized/square_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_pearson_kmedoids1123,labels,paste("geneCluster_figure1_compareSquare_hcluster_normalized.pdf",sep=""))
plot_figure2(data_pearson_kmedoids1123,labels,paste("geneCluster_figure2_compareSquare_hcluster_normalized.pdf",sep=""))


data_pearson_kmedoids1123=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-pam-normalized/square_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_pearson_kmedoids1123,labels,paste("geneCluster_figure1_compareSquare_pam_normalized.pdf",sep=""),1250,1250)
plot_figure2(data_pearson_kmedoids1123,labels,paste("geneCluster_figure2_compareSquare_pam_normalized.pdf",sep=""),-1.7,1.7)

######following added in 2019.5.6#################
input_dir="benchmark_datasets_normalized/"
list_file=paste(input_dir,"list.txt",sep="")
output_dir="./"
fileList=as.matrix(read.table(list_file,header = F))
labels=convertToID(fileList,".csv")


data_pearson_kmedoids1123=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-spearman-pam-normalized-2019.5.6/square_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_pearson_kmedoids1123,labels,paste("geneCluster_figure1_compareSquare_spearman_pam_normalized.pdf",sep=""))
plot_figure2(data_pearson_kmedoids1123,labels,paste("geneCluster_figure2_compareSquare_spearman_pam_normalized.pdf",sep=""),-1.7,1.7)

data_pearson_kmedoids1123=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-uncentered-pam-normalized-2019.5.8/square_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_pearson_kmedoids1123,labels,paste("geneCluster_figure1_compareSquare_uncentered_pam_normalized.pdf",sep=""))
plot_figure2(data_pearson_kmedoids1123,labels,paste("geneCluster_figure2_compareSquare_uncentered_pam_normalized.pdf",sep=""))

data_pearson_kmedoids1123=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-spearman-hcluster-normalized-2019.5.6/square_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_pearson_kmedoids1123,labels,paste("geneCluster_figure1_compareSquare_spearman_hcluster_normalized.pdf",sep=""))
plot_figure2(data_pearson_kmedoids1123,labels,paste("geneCluster_figure2_compareSquare_spearman_hcluster_normalized.pdf",sep=""))

data_pearson_kmedoids1123=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-uncentered-hcluster-normalized-2019.5.6/square_vs_root.csv",header = T,row.names = 1)
plot_figure1(data_pearson_kmedoids1123,labels,paste("geneCluster_figure1_compareSquare_uncentered_hcluster_normalized.pdf",sep=""))
plot_figure2(data_pearson_kmedoids1123,labels,paste("geneCluster_figure2_compareSquare_uncentered_hcluster_normalized.pdf",sep=""))



plot_figure1=function(data,labels,fileName,xend=600,yend=600)
{
  times_r1_smaller=as.numeric(data[2,])
  times_r2_smaller=as.numeric(data[3,])
  ratio=as.numeric(data[1,])
  root_wins=which(ratio<0,arr.ind = T)
  wins=rep("firebrick3",dim(data)[2])
  wins[root_wins]="deepskyblue2"
  
  df=as.data.frame(t(data),labels,wins)

  #r2 is root
  
  pdf(fileName,width = 10,height = 10)
  p1=ggplot(df)+coord_fixed()+geom_segment(aes(x=0,y=0,xend=xend,yend=yend))+
    geom_point(aes(x=times_r1_smaller,y=times_r2_smaller))+
    geom_text_repel(aes(x=times_r1_smaller,y=times_r2_smaller),label=labels,color=wins)+
    xlab("number of GO term r1 have smaller pvalue")+theme_light(base_size = 15)+
    ylab("number of GO term r2 have smaller pvalue")
    
  print(p1)
  dev.off()
}


plot_figure2=function(data,labels,fileName,ylim_left=-1.6,ylim_right=1.6)
{
  times_r1_smaller=as.numeric(data[2,])
  times_r2_smaller=as.numeric(data[3,])
  ratio=as.numeric(data[1,])
  text_y=ratio+0.05
  root_wins=which(ratio<0,arr.ind = T)
  wins=rep("firebrick3",dim(data)[2])
  wins[root_wins]="deepskyblue2"
  text_y[root_wins]=ratio[root_wins]-0.05
  dataset=as.numeric(c(1:length(labels)))
  df=as.data.frame(t(data),labels,wins,dataset,text_y)
  
  pdf(fileName,width = 7,height = 7)
  p1=ggplot(df)+geom_segment(aes(x=dataset,y=0,xend=dataset,yend=ratio),size=10,color=wins)+
    ylim(ylim_left,ylim_right)+theme_light(base_size = 15)+
    ylab("log ratio")+coord_flip()
    #theme(axis.title.y=element_blank(),axis.text.y=element_blank(),
     #                    axis.ticks.y=element_blank())
    
  print(p1)
  dev.off()
}



plot_figure3=function(GO_terms_mat1,GO_terms_mat2)
{
  #compute enrichment quality, calculate fold change?
  list1=c()
  list2=c()
  terms1=GO_terms_mat1[,1]
  pvalues1=GO_terms_mat1[,2]
  unique_terms1=unique(terms1)
  
  
  terms2=GO_terms_mat2[,1]
  pvalues2=GO_terms_mat2[,2]
  unique_terms2=unique(terms2)
  
  overlap_terms=intersect(unique_terms1,unique_terms2)
  
  times_r1_smaller=0
  times_r2_smaller=0
  for(i in c(1:length(overlap_terms)))
  {
    A=which(terms1==overlap_terms[i])
    B=which(terms2==overlap_terms[i])
    list1=c(list1,min(pvalues1[A]))
    
    
    list2=c(list2,min(pvalues2[B]))
    
  }
  

  logp1=log(list1)
  logp2=log(list2)
  tmp2=sort(logp1,index.return=T)
  logp1=tmp2$x
  logp2=logp2[tmp2$ix]
  logp=c(logp1,logp2)
  
  pvalues=as.numeric(c(list1[tmp2$ix],list2[tmp2$ix]))
  
  group=c(rep("ab",length(list1)),rep("root",length(list2)))
  terms=c(c(1:length(list1)),c(1:length(list2)))
  df=as.data.frame(cbind(pvalues,logp,group,terms))
  
  ggplot(df)+geom_point(aes(y=terms,x=logp,color=group))+theme_light(base_size = 15)+
    theme(axis.text.y=element_blank(),
          axis.text.x=element_blank())
}
