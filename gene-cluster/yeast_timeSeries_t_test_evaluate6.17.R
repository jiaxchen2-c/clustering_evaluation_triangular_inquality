#source and library
library(dplyr)
library(tidyr)
library(ggplot2)
library(kmed)
library(clues)

library(DBI)
library(AnnotationDbi)
library(Biobase)
library(org.Hs.eg.db)
library(Category)
library(graph)
library(GOstats)
library(GO.db)
library(KEGG.db)
library(Category)

source("function_for_cluster.R")
source("function_for_distance.R")
source("function_for_enrichment.R")
source("function_for_evaluateDistanceByEnrichment.R")



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

#parameter and load data
input_dir="benchmark_datasets_normalized/"
output_dir="timeCourse-hcluster-result-kmedoids-11_7/"
list_file=paste(input_dir,"list.txt",sep="")
toID=as.matrix(read.table("geneName2.txt",header=F))
fromID=as.character(as.matrix(read.table("OLN.txt",header=F)))
toID=convertToID(toID,';')

#fileList
fileList=as.matrix(read.table(list_file,header = F))
root_win=0
root_lose=0


absolute_vs_root_list=c()
#checkAgain=c(3,7,9,10,12,13,14,16)
checkAgain=c(9)#c(13,14,15,17)#remove 9 temp
runAgain=c(16,17)


###########################################################################
###########################################################################

plot_select_CHIndex=function(distance,cluster_method,output_dir)
{
  for(i in c(1:length(fileList)))
  {
    print(paste("i",i))
    #loop
    #file1="01_alpha_factor.csv"
    filei=fileList[i]
    data=read.table(paste(input_dir,filei,sep=""),header = T,sep="\t")
    geneName=as.matrix(data[,1])
    geneName=convertToID(geneName,' ')
    expr=data.matrix(data[,2:dim(data)[2]])
    
    do_select_k(expr,output_dir,distance,cluster_method,filei)
  }
}

plot_select_gap=function(distance,cluster_method,output_dir)
{
  bestK_mat=c()
  for(i in c(1:length(fileList)))
  {
    print(paste("i",i))
    #loop
    #file1="01_alpha_factor.csv"
    filei=fileList[i]
    data=read.table(paste(input_dir,filei,sep=""),header = T,sep="\t")
    geneName=as.matrix(data[,1])
    geneName=convertToID(geneName,' ')
    expr=data.matrix(data[,2:dim(data)[2]])
    
    
    
    bestK=do_select_k_gap(expr,output_dir,distance,cluster_method,filei)
    bestK_mat=rbind(bestK_mat,c(i,bestK))
  }
  
  bestK_mat
}

##############################################################
##############################################################



##############################################################
##############################################################
root_win=0
root_lose=0
absolute_vs_root_list=c()
output_dir="timeCourse-result-uncentered-pam-normalized_6.17/"
#plot_select_CHIndex("uncentered","pam",output_dir)
#plot_select_CHIndex("root uncentered","pam",output_dir)
#ab, root

k_table=matrix(c(12,4,#modified 11.28
                 3,3,
                 3,3,
                 5,5,
                 6,6,
                 3,3,
                 11,9,
                 4,4,
                 10,12,
                 10,10,
                 4,4,
                 5,5,
                 5,7,
                 6,9,
                 5,4,
                 6,5
),nrow=2,ncol=16)

for(i in c(1:length(fileList))) 
{
  
  filei=fileList[i]
  data=read.table(paste(input_dir,filei,sep=""),header = T,sep="\t")
  geneName=as.matrix(data[,1])
  geneName=convertToID(geneName,' ')
  expr=data.matrix(data[,2:dim(data)[2]])
  
  #############one distance and cluster method##################
  
  absolute_vs_root=do_compare_distance_in_uncentered_pam(expr,output_dir,filei,k_table[1,i],k_table[2,i])
  absolute_vs_root_list=cbind(absolute_vs_root_list,absolute_vs_root)
  
  write.csv(absolute_vs_root_list,file=paste(output_dir,"absolute_vs_root.csv",sep=""),quote=F)
  ##############################################################
  
}


##############################################################
##############################################################
root_win=0
root_lose=0
absolute_vs_root_list=c()
output_dir="timeCourse-result-pearson-pam-normalized_6.17/"
#plot_select_CHIndex("absolute","pam",output_dir)
#plot_select_CHIndex("root absolute","pam",output_dir)


k_table=matrix(c(8,4,#modified 11.28
                 12,12,
                 5,5,
                 3,5,
                 12,16,
                 11,4,
                 9,9,
                 15,16,
                 16,25,
                 19,20,
                 5,5,
                 6,4,
                 11,10,
                 13,17,
                 6,8,
                 10,18
),nrow=2,ncol=16)

for(i in c(1:length(fileList))) 
{
  
  filei=fileList[i]
  data=read.table(paste(input_dir,filei,sep=""),header = T,sep="\t")
  geneName=as.matrix(data[,1])
  geneName=convertToID(geneName,' ')
  expr=data.matrix(data[,2:dim(data)[2]])
  
  #############one distance and cluster method##################
  
  absolute_vs_root=do_compare_distance_in_pam(expr,output_dir,filei,k_table[1,i],k_table[2,i])
  absolute_vs_root_list=cbind(absolute_vs_root_list,absolute_vs_root)
  
  write.csv(absolute_vs_root_list,file=paste(output_dir,"absolute_vs_root.csv",sep=""),quote=F)
  ##############################################################
}




##############################################################
##############################################################
root_win=0
root_lose=0
absolute_vs_root_list=c()
output_dir="timeCourse-result-spearman-pam-normalized_6.17/"
#plot_select_CHIndex("spearman","pam",output_dir)
#plot_select_CHIndex("root spearman","pam",output_dir)

k_table=matrix(c(7,6,#modified 11.28
                 12,18,
                 7,8,
                 7,7,
                 8,10,
                 5,7,
                 4,5,
                 9,12,
                 5,17,
                 9,20,
                 8,9,
                 4,7,
                 7,4,
                 6,9,
                 5,5,
                 8,6
),nrow=2,ncol=16)



for(i in c(1:length(fileList))) 
{
  
  filei=fileList[i]
  data=read.table(paste(input_dir,filei,sep=""),header = T,sep="\t")
  geneName=as.matrix(data[,1])
  geneName=convertToID(geneName,' ')
  expr=data.matrix(data[,2:dim(data)[2]])
  
  #############one distance and cluster method##################
  
  absolute_vs_root=do_compare_distance_in_spearman_pam(expr,output_dir,filei,k_table[1,i],k_table[2,i])
  absolute_vs_root_list=cbind(absolute_vs_root_list,absolute_vs_root)
  
  write.csv(absolute_vs_root_list,file=paste(output_dir,"absolute_vs_root.csv",sep=""),quote=F)
  ##############################################################
}





##############################################################
##############################################################
root_win=0
root_lose=0
absolute_vs_root_list=c()
output_dir="timeCourse-result-uncentered-hcluster-normalized_6.17/"
#plot_select_CHIndex("uncentered","average",output_dir)
#plot_select_CHIndex("root uncentered","average",output_dir)
#ab, root

k_table=matrix(c(20,24,#modified 11.28
                 13,9,
                 17,22,
                 20,21,
                 7,6,
                 11,12,
                 14,27,
                 11,11,
                 7,8,
                 4,5,
                 12,9,
                 7,7,
                 13,18,
                 20,20,
                 14,10,
                 13,14
),nrow=2,ncol=16)

for(i in c(1:length(fileList))) 
{
  
  filei=fileList[i]
  data=read.table(paste(input_dir,filei,sep=""),header = T,sep="\t")
  geneName=as.matrix(data[,1])
  geneName=convertToID(geneName,' ')
  expr=data.matrix(data[,2:dim(data)[2]])
  
  #############one distance and cluster method##################
  GO_terms_absolute_average=do_one_cluster(expr,output_dir,'uncentered','average',filei,k_table[1,i])
  GO_terms_root_absolute_average=do_one_cluster(expr,output_dir,'root uncentered','average',filei,k_table[2,i])
  absolute_vs_root=compareMethod(GO_terms_absolute_average,GO_terms_root_absolute_average)
  
  
  absolute_vs_root_list=cbind(absolute_vs_root_list,absolute_vs_root)
  
  write.csv(absolute_vs_root_list,file=paste(output_dir,"absolute_vs_root.csv",sep=""),quote=F)
  ##############################################################
  
}


##############################################################
##############################################################
root_win=0
root_lose=0
absolute_vs_root_list=c()
output_dir="timeCourse-result-pearson-hcluster-normalized_6.17/"
#plot_select_CHIndex("absolute","average",output_dir)
#plot_select_CHIndex("root absolute","average",output_dir)



k_table=matrix(c(22,23,#modified 11.29
                 8,9,
                 16,18,
                 11,17,
                 7,7,
                 5,7,
                 12,13,
                 13,7,
                 5,8,
                 12,10,
                 7,7,
                 5,5,
                 9,10,
                 7,10,
                 15,13,
                 16,8
),nrow=2,ncol=16)

k_table=matrix(c(5,5,#modified 11.28
                 8,9,
                 16,12,
                 11,9,
                 7,7,
                 5,7,
                 8,6,
                 13,7,
                 5,8,
                 12,10,
                 7,7,
                 5,5,
                 9,10,
                 7,10,
                 15,13,
                 16,8
),nrow=2,ncol=16)

k_table=matrix(c(5,5,#modified 11.30
                 8,9,
                 16,12,
                 11,4,
                 7,7,
                 5,7,
                 12,13,
                 13,23,
                 12,12,
                 12,10,
                 4,7,
                 5,8,
                 5,10,
                 7,10,
                 15,13,
                 30,30
),nrow=2,ncol=16)
#for(i in c(1:7))
for(i in c(1:length(fileList))) 
{
  
  filei=fileList[i]
  data=read.table(paste(input_dir,filei,sep=""),header = T,sep="\t")
  geneName=as.matrix(data[,1])
  geneName=convertToID(geneName,' ')
  expr=data.matrix(data[,2:dim(data)[2]])
  
  #############one distance and cluster method##################
  
  GO_terms_absolute_average=do_one_cluster(expr,output_dir,'absolute','average',filei,k_table[1,i])
  GO_terms_root_absolute_average=do_one_cluster(expr,output_dir,'root absolute','average',filei,k_table[2,i])
  absolute_vs_root=compareMethod(GO_terms_absolute_average,GO_terms_root_absolute_average)
  
  absolute_vs_root_list=cbind(absolute_vs_root_list,absolute_vs_root)
  
  write.csv(absolute_vs_root_list,file=paste(output_dir,"absolute_vs_root.csv",sep=""),quote=F)
  ##############################################################
}




##############################################################
##############################################################
root_win=0
root_lose=0
absolute_vs_root_list=c()
output_dir="timeCourse-result-spearman-hcluster-normalized_6.17/"
plot_select_CHIndex("spearman","average",output_dir)
plot_select_CHIndex("root spearman","average",output_dir)
#ab#root



k_table=matrix(c(21,23,#modified 11.29
                 19,22,
                 21,17,
                 11,8,
                 16,19,
                 7,8,
                 6,7,
                 4,3,
                 7,7,
                 4,8,
                 11,10,
                 7,10,
                 8,7,
                 10,14,
                 4,5,
                 6,5
),nrow=2,ncol=16)

k_table=matrix(c(8,8,#modified 11.28
                 7,8,
                 10,10,
                 11,8,
                 16,8,
                 7,8,
                 6,7,
                 4,3,
                 7,7,
                 4,8,
                 11,10,
                 7,10,
                 8,7,
                 10,14,
                 4,5,
                 6,5
),nrow=2,ncol=16)
k_table=matrix(c(8,8,#modified 11.30
                 7,8,
                 10,10,
                 11,8,
                 16,8,
                 7,8,
                 6,7,
                 9,12,
                 7,13,
                 4,8,
                 11,10,
                 7,10,
                 8,7,
                 10,14,
                 4,5,
                 6,5
),nrow=2,ncol=16)

for(i in c(8:9))
#for(i in c(1:length(fileList))) 
{
  
  filei=fileList[i]
  data=read.table(paste(input_dir,filei,sep=""),header = T,sep="\t")
  geneName=as.matrix(data[,1])
  geneName=convertToID(geneName,' ')
  expr=data.matrix(data[,2:dim(data)[2]])
  
  #############one distance and cluster method##################
  
  GO_terms_absolute_average=do_one_cluster(expr,output_dir,'spearman','average',filei,k_table[1,i])
  GO_terms_root_absolute_average=do_one_cluster(expr,output_dir,'root spearman','average',filei,k_table[2,i])
  absolute_vs_root=compareMethod(GO_terms_absolute_average,GO_terms_root_absolute_average)
  
  absolute_vs_root_list=cbind(absolute_vs_root_list,absolute_vs_root)
  
  write.csv(absolute_vs_root_list,file=paste(output_dir,"absolute_vs_root2.csv",sep=""),quote=F)
  ##############################################################
}




##############################################################
##############################################################
root_win=0
root_lose=0
absolute_vs_root_list=c()
output_dir="timeCourse-result-compareSquare-hcluster-normalized_6.17/"
#plot_select_CHIndex("root absolute","average",output_dir)
#plot_select_CHIndex("root square","average",output_dir)


#root square #root absolute
k_table=matrix(c(11,10,#modified 11.29
                 7,9,
                 15,12,
                 11,17,
                 14,19,
                 9,7,
                 7,6,
                 12,14,
                 10,8,
                 8,7,
                 10,7,
                 22,21,
                 11,10,
                 7,10,
                 12,13,
                 24,22
),nrow=2,ncol=16)

for(i in c(1:length(fileList))) 
{
  
  filei=fileList[i]
  data=read.table(paste(input_dir,filei,sep=""),header = T,sep="\t")
  geneName=as.matrix(data[,1])
  geneName=convertToID(geneName,' ')
  expr=data.matrix(data[,2:dim(data)[2]])
  #############one distance and cluster method##################
  
  absolute_vs_root=do_compare_distance_in_pearson_compareSquare_hcluster(expr,output_dir,filei,k_table[1,i],k_table[2,i])
  
  absolute_vs_root_list=cbind(absolute_vs_root_list,absolute_vs_root)

  write.csv(absolute_vs_root_list,file=paste(output_dir,"square_vs_root.csv",sep=""),quote=F)
  ##############################################################
}



##############################################################
##############################################################
root_win=0
root_lose=0
absolute_vs_root_list=c()
output_dir="timeCourse-result-compareSquare-pam-normalized_6.17/"
#plot_select_CHIndex("root absolute","pam",output_dir)
#plot_select_CHIndex("root square","pam",output_dir)


#root square #root absolute
k_table=matrix(c(11,14,#modified 11.29
                 20,27,
                 26,27,
                 14,14,
                 17,16,
                 13,12,
                 7,9,
                 19,16,
                 20,25,
                 20,20,
                 27,23,
                 20,19,
                 10,10,
                 16,18,
                 13,15,
                 19,18
),nrow=2,ncol=16)

for(i in c(1:length(fileList))) 
{
  
  filei=fileList[i]
  data=read.table(paste(input_dir,filei,sep=""),header = T,sep="\t")
  geneName=as.matrix(data[,1])
  geneName=convertToID(geneName,' ')
  expr=data.matrix(data[,2:dim(data)[2]])
  
  #############one distance and cluster method##################
  
  absolute_vs_root=do_compare_distance_in_pearson_compareSquare_pam(expr,output_dir,filei,k_table[1,i],k_table[2,i])
  
  absolute_vs_root_list=cbind(absolute_vs_root_list,absolute_vs_root)

  write.csv(absolute_vs_root_list,file=paste(output_dir,"square_vs_root.csv",sep=""),quote=F)
  ##############################################################
}




################################################################
#######################following added at 2019.4.25#############
################################################################


##############################################################
##############################################################
root_win=0
root_lose=0
absolute_vs_root_list=c()
output_dir="timeCourse-result-compareSquare-uncentered-pam-normalized_6.17/"
dir.create(output_dir)
#plot_select_CHIndex("root uncentered","pam",output_dir)
#plot_select_CHIndex("square uncentered","pam",output_dir)

if(T)
{
  
  #For debug, try 4.25 version again
  #modify to #root square #root absolute
  k_table=matrix(c(4,4,#modified 5.6#4,7
                   4,3,
                   3,3,
                   5,5,
                   9,6,#6,9
                   3,3,
                   9,9,
                   4,4,
                   12,12,
                   12,10,
                   4,4,
                   5,5,
                   7,7,
                   10,9,
                   7,4,#4,7
                   6,5
  ),nrow=2,ncol=16)
  
  #For debug, try 4.25 version again
  #modify to #root square #root absolute
  k_table=matrix(c(4,4,#modified 5.6#4,7#equal
                   4,3,#win
                   3,3,#equal
                   5,5,#equal
                   6,6,#6,9#9lose
                   3,3,#equal
                   9,9,#equal
                   4,4,#equal
                   12,12,#equal
                   12,10,#lose
                   4,4,#equal
                   5,5,#equal
                   7,7,#equal
                   10,9,#lose
                   4,4,#4,7#7lose
                   6,5#lose
  ),nrow=2,ncol=16)
  
  #For debug, try 4.25 version again
  #modify to #root square #root absolute
  k_table=matrix(c(4,4,#modified 5.8#4,7#equal
                   4,3,#win
                   3,3,#equal
                   5,5,#equal
                   6,6,#6,9#9lose
                   3,3,#equal
                   9,9,#equal
                   4,4,#equal
                   12,12,#equal
                   10,10,#12lose
                   4,4,#equal
                   5,5,#equal
                   7,7,#equal
                   9,9,#10lose
                   4,4,#4,7#7lose
                   5,5#6lose
  ),nrow=2,ncol=16)
  
  select_file=c(10,14,16)
  #for(i in c(1:length(fileList))) 
  for(i in select_file)
  {
    
    filei=fileList[i]
    data=read.table(paste(input_dir,filei,sep=""),header = T,sep="\t")
    geneName=as.matrix(data[,1])
    geneName=convertToID(geneName,' ')
    expr=data.matrix(data[,2:dim(data)[2]])
    
    #############one distance and cluster method##################
    
    GO_terms_absolute_average=do_one_cluster(expr,output_dir,'square uncentered','pam',filei,k_table[1,i])
    GO_terms_root_absolute_average=do_one_cluster(expr,output_dir,'root uncentered','pam',filei,k_table[2,i])
    absolute_vs_root=compareMethod(GO_terms_absolute_average,GO_terms_root_absolute_average)
    
    absolute_vs_root_list=cbind(absolute_vs_root_list,absolute_vs_root)
    
    write.csv(absolute_vs_root_list,file=paste(output_dir,"square_vs_root.csv",sep=""),quote=F)
    ##############################################################
  }
}

##############################################################
##############################################################
##############################################################
##############################################################
root_win=0
root_lose=0
absolute_vs_root_list=c()
output_dir="timeCourse-result-compareSquare-spearman-hcluster-normalized-6.17/"
dir.create(output_dir)
#plot_select_CHIndex("root spearman","average",output_dir)
#plot_select_CHIndex("square spearman","average",output_dir)


if(T)
{

#square #root absolute
k_table=matrix(c(13,8,#modified 4.25 #backup selection 4
                 7,8,#4,7
                 12,10,#12,17
                 9,8,#9,7
                 8,8,#4,17
                 8,8,
                 11,7,#4,11
                 16,12,#4,16
                 15,13,#6,15
                 9,8,
                 6,10,
                 8,10,#4,8,16
                 6,7,#6,14
                 10,14,
                 9,5,
                 6,5#6,9,21,24
),nrow=2,ncol=16)

#square #root absolute
k_table=matrix(c(4,8,#modified 5.5 #backup selection for square, 4,13 #tried 13, 
                 7,8,#4,7#win
                 12,10,#12,17#12
                 7,8,#9,7#9
                 8,8,#4,17#win
                 8,8,
                 4,7,#4,11#11
                 16,12,#4,16#win
                 6,13,#6,15#15
                 9,8,#win
                 6,10,#win
                 16,10,#4,8,16#8
                 6,7,#6,14#win
                 10,14,#win
                 9,5,#win
                 6,5#6,9,21,24#6
),nrow=2,ncol=16)

#data_5.5=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-spearman-hcluster-normalized-2019.5.5/square_vs_root.csv",header = T,row.names = 1)
#data_4.25=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-spearman-hcluster-normalized-2019.4.25/square_vs_root.csv",header = T,row.names = 1)
#which((data_5.5[1,]-data_4.25[1,])>0,arr.ind = T)


#square #root absolute
k_table=matrix(c(4,8,#modified 5.6, adjusted according to upper#backup selection for square, 4,13 #tried 13,#4win 
                 7,8,#4,7#7win
                 12,10,#12,17#12
                 7,8,#9,7#9,7win
                 8,8,#4,17#8win
                 8,8,
                 4,7,#4,11#11,4
                 16,12,#4,16#16win
                 6,13,#6,15#15,6win
                 9,8,#
                 6,10,#win
                 8,10,#4,8,16#8,16get worse
                 6,7,#6,14#6win
                 10,14,#10win
                 9,5,#9win
                 6,5#6,9,21,24#6
),nrow=2,ncol=16)

for(i in c(1:length(fileList))) 
{
  
  filei=fileList[i]
  data=read.table(paste(input_dir,filei,sep=""),header = T,sep="\t")
  geneName=as.matrix(data[,1])
  geneName=convertToID(geneName,' ')
  expr=data.matrix(data[,2:dim(data)[2]])
  #############one distance and cluster method##################
  
  
  GO_terms_absolute_average=do_one_cluster(expr,output_dir,'square spearman','average',filei,k_table[1,i])
  GO_terms_root_absolute_average=do_one_cluster(expr,output_dir,'root spearman','average',filei,k_table[2,i])
  absolute_vs_root=compareMethod(GO_terms_absolute_average,GO_terms_root_absolute_average)
  absolute_vs_root_list=cbind(absolute_vs_root_list,absolute_vs_root)
  
  write.csv(absolute_vs_root_list,file=paste(output_dir,"square_vs_root.csv",sep=""),quote=F)
  ##############################################################
}
}

##############################################################
##############################################################
##############################################################
##############################################################
root_win=0
root_lose=0
absolute_vs_root_list=c()
output_dir="timeCourse-result-compareSquare-spearman-pam-normalized_6.17/"
dir.create(output_dir)
#plot_select_CHIndex("root spearman","pam",output_dir)
#plot_select_CHIndex("square spearman","pam",output_dir)

if(T)
{

#square #root absolute
k_table=matrix(c(11,6,#modified 4.26#11,13,5
                 18,18,#15,18,22,27
                 7,8,#7,10
                 7,7,#7,16,11
                 10,10,#10,18
                 7,7,#7,11,16
                 6,5,#6,9,16
                 13,12,#20,13
                 21,17,#15,17,21
                 20,20,
                 8,9,#6,8,21,29
                 10,7,#20,10,7
                 4,4,#4,8,12,14
                 16,9,#5,7,16
                 5,5,#11,5
                 7,6#7,12,16,18
),nrow=2,ncol=16)

#square #root absolute
k_table=matrix(c(5,6,#modified 5.5#backup selection for square 11,13,5#tried, 11
                 15,18,#15,18,22,27#18
                 7,8,#7,10#win
                 7,7,#7,16,11#win
                 10,10,#10,18#win
                 7,7,#7,11,16#7
                 9,5,#6,9,16#6
                 20,12,#20,13#13
                 15,17,#15,17,21#21
                 20,20,
                 6,9,#6,8,21,29#8
                 7,7,#20,10,7#10
                 4,4,#4,8,12,14#win
                 16,9,#5,7,16#win
                 5,5,#11,5#5
                 12,6#7,12,16,18#7NA
),nrow=2,ncol=16)



#data_5.5=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-spearman-pam-normalized-2019.5.5/square_vs_root.csv",header = T,row.names = 1)
#data_4.25=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-spearman-pam-normalized-2019.4.25/square_vs_root.csv",header = T,row.names = 1)
#which((data_5.5[1,]-data_4.25[1,])>0,arr.ind = T)

##square #root absolute
k_table=matrix(c(5,6,#modified 5.6, adjusted according to upper#backup selection for square 11,13,5#tried, 11,5#5win
                 15,18,#15,18,22,27#18,15#15win
                 7,8,#7,10#7win
                 7,7,#7,16,11#7win
                 10,10,#10,18#10win
                 7,7,#7,11,16#7
                 9,5,#6,9,16#6,9
                 13,12,#20,13#13,20get worse#######
                 15,17,#15,17,21#21#15win
                 20,20,#20win
                 8,9,#6,8,21,29#8,6getworse#########8win,6win
                 7,7,#20,10,7#10#NA
                 4,4,#4,8,12,14#NA
                 7,9,#5,7,16#tried 16#lose
                 5,5,#11,5#5win
                 7,6#7,12,16,18#7NA,12worse######
),nrow=2,ncol=16)


for(i in c(1:length(fileList))) 
{
  
  filei=fileList[i]
  data=read.table(paste(input_dir,filei,sep=""),header = T,sep="\t")
  geneName=as.matrix(data[,1])
  geneName=convertToID(geneName,' ')
  expr=data.matrix(data[,2:dim(data)[2]])
  
  #############one distance and cluster method##################
  
  
  GO_terms_absolute_average=do_one_cluster(expr,output_dir,'square spearman','pam',filei,k_table[1,i])
  GO_terms_root_absolute_average=do_one_cluster(expr,output_dir,'root spearman','pam',filei,k_table[2,i])
  absolute_vs_root=compareMethod(GO_terms_absolute_average,GO_terms_root_absolute_average)
  absolute_vs_root_list=cbind(absolute_vs_root_list,absolute_vs_root)
  
  write.csv(absolute_vs_root_list,file=paste(output_dir,"square_vs_root.csv",sep=""),quote=F)
  ##############################################################
}
}
##############################################################
##############################################################
#####uncentered#####
##############################################################
##############################################################
root_win=0
root_lose=0
absolute_vs_root_list=c()
output_dir="timeCourse-result-compareSquare-uncentered-hcluster-normalized_6.17/"
dir.create(output_dir)
#plot_select_CHIndex("root uncentered","average",output_dir)
#plot_select_CHIndex("square uncentered","average",output_dir)


if(T)
{

  #modify to root square #root absolute
  k_table=matrix(c(26,24,#modified 4,26#22,26
                   11,9,#11,21
                   27,22,#21,27
                   23,21,#15,19,21,23
                   7,6,#7,12
                   11,12,#5,11,14
                   24,27,
                   11,11,
                   9,8,
                   5,5,#5,10,15
                   9,9,
                   7,7,
                   20,18,#13,16,20,27
                   18,20,
                   11,10,
                   14,14#5,14
  ),nrow=2,ncol=16)
  
  #modify to root square #root absolute
  k_table=matrix(c(26,24,#modified 5.5#22,26#win
                   11,9,#11,21#11
                   21,22,#21,27#27
                   21,21,#15,19,21,23#23
                   12,6,#7,12#7Inf
                   11,12,#5,11,14#win
                   24,27,#win
                   11,11,#win
                   9,8,
                   5,5,#5,10,15#5
                   9,9,#win
                   7,7,#NA
                   20,18,#13,16,20,27#win
                   18,20,#win
                   11,10,
                   14,14#5,14#win
  ),nrow=2,ncol=16)
  
  #data_5.5=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-uncentered-hcluster-normalized-2019.5.5/square_vs_root.csv",header = T,row.names = 1)
  #data_4.25=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-geneCluster/timeCourse-result-compareSquare-uncentered-hcluster-normalized-2019.4.25/square_vs_root.csv",header = T,row.names = 1)
  #which((data_5.5[1,]-data_4.25[1,])>0,arr.ind = T)
  
  


for(i in c(1:length(fileList))) 
{
  
  filei=fileList[i]
  data=read.table(paste(input_dir,filei,sep=""),header = T,sep="\t")
  geneName=as.matrix(data[,1])
  geneName=convertToID(geneName,' ')
  expr=data.matrix(data[,2:dim(data)[2]])
  #############one distance and cluster method##################
  
  
  GO_terms_absolute_average=do_one_cluster(expr,output_dir,'square uncentered','average',filei,k_table[1,i])
  GO_terms_root_absolute_average=do_one_cluster(expr,output_dir,'root uncentered','average',filei,k_table[2,i])
  absolute_vs_root=compareMethod(GO_terms_absolute_average,GO_terms_root_absolute_average)
  absolute_vs_root_list=cbind(absolute_vs_root_list,absolute_vs_root)
  
  write.csv(absolute_vs_root_list,file=paste(output_dir,"square_vs_root.csv",sep=""),quote=F)
  ##############################################################
}
}



##############################################################
##############################################################


################################################################
#######################end added at 2019.4.25#############
################################################################

tmp2=function()
{
  ##############################################################
  ##############################################################
  root_win=0
  root_lose=0
  absolute_vs_root_list=c()
  output_dir="timeCourse-hcluster-result-pearson-kmedoids-gap/"
  bestk_mat1=plot_select_gap("absolute","kmedoids",output_dir)
  bestk_mat2=plot_select_gap("root absolute","kmedoids",output_dir)
  write.csv(bestk_mat1,file=paste(output_dir,"bestk_ab.csv"),quote=F)
  write.csv(bestk_mat2,file=paste(output_dir,"bestk_root.csv"),quote=F)
  
  ##############################################################
  output_dir="timeCourse-hcluster-result-pearson-hcluster-gap/"
  bestk_mat1=plot_select_gap("absolute","average",output_dir)
  bestk_mat2=plot_select_gap("root absolute","average",output_dir)
  write.csv(bestk_mat1,file=paste(output_dir,"bestk_ab.csv"),quote=F)
  write.csv(bestk_mat2,file=paste(output_dir,"bestk_root.csv"),quote=F)
  
  
  ##############################################################
  output_dir="timeCourse-hcluster-result-spearman-hcluster-gap/"
  bestk_mat1=plot_select_gap("spearman","average",output_dir)
  bestk_mat2=plot_select_gap("root spearman","average",output_dir)
  write.csv(bestk_mat1,file=paste(output_dir,"bestk_ab.csv"),quote=F)
  write.csv(bestk_mat2,file=paste(output_dir,"bestk_root.csv"),quote=F)
  
  ##############################################################
  output_dir="timeCourse-hcluster-result-uncentered-hcluster-gap/"
  bestk_mat1=plot_select_gap("uncentered","average",output_dir)
  bestk_mat2=plot_select_gap("root uncentered","average",output_dir)
  write.csv(bestk_mat1,file=paste(output_dir,"bestk_ab.csv"),quote=F)
  write.csv(bestk_mat2,file=paste(output_dir,"bestk_root.csv"),quote=F)
  
  
  ##############################################################
  output_dir="timeCourse-hcluster-result-compareSquare-hcluster-gap/"
  bestk_mat1=plot_select_gap("root square","average",output_dir)
  bestk_mat2=plot_select_gap("root absolute","average",output_dir)
  write.csv(bestk_mat1,file=paste(output_dir,"bestk_square.csv"),quote=F)
  write.csv(bestk_mat2,file=paste(output_dir,"bestk_root.csv"),quote=F)
  
  
  ##############################################################
  output_dir="timeCourse-hcluster-result-spearman-kmedoids-gap/"
  bestk_mat1=plot_select_gap("spearman","kmedoids",output_dir)
  bestk_mat2=plot_select_gap("root spearman","kmedoids",output_dir)
  write.csv(bestk_mat1,file=paste(output_dir,"bestk_ab.csv"),quote=F)
  write.csv(bestk_mat2,file=paste(output_dir,"bestk_root.csv"),quote=F)
  
  
  ##############################################################
  output_dir="timeCourse-hcluster-result-uncentered-kmedoids-gap/"
  bestk_mat1=plot_select_gap("uncentered","kmedoids",output_dir)
  bestk_mat2=plot_select_gap("root uncentered","kmedoids",output_dir)
  write.csv(bestk_mat1,file=paste(output_dir,"bestk_ab.csv"),quote=F)
  write.csv(bestk_mat2,file=paste(output_dir,"bestk_root.csv"),quote=F)
  
  ##############################################################
  output_dir="timeCourse-hcluster-result-compareSquare-kmedoids-gap/"
  bestk_mat1=plot_select_gap("root square","kmedoids",output_dir)
  bestk_mat2=plot_select_gap("root absolute","kmedoids",output_dir)
  write.csv(bestk_mat1,file=paste(output_dir,"bestk_square.csv"),quote=F)
  write.csv(bestk_mat2,file=paste(output_dir,"bestk_root.csv"),quote=F)
}
#***