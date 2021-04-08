
library(ggplot2)
library(reshape2)
library(ggrepel)
source("function_for_test_difference.R")
source("function_for_cluster.R")
source("function_for_distance.R")
source("support_function.R")
source("function_for_plot_triangularInequality.R")

#parameter and load data
input_dir="benchmark_datasets_normalized/"
output_dir="timeCourse-hcluster-result-kmedoids-11_7/"
list_file=paste(input_dir,"list.txt",sep="")
toID=as.matrix(read.table("geneName2.txt",header=F))
fromID=as.character(as.matrix(read.table("OLN.txt",header=F)))
toID=convertToID(toID,';')

#fileList
fileList=as.matrix(read.table(list_file,header = F))

############################load data######################################################

trils=20
frequencys_mat=matrix(0,trils,16)
distance_method_absolute="absolute"
distance_method_root="root absolute"

for(filei in c(1:16))
{
  print(filei)
  
 
  fileName=fileList[filei]
  data=read.table(paste(input_dir,fileName,sep=""),header = T,sep="\t")
  geneName=as.matrix(data[,1])
  geneName=convertToID(geneName,' ')
  expr=data.matrix(data[,2:dim(data)[2]])
  
  #
  
  
  mat_distance_absolute1=produceDistanceMatrices(expr,distance_method_absolute)
  
  mat_distance_root1=produceDistanceMatrices(expr,distance_method_root)
  
  
  ################################
  hc_ab=hclust(as.dist(mat_distance_absolute1),method = "average")
  hc_root=hclust(as.dist(mat_distance_root1),method = "average")
  
  total_ranks=c()
  rank_for_plot=NA

  frequencys_array=c()
  for(i in c(1:trils))
  {#for 1:16, because after 16, the order of cluster and a will be different in two result
    #so we can not use rank_mat1-rank_mat2 to compare or to display the difference.
    #print(i)
    all_x_all_ab=hclust_cal_distance_by_merge(mat_distance_absolute1,hc_ab$merge,i)
    all_x_all_root=hclust_cal_distance_by_merge(mat_distance_root1,hc_root$merge,i)
    
    #method 1, frequency by (take the distance[i,j] into account if it form disaobey loop with any two other distance[k,l] in the dataset)
    #record_mat=look_for_disobeyTrianglularInequality(all_x_all_ab$distance_mat,all_x_all_ab$colOrder)
    #frequencys=sum(record_mat)/((dim(all_x_all_ab$distance_mat)[1]*(dim(all_x_all_ab$distance_mat)[2]-1))/2)
    
    #method 2, frequency by (disobey_loop_pair_count/total_loop_pair_count)
    record_count=look_for_disobeyTrianglularInequality_recordByPair(all_x_all_ab$distance_mat,all_x_all_ab$colOrder)
    frequencys=(record_count$disobey_pair_count)/(record_count$total_pair_count)
    print(frequencys)
    if(!frequencys>0)
    {
      break
    }else{
      frequencys_mat[i,filei]=frequencys
    }
    frequencys_array=c(frequencys_array,frequencys)
    #plot the range and frequency
    
    
  }
  #frequencys_mat=cbind(frequencys_mat,frequencys_array)
  
}
write.csv(frequencys_mat,file="disobey_Triangular_inequality_frequence_trend_gene_byPair.csv",quote=F)
plot_curve_mat(frequencys_mat,"disobey_Triangular_inequality_frequence_trend_gene_byPair.pdf")
