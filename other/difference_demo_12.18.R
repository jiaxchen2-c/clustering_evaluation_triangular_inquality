#source and library
source("function_for_test_difference.R")
#input and parameter
input_dir="./CompCancer/"
list_file=paste(input_dir,"list.txt",sep="")

filei=18
distance_method_absolute="absolute"
distance_method_root="root absolute"
############################load data######################################################
#load data and reference cluster
fileList=as.matrix(read.table(list_file,header = F))

filei=fileList[filei]
fileName=paste(input_dir,filei,sep="")

#load data and reference cluster
data=load_data_and_normalize(fileName)
expr=data$expr
geneID=data$geneID
sampleID=data$sampleID
ref_cluster_array=data$ref_cluster_array
ref_cluster=data$ref_cluster

unique_ref_cluster=unique(as.character(ref_cluster))
k=length(unique_ref_cluster)
#print("estimate k")
#print(k)
##############################end load data#################################################




mat_distance_absolute1=produceDistanceMatrices(expr,distance_method_absolute)

mat_distance_root1=produceDistanceMatrices(expr,distance_method_root)

############################perform hclust on two distance#######################

hc_ab=hclust(as.dist(mat_distance_absolute1),method = "average")
hc_root=hclust(as.dist(mat_distance_root1),method = "average")


################################################################################

##################################test by hclust process########################
mat_distance_absolute1=as.matrix(mat_distance_absolute1)

total_ranks=c()
rank_for_plot=NA
trils=53
for(i in c(2:trils))
{#for 1:16, because after 16, the order of cluster and a will be different in two result
  #so we can not use rank_mat1-rank_mat2 to compare or to display the difference.
  print(i)
  all_x_all_ab=hclust_cal_distance_by_merge(mat_distance_absolute1,hc_ab$merge,i)
  all_x_all_root=hclust_cal_distance_by_merge(mat_distance_root1,hc_root$merge,i)
  
  iteration_look_for_disobeyTrianglularInequality(all_x_all_ab$distance_mat,all_x_all_ab$colOrder)
  iteration_look_for_disobeyTrianglularInequality(all_x_all_root$distance_mat,all_x_all_root$colOrder)
  
  
  pair_and_rank=display_rank_diff(all_x_all_ab,all_x_all_root,i,sampleID,T)
  
  #tmp=as.numeric(as.matrix(pair_and_rank[,3]))
  #A=which(tmp<30)
  #print(pair_and_rank[A,])
  
  n=dim(all_x_all_ab$rank_mat)[1]
  total_rank_in_this_round=n*(n+1)/2
  
  total_ranks=c(total_ranks,total_rank_in_this_round)
  rank_for_plot_this_round=as.numeric(as.matrix(pair_and_rank[,3]))
  
  if(is.na(rank_for_plot))
  {
    rank_for_plot=matrix(0,total_rank_in_this_round,trils)
  }
  rank_for_plot[rank_for_plot_this_round,i]=1
  rank_for_plot[total_rank_in_this_round,i]=rank_for_plot[total_rank_in_this_round,i]+2
  
}
ggplot_diff_rank_changes(rank_for_plot = rank_for_plot,T,"rank_vary18.pdf")
ggplot_diff_rank_changes(rank_for_plot = rank_for_plot,F,"rank_vary_all18.pdf")



iteration_look_for_disobeyTrianglularInequality(mat_distance_root1)
iteration_look_for_disobeyTrianglularInequality(mat_distance_absolute1)


##################################end by hclust process########################

##################################begin plot dendrogram difference#############
input_dir="./CompCancer/"
output_dir="./"
list_file=paste(input_dir,"list.txt",sep="")
#load data and reference cluster
fileList=as.matrix(read.table(list_file,header = F))
distance_method_absolute="absolute"
distance_method_root="root absolute"


select_plot_diff_dendrogram=c(1:35)
length_of_sample=rep(0,35)

for(i in select_plot_diff_dendrogram)
{
  
  filei=fileList[i]
  fileName_i=paste(input_dir,filei,sep="")
  
  #
  data=load_data_and_normalize(fileName_i)
  expr=data$expr
  sampleID_i=data$sampleID
  print(paste("i",i,"length of samples",length(data$sampleID)))
  length_of_sample[i]=length(data$sampleID)
  
  mat_distance_absolute_i=produceDistanceMatrices(expr,distance_method_absolute)
  
  mat_distance_root_i=produceDistanceMatrices(expr,distance_method_root)
  
  hc_ab_i=hclust(as.dist(mat_distance_absolute_i),method = "average")
  hc_root_i=hclust(as.dist(mat_distance_root_i),method = "average")
  
  plot_diff_dendrogram(hc_ab_i$merge,hc_root_i$merge,sampleID_i,fileHeader = paste(output_dir,i,sep=""),"absolute","root")
  
  
}
file=c(1:35)
df1=data.frame(file,length_of_sample)

ggplot(df1,aes(x=file,y=length_of_sample))+geom_bar(stat="identity")+theme_light()+ylab("number of samples")

##################################end plot dendrogram difference###############
