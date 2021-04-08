#source and library
source("function_for_test_difference.R")
#input and parameter
input_dir="./CompCancer/"
list_file=paste(input_dir,"list.txt",sep="")

#load data and reference cluster
fileList=as.matrix(read.table(list_file,header = F))


############################load data######################################################
i=28
filei=fileList[i]
fileName=paste(input_dir,filei,sep="")
data=read.table(fileName,header = T,sep="\t")
sampleID=as.matrix(colnames(data)[2:dim(data)[2]])

mat_distance_absolute1=as.matrix(read.csv("result/28_absolute_averagedistanceMatrix.csv",header = F))

mat_distance_root1=as.matrix(read.csv("result/28_root absolute_averagedistanceMatrix.csv",header = F))

mat_distance_absolute2=as.matrix(read.csv("result/29_absolute_averagedistanceMatrix.csv",header = F))

mat_distance_root2=as.matrix(read.csv("result/29_root absolute_averagedistanceMatrix.csv",header = F))


############################perform hclust on two distance#######################

hc_ab=hclust(as.dist(mat_distance_absolute1),method = "average")
hc_root=hclust(as.dist(mat_distance_root1),method = "average")


#for hc_ab$merge
A=which(hc_ab$merge<0,arr.ind = T)
hc_ab$merge[A]=sampleID[-hc_ab$merge[A]]
write.csv(hc_ab$merge,file="hc_ab_merge.csv",quote=F)

A=which(hc_root$merge<0,arr.ind = T)
hc_root$merge[A]=sampleID[-hc_root$merge[A]]
write.csv(hc_root$merge,file="hc_root_merge.csv",quote=F)
################################################################################

##############################begin test rank###################################
mat_rank_absolute1=t(apply(mat_distance_absolute1,1,order))
mat_rank_root1=t(apply(mat_distance_root1,1,order))
mat_rank_absolute2=t(apply(mat_distance_absolute2,1,order))
mat_rank_root2=t(apply(mat_distance_root2,1,order))

diff1=mat_rank_absolute1-mat_rank_root1
diff2=mat_rank_absolute2-mat_rank_root2

################################end test rank###################################
################################begin test D27##################################
#clusterA
D27=which(sampleID=="DLBC27")
#clusterB
D7=which(sampleID=="DLBC7")
#clusterC
D17=which(sampleID=="DLBC17")
F17=which(sampleID=="FSCC17")
D54=which(sampleID=="DLBC54")
F12=which(sampleID=="FSCC12")


#for root
#between A and B
mat_distance_root1[D27,D7]

#between A and C #root
dist_s=c(mat_distance_root1[D27,D17],mat_distance_root1[D27,F17],
         mat_distance_root1[D27,D54],mat_distance_root1[D27,F12])
mean(dist_s)

#for absolute
#between A and B #absolute
mat_distance_absolute1[D27,D7]

#between A and C 
dist_s2=c(mat_distance_absolute1[D27,D17],mat_distance_absolute1[D27,F17],
         mat_distance_absolute1[D27,D54],mat_distance_absolute1[D27,F12])
mean(dist_s2)

##################################end test D27##################################

##################################test by hclust process########################
mat_distance_absolute1=as.matrix(mat_distance_absolute1)
mat_rank_absolute1=t(apply(mat_distance_absolute1,1,order))
mat_rank_root1=t(apply(mat_distance_root1,1,order))

test1=mat_distance_root1[upper.tri(mat_distance_root1)]
test2=mat_distance_absolute1[upper.tri(mat_distance_absolute1)]
order1= which(test1==min(min(test1)),arr.ind = T)
order2=which(test2==min(min(test2)),arr.ind = T)

diff_order=order1-order2
print(diff_order)
#diff_rank=mat_rank_absolute1-mat_rank_root1
#A=which(diff_rank!=0,arr.ind = T)
#print(length(A))
####################
hc_ab=hclust(as.dist(mat_distance_absolute1),method = "average")
hc_root=hclust(as.dist(mat_distance_root1),method = "average")


total_ranks=c()
rank_for_plot=NA
trils=16
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
ggplot_diff_rank_changes(rank_for_plot = rank_for_plot,T,"rank_vary.pdf")
ggplot_diff_rank_changes(rank_for_plot = rank_for_plot,F,"rank_vary_all.pdf")



iteration_look_for_disobeyTrianglularInequality(mat_distance_root1)
iteration_look_for_disobeyTrianglularInequality(mat_distance_absolute1)


##################################end by hclust process########################

##################################begin plot dendrogram difference#############
input_dir="./CompCancer/"
output_dir="./"
list_file=paste(input_dir,"list.txt",sep="")
#load data and reference cluster
fileList=as.matrix(read.table(list_file,header = F))



select_plot_diff_dendrogram=c(28)

for(i in select_plot_diff_dendrogram)
{
  
  filei=fileList[i]
  fileName_i=paste(input_dir,filei,sep="")
  data_i=read.table(fileName_i,header = T,sep="\t")
  sampleID_i=as.matrix(colnames(data_i)[2:dim(data_i)[2]])
  
  mat_distance_absolute_i=as.matrix(read.csv(paste("result/",i,"_absolute_averagedistanceMatrix.csv",sep=""),header = F))
  
  mat_distance_root_i=as.matrix(read.csv(paste("result/",i,"_root absolute_averagedistanceMatrix.csv",sep=""),header = F))
  
  hc_ab_i=hclust(as.dist(mat_distance_absolute_i),method = "average")
  hc_root_i=hclust(as.dist(mat_distance_root_i),method = "average")
  
  plot_diff_dendrogram(hc_ab_i$merge,hc_root_i$merge,sampleID_i,fileHeader = paste(output_dir,i,sep=""),"absolute","root")
  
  
}


##################################end plot dendrogram difference###############
