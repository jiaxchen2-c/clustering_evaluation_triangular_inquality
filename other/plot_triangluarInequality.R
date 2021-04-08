library(ggplot2)
library(reshape2)
library(ggrepel)
source("function_for_distance.R")
source("function_for_test_difference.R")
source("function_for_cluster.R")
#input and parameter
input_dir="./CompCancer/"
list_file=paste(input_dir,"list.txt",sep="")
#load data and reference cluster
fileList=as.matrix(read.table(list_file,header = F))
############################load data######################################################
#i=28
#filei=fileList[i]
#fileName=paste(input_dir,filei,sep="")
#data=read.table(fileName,header = T,sep="\t")
#sampleID=as.matrix(colnames(data)[2:dim(data)[2]])
trils=20
frequencys_mat=matrix(0,trils,35)
distance_method_absolute="absolute"
distance_method_root="root absolute"

for(filei in c(1:35))
{
  print(filei)
  
  fileName_i=paste(input_dir,fileList[filei],sep="")
  
  #
  data=load_data_and_normalize(fileName_i,norm=T)
  expr=data$expr
  sampleID_i=data$sampleID
  
  
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
    
    #
    
  }
  #frequencys_mat=cbind(frequencys_mat,frequencys_array)
  
}
pdf("disobey_Triangular_inequality_heatmap.pdf")
plot_heatmap(record_mat)
dev.off()
plot_curve_mat(frequencys_mat,"disobey_Triangular_inequality_frequence_trend_byPair_norm.pdf")

##################################define functions##############################
plot_heatmap=function(mat)
{
  test=mat
  test=as.data.frame(test)
  test$id=rownames(test)
  mat_melt <- melt(test, id.var="id")
  p1=ggplot(mat_melt,aes(x=as.integer(variable),y=as.integer(id)))+
    geom_tile(aes(fill=mat_melt$value),colour="black")+theme_light()+
    #scale_fill_gradient("ratio of together",low="steelblue",medium="white",high="firebrick")
  scale_fill_gradientn(colours = c("white", "steelblue"),
                       values = scales::rescale(c( 0, 1)))
  print(p1)
  
}

plot_curve=function(frequencys_array)
{
  y=frequencys_array
  x=c(1:length(frequencys_array))
  df=data.frame(x,y)
  ggplot(df)+geom_line(aes(x=x,y=y))
}

plot_curve_mat=function(frequencys_mat,fileName)
{
  y=c()
  x=c()
  type=c()
  for(i in c(1:dim(frequencys_mat)[2]))
  {
    y=c(y,as.numeric(frequencys_mat[,i]))
    x=c(x,c(1:dim(frequencys_mat)[1]))
    type=c(type,rep(i,dim(frequencys_mat)[1]))
  }
  
  type=as.character(type)
  df=data.frame(x,y,type)
  pdf(fileName,width = 12,height = 10)
  p1=ggplot(df)+geom_line(aes(x=x,y=y,color=type))+theme_light(base_size = 13)+
    xlab("step")+ylab("frequency")
  print(p1)
  dev.off()
}


look_for_disobeyTrianglularInequality=function(distance_mat,colOrder=NA)
{
  disobey_record_mat=matrix(0,dim(distance_mat)[1],dim(distance_mat)[2])
  if(is.na(colOrder)){
    colOrder=c(1:dim(distance_mat)[1])
  }
  for(i in c(1:dim(distance_mat)[1]))
  {
    for(j in c(i+1:dim(distance_mat)[2]))
    {
      if(j>dim(distance_mat)[2])
      {
        next
      }
      dist_ij=distance_mat[i,j]
      
      useful_index=setdiff(c(1:dim(distance_mat)[1]),c(i,j))
      
      dist_i=distance_mat[i,useful_index]
      dist_j=distance_mat[j,useful_index]
      
      plus_dist=dist_i+dist_j
      minus_dist=abs(dist_i-dist_j)
      
      A=which(plus_dist<dist_ij)
      B=which(minus_dist>dist_ij)
      
      if(length(A)>0||length(B)>0)
      {
        disobey_record_mat[i,j]=1
        #print(paste("ij",i,j,colOrder[i],colOrder[j],"disobey trianglur inequality"))
        
        C=setdiff(A,B)
        length(C)
        #for(k in useful_index[A])
        #{
          #disobey_record_mat[i,k]=1
          #disobey_record_mat[j,k]=1
          #print(paste("(",colOrder[i],colOrder[j],")",distance_mat[i,j],"(",colOrder[i],colOrder[k],")",distance_mat[i,k],"(",colOrder[j],colOrder[k],")",distance_mat[k,j]))
        #}
        #for(k in useful_index[B])
        #{
          #disobey_record_mat[i,k]=1
          #disobey_record_mat[j,k]=1
          #print(paste("(",colOrder[i],colOrder[j],")",distance_mat[i,j],"(",colOrder[i],colOrder[k],")",distance_mat[i,k],"(",colOrder[j],colOrder[k],")",distance_mat[k,j]))
          
        #}
      }
    }
  }
  disobey_record_mat
}

look_for_disobeyTrianglularInequality_recordByPair=function(distance_mat,colOrder=NA)
{
  #pair with duplicate, by the ratio should be the same
  disobey_pair_count=0
  total_pair_count=0
  
  if(is.na(colOrder)){
    colOrder=c(1:dim(distance_mat)[1])
  }
  for(i in c(1:dim(distance_mat)[1]))
  {
    for(j in c(i+1:dim(distance_mat)[2]))
    {
      if(j>dim(distance_mat)[2])
      {
        next
      }
      dist_ij=distance_mat[i,j]
      
      useful_index=setdiff(c(1:dim(distance_mat)[1]),c(i,j))
      
      dist_i=distance_mat[i,useful_index]
      dist_j=distance_mat[j,useful_index]
      
      plus_dist=dist_i+dist_j
      minus_dist=abs(dist_i-dist_j)
      
      A=which(plus_dist<dist_ij)
      B=which(minus_dist>dist_ij)
      
      C=union(A,B)
      disobey_pair_count=disobey_pair_count+length(C)
      total_pair_count=total_pair_count+length(useful_index)
      
    }
  }
  
  list(disobey_pair_count=disobey_pair_count,
       total_pair_count=total_pair_count)
}