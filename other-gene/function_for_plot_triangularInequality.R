library(ggplot2)
library(ggrepel)
library(reshape2)
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