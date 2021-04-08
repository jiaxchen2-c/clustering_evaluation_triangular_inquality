data1= read.csv("~/Desktop/1_normalized_pcc.csv",header=T)
data1=data1[1:7,1:6]

data2= read.csv("~/Desktop/2.csv",header=T)
data2=data2[1:7,1:6]

data3= read.csv("~/Desktop/3.csv",header=T)
data3=data3[1:7,1:6]

data4= read.csv("~/Desktop/4.csv",header=T)
data4=data4[1:7,1:6]

data5= read.csv("~/Desktop/5.csv",header=T)
data5=data5[1:7,1:6]

data6= read.csv("~/Desktop/6.csv",header=T)
data6=data6[1:7,1:6]

data7= read.csv("~/Desktop/7.csv",header=T)
data7=data7[1:7,1:6]

matrix_list=list(data1,data2,data3,data4,data5,data6,data7)

library(ggplot2)

x1=c()
y1=c()
x2=c()
y2=c()
distance_method=c()
cluster_method=c()

for(k in c(1:7))#matrixs
{
  data=matrix_list[k][[1]]
  for(i in c(1:7))#dataset
  {
    for(j in c(1:6))#method
    {
      
      x1_ij=(1/7)*j+k
      x2_ij=(1/7)*j+k
      y1_ij=i
      y2_ij=i+data[i,j]
      
      if(j<4)
      {
        cluster_method_ij="hclust"
        
      }else{
        cluster_method_ij="pam"
      }
      if(j==1)
      {
        distance_method_ij="1 Hierarchical cluster, dr"
        
      }
      if(j==2){
        distance_method_ij="2 Hierarchical cluster, da"
      }
      if(j==3){
        distance_method_ij="3 Hierarchical cluster, ds"
      }
      if(j==4)
      {
        distance_method_ij="4 PAM, dr"
      }
      if(j==5)
      {
        distance_method_ij="5 PAM, da"
      }
      if(j==6)
      {
        distance_method_ij="6 PAM, ds"
      }
      
      x1=c(x1,x1_ij)
      x2=c(x2,x2_ij)
      y1=c(y1,y1_ij)
      y2=c(y2,y2_ij)
      distance_method=c(distance_method,distance_method_ij)
      cluster_method=c(cluster_method,cluster_method_ij)
    }
  }
}

df=data.frame(x1,x2,y1,y2,distance_method,cluster_method)

p1=ggplot(df)+geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2,colour=distance_method))+theme_minimal()

pdf("/Users/jiaxchen2/Desktop/2018.10.9/experiment-singlecell/with_cluster_partition_goldstandard/line_chart.pdf", width = 7, height=7)
print(p1)
dev.off()

