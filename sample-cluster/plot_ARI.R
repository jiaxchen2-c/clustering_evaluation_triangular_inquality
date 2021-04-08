require(ggplot2)

ggplot_ARI_mat=function(fileName,ARI_mat,select_col,label1="root",label2="absolute")
{
  n=dim(ARI_mat)[1]
  m=dim(ARI_mat)[2]
  ARI_for_plot=matrix(NA,n*2,6)
  colnames(ARI_for_plot)=c("ARI","dataset","mode","distance","group","x")
  index=1
  ARIs=rep(NA,m*n)
  datasets=rep(NA,m*n)
  ARI_for_plot=as.data.frame(ARI_for_plot)
  current_x=1
  
  ARI_method1=ARI_mat[,select_col[1]]
  ARI_method2=ARI_mat[,select_col[2]]
  
  wilcox_test_ARI = wilcox.test(ARI_method1,ARI_method2)
  print(wilcox_test_ARI)
  
  t_test_ARI = t.test(ARI_method1,ARI_method2)
  print(t_test_ARI)
  
  for(i in c(1:n))
  {
    dataset=i
    for(j in select_col)
    {
      ARI_for_plot[index,1]=as.numeric(ARI_mat[i,j])
      ARI_for_plot[index,2]=as.numeric(dataset)
      
      if(j==select_col[1]){
        ARI_for_plot[index,4]=label1
        ARI_for_plot[index,6]=current_x
        current_x=current_x+1
      }
      if(j==select_col[2]){
        ARI_for_plot[index,4]=label2
        ARI_for_plot[index,6]=current_x
        current_x=current_x+1
      }
      index=index+1
      current_x=current_x+2
    }
  }
  
  pdf(fileName,width=12,height = 8)
  p1=ggplot(data=ARI_for_plot, aes(x=dataset, y=ARI, fill=distance)) +
    geom_bar(stat="identity", position=position_dodge())+theme_light()
  scale_x_continuous(breaks= seq(1:n) ,label = seq(1:n))
  print(p1)
  dev.off()
}



####################################################################

ARI_mat_pearson=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-sampleCluster/otherDistance_ARI/ARI_mat_pearson.csv",row.names = 1,header = T)
ggplot_ARI_mat("pearson_hcluster.pdf",ARI_mat_pearson,c(1,2))
ggplot_ARI_mat("pearson_kmedoids.pdf",ARI_mat_pearson,c(3,4))
  
ARI_mat_spearman=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-sampleCluster/otherDistance_ARI/otherDistance_ARIARI_mat_spearman.csv",row.names = 1,header = T)
ggplot_ARI_mat("spearman_hcluster.pdf",ARI_mat_spearman,c(1,2))
ggplot_ARI_mat("spearman_kmedoids.pdf",ARI_mat_spearman,c(3,4))


ARI_mat_uncentered=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-sampleCluster/otherDistance_ARI/otherDistance_ARIARI_mat_uncentered.csv",row.names = 1,header = T)
ggplot_ARI_mat("uncentered_hcluster.pdf",ARI_mat_uncentered,c(1,2))
ggplot_ARI_mat("uncentered_kmedoids.pdf",ARI_mat_uncentered,c(3,4))

ARI_mat_pearson=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-sampleCluster/compareSquare_ARI/ARI_mat_pearson_compareSquare.csv",row.names = 1,header = T)
ggplot_ARI_mat("pearson_compareSquare_hcluster.pdf",ARI_mat_pearson,c(1,2),"root","square")
ggplot_ARI_mat("pearson_compareSquare_kmedoids.pdf",ARI_mat_pearson,c(3,4),"root","square")

ARI_mat_pearson=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-sampleCluster/compareSquare_ARI/ARI_mat_pearson_compareSquare_otherHcluster.csv",row.names = 1,header = T)
ggplot_ARI_mat("pearson_compareSquare_complete.pdf",ARI_mat_pearson,c(1,2),"root","square")
ggplot_ARI_mat("pearson_compareSquare_single.pdf",ARI_mat_pearson,c(3,4),"root","square")


ARI_mat_pearson=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-sampleCluster/compareSquare_ARI/ARI_mat_pearson_compareSquare_hclusterAndPam.csv",row.names = 1,header = T)
ggplot_ARI_mat("pearson_compareSquare_hcluster.pdf",ARI_mat_pearson,c(1,2),"root","square")
ggplot_ARI_mat("pearson_compareSquare_pam.pdf",ARI_mat_pearson,c(3,4),"root","square")


####################################################################

ARI_mat_pearson=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-sampleCluster/otherDistance_ARI/ARI_mat_pearson_hclusterAndPam.csv",row.names = 1,header = T)
ggplot_ARI_mat("pearson_hcluster.pdf",ARI_mat_pearson,c(1,2))
ggplot_ARI_mat("pearson_pam.pdf",ARI_mat_pearson,c(3,4))

ARI_mat_spearman=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-sampleCluster/otherDistance_ARI/ARI_mat_spearman_hclusterAndPam.csv",row.names = 1,header = T)
ggplot_ARI_mat("spearman_hcluster.pdf",ARI_mat_spearman,c(1,2))
ggplot_ARI_mat("spearman_pam.pdf",ARI_mat_spearman,c(3,4))


ARI_mat_uncentered=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-sampleCluster/otherDistance_ARI/ARI_mat_uncentered_hclusterAndPam.csv",row.names = 1,header = T)
ggplot_ARI_mat("uncentered_hcluster.pdf",ARI_mat_uncentered,c(1,2))
ggplot_ARI_mat("uncentered_pam.pdf",ARI_mat_uncentered,c(3,4))

################################result in paper####################################
root_win_times=c()
root_lose_times=c()
ARI_mat_pearson=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-sampleCluster/otherDistance_ARI/ARI_mat_pearson_hclusterAndPam_normalized.csv",row.names = 1,header = T)
ggplot_ARI_mat("pearson_hcluster_normalized.pdf",ARI_mat_pearson,c(1,2))
ggplot_ARI_mat("pearson_pam_normalized.pdf",ARI_mat_pearson,c(3,4))

A=ARI_mat_uncentered[,1]-ARI_mat_uncentered[,2]
length(which(A>0)) #root win
root_win_times=c(root_win_times,length(which(A>0)))
length(which(A<0)) #root lose
root_lose_times=c(root_lose_times,length(which(A<0)) )
B=ARI_mat_uncentered[,3]-ARI_mat_uncentered[,4]
length(which(B>0)) #root win
root_win_times=c(root_win_times,length(which(B>0)))
length(which(B<0)) #root lose
root_lose_times=c(root_lose_times,length(which(B<0)) )

ARI_mat_spearman=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-sampleCluster/otherDistance_ARI/ARI_mat_spearman_hclusterAndPam_normalized.csv",row.names = 1,header = T)
ggplot_ARI_mat("spearman_hcluster_normalized.pdf",ARI_mat_spearman,c(1,2))
ggplot_ARI_mat("spearman_pam_normalized.pdf",ARI_mat_spearman,c(3,4))
A=ARI_mat_uncentered[,1]-ARI_mat_uncentered[,2]
length(which(A>0)) #root win
root_win_times=c(root_win_times,length(which(A>0)))
length(which(A<0)) #root lose
root_lose_times=c(root_lose_times,length(which(A<0)) )
B=ARI_mat_uncentered[,3]-ARI_mat_uncentered[,4]
length(which(B>0)) #root win
root_win_times=c(root_win_times,length(which(B>0)))
length(which(B<0)) #root lose
root_lose_times=c(root_lose_times,length(which(B<0)) )


ARI_mat_uncentered=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-sampleCluster/otherDistance_ARI/ARI_mat_uncentered_hclusterAndPam_normalized.csv",row.names = 1,header = T)
ggplot_ARI_mat("uncentered_hcluster_normalized.pdf",ARI_mat_uncentered,c(1,2))
ggplot_ARI_mat("uncentered_pam_normalized.pdf",ARI_mat_uncentered,c(3,4))

A=ARI_mat_uncentered[,1]-ARI_mat_uncentered[,2]
length(which(A>0)) #root win
root_win_times=c(root_win_times,length(which(A>0)))
length(which(A<0)) #root lose
root_lose_times=c(root_lose_times,length(which(A<0)) )
B=ARI_mat_uncentered[,3]-ARI_mat_uncentered[,4]
length(which(B>0)) #root win
root_win_times=c(root_win_times,length(which(B>0)))
length(which(B<0)) #root lose
root_lose_times=c(root_lose_times,length(which(B<0)) )

ARI_mat_uncentered=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-sampleCluster/otherDistance_ARI/ARI_mat_pearson_compareSquare_hclusterAndPam_normalized.csv",row.names = 1,header = T)
ggplot_ARI_mat("compareSquare_hcluster_normalized.pdf",ARI_mat_uncentered,c(1,2),"root","square")
ggplot_ARI_mat("compareSquare_pam_normalized.pdf",ARI_mat_uncentered,c(3,4),"root","square")



##############################result in paper 2019.4.24######################################
root_win_times=c()
root_lose_times=c()
ARI_mat_uncentered=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-sampleCluster/otherDistance_ARI/ARI_mat_pearson_compareSquare_hclusterAndPam_normalized_checkAgain.csv",row.names = 1,header = T)
ggplot_ARI_mat("compareSquare_hcluster_normalized_2019.4.24.pdf",ARI_mat_uncentered,c(1,2),"2root","1square")
ggplot_ARI_mat("compareSquare_pam_normalized_2019.4.24.pdf",ARI_mat_uncentered,c(3,4),"2root","1square")

A=ARI_mat_uncentered[,1]-ARI_mat_uncentered[,2]
length(which(A>0)) #root win
root_win_times=c(root_win_times,length(which(A>0)))
length(which(A<0)) #root lose
root_lose_times=c(root_lose_times,length(which(A<0)) )
B=ARI_mat_uncentered[,3]-ARI_mat_uncentered[,4]
length(which(B>0)) #root win
root_win_times=c(root_win_times,length(which(B>0)))
length(which(B<0)) #root lose
root_lose_times=c(root_lose_times,length(which(B<0)) )

ARI_mat_uncentered=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-sampleCluster/otherDistance_ARI/ARI_mat_spearman_compareSquare_hclusterAndPam_normalized_2019.4.24.csv",row.names = 1,header = T)
ggplot_ARI_mat("compareSquare_spearman_hcluster_normalized_2019.4.24.pdf",ARI_mat_uncentered,c(1,2),"2root","1square")
ggplot_ARI_mat("compareSquare_spearman_pam_normalized_2019.4.24.pdf",ARI_mat_uncentered,c(3,4),"2root","1square")

A=ARI_mat_uncentered[,1]-ARI_mat_uncentered[,2]
length(which(A>0)) #root win
root_win_times=c(root_win_times,length(which(A>0)))
length(which(A<0)) #root lose
root_lose_times=c(root_lose_times,length(which(A<0)) )
B=ARI_mat_uncentered[,3]-ARI_mat_uncentered[,4]
length(which(B>0)) #root win
root_win_times=c(root_win_times,length(which(B>0)))
length(which(B<0)) #root lose
root_lose_times=c(root_lose_times,length(which(B<0)) )


ARI_mat_uncentered=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-sampleCluster/otherDistance_ARI/ARI_mat_uncentered_compareSquare_hclusterAndPam_normalized_2019.4.24.csv",row.names = 1,header = T)
ggplot_ARI_mat("compareSquare_uncentered_hcluster_normalized_2019.4.24.pdf",ARI_mat_uncentered,c(1,2),"2root","1square")
ggplot_ARI_mat("compareSquare_uncentered_pam_normalized_2019.4.24.pdf",ARI_mat_uncentered,c(3,4),"2root","1square")

A=ARI_mat_uncentered[,1]-ARI_mat_uncentered[,2]
length(which(A>0)) #root win
root_win_times=c(root_win_times,length(which(A>0)))
length(which(A<0)) #root lose
root_lose_times=c(root_lose_times,length(which(A<0)) )
B=ARI_mat_uncentered[,3]-ARI_mat_uncentered[,4]
length(which(B>0)) #root win
root_win_times=c(root_win_times,length(which(B>0)))
length(which(B<0)) #root lose
root_lose_times=c(root_lose_times,length(which(B<0)) )


####################################################################

ARI_mat_pearson=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-sampleCluster/otherDistance_ARI/ARI_mat_pearson_hclusterAndPam_normalized2.csv",row.names = 1,header = T)
ggplot_ARI_mat("pearson_hcluster_normalized2.pdf",ARI_mat_pearson,c(1,2))
ggplot_ARI_mat("pearson_pam_normalized2.pdf",ARI_mat_pearson,c(3,4))

ARI_mat_spearman=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-sampleCluster/otherDistance_ARI/ARI_mat_spearman_hclusterAndPam_normalized2.csv",row.names = 1,header = T)
ggplot_ARI_mat("spearman_hcluster_normalized2.pdf",ARI_mat_spearman,c(1,2))
ggplot_ARI_mat("spearman_pam_normalized2.pdf",ARI_mat_spearman,c(3,4))


ARI_mat_uncentered=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-sampleCluster/otherDistance_ARI/ARI_mat_uncentered_hclusterAndPam_normalized2.csv",row.names = 1,header = T)
ggplot_ARI_mat("uncentered_hcluster_normalized2.pdf",ARI_mat_uncentered,c(1,2))
ggplot_ARI_mat("uncentered_pam_normalized2.pdf",ARI_mat_uncentered,c(3,4))

ARI_mat_uncentered=read.csv("/Users/jiaxchen2/Desktop/2018.10.9/experiment-sampleCluster/otherDistance_ARI/ARI_mat_pearson_compareSquare_hclusterAndPam_normalized2.csv",row.names = 1,header = T)
ggplot_ARI_mat("compareSquare_hcluster_normalized2.pdf",ARI_mat_uncentered,c(1,2),"root","square")
ggplot_ARI_mat("compareSquare_pam_normalized2.pdf",ARI_mat_uncentered,c(3,4),"root","square")

  