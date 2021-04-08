
compareMethod=function(GO_terms_mat1,GO_terms_mat2)
{
  #compute enrichment quality, calculate fold change?
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
    if(min(pvalues1[A])<min(pvalues2[B]))
    {
      times_r1_smaller=times_r1_smaller+1
    }else{
      if(min(pvalues1[A])>min(pvalues2[B]))
      {
        times_r2_smaller=times_r2_smaller+1
      }
      
    }
  }
  comparison_r1_r2=log(times_r1_smaller/times_r2_smaller)
  comparison_r1_r2
}



do_compare_distance_in_average_hcluster=function(expr,output_dir,data_prefix,k1=8,k2=8)
{
  GO_terms_absolute_average=do_one_cluster(expr,output_dir,'absolute','average',data_prefix,k1)
  GO_terms_root_absolute_average=do_one_cluster(expr,output_dir,'root absolute','average',data_prefix,k2)
  average_absolute_compare_root=compareMethod(GO_terms_absolute_average,GO_terms_root_absolute_average)
  average_absolute_compare_root
}


do_select_k=function(expr,output_dir,distance_method,cluster_method,data_prefix)
{
  #distance_method="absolute" or "root absolute"
  #cluster_method="single" or "complete" or "average"
  prefix=paste(data_prefix,distance_method,cluster_method,sep="_")
  distanceMatrices1=produceDistanceMatrices(expr,distance_method)
  selectK(expr,cluster_method,kmax=round(sqrt(dim(expr)[1])),distanceMatrices1,paste(output_dir,prefix,"CHIndex.pdf",sep=""))
}

do_one_cluster=function(expr,output_dir,distance_method,cluster_method,data_prefix,k)
{
  #distance_method="absolute" or "root absolute"
  #cluster_method="single" or "complete" or "average"
  prefix=paste(data_prefix,distance_method,cluster_method,sep="_")
  distanceMatrices1=produceDistanceMatrices(expr,distance_method)
  selectK(expr,cluster_method,kmax=round(sqrt(dim(expr)[1])),distanceMatrices1,paste(output_dir,prefix,"CHIndex.pdf",sep=""))
  
  memb1=doCluster(distanceMatrices1,cluster_method,k,paste(output_dir,prefix,"dendrogram.pdf",sep=""))
  write.table(cbind(as.character(geneName),memb1),file = paste(output_dir,prefix,".csv",sep=""),quote=F,sep="",row.names = F,col.names = F)
  GO_terms=doEnrichment_yeast(geneName,memb1,output_dir,prefix)
  GO_terms
}
