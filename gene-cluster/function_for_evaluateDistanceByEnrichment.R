
compareMethod=function(GO_terms_mat1,GO_terms_mat2)
{
  tmp=function(){
  if(is.na(GO_terms_mat1))
  {
    return(list(comparison_r1_r2=NA,
                 times_r1_smaller=NA,
                 times_r2_smaller=NA))
  }
  if(is.na(GO_terms_mat2))
  {
    return(list(comparison_r1_r2=NA,
                times_r1_smaller=NA,
                times_r2_smaller=NA))
  }
  }
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
  overlap_term_pvalue1=c()
  overlap_term_pvalue2=c()
  for(i in c(1:length(overlap_terms)))
  {
    A=which(terms1==overlap_terms[i])
    B=which(terms2==overlap_terms[i])
    overlap_term_pvalue1=c(overlap_term_pvalue1, min(pvalues1[A]))
    overlap_term_pvalue2=c(overlap_term_pvalue2, min(pvalues2[B]))
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
  
  t_test_GO_terms = t.test(overlap_term_pvalue1,overlap_term_pvalue2)
  t_test_pvalue = t_test_GO_terms$p.value
  comparison_r1_r2=log(times_r1_smaller/times_r2_smaller)
  comparison_r1_r2
  list(comparison_r1_r2=comparison_r1_r2,
       times_r1_smaller=times_r1_smaller,
       times_r2_smaller=times_r2_smaller,
       t_test_pvalue)
}


do_compare_distance_in_pearson_compareSquare_hcluster=function(expr,output_dir,data_prefix,k1=8,k2=8)
{
  GO_terms_absolute_average=do_one_cluster(expr,output_dir,'root square','average',data_prefix,k1)
  GO_terms_root_absolute_average=do_one_cluster(expr,output_dir,'root absolute','average',data_prefix,k2)
  average_absolute_compare_root=compareMethod(GO_terms_absolute_average,GO_terms_root_absolute_average)
  average_absolute_compare_root
}

do_compare_distance_in_pearson_compareSquare_kmedoids=function(expr,output_dir,data_prefix,k1=8,k2=8)
{
  GO_terms_absolute_average=do_one_cluster(expr,output_dir,'root square','kmedoids',data_prefix,k1)
  GO_terms_root_absolute_average=do_one_cluster(expr,output_dir,'root absolute','kmedoids',data_prefix,k2)
  average_absolute_compare_root=compareMethod(GO_terms_absolute_average,GO_terms_root_absolute_average)
  average_absolute_compare_root
}

do_compare_distance_in_pearson_compareSquare_pam=function(expr,output_dir,data_prefix,k1=8,k2=8)
{
  GO_terms_absolute_average=do_one_cluster(expr,output_dir,'root square','pam',data_prefix,k1)
  GO_terms_root_absolute_average=do_one_cluster(expr,output_dir,'root absolute','pam',data_prefix,k2)
  average_absolute_compare_root=compareMethod(GO_terms_absolute_average,GO_terms_root_absolute_average)
  average_absolute_compare_root
}

do_compare_distance_in_average_hcluster=function(expr,output_dir,data_prefix,k1=8,k2=8)
{
  GO_terms_absolute_average=do_one_cluster(expr,output_dir,'absolute','average',data_prefix,k1)
  GO_terms_root_absolute_average=do_one_cluster(expr,output_dir,'root absolute','average',data_prefix,k2)
  average_absolute_compare_root=compareMethod(GO_terms_absolute_average,GO_terms_root_absolute_average)
  average_absolute_compare_root
}
do_compare_distance_in_spearman_hcluster=function(expr,output_dir,data_prefix,k1=8,k2=8)
{
  GO_terms_absolute_average=do_one_cluster(expr,output_dir,'spearman','average',data_prefix,k1)
  GO_terms_root_absolute_average=do_one_cluster(expr,output_dir,'root spearman','average',data_prefix,k2)
  average_absolute_compare_root=compareMethod(GO_terms_absolute_average,GO_terms_root_absolute_average)
  average_absolute_compare_root
}
do_compare_distance_in_uncentered_hcluster=function(expr,output_dir,data_prefix,k1=8,k2=8)
{
  GO_terms_absolute_average=do_one_cluster(expr,output_dir,'uncentered','average',data_prefix,k1)
  GO_terms_root_absolute_average=do_one_cluster(expr,output_dir,'root uncentered','average',data_prefix,k2)
  average_absolute_compare_root=compareMethod(GO_terms_absolute_average,GO_terms_root_absolute_average)
  average_absolute_compare_root
}

do_compare_distance_in_kmedoids=function(expr,output_dir,data_prefix,k1=8,k2=8)
{
  GO_terms_absolute=do_one_cluster(expr,output_dir,'absolute','kmedoids',data_prefix,k1)
  GO_terms_root_absolute=do_one_cluster(expr,output_dir,'root absolute','kmedoids',data_prefix,k2)
  average_absolute_compare_root=compareMethod(GO_terms_absolute,GO_terms_root_absolute)
  average_absolute_compare_root
}

do_compare_distance_in_spearman_kmedoids=function(expr,output_dir,data_prefix,k1=8,k2=8)
{
  GO_terms_absolute=do_one_cluster(expr,output_dir,'spearman','kmedoids',data_prefix,k1)
  GO_terms_root_absolute=do_one_cluster(expr,output_dir,'root spearman','kmedoids',data_prefix,k2)
  average_absolute_compare_root=compareMethod(GO_terms_absolute,GO_terms_root_absolute)
  average_absolute_compare_root
}

do_compare_distance_in_uncentered_kmedoids=function(expr,output_dir,data_prefix,k1=8,k2=8)
{
  GO_terms_absolute=do_one_cluster(expr,output_dir,'uncentered','kmedoids',data_prefix,k1)
  GO_terms_root_absolute=do_one_cluster(expr,output_dir,'root uncentered','kmedoids',data_prefix,k2)
  average_absolute_compare_root=compareMethod(GO_terms_absolute,GO_terms_root_absolute)
  average_absolute_compare_root
}

do_compare_distance_in_pam=function(expr,output_dir,data_prefix,k1=8,k2=8)
{
  GO_terms_absolute=do_one_cluster(expr,output_dir,'absolute','pam',data_prefix,k1)
  GO_terms_root_absolute=do_one_cluster(expr,output_dir,'root absolute','pam',data_prefix,k2)
  average_absolute_compare_root=compareMethod(GO_terms_absolute,GO_terms_root_absolute)
  average_absolute_compare_root
}

do_compare_distance_in_spearman_pam=function(expr,output_dir,data_prefix,k1=8,k2=8)
{
  GO_terms_absolute=do_one_cluster(expr,output_dir,'spearman','pam',data_prefix,k1)
  GO_terms_root_absolute=do_one_cluster(expr,output_dir,'root spearman','pam',data_prefix,k2)
  average_absolute_compare_root=compareMethod(GO_terms_absolute,GO_terms_root_absolute)
  average_absolute_compare_root
}

do_compare_distance_in_uncentered_pam=function(expr,output_dir,data_prefix,k1=8,k2=8)
{
  GO_terms_absolute=do_one_cluster(expr,output_dir,'uncentered','pam',data_prefix,k1)
  GO_terms_root_absolute=do_one_cluster(expr,output_dir,'root uncentered','pam',data_prefix,k2)
  average_absolute_compare_root=compareMethod(GO_terms_absolute,GO_terms_root_absolute)
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

do_select_k_gap=function(expr,output_dir,distance_method,cluster_method,data_prefix)
{
  #distance_method="absolute" or "root absolute"
  #cluster_method="single" or "complete" or "average"
  prefix=paste(data_prefix,distance_method,cluster_method,sep="_")
  distanceMatrices1=produceDistanceMatrices(expr,distance_method)
  bestK=selectK_gap(expr,cluster_method,kmax=round(sqrt(dim(expr)[1])),distanceMatrices1,paste(output_dir,prefix,"CHIndex.pdf",sep=""),distance_method)
  bestK
}

do_one_cluster=function(expr,output_dir,distance_method,cluster_method,data_prefix,k)
{
  #distance_method="absolute" or "root absolute"
  #cluster_method="single" or "complete" or "average"
  prefix=paste(data_prefix,distance_method,cluster_method,sep="_")
  distanceMatrices1=produceDistanceMatrices(expr,distance_method)
  #selectK(expr,cluster_method,kmax=round(sqrt(dim(expr)[1])),distanceMatrices1,paste(output_dir,prefix,"CHIndex.pdf",sep=""))
  
  memb1=doCluster(distanceMatrices1,cluster_method,k,paste(output_dir,prefix,"dendrogram.pdf",sep=""))
  write.table(cbind(as.character(geneName),memb1),file = paste(output_dir,prefix,".csv",sep=""),quote=F,sep="",row.names = F,col.names = F)
  print(memb1)
  GO_terms=doEnrichment_yeast(geneName,memb1,output_dir,prefix)
  print(GO_terms)
  write.table(GO_terms,file=paste(output_dir,prefix,"GO_allclass_mat.csv",sep=""),quote=F)
  GO_terms
}
