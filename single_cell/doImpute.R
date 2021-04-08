library(scImpute)

expr_file_list="filtered_expr_list_noNormalized.txt"
label_file_list="filtered_partition_list.txt"
out_dir_list="out_dir_list.txt"


fileList=as.matrix(read.table(label_file_list,header = F))
expr_fileList=as.matrix(read.table(expr_file_list,header = F))
out_list=as.matrix(read.table(out_dir_list,header = F))


for(i in c(2:length(fileList)))
{
  partition_file=fileList[i]
  data=read.table(partition_file,header = F)
  ref_cluster=data[,2]
  unique_ref_cluster=unique(ref_cluster)
 
  Kcluster=length(unique_ref_cluster)
  expr_filei=expr_fileList[i]
  out_diri=out_list[i]
  print(paste("Kcluster",Kcluster))
  print(paste("expr_filei",expr_filei))
  print(paste("out_diri",out_diri))
  
  scimpute(count_path=expr_filei,infile="txt",outfile="txt",out_dir=out_diri,labeled=F,drop_thre=0.5,Kcluster=Kcluster,ncores=10)
}

