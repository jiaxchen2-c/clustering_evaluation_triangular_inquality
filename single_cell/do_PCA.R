
expr_file_list="filtered_expr_list.txt"
out_dir_list="out_dir_list.txt"


expr_fileList=as.matrix(read.table(expr_file_list,header = F))
out_list=as.matrix(read.table(out_dir_list,header = F))


for(i in c(1:length(expr_fileList)))
{

  expr_filei=expr_fileList[i]
  out_diri=out_list[i]
  print(paste("expr_filei",expr_filei))
  print(paste("out_diri",out_diri))
  
  expr=read.table(expr_filei,header = T,sep="\t",row.names = 1)
  
  t_expr=t(expr)
  
  expr.pca <- prcomp(t_expr, center = TRUE,scale. = TRUE)
  
  out_pca=expr.pca$x[,1:50]
  tout_pca=t(out_pca)
  write.table(tout_pca,file=paste(out_diri,"_PCA_X.txt",sep=""),sep="\t",quote=F)
}





