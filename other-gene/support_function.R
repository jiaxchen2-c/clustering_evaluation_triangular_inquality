







do_geneCluster_oneTest_bootstrapTest=function(fileName,distance_method,cluster_method,output_dir,ith,k)
{
  data=read.table(fileName,header = T,sep="\t")
  geneName=as.matrix(data[,2])
  geneName=convertToID(geneName,' ')
  expr=data[,3:dim(data)[2]]
  
  #compute distance
  prefix=paste(ith,distance_method,cluster_method,sep="_")
  distanceMatrices1=produceDistanceMatrices(expr,distance_method)
  write.table(distanceMatrices1,file=paste(output_dir,prefix,"distanceMatrix.csv",sep=""),quote=F,row.names = F,col.names = F,sep=",")
  
  boot_result=ClusterBootstrap_fromDist(distanceMatrices1, k, noise.cut = 0, bootstrap = 100, 
                                        dissolve = .5, cluster_method)
  boot_result$partition
  
  write.table(cbind(as.character(geneName),boot_result$partition),file = paste(output_dir,prefix,".csv",sep=""),quote=F,sep=",",row.names = F,col.names = F)
  
  
  print(prefix)
  print(paste("bootmean",boot_result$bootmean,"bootdissolved", boot_result$bootdissolved))
  
}

convertToID=function(toID,sep)
{
  resultID=c()
  for(i in c(1:dim(toID)[1]))
  {
    resultID[i]=unlist(strsplit(toID[i],sep))[1]
  }
  resultID
}