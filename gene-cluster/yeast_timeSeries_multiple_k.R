#source and library
library(dplyr)
library(tidyr)
library(ggplot2)
library(kmed)
library(clues)

library(DBI)
library(AnnotationDbi)
library(Biobase)
library(GO.db)


library(RSQLite)#version2.1.0
library(org.Hs.eg.db)
library(graph)
library(Category)
library(GOstats)
library(KEGG.db)
library(org.Sc.sgd.db)




source("function_for_cluster.R")
source("function_for_distance.R")
source("function_for_enrichment.R")
source("function_for_evaluateDistanceByEnrichment.R")



###################################################################################
##################################define functions#################################
convertToID=function(toID,sep)
{
  resultID=c()
  for(i in c(1:dim(toID)[1]))
  {
    resultID[i]=unlist(strsplit(toID[i],sep))[1]
  }
  resultID
}

plot_select_CHIndex=function(distance,cluster_method,output_dir)
{
  for(i in c(1:length(fileList)))
  {
    print(paste("i",i))
    #loop
    #file1="01_alpha_factor.csv"
    filei=fileList[i]
    data=read.table(paste(input_dir,filei,sep=""),header = T,sep="\t")
    geneName=as.matrix(data[,1])
    geneName=convertToID(geneName,' ')
    expr=data.matrix(data[,2:dim(data)[2]])
    
    do_select_k(expr,output_dir,distance,cluster_method,filei)
  }
}

plot_select_gap=function(distance,cluster_method,output_dir)
{
  bestK_mat=c()
  for(i in c(1:length(fileList)))
  {
    print(paste("i",i))
    #loop
    #file1="01_alpha_factor.csv"
    filei=fileList[i]
    data=read.table(paste(input_dir,filei,sep=""),header = T,sep="\t")
    geneName=as.matrix(data[,1])
    geneName=convertToID(geneName,' ')
    expr=data.matrix(data[,2:dim(data)[2]])
    
    
    
    bestK=do_select_k_gap(expr,output_dir,distance,cluster_method,filei)
    bestK_mat=rbind(bestK_mat,c(i,bestK))
  }
  
  bestK_mat
}

##################################end define functions#############################
###################################################################################

#parameter and load data
input_dir="benchmark_datasets_normalized/"
output_dir="timeCourse-hcluster-result-kmedoids-11_7/"
list_file=paste(input_dir,"list.txt",sep="")
toID=as.matrix(read.table("geneName2.txt",header=F))
fromID=as.character(as.matrix(read.table("OLN.txt",header=F)))
toID=convertToID(toID,';')

#fileList
fileList=as.matrix(read.table(list_file,header = F))
root_win=0
root_lose=0


absolute_vs_root_list=c()
#checkAgain=c(3,7,9,10,12,13,14,16)
checkAgain=c(9)#c(13,14,15,17)#remove 9 temp
runAgain=c(16,17)


###########################################################################
###########################################################################
##############################################################



work_step1=function()
{
  
  list_file=paste(input_dir,"list.txt",sep="")
  toID=as.matrix(read.table("geneName2.txt",header=F))
  fromID=as.character(as.matrix(read.table("OLN.txt",header=F)))
  toID=convertToID(toID,';')
  
  input_dir="benchmark_datasets_normalized/"
  list_file=paste(input_dir,"list.txt",sep="")
  toID=as.matrix(read.table("geneName2.txt",header=F))
  fromID=as.character(as.matrix(read.table("OLN.txt",header=F)))
  toID=convertToID(toID,';')
  
  
  #fileList
  fileList=as.matrix(read.table(list_file,header = F))
  
  
  #plot, then for each method, give candidate k (at least peak)
  cluster_methods=c("pam","average")
  distance_methods=c("root absolute","absolute","root square", "spearman","root spearman", "square spearman", 
                     "uncentered", "root uncentered", "square uncentered")
  
  for(i in 1:length(cluster_methods))
  {
    for(j in 1:length(distance_methods))
    {
      cluster_method=cluster_methods[i]
      distance_method=distance_methods[j]
      one_run_combine_method_dist_step1(cluster_method,distance_method)

    }
  }
  
  #read result, then store the candidate k in file, each dataset a line, split by ,
  #store in timeCourse_normalized_multiK_2019sep_step1_result
  #mkdir "timeCourse_normalized_multiK_2019sep_step1_result/"
  #file = paste("timeCourse_normalized_multiK_2019sep_step1_result/",cluster_method,"_",distance_method,sep="")
  
}

work_step2=function(k_all=T)
{
  #run GO enrichment for each runn (cluster_method,distance_method,dataset,k), store complete GO table
  
  
  #parameter and load data
  output_dir="timeCourse_normalized_multiK_2019sep_step2_rerun/"
  dir.create(output_dir)
  
  
  input_dir="benchmark_datasets_normalized/"
  list_file=paste(input_dir,"list.txt",sep="")
  toID=as.matrix(read.table("geneName2.txt",header=F))
  fromID=as.character(as.matrix(read.table("OLN.txt",header=F)))
  toID=convertToID(toID,';')
  
  
  #fileList
  fileList=as.matrix(read.table(list_file,header = F))
  
  #read the k table, for each data a k_array
  #end in 16_square spearman_pam_12
  
  #cluster_methods=c("average") 
  cluster_methods=c("pam")
  distance_methods=c("root absolute", "absolute","root square", "spearman","root spearman", "square spearman")
  #distance_methods=c("square spearman")
  
  for(cluster_method in cluster_methods)
  {
    for(distance_method in distance_methods)
    {
      
      k_array_file = paste("timeCourse_normalized_multiK_2019sep_step1_result/",cluster_method,"_",distance_method,sep="")
      k_table=as.matrix(read.table(k_array_file,header=F,sep="\n") )
      
      #deal with format for k_array!!!!!!!!!!!!!!!!!!
      #for(i in c(1:length(fileList))) 
      for(i in c(8:length(fileList))) 
      {
        
        if(k_all)
        {
          k_array=seq(5,26,by=3)
        }else{
          k_array=as.numeric(unlist(strsplit(as.character(k_table[i]),split=",")))
        }
        
        one_run_enrichment_step2(fileList,cluster_method,distance_method,i,k_array) 
        
      }
      
    }
  }
  
  
  
  
}




one_run_combine_method_dist_step1=function(cluster_method,distance_method)
{
  print(cluster_method)
  print(distance_method)
  
  output_dir=paste("timeCourse_normalized_multiK_2019sep_step1_",cluster_method,"_",distance_method,"/",sep="")
  print(output_dir)
  dir.create(output_dir)
  plot_select_CHIndex(distance_method,cluster_method,output_dir)
  plot_select_CHIndex(distance_method,cluster_method,output_dir)
  
  #plot, then for each method, give candidate k (at least peak)
  

}


one_run_enrichment_step2=function(fileList,cluster_method,distance_method,dataseti,k_array)
{
  output_dir=paste("timeCourse_normalized_multiK_2019sep_step2_rerun_",dataseti,"/",sep="")
  dir.create(output_dir)
  
  
  
  input_dir="benchmark_datasets_normalized/"
  
  filei=fileList[dataseti]
  data=read.table(paste(input_dir,filei,sep=""),header = T,sep="\t")
  geneName=as.matrix(data[,1])
  geneName=convertToID(geneName,' ')
  expr=data.matrix(data[,2:dim(data)[2]])

  print(k_array)
  for(k in k_array)
  {
    data_prefix=paste(dataseti,distance_method,cluster_method,k,sep="_")
    print(data_prefix)
    do_one_cluster(expr,output_dir,distance_method,cluster_method,data_prefix,k,geneName)
  }
}



one_run_compare_step3=function(enrichment_dir_header,cluster_method1, distance_method1, k_array1, cluster_method2, distance_method2, k_array2, dataseti)
{
  #enrichment_dir="timeCourse_normalized_multiK_2019sep_step2_rerun_"
  #consider how to deal with multiple k  
  result=c()
  input_dir=paste(enrichment_dir_header,dataseti,"/",sep="")
  compare_pair_list=c()
  for(k1 in k_array1)
  {
    for(k2 in k_array2)
    {
      data_prefix=paste(dataseti,distance_method1,cluster_method1,k1,sep="_")
      prefix=paste(data_prefix,distance_method1,cluster_method1,sep="_")
      file=paste(input_dir,prefix,"GO_allclass_mat.csv",sep="")
      GO_term1 = read.table(file,header=T)
      
      
      data_prefix=paste(dataseti,distance_method2,cluster_method2,k2,sep="_")
      prefix=paste(data_prefix,distance_method2,cluster_method2,sep="_")
      file=paste(input_dir,prefix,"GO_allclass_mat.csv",sep="")
      GO_term2 = read.table(file,header=T)
      
      result1=compareMethod_t_test(GO_term1,GO_term2)
      tmp=c(result1$comparison_r1_r2,result1$diff_pvalue)
      #print(paste(k1,k2,dataseti,cluster_method1,distance_method1,distance_method2))
      #print(result1)
      compare_pair_list=c(compare_pair_list,paste(dataseti,k1,k2,sep="_"))
      result=cbind(result,result1)
      #print(result)
    }
  }
  colnames(result)=compare_pair_list
  result
}


compareMethod_t_test=function(GO_terms_mat1,GO_terms_mat2)
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
       diff_pvalue=t_test_pvalue)
}


work_step3=function(output_file="multiple_all_k_compare_result.txt",k_all=T,enrichment_dir_header)
{
  #run GO enrichment for each runn (cluster_method,distance_method,dataset,k), store complete GO table
  
  
  #parameter and load data
  output_dir="timeCourse_normalized_multiK_2019sep_step3/"
  dir.create(output_dir)

  #read the k table, for each data a k_array
  #end in 16_square spearman_pam_12
  
  
  cluster_methods=c("average")
  distance_methods=c("root absolute", "absolute","root square", "spearman","root spearman", "square spearman")
  for(cluster_method in cluster_methods)
  {
    for(d1 in 1:length(distance_methods))
    {
      for(d2 in 1:length(distance_methods))
      {
        if(d2>d1)
        {
          
          distance_method1=distance_methods[d1]
          distance_method2=distance_methods[d2]
          
          
          
          k_array_file1 = paste("timeCourse_normalized_multiK_2019sep_step1_result/",cluster_method,"_",distance_method1,sep="")
          k_table1=as.matrix(read.table(k_array_file1,header=F,sep="\n") )
            
          k_array_file2 = paste("timeCourse_normalized_multiK_2019sep_step1_result/",cluster_method,"_",distance_method2,sep="")
          k_table2=as.matrix(read.table(k_array_file2,header=F,sep="\n") )
            
          
          
          #deal with format for k_array!!!!!!!!!!!!
          for(i in c(1:length(fileList))) 
          {
            if(k_all)
            {
              k_array1=seq(5,26,by=3)
              k_array2=seq(5,26,by=3)
            }else{
              k_array1=as.numeric(unlist(strsplit(as.character(k_table1[i]),split=",")))
              k_array2=as.numeric(unlist(strsplit(as.character(k_table2[i]),split=",")))
              
            }
            
            
            cat(paste("Data",i,"," ,cluster_method,",",distance_method1,",",distance_method2,"\n"),file=output_file,append=TRUE)
            #print(paste(distance_method1,distance_method2,i))
            result=one_run_compare_step3(enrichment_dir_header=enrichment_dir_header,cluster_method,distance_method1,k_array1,cluster_method,distance_method2,k_array2,i) 
            cat(paste(colnames(result),","),file=output_file,append=TRUE)
            cat("\n",file=output_file,append=TRUE)
            for(linei in c(1:dim(result)[1]))
            {
              for(colj in c(1:dim(result)[2]))
              {
                cat(paste(result[linei,colj],","),file=output_file,append=TRUE)
              }
              cat("\n",file=output_file,append=TRUE)
              
            }
            print(paste("Data",i,"," ,cluster_method,",",distance_method1,",",distance_method2,"\n"))
            
            print(result)
            
          } 
        }
      }
      
    }
  }
  
}

work_step4=function(data_file,cluster_methods=c(" pam "))
{
  #cluster_methods=c(" pam ")
  #cluster_methods=c(" average ")
  distance_methods=c(" root absolute ", " absolute "," root square ", " spearman "," root spearman ", " square spearman ")
  mat_for_plot=c()
  for(cluster_method in cluster_methods)
  {
    for(d1 in 1:length(distance_methods))
    {
      for(d2 in 1:length(distance_methods))
      {
        if(d2>d1)
        {
          print(paste(cluster_method, distance_methods[d1], distance_methods[d2]))
          result=for_one_compare(data_file,cluster_method, distance_methods[d1], distance_methods[d2], 3)
          print(result)
          
          line_for_plot=c(distance_methods[d1], distance_methods[d2],d1,d2,result$total_hits, result$total_not_win_hits, result$total_win, result$total_lose)
          mat_for_plot=rbind(mat_for_plot,line_for_plot)
        }
      }
    }
  }
  plot_compare_result_mat(mat_for_plot,cluster_methods)
}


#heatmap plot test
#ggplot_heatmap_continues(noInf_mapping_mat_T_L,timePoints_L,timePoints_L,"mapping_score",color_low="white",color_high="steelblue",is_coord_equal=T)

ggplot_heatmap_continues=function(mat,xlabel=NA,ylabel=NA,fill_legend=NA,color_low="red",color_median="while",color_high="steelblue",is_coord_equal=F)
{
  #mat=noInf_mapping_mat_T_L
  colnames(mat)=c(1:dim(mat)[2])
  rownames(mat)=c(1:dim(mat)[1])
  
  mat <- data.frame(mat)
  mat$id<-rownames(mat)
  melt_mat <- melt(mat, id.var="id")
  
  p1=ggplot(melt_mat,aes(as.integer(variable), as.integer(id)))+
    geom_tile(aes(fill=melt_mat$value),colour="white")+
    scale_fill_gradientn(colours = c("steelblue","white", "firebrick"),
                         values = scales::rescale(c(-1, 0, 1)),lim=c(-10,10))
  
  p1=p1 + scale_x_continuous(breaks=seq(1:length(xlabel)),labels=as.array(xlabel))
  p1=p1 + scale_y_continuous(breaks=seq(1:length(ylabel)),labels=as.array(ylabel))
  p1=p1+theme(axis.text.x = element_text(angle=75, hjust=1)
              ,axis.title.x=element_blank(),axis.title.y = element_blank())
  if(is_coord_equal)
  {
    p1=p1+coord_equal()
  }
  print(p1)
}

plot_compare_result_mat=function(mat_for_plot,fileHeader)
{
  
  df=as.data.frame(mat_for_plot)
  colnames(df)=c("distance_methods1","distance_methods2","d1","d2","total_hits","total_not_win_hits","total_win","total_lose")
  distance_methods=c(" root absolute ", " absolute "," root square ", " spearman "," root spearman ", " square spearman ")
  value_mat=matrix(0,6,6)
  value_mat2=matrix(0,6,6)
  for(i in c(1:dim(mat_for_plot)[1]))
  {
    indexi=as.numeric(as.character(df$d1[i]))
    indexj=as.numeric(as.character(df$d2[i]))
    
    
    value_mat[indexi,indexj]=log2(as.numeric(mat_for_plot[i,5])/as.numeric(mat_for_plot[i,6]))
    #value_mat[indexj,indexi]=log2(as.numeric(mat_for_plot[i,6])/as.numeric(mat_for_plot[i,5]))
    
    value_mat2[indexi,indexj]=log2(as.numeric(mat_for_plot[i,7])/as.numeric(mat_for_plot[i,8]))
    #value_mat2[indexj,indexi]=log2(as.numeric(mat_for_plot[i,8])/as.numeric(mat_for_plot[i,7]))
    
  }
  require(ggplot2)
  require(reshape)
  pdf(paste(fileHeader,"_mat2_geneCluster_pair_heatmap.pdf",sep=""))
  #plot_heatmap(value_mat)
  ggplot_heatmap_continues(value_mat2,xlabel=distance_methods,ylabel=distance_methods,fill_legend=NA,color_low="red",color_median = "white",color_high="steelblue",is_coord_equal=T)
  dev.off()
  
  pdf(paste(fileHeader,"_mat1_geneCluster_pair_heatmap.pdf",sep=""))
  #plot_heatmap(value_mat)
  ggplot_heatmap_continues(value_mat,xlabel=distance_methods,ylabel=distance_methods,fill_legend=NA,color_low="red",color_median = "white",color_high="steelblue",is_coord_equal=T)
  dev.off()
  
  if(F)
  {
    require(ggplot2)
    pdf(paste(fileHeader,"_geneCluster_pair_heatmap.pdf",sep=""))
    
    print(df$value1)
    p1=ggplot(df,aes(x=as.Date(d1),y=as.Date(d2)))+geom_tile(aes(fill=df$value1),colour="white")+scale_fill_gradient("compare score",low="white",high="steelblue")
    p1=p1 + scale_x_continuous(breaks=seq(1:6),labels=as.array(distance_methods))                                                    
    p1=p1 + scale_y_continuous(breaks=seq(1:6),labels=as.array(distance_methods)) 
    p1=p1+theme(axis.text.x = element_text(angle=75, hjust=1)
                ,axis.title.x=element_blank(),axis.title.y = element_blank())
    p1=p1+coord_equal()
    print(p1)
    dev.off()
  }
  
}

get_best_k_of_GO_for_one=function(data_file,cluster_method, distance_method1, data_i)
{
  data=read.csv(data_file,sep="\n",header = F)
  compare_pair=seq(1,dim(data)[1],by=6)
  k_pair_list=seq(2,dim(data)[1],by=6)
  sigh_list=seq(3,dim(data)[1],by=6)
  pvalue_list=seq(6,dim(data)[1],by=6)
  
  total_win=0
  total_lose=0
  total_hits=0
  total_not_win_hits=0
  
  
  return_bestk=0
  for(i in c(1:length(sigh_list)))
  {
    win_for_each_k=rep(0,35)
    pair_this_compare=unlist(strsplit(as.character(data[compare_pair[i],1]),split=","))
    if(pair_this_compare[1]==data_i)
    {
      if(pair_this_compare[2]==cluster_method)
      {
        if(pair_this_compare[3]==distance_method1)
        {
          if(pair_this_compare[4]==distance_method1)
          {
            index=sigh_list[i]
            tmp=unlist(strsplit(as.character(data[index,1]),split=","))
            sigh_this_compare=as.numeric(tmp)
            
            tmp2=unlist(strsplit(as.character(data[pvalue_list[i],1]),split=","))
            pvalue_this_compare=as.numeric(tmp2)
            
            tmp3=unlist(strsplit(as.character(data[k_pair_list[i],1]),split=","))
            
            #select_significant_compare=which(pvalue_this_compare<0.05,arr.ind = T)
            #sigh_this_compare=sigh_this_compare[select_significant_compare]
            
            
            k1_k2_diff_all=c()
            k1_all=c()
            k2_all=c()
            for(j in c(1:length(tmp3)))
            {
              tmp4=unlist(strsplit(tmp3[j],split="_"))
              k1=tmp4[2]
              k2=tmp4[3]
              k1_k2_diff=as.numeric(k1)-as.numeric(k2)
              k1_all=c(k1_all,as.numeric(k1))
              k2_all=c(k2_all,as.numeric(k2))
              k1_k2_diff_all=c(k1_k2_diff_all,k1_k2_diff)
            }
            
            for(k_i in c(1:length(sigh_this_compare)))
            {
              if(!is.na(pvalue_this_compare[k_i]))
              {
                if(pvalue_this_compare[k_i]<0.05)
                {
                  if(sigh_this_compare[k_i]>0)
                  {
                    wins_k=k1_all[k_i]
                  }else{
                    wins_k=k2_all[k_i]
                  }
                  win_for_each_k[wins_k]=win_for_each_k[wins_k]+1
                }
              }else{
                return_bestk=k1_all[k_i]
                
              }
              
            }
            
            bestk=which(win_for_each_k==max(win_for_each_k),arr.ind = T)
            #one_row=c(compare_pair[i],bestk)
            #best_k_for_cluster_distance_method=rbind(best_k_for_cluster_distance_method, one_row)
            if(!is.na(bestk[1]))
            {
              return_bestk=bestk[1]
            }
            
          }
        }
      }
    }
  }
  return_bestk
}

for_one_compare_best_k=function(data_file,data_file_bestk,cluster_method, distance_method1, distance_method2,diff_range=5)
{
  data=read.csv(data_file,sep="\n",header = F)
  compare_pair=seq(1,dim(data)[1],by=6)
  k_pair_list=seq(2,dim(data)[1],by=6)
  sigh_list=seq(3,dim(data)[1],by=6)
  pvalue_list=seq(6,dim(data)[1],by=6)
  
  total_win=0
  total_lose=0
  total_hits=0
  total_not_win_hits=0
  bestk1=0
  bestk2=0
  for(i in c(1:length(sigh_list)))
  {
    pair_this_compare=unlist(strsplit(as.character(data[compare_pair[i],1]),split=","))
    if(pair_this_compare[2]==cluster_method)
    {
      if(pair_this_compare[3]==distance_method1)
      {
        if(pair_this_compare[4]==distance_method2)
        {
          index=sigh_list[i]
          tmp=unlist(strsplit(as.character(data[index,1]),split=","))
          sigh_this_compare=as.numeric(tmp)
          
          tmp2=unlist(strsplit(as.character(data[pvalue_list[i],1]),split=","))
          pvalue_this_compare=as.numeric(tmp2)
          
          tmp3=unlist(strsplit(as.character(data[k_pair_list[i],1]),split=","))
          
          
          bestk1=get_best_k_of_GO_for_one(data_file_bestk,cluster_method, distance_method1, pair_this_compare[1])
          bestk2=get_best_k_of_GO_for_one(data_file_bestk,cluster_method, distance_method2, pair_this_compare[1])
          
          print(paste(pair_this_compare[1],cluster_method, distance_method1, distance_method2))
          print(paste("bestk1",bestk1,"bestk2",bestk2))
          
          if(bestk1!=0)
          {
            if(bestk2!=0)
            {
              k1_k2_diff_all=c()
              k1_all=c()
              k2_all=c()
              for(j in c(1:length(tmp3)))
              {
                tmp4=unlist(strsplit(tmp3[j],split="_"))
                k1=tmp4[2]
                k2=tmp4[3]
                print(paste("kpair list[i]",k_pair_list[i]))
                print("tmp4")
                print(tmp4)
                print(paste("k1",k1,"bestk1",bestk1))
                if(as.numeric(k1)==bestk1){
                  if(as.numeric(k2)==bestk2)
                  {
                    if(pvalue_this_compare[j]<0.05)
                    {
                      if(sigh_this_compare[j]>0)
                      {
                        total_win=total_win+1
                        total_hits=total_hits+1
                      }else{
                        total_lose=total_lose+1
                        total_not_win_hits=total_not_win_hits+1
                      }
                    }
                  }
                }
            }
          }
          
          }
        }
      }
    }
  }
  
  #print(paste("total_hist",total_hits,"total_not_win_hits",total_not_win_hits))
  #print(paste("total_win",total_win,"total_lose",total_lose))
  list(total_hits=total_hits,total_not_win_hits=total_not_win_hits,total_win=total_win,total_lose=total_lose)
}

work_step4_get_best_k_of_GO_for_one=function(data_file,data_file_bestk,cluster_methods=c(" average "))
{
  #cluster_methods=c(" pam ")
  #cluster_methods=c(" average ")
  distance_methods=c(" root absolute ", " absolute "," root square ", " spearman "," root spearman ", " square spearman ")
  mat_for_plot=c()
  for(cluster_method in cluster_methods)
  {
    for(d1 in 1:length(distance_methods))
    {
      for(d2 in 1:length(distance_methods))
      {
        if(d2>d1)
        {
          print(paste(cluster_method, distance_methods[d1], distance_methods[d2]))
          result=for_one_compare_best_k(data_file,data_file_bestk,cluster_method, distance_methods[d1], distance_methods[d2], 3)
          print(result)
          
          line_for_plot=c(distance_methods[d1], distance_methods[d2],d1,d2,result$total_hits, result$total_not_win_hits, result$total_win, result$total_lose)
          mat_for_plot=rbind(mat_for_plot,line_for_plot)
        }
      }
    }
  }
  plot_compare_result_mat(mat_for_plot,cluster_methods)
}

work_step3_get_best_k_of_GO_for_one=function(enrichment_dir_header,output_file="multiple_all_k_compare_result.txt",k_all=F)
{
  #run GO enrichment for each runn (cluster_method,distance_method,dataset,k), store complete GO table
  
  
  #parameter and load data
  output_dir="timeCourse_normalized_multiK_2019sep_step3/"
  dir.create(output_dir)
  
  #read the k table, for each data a k_array
  #end in 16_square spearman_pam_12
  
  
  cluster_methods=c("average")
  distance_methods=c("root absolute", "absolute","root square", "spearman","root spearman", "square spearman")
  for(cluster_method in cluster_methods)
  {
    for(d1 in 1:length(distance_methods))
    {
          distance_method1=distance_methods[d1]
          
          k_array_file1 = paste("timeCourse_normalized_multiK_2019sep_step1_result/",cluster_method,"_",distance_method1,sep="")
          k_table1=as.matrix(read.table(k_array_file1,header=F,sep="\n") )
          
          
          #deal with format for k_array!!!!!!!!!!!!
          for(i in c(1:length(fileList))) 
          {
            if(k_all)
            {
              k_array1=seq(5,26,by=3)
              
            }else{
              k_array1=as.numeric(unlist(strsplit(as.character(k_table1[i]),split=",")))
              
            }
            cat(paste("Data",i,"," ,cluster_method,",",distance_method1,",",distance_method1,"\n"),file=output_file,append=TRUE)
            #print(paste(distance_method1,distance_method2,i))
            result=one_run_compare_step3(enrichment_dir_header,cluster_method,distance_method1,k_array1,cluster_method,distance_method1,k_array1,i)
            cat(paste(colnames(result),","),file=output_file,append=TRUE)
            cat("\n",file=output_file,append=TRUE)
            for(linei in c(1:dim(result)[1]))
            {
              for(colj in c(1:dim(result)[2]))
              {
                cat(paste(result[linei,colj],","),file=output_file,append=TRUE)
              }
              cat("\n",file=output_file,append=TRUE)
              
            }
            print(paste("Data",i,"," ,cluster_method,",",distance_method1,",",distance_method1,"\n"))
            
            print(result)
            }
      
    }
  }
  
}



for_one_compare=function(data_file,cluster_method, distance_method1, distance_method2,diff_range=5)
{
  data=read.csv(data_file,sep="\n",header = F)
  compare_pair=seq(1,dim(data)[1],by=6)
  k_pair_list=seq(2,dim(data)[1],by=6)
  sigh_list=seq(3,dim(data)[1],by=6)
  pvalue_list=seq(6,dim(data)[1],by=6)
  
  total_win=0
  total_lose=0
  total_hits=0
  total_not_win_hits=0
  for(i in c(1:length(sigh_list)))
  {
    pair_this_compare=unlist(strsplit(as.character(data[compare_pair[i],1]),split=","))
    if(pair_this_compare[2]==cluster_method)
    {
      if(pair_this_compare[3]==distance_method1)
      {
        if(pair_this_compare[4]==distance_method2)
        {
          index=sigh_list[i]
          tmp=unlist(strsplit(as.character(data[index,1]),split=","))
          sigh_this_compare=as.numeric(tmp)
          
          tmp2=unlist(strsplit(as.character(data[pvalue_list[i],1]),split=","))
          pvalue_this_compare=as.numeric(tmp2)
          
          tmp3=unlist(strsplit(as.character(data[k_pair_list[i],1]),split=","))
          
          select_significant_compare=which(pvalue_this_compare<0.05,arr.ind = T)
          sigh_this_compare=sigh_this_compare[select_significant_compare]
          
          
          k1_k2_diff_all=c()
          k1_all=c()
          k2_all=c()
          for(j in c(1:length(tmp3)))
          {
            tmp4=unlist(strsplit(tmp3[j],split="_"))
            k1=tmp4[2]
            k2=tmp4[3]
            k1_k2_diff=as.numeric(k1)-as.numeric(k2)
            k1_all=c(k1_all,as.numeric(k1))
            k2_all=c(k2_all,as.numeric(k2))
            k1_k2_diff_all=c(k1_k2_diff_all,k1_k2_diff)
          }
          
          select_k_diff=which(abs(k1_k2_diff)<20,arr.ind=T)
          sigh_this_compare=sigh_this_compare[select_k_diff]
          
          A=which(sigh_this_compare>0,arr.ind = T)
          total_win=total_win+length(A)
          B=which(sigh_this_compare<0,arr.ind = T)
          total_lose=total_lose+length(B)
          if(length(A)>length(B))
          {
            total_hits=total_hits+1
          }
          
          if(length(A)<length(B)){
            total_not_win_hits=total_not_win_hits+1
          }
        }
      }
    }
  }
  print(paste("total_hist",total_hits,"total_not_win_hits",total_not_win_hits))
  print(paste("total_win",total_win,"total_lose",total_lose))
  list(total_hits=total_hits,total_not_win_hits=total_not_win_hits,total_win=total_win,total_lose=total_lose)
}

################# end function definition #########################

################# begin work ######################################

#note: modify the input and outputdir if want to change k_all
#work_step1() #estimate k
#work_step2(k_all=T) # do enrichment
#for all k
work_step3(output_file="multiple_all_k_compare_result.txt",k_all=T,enrichment_dir_header="timeCourse_normalized_multiK_2019sep_step2_rerun_") # collect enrichment
work_step4(data_file="multiple_all_k_compare_result.txt", " average ")

#for compare by select their best k(from all k) by GO
work_step3(output_file="multiple_all_k_compare_result.txt",k_all=T,enrichment_dir_header="timeCourse_normalized_multiK_2019sep_step2_rerun_") # collect enrichment
work_step3_get_best_k_of_GO_for_one(enrichment_dir_header="timeCourse_normalized_multiK_2019sep_step2_rerun_", output_file="multiple_k_compare_result4.txt",k_all=T)
work_step4_get_best_k_of_GO_for_one(data_file="multiple_all_k_compare_result.txt", data_file_bestk="multiple_k_compare_result4.txt",cluster_methods=c(" average "))

#for compare by select their best k(from chindex k) by GO
work_step3(output_file="multiple_CHIndex_k_compare_result.txt",k_all=F,enrichment_dir_header="timeCourse_normalized_multiK_2019sep_step2_") # collect enrichment
work_step3_get_best_k_of_GO_for_one(enrichment_dir_header="timeCourse_normalized_multiK_2019sep_step2_", output_file="multiple_k_compare_result5.txt",k_all=F)
work_step4_get_best_k_of_GO_for_one(data_file="multiple_CHIndex_k_compare_result.txt", data_file_bestk="multiple_k_compare_result5.txt",cluster_methods=c(" average "))
