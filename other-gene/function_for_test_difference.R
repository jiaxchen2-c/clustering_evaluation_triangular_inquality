library(reshape2)

##################################define functions##############################

load_data_and_normalize=function(fileName,norm=T)
{
  #load data and reference cluster
  data=read.table(fileName,header = T,sep="\t")
  
  
  geneID=as.matrix(data[2:dim(data)[1],1])
  ref_cluster=c(as.matrix(data[1,c(2:dim(data)[2])]))
  ref_cluster_array=rep(0,length(ref_cluster))
  unique_ref_cluster=unique(as.character(ref_cluster))
  for(j in c(1:length(unique_ref_cluster)))
  {
    A=which(ref_cluster==unique_ref_cluster[j],arr.ind = T)
    ref_cluster_array[A]=j
  }
  
  expr=as.matrix(data[c(2:dim(data)[1]),c(2:dim(data)[2])])
  class(expr)="numeric"
  sampleID=as.matrix(colnames(data)[2:dim(data)[2]])
  B=which(is.na(expr),arr.ind = T)
  if(dim(B)[1]>0)
  {
    selectSample=setdiff(c(1:length(sampleID)),B[,2])
    selectGene=setdiff(c(1:length(geneID)),B[,1])
    ref_cluster_array=ref_cluster_array[selectSample]
    geneID=geneID[selectGene]
    sampleID=sampleID[selectSample]
    expr=expr[selectGene,selectSample]
  }
  expr=t(expr)
  if(norm)
  {
    expr=scale(expr)
  }
  
  
  list(expr=expr,
       sampleID=sampleID,
       geneID=geneID,
       ref_cluster_array=ref_cluster_array,
       ref_cluster=ref_cluster
  )
}

iteration_look_for_disobeyTrianglularInequality=function(distance_mat,colOrder=NA)
{
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
        
        print(paste("ij",i,j,colOrder[i],colOrder[j],"disobey trianglur inequality"))
        
        for(k in useful_index[A])
        {
          
          print(paste("(",colOrder[i],colOrder[j],")",distance_mat[i,j],"(",colOrder[i],colOrder[k],")",distance_mat[i,k],"(",colOrder[j],colOrder[k],")",distance_mat[k,j]))
        }
        for(k in useful_index[B])
        {
          print(paste("(",colOrder[i],colOrder[j],")",distance_mat[i,j],"(",colOrder[i],colOrder[k],")",distance_mat[i,k],"(",colOrder[j],colOrder[k],")",distance_mat[k,j]))
          
        }
      }
    }
  }
}


ggplot_diff_rank_changes=function(rank_for_plot,flag_top100=F,fileName)
{
  if(flag_top100)
  {
    rank_for_plot=rank_for_plot[1:100,]
    
  }
  line_length=2
  gap=0.2
  plot_ranks=which(rank_for_plot>0,arr.ind = T)
  colnames(plot_ranks)=c("rank","iteration")
  as.data.frame(plot_ranks)
  
  x1=(plot_ranks[,2]-1)*line_length+(plot_ranks[,2])*gap
  y1=plot_ranks[,1]
  x2=x1+line_length
  y2=y1
  group=as.character(rank_for_plot[plot_ranks])
  
  xbreaks=seq(1:dim(rank_for_plot)[2])
  xticks=(xbreaks*2*gap+(xbreaks-1)*2*line_length+line_length)/2
  
  df=data.frame(x1,y1,x2,y2,group)
  p1=ggplot()+ 
    geom_line(aes(x=c(1,max(x2)),y=c(1,1)),linetype = "dashed")+
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = group), data = df)+
    scale_color_manual(values=(c("firebrick3","deepskyblue3")))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())+ylab("rank")+
    scale_x_continuous(breaks= xticks ,label = seq(1:dim(rank_for_plot)[2]))+
    xlab("step")+theme(axis.ticks.x=element_blank())
  
  pdf(fileName,width=10,height=8)
  print(p1)
  dev.off()
}




hclust_cal_distance_by_merge=function(distanceMatrices, merge,selectedi=16)
{
  rank_all_x_all=NA
  distanceMatrices=as.matrix(distanceMatrices)
  A=which(merge<0,arr.ind = T)
  cluster_elements=matrix(0,dim(merge)[1],dim(distanceMatrices)[1])
  alone_elements=matrix(1,dim(merge)[1],dim(distanceMatrices)[1])
  clusters_current_iteration=c()
  alone_current_iteration=c(1:dim(distanceMatrices)[1])
  
  
  return_cluster_elements=NA
  return_clusters_current_iteration=NA
  return_alone_current_iteration=NA
  
  for(i in c(1:(dim(merge)[1])))
  {
    #print(i)
    cluster_n=length(clusters_current_iteration)
    alone_n=sum(alone_elements[i,])
    
    cluster_x_alone=matrix(0,cluster_n,alone_n)
    alone_x_alone=matrix(0,alone_n,alone_n)
    cluster_x_cluster=matrix(0,cluster_n,cluster_n)
    all_x_all=matrix(0,(cluster_n+alone_n),(cluster_n+alone_n))
    if(cluster_n>0)
    {
      for(j in c(1:cluster_n))
      {
        if(alone_n>0)
        {
          for(k in c(1:alone_n))
          {
            clusterj=clusters_current_iteration[j]
            alonek=alone_current_iteration[k]
            cluster_elements_j=which(cluster_elements[clusterj,]>0,arr.ind = T)
            dist_jk=mean(distanceMatrices[alonek,cluster_elements_j])
            cluster_x_alone[j,k]=dist_jk
            
          }
        }
        
      }
    }
    if(cluster_n>0)
    {
      if(alone_n>0)
      {
        rownames(cluster_x_alone)=clusters_current_iteration
        colnames(cluster_x_alone)=alone_current_iteration
      }else{
        cluster_x_alone=NA
      }
    }else{
      cluster_x_alone=NA
    }
    
    
    if(cluster_n>0)
    {
      for(j in c(1:cluster_n))
      {
        for(k in c(1:cluster_n))
        {
          clusterj=clusters_current_iteration[j]
          clusterk=clusters_current_iteration[k]
          cluster_elements_j=which(cluster_elements[clusterj,]>0,arr.ind = T)
          cluster_elements_k=which(cluster_elements[clusterk,]>0,arr.ind = T)
          
          dist_jk=mean(distanceMatrices[cluster_elements_k,cluster_elements_j])
          cluster_x_cluster[j,k]=dist_jk
        }
      }
      rownames(cluster_x_cluster)=clusters_current_iteration
      colnames(cluster_x_cluster)=clusters_current_iteration
    }else{
      cluster_x_cluster=NA
    }
    
    
    if(alone_n>0)
    {
      alone_x_alone=distanceMatrices[alone_current_iteration,alone_current_iteration]
    }else{
      alone_x_alone=NA
      rownames(alone_x_alone)=alone_current_iteration
      colnames(alone_x_alone)=alone_current_iteration
    }
    
    all_x_all[1:cluster_n,1:cluster_n]=cluster_x_cluster
    all_x_all[1:cluster_n,(cluster_n+1):(cluster_n+alone_n)]=cluster_x_alone
    all_x_all[(cluster_n+1):(cluster_n+alone_n),1:cluster_n]=t(cluster_x_alone)
    all_x_all[(cluster_n+1):(cluster_n+alone_n),(cluster_n+1):(cluster_n+alone_n)]=alone_x_alone
    #paste colnames and rownames
    
    ###############end test min#################
    
    if(i==selectedi)
    {
      #######################################################
      #cluster 5,8
      #31 51
      #cluster5=which(clusters_current_iteration==5,arr.ind = T)
      #cluster8=which(clusters_current_iteration==8,arr.ind = T)
      
      #######################################################
      write.csv(cluster_x_cluster,file = "cluster_x_cluster.csv",quote=F)
      write.csv(alone_x_alone,file = "alone_x_alone.csv",quote=F)
      write.csv(cluster_x_alone,file = "cluster_x_alone.csv",quote=F)
      write.csv(all_x_all,file="all_x_all.csv",quote=F)
      
      #rank_all_x_all=t(apply(all_x_all,1,order))
      rank_all_x_all=rank_of_matrix(all_x_all)
      return_cluster_elements=cluster_elements
      return_clusters_current_iteration=clusters_current_iteration
      return_alone_current_iteration=alone_current_iteration
      break
    }
    
    for(j in c(1:(dim(merge)[2])))
    {
      one=merge[i,j]
      if(one>0)
      {

        cluster_elements[i,]=cluster_elements[i,]+cluster_elements[one,]
        
        remove=c(one)
        clusters_current_iteration=clusters_current_iteration[!clusters_current_iteration %in% remove]
        
      }else{
        cluster_elements[i,(-one)]=1
        alone_elements[(i:dim(merge)[1]),(-one)]=0
        
        remove=c(-one)
        alone_current_iteration=alone_current_iteration[!alone_current_iteration %in% remove]
        
      }
      
    }
    clusters_current_iteration=c(clusters_current_iteration,i)
    
  }
  list(return_cluster_elements=return_cluster_elements,
       return_clusters_current_iteration=return_clusters_current_iteration,
       return_alone_current_iteration=return_alone_current_iteration,
       rank_mat = rank_all_x_all,
       distance_mat = all_x_all,
       colOrder=c(paste('c',return_clusters_current_iteration,sep="_"),paste('a',return_alone_current_iteration,sep="_")))
  
  
}


rank_of_matrix=function(mat)#return, 1 means minimum
{
  
  tmp=as.numeric(mat[upper.tri(mat)])
  sorted=sort(tmp,index.return=T)
  sorted$ix
  ranks=sort(sorted$ix,index.return=T)
  ranks$ix
  
  rank_mat=matrix(0,dim(mat)[1],dim(mat)[2])
  ut <- upper.tri(rank_mat)
  lt <- lower.tri(rank_mat)
  rank_mat[ut] <- ranks$ix
  
  rank_mat[lt] <- t(rank_mat)[lt]
  rank_mat
}

display_rank_diff=function(all_x_all1, all_x_all2,iteration,sampleID,show_loops=F)
{
  #all_x_all1=all_x_all_ab
  #all_x_all2=all_x_all_root
  
  rank_mat1=all_x_all1$rank_mat
  rank_mat2=all_x_all2$rank_mat
  
  diff_rank=rank_mat1-rank_mat2
  diff_rank[lower.tri(diff_rank)]=0
  A=which(diff_rank!=0,arr.ind = T)
  
  if(length(A)>0)
  {
    print(paste("num",length(A),"min",min(rank_mat1[A]),"max",max(rank_mat2[A])))
  }else{
    print("num",length(A))
  }
  
  
  unique_elements=sort(unique(as.numeric(A)))
  names=all_x_all1$colOrder[unique_elements]
  
  all_names=all_x_all1$colOrder
  tmp_source=all_names[A[,1]]
  tmp_target=all_names[A[,2]]
  tmp_rank1=rank_mat1[A]
  tmp_rank2=rank_mat2[A]
  pair_and_rank=data.frame(tmp_source,tmp_target,tmp_rank1,tmp_rank2)
  colors=check_loops(pair_and_rank)
  display_loops(pair_and_rank,colors,show_loops)
  
  
  
  clusters=c()
  alones=c()
  for(i in c(1:length(names)))
  {
    
    tmp=unlist(strsplit(names[i],'_'))
    
    if(tmp[1]=='c')
    {
      clusters=c(clusters,as.numeric(tmp[2]))
    }else{
      alones=c(alones,as.numeric(tmp[2]))
    }
  }
  
  node_names1=get_node_name(all_x_all1$colOrder,sampleID)
  node_names2=get_node_name(all_x_all2$colOrder,sampleID)
  plot_circle_network(paste(iteration,"_1_network.pdf",sep=""),c(1:length(all_x_all1$colOrder)),node_names1,A,colors,pair_and_rank)
  plot_circle_network(paste(iteration,"_2_network.pdf",sep=""),c(1:length(all_x_all2$colOrder)),node_names2,A,colors,pair_and_rank)
  
  pair_and_rank
}

get_node_name=function(colOrder,sampleID)
{
  node_names=colOrder
  for(i in c(1:length(colOrder)))
  {
    tmp=unlist(strsplit(colOrder[i],'_'))
    if(tmp[1]=='a')
    {
      index=as.numeric(tmp[2])
      node_names[i]=sampleID[index]
    }
  }
  node_names
}

display_loops=function(pair_and_rank, cluster_record,print_each_loop=T)
{
  unique_colors=unique(cluster_record)
  all_flag=T
  for(j in c(1:length(unique_colors)))
  {
    this_color=unique_colors[j]
    A=which(cluster_record==this_color,arr.ind = T)
    if(print_each_loop)
    {
      print(pair_and_rank[A,])
    }
    
    
    #check all contain c~
    flag=F
    selected_pairs=pair_and_rank[A,]
    names=as.matrix(selected_pairs[,1])
    
    
    for(k in c(1:length(names)))
    {
      
      tmp=unlist(strsplit(names[k],'_'))
      
      if(tmp[1]=='c')
      {
        flag=T
        break
      }
    }
    
    #print(paste("flag",flag))
    if(!flag)
    {
      all_flag=F
      print(pair_and_rank[A,])
    }
    
    
  }
  print(all_flag)
}

check_loops=function(pair_and_rank)
{
  #pair_and_rank=tmp
  pair_and_rank=as.matrix(pair_and_rank)
  unique_names=unique(c(pair_and_rank[,1],pair_and_rank[,2]))
  
  ranks1=pair_and_rank[,3]
  ranks2=pair_and_rank[,4]
  
  pair_used=rep(0,length(ranks1))
  
  cluster_record=rep(0,length(ranks1))
  
  
  for(j in c(1:length(ranks1)))
  {
    
    #if should be assign to existed class
    
    number1=ranks1[j]
    number2=ranks2[j]
    A=which(ranks2==number1,arr.ind = T)
    index1=cluster_record[A]
    B=which(ranks1==number2,arr.ind = T)
    index2=cluster_record[B]
    
    if(index1+index2==0)
    {
      #a new one
      this_index=max(cluster_record)+1
      cluster_record[j]=this_index
    }else{
      if(A==B)
      {
        #merge into cluster[A]
        this_index=cluster_record[A]
        cluster_record[j]=this_index
        
      }else{
        #merge cluster[A], cluster[B], and this new one
        this_index=max(cluster_record)+1
        cluster_record[j]=this_index
        
        if(index1>0)
        {
          cluster_A=which(cluster_record==index1,arr.ind = T)
          cluster_record[cluster_A]=this_index
        }
        if(index2>0)
        {
          cluster_B=which(cluster_record==index2,arr.ind = T)
          cluster_record[cluster_B]=this_index
        }
        
      }
    }
    
  }
  cluster_record
}

plot_circle_network=function(fileName,nodes_index,nodes_label,edges,loop_colors,pair_and_rank,r=10,plot_compare_rank=F)
{
  #nodes_index=c(1:length(all_x_all1$colOrder))
  #nodes_label=all_x_all1$colOrder
  #edges=A
  
  n=length(nodes_index)
  
  each_angle=360/n
  
  nodes_angle=c(0)
  for(i in c(2:n))
  {
    nodes_angle=c(nodes_angle,(max(nodes_angle)+each_angle))
  }
  
  x=rep(0,n)
  y=rep(0,n)
  x=r*cospi(nodes_angle/180)
  y=r*sinpi(nodes_angle/180)
  
  x1=x[edges[,1]]
  y1=y[edges[,1]]
  x2=x[edges[,2]]
  y2=y[edges[,2]]
  
  loop_colors=as.character(loop_colors)
  if(!plot_compare_rank)
  {
    #loop_colors=rep(2,length(loop_colors))
  }else{
    loop_colors=as.character(loop_colors)
  }
  
  edge_position=data.frame(x1,y1,x2,y2,loop_colors)
  df=data.frame(nodes_index,nodes_label,nodes_angle,x,y)
  
  rank1=as.matrix(pair_and_rank[,3])
  rank2=as.matrix(pair_and_rank[,4])
  compare_rank_x=(x1+x2)/2
  compare_rank_y=(y1+y2)/2-0.2
  compare_rank_text=rep(NA,length(compare_rank_x))
  for(i in c(1:length(compare_rank_x)))
  {
    compare_rank_text[i]=paste(rank1[i],rank2[i],sep=",")
  }
  
  compare_rank=data.frame(compare_rank_x,compare_rank_y,compare_rank_text)
  
  
  pdf(fileName,width=12,height = 12)
  p1=ggplot() +
    coord_fixed()+geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour = loop_colors), data = edge_position)+
    geom_text(data=df, aes(x=x, y=y, group=1),label=nodes_label,angle=nodes_angle,check_overlap = TRUE)+
    theme_void()+theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())
  if(plot_compare_rank)
  {
    p1=p1+geom_text(data=compare_rank,aes(x=compare_rank_x,y=compare_rank_y),label=compare_rank_text,check_overlap = TRUE)
    
  }
  print(p1)
  dev.off()
}


plot_diff_dendrogram=function(merge1,merge2,sampleID,fileHeader,distance_one,distance_two)
{
  record1=generate_one_tree(merge1,sampleID)
  record2=generate_one_tree(merge2,sampleID)
  plot_one_tree(record1,record2,paste(fileHeader,distance_one,"dendrogram_compare.pdf",sep="_"))
  plot_one_tree(record2,record1,paste(fileHeader,distance_two,"dendrogram_compare.pdf",sep="_"))
  
}

plot_one_tree=function(record1,record_control,fileName)
{
  
  line_record=record1$line_record
  leaf_record=record1$leaf_record
  
  link_with_left=record1$record[,1]
  link_with_right=record1$record[,2]
  A=which(!is.na(link_with_left),arr.ind = T)
  B=which(!is.na(link_with_right),arr.ind = T)
  link_with_left=link_with_left[A]
  link_with_right=link_with_right[B]
  
  link_with_left_control=record_control$record[,1]
  link_with_right_control=record_control$record[,2]
  link_with_left_control=link_with_left_control[A]
  link_with_right_control=link_with_right_control[B]
  
  for(i in c(1:length(A)))
  {
    if(link_with_left[i]!=link_with_left_control[i])
    {
      if(link_with_left[i]!=link_with_right_control[i])
      {
        line_record[i,5]=1
        line_record[(length(A)+i),5]=1
      }
    }
  }
  for(i in c(1:length(B)))
  {
    if(link_with_right[i]!=link_with_left_control[i])
    {
      if(link_with_right[i]!=link_with_right_control[i])
      {
        line_record[(length(A)+i),5]=1
        line_record[i,5]=1
      }
    }
  }
  group=as.character(line_record[,5])
  line_record=as.data.frame(line_record,group)
  pdf(fileName,width = 12,height=10)
  #ggplot()+geom_text(data=record1$record,aes(x=(x),y=(y+0.2),label=V7,angle = 60))+geom_segment(aes(x = fromx, y = fromy, xend = tox, yend = toy, colour = group), data = line_record)+scale_color_manual(values=(c( "deepskyblue1","brown1")))
  p1=ggplot()+
    geom_segment(aes(x = fromx, y = fromy, xend = tox, yend = toy, colour = group), data = line_record)+
    geom_text(data=leaf_record,aes(x=(x),y=(y-1),label=V7,angle = 60))+
    scale_color_manual(values=(c( "deepskyblue1","brown1")))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())+ylab("step")+
    theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.text=element_blank())
  print(p1)
  dev.off()
  
}

generate_one_tree=function(merge1,sampleID)
{
  #initialize left,right, partent
  #for each individual and cluster from 1 to 76(minus and postive in merge)
  unique_element=sort(unique(as.numeric(as.matrix(merge1))))
  unique_element=c(unique_element,(max(unique_element)+1))
  record=matrix(NA,length(unique_element),6)
  rownames(record)=unique_element
  colnames(record)=c("left","right","parent","clust_num","x","y")
  for(i in c(1:length(unique_element)))
  {
    element_i=unique_element[i]
    if(element_i>0)
    {
      record[i,1]=merge1[element_i,1]
      record[i,2]=merge1[element_i,2]
    }
    A=which(merge1==element_i,arr.ind = T)
    record[i,3]=A[1]
    
  }
  #find parent most, whose parent is NA, then
  parent_most_index=which(is.na(record[,3]),arr.ind = T)
  parent_most_element=unique_element[parent_most_index]
  
  #
  #record cluster_num(number of nodes beneath this, not include itself)
  #！
  cluster_size=rep(0,length(unique_element))
  sub_size=get_subTree_size(parent_most_element,unique_element=unique_element,record,cluster_size)
  record[,4]=sub_size$cluster_size
  
  #from parent most,given range min_this, max_this, (1:number of elements)
  #call visit node function, given range and node index
  
  min_parent_most=1
  max_parent_most=length(unique_element)
  x_array=rep(0,length(unique_element))
  y_array=rep(0,length(unique_element))
  current_y=0
  x_y_array=visit_node(min_parent_most,max_parent_most,parent_most_element,unique_element,record,x_array,y_array,current_y)
  record[,5]=x_y_array$x_array
  record[,6]=x_y_array$y_array
  
  for(i in c(1:length(unique_element)))
  {
    if(unique_element[i]<0)
    {
      record[i,6]=0
    }else{
      record[i,6]=unique_element[i]
    }
  }
  
  colnames(record)[5]="x"
  record=cbind(record,rownames(record))
  record=as.matrix(record)
  class(record)="numeric"
  record=as.data.frame(record)
  
  #########generate link_record, leaf_record, and plot
  link_with_left=record[,1]
  link_with_right=record[,2]
  A=which(!is.na(link_with_left),arr.ind = T)
  B=which(!is.na(link_with_right),arr.ind = T)
  line_record=matrix(0,(length(A)+length(B)),5)
  colnames(line_record)=c("fromx","fromy","tox","toy","different")
  line_record[1:length(A),1]=record[A,5]
  line_record[1:length(A),2]=record[A,6]
  line_record[(length(A)+1):(length(A)+length(B)),1]=record[B,5]
  line_record[(length(A)+1):(length(A)+length(B)),2]=record[B,6]
  link_with_left=link_with_left[A]
  link_with_right=link_with_right[B]
  
  
  for(i in c(1:length(A)))
  {
    index_i=which(unique_element==link_with_left[i],arr.ind = T)
    line_record[i,3]=record[index_i,5]
    line_record[i,4]=record[index_i,6]
  }
  for(i in c(1:length(B)))
  {
    index_i=which(unique_element==link_with_right[i],arr.ind = T)
    line_record[(length(A)+i),3]=record[index_i,5]
    line_record[(length(A)+i),4]=record[index_i,6]
  }
  
  line_record=as.data.frame(line_record)
  C=which(unique_element<0,arr.ind = T)
  leaf_record=record
  leaf_record[C,7]=sampleID[-leaf_record[C,7]]
  leaf_record=as.data.frame(leaf_record)
  #ggplot()+geom_text(data=record,aes(x=(x),y=(y+0.2),label=V7,angle = 60))+geom_segment(aes(x = fromx, y = fromy, xend = tox, yend = toy, colour = "segment"), data = line_record)
  #ggplot()+geom_text(data=leaf_record,aes(x=(x),y=(y+0.2),label=V7,angle = 60))+geom_segment(aes(x = fromx, y = fromy, xend = tox, yend = toy, colour = "segment"), data = line_record)
  
  list(record=record,
       line_record=line_record,
       leaf_record=leaf_record)
}

get_subTree_size=function(this_element,unique_element,record,cluster_size) #record global
{
  this_index=which(unique_element==this_element,arr.ind = T)
  this_size=NA
  if(is.na(record[this_index,1]))
  {
    this_size=0
  }else{
    #get left subTree_size
    left_element=record[this_index,1]
    left=get_subTree_size(left_element,unique_element,record,cluster_size)
    left_size=left$this_size
    cluster_size=pmax(cluster_size,left$cluster_size)
    #get right subTree size
    right_element=record[this_index,2]
    
    right=get_subTree_size(right_element,unique_element,record,cluster_size)
    right_size=right$this_size
    
    cluster_size=pmax(cluster_size,right$cluster_size)
    #size =left+right+2
    this_size=left_size+right_size+2
    
    #print(paste("element",this_element,"cluster_size",this_size,"left",left_size,"right",right_size,"right_element",right_element,"this_index",this_index,"left_element",left_element))
    
  }
  
  #record size
  cluster_size[this_index]=this_size
  #return size
  list(this_size=this_size,
       cluster_size=cluster_size)
  
}
plot_heatmap=function(mat)
{
  test=mat
  test=as.data.frame(test)
  test$id=rownames(test)
  mat_melt <- melt(test, id.var="id")
  p1=ggplot(mat_melt,aes(x=as.integer(variable),y=as.integer(id)))+
    geom_tile(aes(fill=mat_melt$value),colour="white")+
    #scale_fill_gradient("ratio of together",low="steelblue",medium="white",high="firebrick")
    scale_fill_gradientn(colours = c("cyan", "white", "red"),
                         values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)))
  print(p1)
  
}
visit_node=function(min_this, max_this, this_element,unique_element,record,x_array,y_array,current_y)#range, node, record global 
{
  #the range= min_this, max_this
  this_index=which(unique_element==this_element,arr.ind = T)
  #the range of left, min=min_this, max=min_this+clust_num, 
  min_left=min_this
  left_element=record[this_index,1]
  if(is.na(left_element))
  {
    x_this=min_this
    x_array[this_index]=x_this
    y_array[this_index]=current_y
  }else{
    
    left_index=which(unique_element==left_element,arr.ind = T)
    left_clust_num=record[left_index,4]
    
    max_left=min_this+left_clust_num
    
    #x of this, max(left)+1
    #record this x
    x_this=max_left+1
    x_array[this_index]=x_this
    y_array[this_index]=current_y
    
    current_y=current_y+1
    #visit left, function
    #print(paste("left",min_left,max_left,left_element))
    left_x_y_array=visit_node(min_left,max_left,left_element,unique_element,record,x_array,y_array,current_y)
    
    x_array=pmax(x_array,left_x_y_array$x_array)
    y_array=pmax(y_array,left_x_y_array$y_array)
    
    
    #the range of right, min =max(left)+2, max=min(right)+clust_num of max_this
    #visit right, function
    min_right=max_left+2
    max_right=max_this
    right_element=record[this_index,2]
    #print(paste("right",min_right,max_right,right_element))
    right_x_y_array=visit_node(min_right,max_right,right_element,unique_element,record,x_array,y_array,current_y)
    x_array=pmax(x_array,right_x_y_array$x_array)
    y_array=pmax(y_array,right_x_y_array$y_array)
  }
  
  
  list(x_array=x_array,
       y_array=y_array)
  
}
###############################################################################
