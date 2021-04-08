data_G1=as.matrix(read.csv("data_mesc/G1_singlecells_magic.csv",header = T))
mat_G1=as.matrix(read.table("data_mesc/G1_singlecells_counts.txt",header = T,row.names = 1))
mat_G1=mat_G1[,4:dim(mat_G1)[2]]


data_G2M=as.matrix(read.csv("data_mesc/G2M_singlecells_magic.csv",header = T))
mat_G2M=as.matrix(read.table("data_mesc/G2M_singlecells_counts.txt",header = T,row.names = 1))
mat_G2M=mat_G2M[,4:dim(mat_G2M)[2]]


data_S=as.matrix(read.csv("data_mesc/S_singlecells_magic.csv",header = T))
mat_S=as.matrix(read.table("data_mesc/S_singlecells_counts.txt",header = T,row.names = 1))
mat_S=mat_S[,4:dim(mat_S)[2]]

data_G1=t(data_G1)
#rownames(data_G1)=rownames(mat_G1)
colnames(data_G1)=colnames(mat_G1)

data_G2M=t(data_G2M)
#rownames(data_G2M)=rownames(mat_G2M)
colnames(data_G2M)=colnames(mat_G2M)

data_S=t(data_S)
rownames(data_S)=rownames(mat_S)
colnames(data_S)=colnames(mat_S)


expr=cbind(data_G1,data_G2M,data_S)
cate=c(rep('G1',dim(data_G1)[2]),
       rep('G2M',dim(data_G2M)[2]),
       rep('S',dim(data_S)[2]))


result=rbind(cate,expr)
colnames(result)=c(colnames(data_G1),colnames(data_G2M),colnames(data_S))
rownames(result)=c('',rownames(data_G1))

sampleID=colnames(result)

result2=rbind(sampleID,result)
rownames(result2)=c('',rownames(result))

write.table(result2,file="data_mesc/concated_magic_count_mesc.txt",quote=F, row.names = T, col.names = F,sep="\t")

A=which(is.nan(result2),arr.ind=T)
print(A)

##################collect for liver data###################
data_magic=read.csv("data_liver/gse124395_normalhumanlivercellatlasdata_magic.csv",header=F)
data=read.table("data_liver/gse124395_normalhumanlivercellatlasdata.txt",header = T, row.names = 1)
data_magic=t(data_magic)
rownames(data_magic)=rownames(data)
colnames(data_magic)=colnames(data)
write.table(data_magic,file="data_liver/collect_magic.txt",quote=F,row.names = T, col.names = T)

