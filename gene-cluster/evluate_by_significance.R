evaluate_by_t_test = function(result_file)
{
  data=as.matrix(read.csv(result_file,row.names = 1))
  pvalues = data[4,]
  compare_value = data[1,]
  
  A=which(pvalues<0.05,arr.ind = T)
  compare_value[A]
  lose=length(which(compare_value[A]>0))
  win=length(which(compare_value[A]<0))
  
  
  dr_not_better=data[2,A]
  dr_better=data[3,A]
  
  print(t.test(dr_better,dr_not_better))#over test on count for significant diff pair
  print(compare_value[A])
  print(t.test(data[2,],data[3,]))#over test on count for all pair
  list(win,lose,dr_better,dr_not_better)
  
  
}

wins=c()
loses=c()
dr_better=c()
dr_not_better=c()
methods=c()

tmp=evaluate_by_t_test("timeCourse-result-pearson-hcluster-normalized_6.17/absolute_vs_root.csv")
wins=c(wins,unlist(tmp[1]))
loses=c(loses,unlist(tmp[2]))
dr_better=c(dr_better,unlist(tmp[3]))
dr_not_better=c(dr_not_better,unlist(tmp[4]))
methods=c(methods,rep("pearson_hclust",length(unlist(tmp[3]))))
tmp=evaluate_by_t_test("timeCourse-result-spearman-hcluster-normalized_6.17/absolute_vs_root2.csv")
wins=c(wins,unlist(tmp[1]))
loses=c(loses,unlist(tmp[2]))
dr_better=c(dr_better,unlist(tmp[3]))
dr_not_better=c(dr_not_better,unlist(tmp[4]))
methods=c(methods,rep("spearman_hclust",length(unlist(tmp[3]))))
tmp=evaluate_by_t_test("timeCourse-result-uncentered-hcluster-normalized_6.17/absolute_vs_root.csv")
wins=c(wins,unlist(tmp[1]))
loses=c(loses,unlist(tmp[2]))
dr_better=c(dr_better,unlist(tmp[3]))
dr_not_better=c(dr_not_better,unlist(tmp[4]))
methods=c(methods,rep("uncentered_hclust",length(unlist(tmp[3]))))

tmp=evaluate_by_t_test("timeCourse-result-pearson-pam-normalized_6.17/absolute_vs_root.csv")
wins=c(wins,unlist(tmp[1]))
loses=c(loses,unlist(tmp[2]))
dr_better=c(dr_better,unlist(tmp[3]))
dr_not_better=c(dr_not_better,unlist(tmp[4]))
methods=c(methods,rep("pearson_pam",length(unlist(tmp[3]))))
tmp=evaluate_by_t_test("timeCourse-result-spearman-pam-normalized_6.17/absolute_vs_root.csv")
wins=c(wins,unlist(tmp[1]))
loses=c(loses,unlist(tmp[2]))
dr_better=c(dr_better,unlist(tmp[3]))
dr_not_better=c(dr_not_better,unlist(tmp[4]))
methods=c(methods,rep("spearman_pam",length(unlist(tmp[3]))))


wilcox.test(wins,loses)
wilcox.test(dr_better,dr_not_better)

tmp=evaluate_by_t_test("timeCourse-result-uncentered-pam-normalized_6.17/absolute_vs_root.csv")
wins=c(wins,unlist(tmp[1]))
loses=c(loses,unlist(tmp[2]))
dr_better=c(dr_better,unlist(tmp[3]))
dr_not_better=c(dr_not_better,unlist(tmp[4]))
methods=c(methods,rep("uncentered_pam",length(unlist(tmp[3]))))

t.test(wins,loses)#over test on win/lose event(with significant diff) for all combination~
wilcox.test(wins,loses)
wilcox.test(dr_better,dr_not_better)


x=c(c(1:length(dr_better)),c(1:length(dr_not_better)))
y=c(dr_better,dr_not_better)
cate=c(rep("dr_better",length(dr_better)),rep("da_better",length(dr_not_better)))
df=data.frame(x,y,cate,methods)
pdf("GO_evaluate_only_significant_dr_da.pdf",width=12,height = 8)
p1=ggplot(data=df, aes(x=x, y=y, fill=cate)) +
  geom_bar(stat="identity", position=position_dodge())+theme_light()
print(p1)
dev.off()
pdf("GO_evaluate_only_significant_dr_da2.pdf",width=12,height = 8)
p1=ggplot(data=df, aes(x=x, y=y, color=cate,fill=methods)) +
  geom_bar(stat="identity", position=position_dodge())+theme_light()
print(p1)
dev.off()

###################start compare between dr and ds########################
############################################################################
wins=c()
loses=c()
dr_better=c()
dr_not_better=c()
methods=c()

tmp=evaluate_by_t_test("timeCourse-result-compareSquare-hcluster-normalized_6.17/square_vs_root.csv")
wins=c(wins,unlist(tmp[1]))
loses=c(loses,unlist(tmp[2]))
dr_better=c(dr_better,unlist(tmp[3]))
dr_not_better=c(dr_not_better,unlist(tmp[4]))
methods=c(methods,rep("pearson_hclust",length(unlist(tmp[3]))))
tmp=evaluate_by_t_test("timeCourse-result-compareSquare-spearman-hcluster-normalized-6.17/square_vs_root.csv")
wins=c(wins,unlist(tmp[1]))
loses=c(loses,unlist(tmp[2]))
dr_better=c(dr_better,unlist(tmp[3]))
dr_not_better=c(dr_not_better,unlist(tmp[4]))
methods=c(methods,rep("spearman_hclust",length(unlist(tmp[3]))))
tmp=evaluate_by_t_test("timeCourse-result-compareSquare-uncentered-hcluster-normalized_6.17/square_vs_root.csv")
wins=c(wins,unlist(tmp[1]))
loses=c(loses,unlist(tmp[2]))
dr_better=c(dr_better,unlist(tmp[3]))
dr_not_better=c(dr_not_better,unlist(tmp[4]))
methods=c(methods,rep("uncentered_hclust",length(unlist(tmp[3]))))
tmp=evaluate_by_t_test("timeCourse-result-compareSquare-pam-normalized_6.17/square_vs_root.csv")
wins=c(wins,unlist(tmp[1]))
loses=c(loses,unlist(tmp[2]))
dr_better=c(dr_better,unlist(tmp[3]))
dr_not_better=c(dr_not_better,unlist(tmp[4]))
methods=c(methods,rep("pearson_pam",length(unlist(tmp[3]))))
tmp=evaluate_by_t_test("timeCourse-result-compareSquare-spearman-pam-normalized_6.17/square_vs_root.csv")
wins=c(wins,unlist(tmp[1]))
loses=c(loses,unlist(tmp[2]))
dr_better=c(dr_better,unlist(tmp[3]))
dr_not_better=c(dr_not_better,unlist(tmp[4]))
methods=c(methods,rep("spearman_pam",length(unlist(tmp[3]))))
tmp=evaluate_by_t_test("timeCourse-result-compareSquare-uncentered-pam-normalized_6.17/square_vs_root.csv")
wins=c(wins,unlist(tmp[1]))
loses=c(loses,unlist(tmp[2]))
dr_better=c(dr_better,unlist(tmp[3]))
dr_not_better=c(dr_not_better,unlist(tmp[4]))
methods=c(methods,rep("uncentered_pam",length(unlist(tmp[3]))))




wilcox.test(wins,loses)
wilcox.test(dr_better,dr_not_better)

x=c(c(1:length(dr_better)),c(1:length(dr_not_better)))
y=c(dr_better,dr_not_better)
cate=c(rep("2dr_better",length(dr_better)),rep("1ds_better",length(dr_not_better)))
df=data.frame(x,y,cate,methods)
pdf("GO_evaluate_only_significant_dr_ds.pdf",width=12,height = 8)
p1=ggplot(data=df, aes(x=x, y=y, fill=cate)) +
  geom_bar(stat="identity", position=position_dodge())+theme_light()
print(p1)
dev.off()
pdf("GO_evaluate_only_significant_dr_ds2.pdf",width=12,height = 8)
p1=ggplot(data=df, aes(x=x, y=y, color=cate,fill=methods)) +
  geom_bar(stat="identity", position=position_dodge())+theme_light()
print(p1)
dev.off()
