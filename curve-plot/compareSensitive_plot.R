library(ggplot2)

x=c(-1,0.1,1)

#r
x0=-1
gap=0.0001
x=c()
y_abs=c()
y_root=c()
y_square=c()
y_half=c()

xi=x0
while(T)
{
  yi_r=sqrt(1-abs(xi))
  yi_s=sqrt(1-xi^2)
  yi_a=1-abs(xi)
  yi_h=sqrt((1-xi)/2)
  x=c(x,xi)
  y_root=c(y_root,yi_r)
  y_square=c(y_square,yi_s)
  y_abs=c(y_abs,yi_a)
  y_half=c(y_half,yi_h)
  xi=xi+gap
  if(xi>1)
  {
    break
  }
}

####form df for ggplot######
x_all=c(as.numeric(x),as.numeric(x),as.numeric(x),as.numeric(x))
y_all=c(as.numeric(y_abs),as.numeric(y_root),as.numeric(y_square),as.numeric(y_half))
category_all=c(rep("1abs",length(y_abs)),rep("2root",length(y_root)),rep("3square",length(y_square)),rep("4half",length(y_half)))
  


df=data.frame(x_all,y_all,category_all)


p1=ggplot(data=df,aes(x=c(df$x_all),y=c(df$y_all),color=as.character(c(category_all)),linetype=as.character(c(category_all))))+geom_path(size=2)+
  theme_classic()+theme(text=element_text(size=25))

pdf("test1_compareSensitive.pdf",width = 12,height = 10)
#plot(x,y_root)
print(p1)

dev.off()

#######traditional plot#####
pdf("test_root_compareSensitive.pdf")
plot(x,y_root)
dev.off()

pdf("test1_abs_compareSensitive.pdf")
plot(x,y_abs)
dev.off()

pdf("test1_square_compareSensitive.pdf")
plot(x,y_square)
dev.off()


