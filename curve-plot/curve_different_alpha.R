library(ggplot2)

x=c(-1,0.1,1)

#r
x0=-1
gap=0.0001

x_all=c()
y_all=c()
alpha_all=c()
for(alpha in c(1:9))
{
  y_alpha=c()
  x=c()
  xi=x0
  while(T)
  {
    yi_alpha_i=sqrt(1-(abs(xi))^alpha)
    y_alpha=c(y_alpha,yi_alpha_i)
    x=c(x,xi)
    xi=xi+gap
    if(xi>1)
    {
      break
    }
  }
  x_all=c(x_all, x)
  y_all=c(y_all, y_alpha)
  alpha_all=c(alpha_all,rep(alpha,length(x)))
}


####form df for ggplot######



df=data.frame(x_all,y_all,alpha_all)


p1=ggplot(data=df,aes(x=c(df$x_all),y=c(df$y_all),color=as.character(c(alpha_all)),linetype=as.character(c(alpha_all))))+geom_path(size=1)+
  theme_classic()+theme(text=element_text(size=25))

pdf("curve_different_alpha.pdf",width = 12,height = 10)
#plot(x,y_root)
print(p1)

dev.off()



