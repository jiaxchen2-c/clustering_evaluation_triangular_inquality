for(alpha in c(1:10))
{
  
i=1

while(i<=100000)
{
  

x=sample(1:1000,10)
y=sample(1:1000,10)
z=sample(1:1000,10)

#dxy=1-abs(cor(x,y))
#dxz=1-abs(cor(x,z))
#dyz=1-abs(cor(y,z))

corxy=cor(x,y,method='pearson')
corxz=cor(x,z,method='pearson')
coryz=cor(y,z,method='pearson')



dxy=1-(abs(corxy))^alpha
dxz=1-(abs(corxz))^alpha
dyz=1-(abs(coryz))^alpha

sqrt_dxy=sqrt(dxy)
sqrt_dxz=sqrt(dxz)
sqrt_dyz=sqrt(dyz)

#if(dxy+dyz<dxz)
if(sqrt_dxy+sqrt_dxz<sqrt_dyz)
{
  print(sqrt_dxy+sqrt_dxz-sqrt_dyz)
  print(x)
  print(y)
  print(z)
  break
}

i=i+1

}

x
y
z
#x=c(79,8)
#y=c(39,41)
#z=c(85,16)
cor(x,y,method='pearson')
cor(x,z,method='pearson')
cor(y,z,method='pearson')

dxy
dxz
dyz
sqrt_dxy
sqrt_dxz
sqrt_dyz
############plot#########
x
y
z
#mat=rbind(x,y,z)
#mat=t(mat)
x0=c(1,2,3)
#plot(x0,x)
#matplot(mat)
#matlines(mat)
i
}