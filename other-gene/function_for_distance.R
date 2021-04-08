require("lsa")

produceDistanceMatrices=function(expr,method)
{

  texpr=t(expr)
  distanceMatrices=NULL
  if(method=="absolute")#absolute Pearson correlation distance
  {
    corrMatrix=cor(texpr)
    distanceMatrices=1-abs(corrMatrix)
  }
  if(method=="root absolute")#square root absolute correlation distance
  {
    corrMatrix=cor(texpr)
    distanceMatrices=1-abs(corrMatrix)
    distanceMatrices=sqrt(distanceMatrices)
  }
  if(method=="root square")#square root absolute correlation distance
  {
    corrMatrix=cor(texpr)
    distanceMatrices=1-corrMatrix^2
    distanceMatrices=sqrt(distanceMatrices)
  }
  
  if(method=="spearman")#absolute Pearson correlation distance
  {
    corrMatrix=cor(texpr,method="spearman")
    distanceMatrices=1-abs(corrMatrix)
  }
  if(method=="root spearman")#square root absolute correlation distance
  {
    corrMatrix=cor(texpr,method="spearman")
    distanceMatrices=1-abs(corrMatrix)
    distanceMatrices=sqrt(distanceMatrices)
  }
  if(method=="square spearman")#square root absolute correlation distance
  {
    corrMatrix=cor(texpr,method="spearman")
    distanceMatrices=1-corrMatrix^2
    distanceMatrices=sqrt(distanceMatrices)
  }
  
  if(method=="kendall")#absolute Pearson correlation distance
  {
    corrMatrix=cor(texpr,method="kendall")
    distanceMatrices=1-abs(corrMatrix)
  }
  
  if(method=="uncentered")#absolute Pearson correlation distance
  {
    corrMatrix=cosine(texpr)
    distanceMatrices=1-abs(corrMatrix)
  }
  if(method=="root uncentered")#square root absolute correlation distance
  {
    corrMatrix=cosine(texpr)
    distanceMatrices=1-abs(corrMatrix)
    distanceMatrices=sqrt(distanceMatrices)
  }
  if(method=="square uncentered")#square root absolute correlation distance
  {
    corrMatrix=cosine(texpr)
    distanceMatrices=1-corrMatrix^2
    distanceMatrices=sqrt(distanceMatrices)
  }
  
  #return distance matrices
  distanceMatrices
}
