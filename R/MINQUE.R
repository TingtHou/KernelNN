MINQUE0<-function(KList,y)
{
  VCs<-length(KList)
  C.mat<-matrix(NA,nrow=VCs,ncol=VCs)
  RightY<-c()
  for (i in 1:VCs)
  {
    for (j in i:VCs)
    {
      C.mat[j,i]<-C.mat[i,j]<-sum(t( KList[[i]])* KList[[j]])
    }
    right<- t(y)%*%KList[[i]]%*%y
    RightY<-c( RightY,right)
  }
  C.mat<-solve(C.mat)
  est.vcs<-C.mat %*% RightY
  return(list(vcs=est.vcs))
}