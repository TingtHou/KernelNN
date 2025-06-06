
library(MASS)
library(ggplot2)
n<-1000
p<-500
it<-100
pvalue<-matrix(NA, nrow = it, ncol = 3)
for (i in 1:it)
{
  G1<-matrix(rbinom(n*p,1,0.02), nrow = n, ncol = p)
  #G2<-matrix(rbinom(n*p), nrow = n, ncol = p)
  I<-diag(n)
  K1<-G1 %*%t(G1)
  # K2<-K1^2

  y<-mvrnorm(n = 1, mu = rep(0,n), Sigma =2*I)
  #J<-matrix(1,nrow = n,ncol = n)

  KernelList<-list(K1,I)
  vcs.result<-MINQUE0(KList = KernelList,y = y)
  zTest<-MNQTest0_Chi(KList = KernelList,vcs = vcs.result$vcs,y=y,TestingID = c(1))
  # vcs.result<-IMINQUE(KList = KernelList,y = y,prior = c(1,1),epoch=100,threshold = 1e-3)
  # if ("err" %in% names(vcs.result))
  # {
  #   next
  # }
  #  zTest<-IMNQTest_Normal(KList = KernelList, est_vcs = vcs.result$vcs,index_interest  = c(1))
  # if ("err" %in% names(zTest))
  # {
  #   next;
  # }
  pvalue[i,1:(length(KernelList)-1)]<-zTest$components[1:(length(KernelList)-1),2]
  pvalue[i,2]<-zTest$overall
  pvalue[i,3]<-1
}

pvalue<-pvalue[which(pvalue[,3]<100),]
xseq<-seq(0,1,by = 0.001)
Power<-data.frame(nrow = length(xseq),ncol = 2)
colnames(Power)<-c("par1","overall")
for(i in seq_along(xseq))
{
  Power[i,]<-colMeans(pvalue[,-3]<xseq[i])
}


Data<-Power

Data$xseq<-xseq
pvalue_MNQ1<-ggplot(Data)+aes(x=xseq, y=par1)+geom_point()+
  geom_abline(aes(slope=1,intercept=0))+xlim(0,1)+ylim(0,1)+
  ylab("Type I error")+xlab("Significance level")

pvalue_Chi_1sided<-ggplot(Data)+aes(x=xseq, y=overall)+geom_point()+
  geom_abline(aes(slope=1,intercept=0))+xlim(0,0.5)+ylim(0,0.5)+
  ylab("Type I error" )+xlab("Significance level")
