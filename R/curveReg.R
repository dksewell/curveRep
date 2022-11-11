curveReg <-
function(Y,X,type=c('fourier','splines'),
                    nComponents=min(dim(dat),ncol(X)),
                    COLSX=NULL,
                    alphaX=c(0.1,0.9),
                    linesX=1){
  if(is.data.frame(X)) X = model.matrix(~X)[,-1]
  X = scale(X); Y = scale(Y)
  
  fun1 = function(x,vec){
    if(NROW(vec)%%2==0){
      a0 = vec[1]*2; vec = c(vec[-1],0)*sqrt(2)
      N = NROW(vec);
      aa = vec[seq(2,N,2)-1]; bb = vec[seq(2,N,2)]
    }else{
      a0 = vec[1]*2; vec = vec[-1]*sqrt(2)
      N = NROW(vec); 
      aa = vec[seq(2,N,2)-1]; bb = vec[seq(2,N,2)]
    }
    ret = a0;
    for(i in 1:(N/2)){
      ret = ret + aa[i]*sin(2*pi*i*x) + bb[i]*cos(2*pi*i*x)
    }
    return(ret)
  }
  
  n = nrow(X);
  pp=ncol(X)
  XY = drop(crossprod(X,Y))
  DD = diag(XY/apply(X,2,function(x)crossprod(x)))
  dat = X%*%DD
  
  dat.svd = svds(dat,k=max(nComponents))
  VV = dat.svd$v%*%diag(dat.svd$d)
  
  xseq = seq(0,1,l=500)
  
  if(length(nComponents)>1){
    temp = ceiling(sqrt(length(nComponents)))
    if(temp*(temp-1)>=length(nComponents)){
      par(mfrow=c(temp,temp-1),mar=rep(0.1,4))
    }else{
      par(mfrow=rep(temp,2),mar=rep(0.1,4))
    }
  }
  
  for(i in nComponents){
    
    if(type[1] == 'fourier'){
      funValsVOrig = t(sapply(1:nrow(VV[,1:i]),
                              function(x)fun1(xseq,VV[x,1:i])))
    }else{
      spl <- bs(xseq,df=i,degree=min(3,i-1),intercept=T)
      Lambda = 
        diff(xseq[1:2])*crossprod(spl,diag(rep(c(1/2,1,1/2),c(1,length(xseq)-2,1)))%*%spl)
      L.eig = eigen(Lambda)
      L.invSqrt = tcrossprod(L.eig$vectors%*%diag(1/sqrt(L.eig$values)),
                             L.eig$vectors)
      PsiOrig = VV[,1:i]%*%L.invSqrt
      funValsVOrig = tcrossprod(PsiOrig,spl)
    }
    
    XYAbs = abs(XY)
    alphaX1 = round(1 + (XYAbs-min(XYAbs))/diff(range(XYAbs))*499)
    if(is.null(COLSX)){
      COLSX = rainbow(500,s = 0.8, v = 0.85,start=0,end=1/6)[500:1][alphaX1]
    }
    alphaX1 = seq(alphaX[1],alphaX[2],l=500)[alphaX1]
    
    if(length(linesX)==1)linesX = rep(linesX,pp)
    if(length(COLSX)==1)COLSX = rep(COLSX,pp)
    
    
    plot(0,xlim=c(0,1.25),ylim=range(funValsVOrig), 
       type='n', lwd=2, col=1,
       xlab="", ylab="",cex.axis=2,cex.lab=2,xaxt='n',yaxt='n',bty='n')
    for(ell in 1:pp)lines(funValsVOrig[ell,]~xseq,
                          lty=linesX[ell],lwd=3,
                          col=adjustcolor(COLSX[ell],alpha.f=alphaX1[ell]))
    segments(0,0,1,0,lty='22',lwd=2)
    
  }
  
}
