curveManova <-
function(dat,X,contr,
                       nPerm=100,
                       type=c("fourier","splines"),
                       nComponents=min(dim(dat),ncol(X)),
                       COLS.Contr='black',
                       COLS.Perm=brewer.pal(9,'Set1'),
                       alpha.Perm1=0.2,
                       alpha.Perm2=0.8,
                       alphaLv=0.05,
                       lines.Contr=1,
                       lines.Perm=1,
                       border.Perm=brewer.pal(9,'Set1')){
  if(is.data.frame(X)) X = model.matrix(~X)
  if(is.data.frame(contr)) contr = model.matrix(~contr)
  if(!is.null(contr))if(is.null(dim(contr))) contr = matrix(contr,nrow=1)
  
  if(length(intersect(class(dat),c('matrix','dgeMatrix','dgCMatrix',
                                   'dgRMatrix','dsyMatrix'))) == 0) dat=as.matrix(dat)
  
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
  
  dat.svd = svds(dat,k=max(nComponents))
  UU = dat.svd$u%*%diag(sqrt(dat.svd$d))
  VV = dat.svd$v%*%diag(sqrt(dat.svd$d))
  XtXInv = qr.solve(crossprod(X))
  HH = tcrossprod(X%*%XtXInv,X)
  ResHat = (diag(nrow(dat))-HH)%*%UU
  SigHat = crossprod(ResHat)/(nrow(UU)-1)
  if(!is.null(contr)){
    HNew = tcrossprod(contr%*%XtXInv,X)
    UHat = HNew%*%UU
  }
  xseq = seq(0,1,l=500)

  if(nrow(contr)>1){
    temp = ceiling(sqrt(nrow(contr)))
    if(temp*(temp-1)>=nrow(contr)){
      par(mfrow=c(temp,temp-1),mar=rep(0.1,4))
    }else{
      par(mfrow=rep(temp,2),mar=rep(0.1,4))
    }
  }else{
    par(mfrow=c(1,1),mar=rep(0.1,4))
  }
  
  if(type[1] == 'fourier'){
    
    funValsUOrig = t(sapply(1:nrow(UHat),
                        function(x)fun1(xseq,UHat[x,1:nComponents])))
    funValsUPerm=list()
    for(cc in 1:nrow(contr)){
      funValsUPerm[[cc]] = matrix(0,nPerm,length(xseq))
      for(it in 1:nPerm){
        UPerm = HNew[cc,]%*%UU[sample(nrow(UU)),]
        funValsUPerm[[cc]][it,] = 
          fun1(xseq,UPerm[1:nComponents])
      }
    }
    
  }else{
    
    spl <- bs(xseq,df=nComponents,degree=min(3,nComponents),intercept=F)
    Lambda = 
      diff(xseq[1:2])*crossprod(spl,diag(rep(c(1/2,1,1/2),c(1,length(xseq)-2,1)))%*%spl)
    L.eig = eigen(Lambda)
    L.invSqrt = tcrossprod(L.eig$vectors%*%diag(1/sqrt(L.eig$values)),
                           L.eig$vectors)
    PsiOrig = UHat[,1:nComponents]%*%L.invSqrt
    funValsUOrig = tcrossprod(PsiOrig,spl)
    funValsUPerm = list()
    for(cc in 1:nrow(contr)){
      funValsUPerm[[cc]] = matrix(0,nPerm,length(xseq))
      for(it in 1:nPerm){
        UPerm = HNew[cc,]%*%UU[sample(nrow(UU)),1:nComponents]%*%L.invSqrt
        funValsUPerm[[cc]][it,] = tcrossprod(UPerm,spl)
      }
    }
  }
  
  for(cc in 1:nrow(contr)){
    
    plot(0,xlim=c(0,1),ylim=range(c(funValsUOrig[cc,],c(funValsUPerm[[cc]]))), 
         type='n', lwd=2, col=1,
         xlab="", ylab="",cex.axis=2,cex.lab=2,xaxt='n',yaxt='n',bty='n')
    
    yU = sapply(1:length(xseq),function(x) 
      quantile(funValsUPerm[[cc]][,x],1-alphaLv/2))
    yL = sapply(1:length(xseq),function(x) 
      quantile(funValsUPerm[[cc]][,x],alphaLv/2))
    polygon(x=c(xseq[1],xseq,xseq[length(xseq)],xseq[length(xseq):1]),
            y=c(yL[1],yU,yL[length(xseq)],yL[length(xseq):1]),
            border=border.Perm,
            col=adjustcolor(COLS.Perm,alpha.f=alpha.Perm2))
    if(!is.null(lines.Perm)){
      matlines(xseq,t(funValsUPerm[[cc]]),lty=lines.Perm,lwd=3,
               col=adjustcolor(COLS.Perm,alpha.f=alpha.Perm1))
    }
    lines(drop(funValsUOrig[cc,])~xseq,lty=lines.Contr,lwd=3,col=COLS.Contr)
    segments(0,0,1,0,lty='22',lwd=2)
    
  }

}
