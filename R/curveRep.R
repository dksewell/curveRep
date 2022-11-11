curveRep <-
function(dat,Y.Obs=NULL,Y.Vars=NULL,
                  type=c('splines','fourier'),
                  nComponents=min(dim(dat)),
                  COLS.Obs=brewer.pal(12,'Set3'),
                  COLS.Vars=brewer.pal(12,'Set3')[12:1],
                  whichPlot=c('observations','variables'),
                  nClust.Obs=10,nClust.Vars=5,
                  alpha.Obs=0.1,alpha.Vars=0.8,
                  lines.Obs=1,lines.Vars=2,
                  poly.Obs=FALSE,poly.Vars=FALSE,
                  border.Obs=NA,border.Vars=NA,
                  splatter.Obs=NULL,splatter.Vars=NULL,
                  splatter.Thresh=0.05,
                  splatter.Weight=1/dim(dat),
                  splatter.Ask=TRUE,
                  splatter.DistU=NULL,
                  splatter.DistV=NULL,
                  metric=TRUE,
                  ylimOverride=c(1,1)){
  
  if(!class(dat)%in%c('matrix','dgeMatrix','dgCMatrix',
                      'dgRMatrix','dsyMatrix')) dat=as.matrix(dat)
  
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
  if(metric){
    UU = dat.svd$u%*%diag(dat.svd$d)
    VV = dat.svd$v%*%diag(dat.svd$d)
  }else{
    UU = dat.svd$u%*%diag(sqrt(dat.svd$d))
    VV = dat.svd$v%*%diag(sqrt(dat.svd$d))
  }
  
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
    if(is.null(Y.Obs)){
      if(nClust.Obs>1){
        Y.Obs = Mclust(UU[,1:i],G=2:nClust.Obs)$class
        }else{Y.Obs = rep(1,nrow(dat)) }
    }
    if(is.null(Y.Vars) & 'variables' %in% whichPlot){
      if(nClust.Vars>1){
        Y.Vars = Mclust(VV[,1:i],G=2:nClust.Vars)$class
      }else{Y.Vars = rep(1,ncol(dat)) }
    }
    
    if('observations'%in%whichPlot){
      if(type[1] == 'fourier'){
        funValsU = t(sapply(1:nrow(UU),
                          function(x)fun1(xseq,UU[x,1:i])))
      }else{
        spl <- bs(xseq,df=i,degree=min(3,i-1),intercept=T)
        Lambda = 
          diff(xseq[1:2])*crossprod(spl,diag(rep(c(1/2,1,1/2),c(1,length(xseq)-2,1)))%*%spl)
        L.eig = eigen(Lambda)
        L.invSqrt = tcrossprod(L.eig$vectors%*%diag(1/sqrt(L.eig$values)),
                               L.eig$vectors)
        Psi = UU[,1:i]%*%L.invSqrt
        funValsU = tcrossprod(Psi,spl)
      }
    }else{funValsU = NULL}
    if('variables'%in% whichPlot){
      if(type[1] == 'fourier'){
        funValsV = t(sapply(1:nrow(VV),
                          function(x)fun1(xseq,VV[x,1:i])))
      }else{
        spl <- bs(xseq,df=i,degree=min(3,i-1),intercept=T)
        Lambda = 
          diff(xseq[1:2])*crossprod(spl,diag(rep(c(1/2,1,1/2),c(1,length(xseq)-2,1)))%*%spl)
        L.eig = eigen(Lambda)
        L.invSqrt = tcrossprod(L.eig$vectors%*%diag(1/sqrt(L.eig$values)),
                               L.eig$vectors)
        Psi = VV[,1:i]%*%L.invSqrt
        funValsV = tcrossprod(Psi,spl)
      }
    }else{funValsV = NULL}
    plot(0,xlim=c(0,1),ylim=range(c(c(funValsU),c(funValsV)))*ylimOverride, type='n', lwd=2, col=1,
         xlab="", ylab="",cex.axis=2,cex.lab=2,xaxt='n',yaxt='n',bty='n')
    
    if('observations' %in% whichPlot){
      if(poly.Obs){
        yU = yL = list()
        for(ii in 1:length(unique(Y.Obs))){
          indY = which(Y.Obs==ii)
          yU[[ii]] = sapply(1:length(xseq),function(x) max(funValsU[indY,x]))
          yL[[ii]] = sapply(1:length(xseq),function(x) min(funValsU[indY,x]))
        }
        for(ii in 1:length(unique(Y.Obs))){
          polygon(x=c(xseq[1],xseq,xseq[length(xseq)],xseq[length(xseq):1]),
                  y=c(yL[[ii]][1],yU[[ii]],yL[[ii]][length(xseq)],yL[[ii]][length(xseq):1]),
                  border=border.Obs[ii],
                  col=adjustcolor(COLS.Obs[ii],alpha.f=alpha.Obs))
        }
      }
      if(!is.null(lines.Obs)){
        if(is.null(splatter.Obs)){
          matlines(xseq,t(funValsU),lty=lines.Obs,lwd=3,
                   col=adjustcolor(COLS.Obs[Y.Obs],alpha.f=alpha.Obs))#alpha.Obs))
        }else{
          if(is.null(splatter.DistU)){
            DD = as.matrix(dist(funValsU))
          }else{
            if(class(splatter.DistU)=='dist'){
              DD = as.matrix(splatter.DistU)
            }else{
              DD= splatter.DistU
            }
            diag(DD) = Inf
          }
          DDThresh = quantile(DD[upper.tri(DD)],probs=splatter.Thresh)
          Opac = rep(1,nrow(dat))
          
          userSaysKeepGoing = TRUE
          while(userSaysKeepGoing){
            for(tt in 1:splatter.Obs){
              for(ss in 1:nrow(dat)){
                neigh = which(DD[ss,]<DDThresh)
                Opac[neigh] = Opac[neigh]*(1 + splatter.Weight[1]*exp(-2*DD[ss,neigh]/DDThresh^2))
              }
              Opac = Opac/max(Opac)
            }
            for(ss in 1:nrow(dat)){
              lines(funValsU[ss,]~xseq,lty=lines.Obs,lwd=3,
                       col=adjustcolor(COLS.Obs[Y.Obs[ss]],alpha.f=Opac[ss]))
            }
            if(!splatter.Ask){userSaysKeepGoing='n'}
            while(!userSaysKeepGoing %in% c('y','n','Y','N')){
              cat("Keep splatting?  Enter 'y' or 'n':")
              userSaysKeepGoing = readline()
            }
            if(userSaysKeepGoing %in% c('y','Y')){
              cat(paste('Just finished',splatter.Obs,'iterations \n '))
              splatter.Obs = NA
              options(warn=-1)
              while(is.na(as.integer(splatter.Obs))){
                cat("How many more iterations?  Enter an integer: \n")
                splatter.Obs = readline()
                splatter.Obs = as.integer(splatter.Obs)
              }
              options(warn=0)
              plot(0,xlim=c(0,1),ylim=range(c(c(funValsU),c(funValsV)))*ylimOverride, 
                   type='n', lwd=2, col=1,
                   xlab="", ylab="",cex.axis=2,cex.lab=2,xaxt='n',yaxt='n',bty='n')
              userSaysKeepGoing = TRUE
            }else{
              userSaysKeepGoing = FALSE
            }
          }
        }
      }
    }
    
    if('variables' %in% whichPlot){
      if(poly.Vars){
        yU = yL = list()
        for(ii in 1:length(unique(Y.Vars))){
          indY = which(Y.Vars==ii)
          yU[[ii]] = sapply(1:length(xseq),function(x) max(funValsV[indY,x]))
          yL[[ii]] = sapply(1:length(xseq),function(x) min(funValsV[indY,x]))
        }
        for(ii in 1:length(unique(Y.Vars))){
          polygon(x=c(xseq[1],xseq,xseq[length(xseq)],xseq[length(xseq):1]),
                  y=c(yL[[ii]][1],yU[[ii]],yL[[ii]][length(xseq)],yL[[ii]][length(xseq):1]),
                  border=border.Vars[ii],
                  col=adjustcolor(COLS.Vars[ii],alpha.f=alpha.Vars))
        }
      }
      if(!is.null(lines.Vars)){
        if(is.null(splatter.Vars)){
          matlines(xseq,t(funValsV),lty=lines.Vars,
                   lwd=3,col=adjustcolor(COLS.Vars[Y.Vars],alpha.f=alpha.Vars))
        }else{
          if(is.null(splatter.DistV)){
            DD = as.matrix(dist(funValsV))
          }else{
            if(class(splatter.DistV)=='dist'){
              DD = as.matrix(splatter.DistV)
            }else{
              DD= splatter.DistV
            }
            diag(DD) = Inf
          }
          DDThresh = quantile(DD[upper.tri(DD)],probs=splatter.Thresh)
          Opac = rep(1,ncol(dat))
          
          userSaysKeepGoing = TRUE
          while(userSaysKeepGoing){
            for(tt in 1:splatter.Vars){
              for(ss in 1:ncol(dat)){
                neigh = which(DD[ss,]<DDThresh)
                Opac[neigh] = Opac[neigh]*(1 + splatter.Weight[2]*exp(-2*DD[ss,neigh]/DDThresh^2))
              }
              Opac = Opac/max(Opac)
            }
            for(ss in 1:ncol(dat)){
              lines(funValsV[ss,]~xseq,lty=lines.Vars,lwd=3,
                    col=adjustcolor(COLS.Vars[Y.Vars[ss]],alpha.f=Opac[ss]))
            }
            if(!splatter.Ask){userSaysKeepGoing='n'}
            while(!userSaysKeepGoing %in% c('y','n','Y','N')){
              cat("Keep splatting?  Enter 'y' or 'n':")
              userSaysKeepGoing = readline()
            }
            if(userSaysKeepGoing %in% c('y','Y')){
              cat(paste('Just finished',splatter.Vars,'iterations \n '))
              splatter.Vars = NA
              options(warn=-1)
              while(is.na(as.integer(splatter.Vars))){
                cat("How many more iterations?  Enter an integer: \n")
                splatter.Vars = readline()
                splatter.Vars = as.integer(splatter.Vars)
              }
              options(warn=0)
              plot(0,xlim=c(0,1),ylim=range(c(c(funValsU),c(funValsV)))*ylimOverride, 
                   type='n', lwd=2, col=1,
                   xlab="", ylab="",cex.axis=2,cex.lab=2,xaxt='n',yaxt='n',bty='n')
              userSaysKeepGoing = TRUE
            }else{
              userSaysKeepGoing = FALSE
            }
          }
        }
      }
    }
    segments(0,0,1,0,lty='22',lwd=2)
  }
  
  par(mfrow=c(1,1),mar=c(5, 4, 4, 2) + 0.1)
}
