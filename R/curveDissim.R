curveDissim <-
  function(dat,Y=NULL,
           nComponents=nrow(dat),
           COLS.Obs=brewer.pal(12,'Set3'),
           nClust=10,
           alpha.Obs=0.1,
           lines.Obs=1,
           poly.Obs=FALSE,
           border.Obs=NA,
           splatter.Obs=NULL,
           splatter.Thresh=0.05,
           splatter.Weight=1/nrow(dat),
           splatter.Ask=TRUE,
           flipV = FALSE,
           ylimOverride=c(1,1),
           xlimOverride=c(1,1)){
    
    if(!class(dat)%in%c('matrix','dgeMatrix','dgCMatrix',
                        'dgRMatrix','dsyMatrix')) dat=as.matrix(dat)
    
    dat.svd = svds(dat,k=max(nComponents))
    UU = dat.svd$u%*%diag(dat.svd$d)
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
      if(is.null(Y)){
        if(nClust>1){
          Y = Mclust(UU[,1:i],G=2:nClust)$class
        }else{Y = rep(1,nrow(dat)) }
      }
      
      spl <- bs(xseq,df=i,degree=min(3,i),intercept=F)
      Lambda = 
        diff(xseq[1:2])*crossprod(spl,diag(rep(c(1/2,1,1/2),c(1,length(xseq)-2,1)))%*%spl)
      L.eig = eigen(Lambda)
      L.invSqrt = tcrossprod(L.eig$vectors%*%diag(1/sqrt(L.eig$values)),
                             L.eig$vectors)
      Psi = UU[,1:i]%*%L.invSqrt
      funValsU = tcrossprod(Psi,spl)
      Psi = VV[,1:i]%*%L.invSqrt
      funValsV = tcrossprod(Psi,spl)*(-1)^flipV
      
      plot(0,xlim=c(-1,1)*xlimOverride,ylim=range(c(c(funValsU),c(funValsV)))*ylimOverride, 
           type='n', lwd=2, col=1,xlab="", ylab="",cex.axis=2,cex.lab=2,
           xaxt='n',yaxt='n',bty='n')
      
      if(poly.Obs){
        yU = yL = list()
        for(ii in 1:length(unique(Y))){
          indY = which(Y==ii)
          yU[[ii]] = sapply(1:length(xseq),function(x) max(funValsU[indY,x]))
          yL[[ii]] = sapply(1:length(xseq),function(x) min(funValsU[indY,x]))
        }
        for(ii in 1:length(unique(Y))){
          polygon(x=c(xseq[1],xseq,xseq[length(xseq)],xseq[length(xseq):1]),
                  y=c(yL[[ii]][1],yU[[ii]],yL[[ii]][length(xseq)],yL[[ii]][length(xseq):1]),
                  border=border.Obs[ii],
                  col=adjustcolor(COLS.Obs[ii],alpha.f=alpha.Obs))
        }
      }
      if(!is.null(lines.Obs)){
        if(is.null(splatter.Obs)){
          matlines(xseq,t(funValsU),lty=lines.Obs,lwd=3,
                   col=adjustcolor(COLS.Obs[Y],alpha.f=alpha.Obs))
          matlines(-xseq,t(funValsV),lty=lines.Obs,lwd=3,
                   col=adjustcolor(COLS.Obs[Y],alpha.f=alpha.Obs))
        }else{
          DD = as.matrix(dist(funValsU)); diag(DD) = Inf
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
                    col=adjustcolor(COLS.Obs[Y[ss]],alpha.f=Opac[ss]))
              lines(funValsV[ss,]~I(-1*xseq),lty=lines.Obs,lwd=3,
                    col=adjustcolor(COLS.Obs[Y[ss]],alpha.f=Opac[ss]))
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
              plot(0,xlim=c(-1,1)*xlimOverride,ylim=range(c(c(funValsU),c(funValsV)))*ylimOverride, 
                   type='n', lwd=2, col=1,xlab="", ylab="",cex.axis=2,cex.lab=2,
                   xaxt='n',yaxt='n',bty='n')
              userSaysKeepGoing = TRUE
            }else{
              userSaysKeepGoing = FALSE
            }
          }
        }
      }
      
      segments(-1,0,1,0,lty='22',lwd=2)
    }
    print(range(c(c(funValsU),c(funValsV)))*ylimOverride)
    par(mfrow=c(1,1),mar=c(5, 4, 4, 2) + 0.1)
  }
