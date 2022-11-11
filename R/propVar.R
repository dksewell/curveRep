propVar <-
function(dat,k=min(dim(dat))){
  dat.svd = svds(dat,k=min(c(k,dim(dat))))
  barplot(cumsum(dat.svd$d^2)/sum(dat.svd$d^2),ylim=c(0,1),
       xlab='Number of Components',
       ylab='Proportion of Variance',col='steelblue',
       cex.lab=2,cex.axis=2)
}
