\name{curveRep-package}
\alias{curveRep-package}
\docType{package}
\title{
Visualizing Data Through Curvilinear Representations of Matrices
}
\description{
This package visualizes any arbitrary matrix by first mapping the matrix to a set of functions. These functions can exactly represent the matrix, distances between objects, or the correlation structure in the data.
}
\details{
\tabular{ll}{
Package: \tab curveRep\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2017-08-23\cr
License: \tab GPL (>= 2)\cr
}

The primary function is the eponymous \code{curveRep()}.  This function can be used on any arbitrary data matrix.  For dissimilarity or distance matrices, you may also use \code{curveDissim()}.  To visualize MANOVA contrasts or compare covariates' relationships with a response variable of interest use \code{curveManova()} and \code{curveReg()} respectively.  See Sewell (2017+) for more details. 
}
\author{
Daniel K. Sewell

Maintainer: Daniel K. Sewell
}
\references{
Sewell, Daniel K. 2017+.  ``Visualizing data through curvilinear representations of matrices.''
}
\keyword{Andrews curves; Basis splines; Fourier series; Singular value decomposition; Visualization}
\examples{
data(wine)
dat = scale(as.matrix(wine[,-1]))
curveRep(dat,Y.Obs=wine[,1],nComponents = 2,
         type='fourier',
         COLS.Obs=brewer.pal(3,'Set1')[c(2,1,3)],
         whichPlot='observations',
         splatter.Obs=35,
         splatter.Thresh=0.01,splatter.Weight=1/nrow(dat),
         splatter.Ask=FALSE,metric=FALSE)
}
