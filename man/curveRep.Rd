\name{curveRep}
\alias{curveRep}
\title{
Curvilinear Representation of a Matrix
}
\description{
Represent the rows or columns or both of any arbitrary matrix.
}
\usage{
curveRep(dat, Y.Obs = NULL, Y.Vars = NULL, 
         type = c("splines", "fourier"), 
         nComponents = min(dim(dat)), 
         COLS.Obs = brewer.pal(12, "Set3"), 
         COLS.Vars = brewer.pal(12, "Set3")[12:1], 
         whichPlot = c("observations", "variables"), 
         nClust.Obs = 10, nClust.Vars = 5, 
         alpha.Obs = 0.1, alpha.Vars = 0.8, 
         lines.Obs = 1, lines.Vars = 2, 
         poly.Obs = FALSE, poly.Vars = FALSE, 
         border.Obs = NA, border.Vars = NA, 
         splatter.Obs = NULL, splatter.Vars = NULL, 
         splatter.Thresh = 0.05, splatter.Weight = 1/dim(dat), 
         splatter.Ask = TRUE, 
         splatter.DistU = NULL, splatter.DistV = NULL, 
         metric = TRUE, ylimOverride = c(1, 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
A matrix object.
}
  \item{Y.Obs}{
Vector of length \code{nrow(dat)}.  Provides grouping information on the objects of the data matrix.
}
  \item{Y.Vars}{
Vector of length \code{ncol(dat)}.  Provides grouping information on the variables of the data matrix.
}
  \item{type}{
Character.  Either "splines" or "fourier", depending on the space of functions desired.
}
  \item{nComponents}{
Number of components for a low rank representation of the data.  Default is to use all components and thus exactly represent the matrix.  It is recommended to incrementally increase \code{nComponents} until no further structure in the data is revealed.
}
  \item{COLS.Obs}{
A vector of colors of length \code{unique(Y.Obs)}.
}
  \item{COLS.Vars}{
A vector of colors of length \code{unique(Y.Vars)}.
}
  \item{whichPlot}{
A character vector of length 1 or 2.  Elements of this vector must be either "observations" or "variables".
}
  \item{nClust.Obs}{
If \code{Y.Obs} is not supplied, the maximum number of clusters to consider using the functionality of the \code{mclust} package (clustering is performed on the row curves' coefficients).  If \code{nClust.Obs=1} or if \code{Y.Obs} is not \code{NULL}, no clustering is performed.
}
  \item{nClust.Vars}{
If \code{Y.Vars} is not supplied, the maximum number of clusters to consider using the functionality of the \code{mclust} package (clustering is performed on the column curves' coefficients).  If \code{nClust.Vars=1} or if \code{Y.Vars} is not \code{NULL}, no clustering is performed.
}
  \item{alpha.Obs}{
Transparency level for the row curves (ranging from 0 to 1).
}
  \item{alpha.Vars}{
Transparency level for the column curves (ranging from 0 to 1).
}
  \item{lines.Obs}{
Line type for row curves.
}
  \item{lines.Vars}{
Line type for column curves.
}
  \item{poly.Obs}{
Logical. Should, for each group of row curves, the region containing the curves be highlighted?
}
  \item{poly.Vars}{
Logical. Should, for each group of column curves, the region containing the curves be highlighted?
}
  \item{border.Obs}{
If \code{poly.Obs=TRUE}, the color to draw the border of the highlighted region.  See \code{?polygon}.
}
  \item{border.Vars}{
If \code{poly.Vars=TRUE}, the color to draw the border of the highlighted region.  See \code{?polygon}.
}
  \item{splatter.Obs}{
Number of splatting iterations to run on the row curves.  If \code{NULL}, no splatting is done.
}
  \item{splatter.Vars}{
Number of splatting iterations to run on the column curves.  If \code{NULL}, no splatting is done.
}
  \item{splatter.Thresh}{
Proportion of neighbors used in the splatting algorithm.  See Zhou et al. 2009.
}
  \item{splatter.Weight}{
Weight used in the splatting algorithm.  See Zhou et al. 2009.
}
  \item{splatter.Ask}{
Logical.  If \code{TRUE}, the splatting algorithm will stop after \code{splatter.Obs} (or \code{splatter.Vars}) and ask the user if/how many more iterations ought to be done.
}
  \item{splatter.DistU}{
Optional distance matrix to use in the splatting algorithm for the row curves.
}
  \item{splatter.DistV}{
Optional distance matrix to use in the splatting algorithm for the column curves.
}
  \item{metric}{
**IMPORTANT!!!**  If \code{TRUE}, the curves must be interpreted according to the integrated squared differences.  These integrals correspond to the Euclidean distance between the rows or the covariance structure (on centered data) for the columns.  Note that the row and the column curves cannot be directly compared.  If \code{FALSE}, the curves must be interpreted according to the integrated product of the row and column functions.  These integrals correspond exactly to the elements of the original data matrix.
}
  \item{ylimOverride}{
Vector of length 2.  Provides the user the opportunity to manually override the default ylim values by scaling them elementwise- hence the default \code{ylimOverride=c(1,1)} does not alter the default ylim values.
}
}

\details{
Any matrix \eqn{M} can be mapped to space of functions such that the dot product, defined as \eqn{\int_0^1 f_i(s)g_j(s)ds} equals \eqn{M_{ij}}.  There also exists a mapping such that the Euclidean distance between rows i and j equals \eqn{\int_0^1(f_i(s)-f_j(s)ds)^2ds}.  This latter mapping is plotted if \code{metric=TRUE}.
}
\references{
Sewell, Daniel K. 2017+.  ``Visualizing data through curvilinear representations of matrices.''
}
\author{
Daniel K. Sewell
}

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
