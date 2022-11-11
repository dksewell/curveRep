\name{curveDissim}
\alias{curveDissim}
\title{
Curvelinear Representation of Dissimilarity Matrices
}
\description{
This function ties the row splines with the column splines to represent the whole dissimilarity matrix.
}
\usage{
curveDissim(dat, Y = NULL, nComponents = nrow(dat), 
            COLS.Obs = brewer.pal(12, "Set3"), nClust = 10, 
            alpha.Obs = 0.1, lines.Obs = 1, 
            poly.Obs = FALSE, border.Obs = NA, 
            splatter.Obs = NULL, splatter.Thresh = 0.05, 
            splatter.Weight = 1/nrow(dat), splatter.Ask = TRUE, 
            flipV = FALSE, ylimOverride = c(1, 1), xlimOverride = c(1, 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
A matrix object.
}
  \item{Y}{
Vector of length \code{nrow(dat)}.  Provides grouping information on the objects of the dissimilarity matrix.
}
  \item{nComponents}{
Number of components for a low rank representation of the data.  Default is to use all components and thus exactly represent the matrix.  It is recommended to incrementally increase \code{nComponents} until no further structure in the data is revealed.
}
  \item{COLS.Obs}{
A vector of colors of length \code{unique(Y.Obs)}.
}
  \item{nClust}{
If \code{Y} is not supplied, the maximum number of clusters to consider using the functionality of the \code{mclust} package (clustering is performed on the row curves' coefficients).  If \code{nClust=1} or if \code{Y.Obs} is not \code{NULL}, no clustering is performed.
}
  \item{alpha.Obs}{
Transparency level for the row curves (ranging from 0 to 1).
}
  \item{lines.Obs}{
Line type for row curves.
}
  \item{poly.Obs}{
Logical. Should, for each group of row curves, the region containing the curves be highlighted?
}
  \item{border.Obs}{
If \code{poly.Obs=TRUE}, the color to draw the border of the highlighted region.  See \code{?polygon}.
}
  \item{splatter.Obs}{
Number of splatting iterations to run on the row curves.  If \code{NULL}, no splatting is done.
}
  \item{splatter.Thresh}{
Proportion of neighbors used in the splatting algorithm.  See Zhou et al. 2009.
}
  \item{splatter.Weight}{
Weight used in the splatting algorithm.  See Zhou et al. 2009.
}
  \item{splatter.Ask}{
Logical.  If \code{TRUE}, the splatting algorithm will stop after \code{splatter.Obs} and ask the user if/how many more iterations ought to be done.
}
  \item{flipV}{
Logical.  If \code{TRUE}, then for \eqn{s<0}, \eqn{-f(s)} is plotted instead.  Often useful for making the visualization appear more smooth.
}
  \item{ylimOverride}{
Vector of length 2.  Provides the user the opportunity to manually override the default ylim values by scaling them elementwise- hence the default \code{ylimOverride=c(1,1)} does not alter the default ylim values.
}
  \item{xlimOverride}{
Vector of length 2.  Provides the user the opportunity to manually override the default xlim values by scaling them elementwise- hence the default \code{xlimOverride=c(1,1)} does not alter the default xlim values.
}
}
\details{
This function essentially performs \code{curveRep} with spline functions for both the rows and columns with no intercept given.  Since these correspond to the same objects in dissimilarity/distance matrices, this function ties the row and column functions together by plotting the column curves \eqn{f(-s)} for \eqn{1<s<0}.
}
\references{
Sewell, Daniel K. 2017+.  ``Visualizing data through curvilinear representations of matrices.''
}
\author{
Daniel K. Sewell
}
\seealso{
curveRep
}
\examples{
data(animalMutation)
Y = c(1,1,2,2,2,2,2,2,3,3,3,3,4,5,6,7,7,8,8,8)
COLS = brewer.pal(9,'Set1')[-6]
COLS[-c(1:2,5:6)] = gray(1,0)
par(mar=rep(0.1,4))
curveDissim(animalMutation,Y=Y,
            nComponents=20,#8,
            COLS.Obs=COLS,
            nClust=10,
            alpha.Obs=0.9,
            lines.Obs=1,
            poly.Obs=FALSE,
            border.Obs=NA,
            splatter.Obs=NULL,
            splatter.Thresh=0.05,
            splatter.Weight=1/nrow(dat),
            splatter.Ask=TRUE,
            flipV = TRUE,
            ylimOverride=c(0.575,0.625))
legend('bottomright',col=adjustcolor(COLS[c(1:2,5:6)]),
       pch=16,cex=1,
       legend=c('Primates','Mammals','Snake','Fish'))
}
