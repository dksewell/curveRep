\name{curveReg}
\alias{curveReg}
\title{
Linear Regression Association Visualization
}
\description{
Visualize and compare the relationships between multiple covariates and a response variable.
}
\usage{
curveReg(Y, X, 
         type = c("fourier", "splines"), 
         nComponents = min(dim(dat), ncol(X)), 
         COLSX = NULL, alphaX = c(0.1, 0.9), 
         linesX = 1)
}
\arguments{
  \item{Y}{
Numeric vector giving the response variable. 
}
  \item{X}{
Matrix of covariates/design matrix.  Do not include the intercept.
}
  \item{type}{
Character.  Either "splines" or "fourier", depending on the space of functions desired.
}
  \item{nComponents}{
Number of components for a low rank representation of the data.  Default is to use all components and thus exactly represent the matrix.  It is recommended to incrementally increase \code{nComponents} until no further structure in the data is revealed.
}
  \item{COLSX}{
Colors for the p covariates.  Either a single color or a vector of length p of colors.  If \code{NULL}, \code{rainbow()} will be used, graded on absolute magnitude of correlations with Y. 
}
  \item{alphaX}{
Vector of length 2 giving lower and upper bound on transparency of the curves.  Those variables with stronger correlation with Y will have more opacity.
}
  \item{linesX}{
Line type for the p covariates.  Either a single line type or a vector of length p.
}
}
\details{
Consider the matrix whose \eqn{j^{th}} row is the one-dimensional projection of the centered response variable on the \eqn{j^{th}} covariate.  We can then visualize the relationships between these projections.  Specifically, we can map this matrix to a set of functions such that the integrated squared distance between two covariates' functions equals the variance in their projections.  Hence we may find a parsimonious model with good explanatory power by finding curves with large magnitudes that explore different parts of the plot. 
}
\references{
Sewell, Daniel K. 2017+.  ``Visualizing data through curvilinear representations of matrices.''
}
\author{
Daniel K. Sewell
}
\examples{
set.seed(1)
n = 200
pp = 96
rho = 0.75
X1 = svd(matrix(rnorm(n*2),n,2))$u
X1 = scale(X1)
X2 = rho*X1 + sqrt(1-rho^2)*matrix(rnorm(n*2),n,2)
X = cbind(X1,X2,matrix(rnorm(n*(pp-4)),n,pp-4))
Y = X[,1] + X[,2] - X[,3] - X[,4] + rnorm(n)/10

curveReg(Y,X,
         type='splines',
         nComponents=4,
         COLSX=rep(c(gray(0.5,1),gray(0.5,0.15)),c(4,96)),
         alphaX=c(0.05,0.9),
         linesX=rep(c('F1','44','13','1343','F1'),c(1,1,1,1,pp-4)))
}
