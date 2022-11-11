\name{curveManova}
\alias{curveManova}
\title{
Visualize MANOVA Contrasts
}
\description{
Visualize MANOVA contrasts and compare the corresponding curves with those obtained from permuted data. 
}
\usage{
curveManova(dat, X, contr, nPerm = 100, 
            type = c("fourier","splines"), 
            nComponents = min(dim(dat), ncol(X)), 
            COLS.Contr = "black", COLS.Perm = brewer.pal(9, "Set1"), 
            alpha.Perm1 = 0.2, alpha.Perm2 = 0.8, 
            alphaLv = 0.05, 
            lines.Contr = 1, lines.Perm = 1, 
            border.Perm = brewer.pal(9, "Set1"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
A matrix object.
}
  \item{X}{
A \code{nrow(dat)}x\eqn{p} matrix object representing the design matrix
}
  \item{contr}{
A \eqn{m}x\eqn{p} matrix object giving the contrast coefficients, where \eqn{m} is the number of contrasts to be evaluated.
}
  \item{nPerm}{
Integer.  Number of permutations with which to compare the original contrast curve.
}
  \item{type}{
Character.  One of either "fourier" or "splines".
}
  \item{nComponents}{
Number of components for a low rank representation of the data.  Default is to use all components and thus exactly represent the matrix.  It is recommended to incrementally increase \code{nComponents} until no further structure in the data is revealed.
}
  \item{COLS.Contr}{
Color for the original contrast curves.
}
  \item{COLS.Perm}{
Color for the contrast curves of the permuted data.
}
  \item{alpha.Perm1}{
Transparency level for the permutation curves (ranging from 0 to 1).
}
  \item{alpha.Perm2}{
Transparency level for the shaded region indicating the (1-\code{alphaLv})100\% pointwise containment interval.
}
  \item{alphaLv}{
Numeric ranging from 0 to 1 giving the (1-\code{alphaLv})100\% pointwise containment interval.
}
  \item{lines.Contr}{
Line type for original contrast curves.
}
  \item{lines.Perm}{
Line type for permutation contrast curves.
}
  \item{border.Perm}{
The color to draw the border of the highlighted region.  See \code{?polygon}.
}
}
\details{
Consider the MANOVA model:
\deqn{
Y = XB + E.
}
For some matrix of contrasts \eqn{C}, we estimate the contrasts with \eqn{C(X'X)^{-1}X'Y}.  A linear combination of Y can be functionally represented simply using appropriately weighted coefficients of the row curves of Y.  By permuting Y and constructing similar contrast curves, we may compare the observed contrasts with what we might expect by chance.  Note that (1) no distributional assumptions are made on E, and (2) this is still not a formal testing procedure- only a way to accrue visual evidence against a null hypothesis.
}
\references{
Sewell, Daniel K. 2017+.  ``Visualizing data through curvilinear representations of matrices.''
}
\author{
Daniel K. Sewell
}
\examples{
data(lymphoma)
dat = scale(lymphoma$x)
X = as.matrix(model.matrix(~factor(lymphoma$y))[,-1])
X[rowSums(X)==0,] = -1
CONTR = matrix(c(-2,-1,
                 -1,-2,
                 1,-1),3,2,byrow = TRUE)
  curveManova(dat,X,
              contr=CONTR,
              nPerm=1000,
              type='fourier',
              nComponents=20,
              COLS.Contr='black',
              COLS.Perm='mediumseagreen',
              alpha.Perm1=0.1,
              alpha.Perm2=0.5,
              lines.Contr=1,
              lines.Perm=1,
              border.Perm='darkseagreen',
              alphaLv=0.05)
}
