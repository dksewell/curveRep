\name{propVar}
\alias{propVar}
\title{
Plot Proportion of Variance
}
\description{
Plots proportion of variance explained by increasing number of components.
}
\usage{
propVar(dat, k = min(dim(dat)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{a matrix object}
  \item{k}{Maximum number of components to evaluate}
}
\details{
Plot of the cumulative sum of the normalized squared singular values (i.e., the eigenvalues of \code{dat}'\code{dat}.
}
\author{
Daniel K. Sewell
}

\seealso{
rARPACK::svds
}
\examples{
data(lymphoma)
propVar(lymphoma$x,25)
}
