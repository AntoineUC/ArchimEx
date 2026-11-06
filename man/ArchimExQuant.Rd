\name{ArchimExQuant}
\alias{ArchimExQuant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extreme conditional quantile estimator for Archimedean copulas
}
\description{
Computes an extreme quantile estimator of Y given X=x, when Y and X are linked by an Archimedean copula.
}
\usage{
ArchimExQuant(tau, X, Y, x, k, copula, ci.level = 0.95,
                          method = "Gaussian", smoothing = FALSE,nboot=1000)
}
\arguments{
  \item{tau}{
Quantile level. Must be between 0 and 1.
}
  \item{X}{
Covariate data matrix.
}
  \item{Y}{
Data vector.
}
  \item{x}{
Prediction point.
}
  \item{k}{
A vector of number of upper order statistics.
}
  \item{copula}{
(Archimedean) Copula model. May be Clayton, Frank, Gumbel or Joe.
}
  \item{ci.level}{
Confidence level. Default: 0.95.
}
  \item{method}{
Confidence interval returned. May be either Gaussian or Bootstrap.
}
  \item{smoothing}{
If TRUE, the c.d.f. estimator of the marginals is smoothed with a Gaussian kernel. Default: FALSE.
}
  \item{nboot}{
Number of bootstrap replication (if method="Boostrap"). Default: 1000.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \item{Point_estimate }{Point estimate of the conditional quantile}
  \item{Lower_bound }{Lower bound of the associated confidence interval}
  \item{Upper_bound }{Upper bound of the associated confidence interval}
}
\references{
Girard, S, Stupfler, G and Usseglio-Carleve, A. Extremal regression in Archimedean copula models.
}
\author{
Antoine Usseglio-Carleve
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%
%}

%\examples{
%
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
