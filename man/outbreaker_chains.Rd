% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.outbreaker_chains.R
\name{print.outbreaker_chains}
\alias{print.outbreaker_chains}
\alias{outbreaker_chains}
\alias{plot.outbreaker_chains}
\alias{summary.outbreaker_chains}
\title{Basic methods for processing outbreaker results}
\usage{
\method{print}{outbreaker_chains}(x, n_row = 3, n_col = 8, type = "chain", ...)

\method{plot}{outbreaker_chains}(
  x,
  y = "post",
  type = c("trace", "hist", "density", "cluster", "alpha", "t_inf", "kappa", "network"),
  burnin = 0,
  min_support = 0.1,
  labels = NULL,
  group_cluster = NULL,
  ...
)

\method{summary}{outbreaker_chains}(object, burnin = 0, group_cluster = NULL, ...)
}
\arguments{
\item{x}{an \code{outbreaker_chains} object as returned by \code{outbreaker}.}

\item{n_row}{the number of rows to display in head and tail; defaults to 3.}

\item{n_col}{the number of columns to display; defaults to 8.}

\item{type}{a character string indicating the kind of plot to be used (see details)}

\item{...}{further arguments to be passed to other methods}

\item{y}{a character string indicating which element of an
\code{outbreaker_chains} object to plot}

\item{burnin}{the number of iterations to be discarded as burnin}

\item{min_support}{a number between 0 and 1 indicating the minimum support of
ancestries to be plotted; only used if 'type' is 'network'}

\item{labels}{a vector of length N indicating the case labels (must be
provided in the same order used for dates of symptom onset)}

\item{group_cluster}{a numeric \code{vector} indicating the breaks to 
aggregate the cluster size distribution.}

\item{object}{an \code{outbreaker_chains} object as returned by \code{outbreaker}.}
}
\value{
The form of the value returned by \code{plot} depends on the \code{type}. 
If the type is set as \code{network}, plot returns a visNetwork object 
containing  the details of the inferred transmission trees. Otherwise, it 
returns a ggplot object containing the elements of the plot.

The function \code{summary} returns a list containing 9 elements:

\itemize{

\item{\code{step}}: contains the first and last values of the iteration 
number; the interval between each iteration retained for the output 
(defined by the parameter \code{sample_every} in \code{create_config}),
 and the number of iterations in the output,

\item{\code{post}}: contains the minimum, maximum, mean, median and 
quartiles of the posterior distribution.

\item{\code{like}}: contains the minimum, maximum, mean, median and 
quartiles of the likelihood distribution.

\item{\code{prior}}: contains the minimum, maximum, mean, median and 
quartiles of the prior distribution.

\item{\code{pi}}: contains the minimum, maximum, mean, median and quartiles 
of the conditional report ratio.

\item{\code{a}}: contains the minimum, maximum, mean, median and quartiles 
of the spatial parameter \code{a}.

\item{\code{b}}: contains the minimum, maximum, mean, median and quartiles 
of the spatial parameter \code{b}.

\item{\code{tree}}: a \code{data.frame} that contains the most likely 
infector, the infection date, and the number of missing generations of each
case. It also contains the \code{support} of the most likely branch (i.e. 
the proportion of iterations where the infector of a case is its most 
likely infector), and \code{import}, the proportion of iteration where the 
case was classified as an importation.

\item{\code{cluster}}: a data frame listing the minimum, maximum, median,
mean and quartile of the cluster size distribution.

}
}
\description{
Several methods are defined for instances of the class
\code{outbreaker_chains}, returned by \code{\link{outbreaker}}, including:
\code{print}, \code{plot}
}
\details{
\code{type} indicates the type of graphic to plot:

\itemize{

\item \code{trace} to visualise MCMC traces for parameters or augmented data (plots the
log-likelihood by default)

\item \code{hist} to plot histograms of quantitative values

\item \code{density} to plot kernel density estimations of quantitative values

\item \code{alpha} to visualise the posterior frequency of ancestries

\item \code{network} to visualise the transmission tree; note that
 this opens up an interactive plot and requires a web browser with
 Javascript enabled; the argument `min_support` is useful to select only the
 most supported ancestries and avoid displaying too many links

\item \code{kappa} to visualise the distributions generations between cases and their
ancestor/infector

\item \code{cluster} to visualise the cluster size distribution, grouped by
the value in group_cluster

}
}
\author{
Initial version by Thibaut Jombart, rewritten by Alexis Robert (\email{alexis.robert@lshtm.ac.uk})
}
