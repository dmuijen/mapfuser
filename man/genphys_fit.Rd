% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genphys_fit.R
\name{genphys_fit}
\alias{genphys_fit}
\title{Model the relationship between genetical and physical genome positions}
\usage{
genphys_fit(MF.obj, type = c("consensus", "map"), z = 5,
  chromosomes = NULL, map = NULL)
}
\arguments{
\item{MF.obj}{A mapfuser object with a reference map loaded and either a consensus map created with mapfuser or a genetic loaded with read.maps}

\item{type}{Fit physical genome positions vs. consensus genetic map made with mapfuser or an individual genetic map}

\item{z}{discard invididual data points based on z-score threshold for scaled pearson residuals}

\item{chromosomes}{The chromosomes to fit a P-spline, default to all chromosomes}

\item{map}{Name of the genetic map to use when the consensus map is not used for fitting a P-spline and recombination rate calculation}
}
\value{
The input object is returned with added components recombination rate at a 0.1 Mbp interval,
the general additive model fit (gam) with penalized spline fits per chromosome,
and predictions of centiMorgan positions used for fitting. Markers that have been removed due to z-threshold are saved to the config slot.
}
\description{
Fit a penalized thin plate regression spline (Wood, 2003) through genetic position in centiMorgan and physical genome positions in (Mega) base pairs.
}
\examples{
\dontshow{
fpath <- system.file("extdata", package="mapfuser")
maps <- list.files(fpath, pattern = "Bur", full.names = TRUE)
MF.obj <- read_maps(mapfiles = maps, sep = ",", header = TRUE, type = "delim")
ref_file <- list.files(fpath, pattern = "reference", full.names = TRUE)
MF.obj <- read_ref(MF.obj = MF.obj, ref_file = ref_file, sep = ",", header = TRUE, type = "delim")
MF.obj <- genphys_fit(MF.obj,chromosomes = 1, map = "Col-0_Bur-0.csv", type = "map")
}
\dontrun{
MF.obj <- genphys_fit(MF.obj, type = "consensus", z = 5, chromosomes = 1:5, map = NULL)
MF.obj <- genphys_fit(MF.obj, type = "map", z = 5, chromosomes = 1:5, map = "Col-0_Blh-1.csv")
# Plot the result
plot(MF.obj, which = "mareymap", maps = "consensus", chr = 1:5)
}
}
\author{
Dennis van Muijen
}
