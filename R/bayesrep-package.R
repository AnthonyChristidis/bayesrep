#' bayesrep: Bayesian Repulsive Ensembles via Ising Priors
#'
#' Fits ensembles of sparse linear and logistic sub-models whose feature
#' allocations repel one another under an Ising prior, so that collinear
#' predictive pathways are separated across the ensemble instead of being
#' arbitrarily collapsed onto a single representative.
#'
#' The entry point is \code{\link{bre}}. See \code{\link{coallocation}} for the
#' label-invariant allocation summary and \code{\link{waic}} for choosing the
#' number of sub-models.
#'
#' @keywords internal
#'
#' @importFrom stats dnorm plogis predict quantile rnorm sd var
#' @importFrom graphics abline axis box image matplot par plot segments title
#' @importFrom grDevices colorRampPalette
#' @importFrom utils head
"_PACKAGE"
