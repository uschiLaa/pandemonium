#' Write WC, coordinates and cluter assingment to a CSV file
#'
#' Settings include: coord, useCov, metric, linkage, k, plotType
#'
#' @param pred prediction matrix
#' @param covInv inverse covariance matrix
#' @param wc matrix specifying the model parameters (on a grid)
#' @param exp observable reference value (e.g. experimental measurement)
#' @param settings list specifying parameters usually selected in the app
#' @param filename path to write the results file to
#' @param user_coord input coordinate matrix (optional)
#' @param user_dist input distance matrix (optional)
#' @export
writeResults <- function(pred, covInv, wc, exp, settings, filename,
                      user_coord=NULL, user_dist=NULL){
  n <- nrow(pred)
  chi2 <- computeChi2(pred, covInv, exp)
  sig <- computeSigma(chi2, 2)
  bf <- which.min(chi2)
  sm <- which.min(rowSums(abs(wc)))
  x <- colnames(wc)[1]
  y <- colnames(wc)[2]
  cond <- 1:nrow(wc)
  coord <- getCoords(settings$coord, settings$useCov, pred, covInv, exp, user_coord)
  dists <- getDists(coord, settings$metric, user_dist)
  fit <- stats::hclust(dists, settings$linkage)
  groups <- stats::cutree(fit, k=settings$k)
  lvl <- unique(groups[stats::order.dendrogram(stats::as.dendrogram(fit))])
  cluster <- as.numeric(factor(groups, levels= lvl))
  benchmarks <- getBenchmarkInformation(as.matrix(dists), groups)
  isBenchmark <- rep(0, nrow(wc))
  isBenchmark[benchmarks$id] <- 1
  coord <- as.data.frame(coord)
  colnames(coord) <- paste0("O", 1:ncol(coord))
  wc %>%
    cbind(coord) %>%
    cbind(cluster) %>%
    cbind(isBenchmark) %>%
    write.csv(filename, row.names = FALSE, quote = FALSE)

}
