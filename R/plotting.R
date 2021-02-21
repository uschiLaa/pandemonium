#' Make parallel coordinate plot
#'
#' @param dat coordinate representation of points
#' @param gr grouping from clustering
#' @param benchmarkIds index values of benchmarks
#' @param s rescale (default=FALSE)
#' @param a alpha transarancy for drawing non-benchmark points (default=0.2)
#' @return ggplot
#' @export
plotPC <- function(dat, gr, benchmarkIds, s=FALSE, a=0.2){

  colnames(dat) <- paste0("O",1:ncol(dat))

  alphalvl <- rep(a, nrow(dat))
  alphalvl[benchmarkIds] <- 1

  dat %>%
    scale(center = T, scale = s) %>%
    tibble::as_tibble() %>%
    tibble::add_column(gr = factor(gr)) %>%
    tibble::add_column(alphalvl = alphalvl) %>%
    GGally::ggparcoord(columns=1:ncol(dat), groupColumn = "gr",
                       scale = "globalminmax", alphaLines = "alphalvl") +
    ggplot2::scale_color_brewer(palette="Dark2") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.y=ggplot2::element_blank(),
          axis.text.y=ggplot2::element_blank(),
          axis.ticks.y=ggplot2::element_blank(),
          legend.position = "none")
}

#' Show clusters in parameter space
#'
#' Parameter values need to be on a regular grid for correct appearance of
#' this plot.
#'
#' @param wc parameter values as matrix
#' @param x,y variables names (as string) to map to x and y axis
#' @param sm,bf index values for the SM and BF point
#' @param benchmarkIds index values of benchmarks
#' @param col color vector according to cluster assignment
#' @param cond row numbers of points used for conditioning
#' @return ggplot
#' @export
plotWC <- function(wc, x, y, sm, bf, benchmarkIds, col, cond=NULL){
  if(is.null(cond)) cond <- 1:nrow(wc)
  x_id <- which(colnames(wc)==x)
  y_id <- which(colnames(wc)==y)
  ggplot2::ggplot(wc[cond,], ggplot2::aes_string(x, y)) +
    ggplot2::geom_tile(fill= col[cond]) +
    ggplot2::geom_point(ggplot2::aes(x = as.numeric(wc[sm,x_id]),
                                     y = as.numeric(wc[sm,y_id])),
                        shape=1, size=3) +
    ggplot2::geom_point(ggplot2::aes(x = as.numeric(wc[bf,x_id]),
                                     y = as.numeric(wc[bf,y_id])),
                        shape=8, size=3) +
    ggplot2::geom_point(data = wc[benchmarkIds,], shape=5, size=3) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Cluster assignment in parameter space") +
    ggplot2::theme(aspect.ratio = 1, legend.position = "none")
}

#' Plot heatmap with dendrogram
#'
#' @param dat coordinate representation of points
#' @param fit result from hclust
#' @param k number of clusters
#' @param pal color palette
#' @return plot
#' @export
plotHeatmap <- function(dat, fit, k, pal){
  dendo <- stats::as.dendrogram(fit) %>%
    dendextend::set("branches_lwd", 3) %>%
    dendextend::color_branches(k = k, col = pal)

  stats::heatmap(dat, Rowv = dendo, Colv = rev(dendo), scale = "none")
}

#' Plot selected cluster statistics
#'
#' @param dist distances
#' @param fit result from hclust
#' @param chivals vector of chi2 values
#' @param stat cluster statistic to draw
#' @param kmax maximum number of clusters to appear in the plot
#' @return ggplot
#' @export
plotCstat <- function(dist, fit, chivals, stat, kmax=8){
  cstats <- getClusterStats(dist, fit, chivals, kmax)
  ggplot2::ggplot(cstats, ggplot2::aes_string("k", stat)) +
    ggplot2::geom_line() +
    ggplot2::xlab("# clusters") +
    ggplot2::ylab(cstat_names[[stat]]) +
    ggplot2::theme_bw()
}

#' Make coordinate plot
#'
#' Parameter values need to be on a regular grid for correct appearance of
#' this plot.
#'
#' @param coord coordinate representation of points
#' @param x,y variables names (as string) to map to x and y axis
#' @param wc parameter values as matrix
#' @param obs observable to plot
#' @param cond row numbers of points used for conditioning
#' @return ggplot
#' @export
plotObs <- function(coord, x, y, wc, obs, cond){
  dat <- coord[cond,]
  colnames(dat) <- paste0("O",1:ncol(coord))
  dat %>%
    tourr::rescale() %>%
    tibble::as_tibble() %>%
    cbind(wc[cond,]) %>%
    ggplot2::ggplot(ggplot2::aes_string(x, y, fill=obs)) +
    ggplot2::geom_tile() +
    ggplot2::guides(fill = FALSE) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(paste0("Centered coordinate values for ", obs)) +
    ggplot2::theme(aspect.ratio = 1)
}

#' Plot chi2
#'
#' Parameter values need to be on a regular grid for correct appearance of
#' this plot.
#'
#' @param wc parameter values as matrix
#' @param chi2 vector with chi2 values
#' @param x,y variables names (as string) to map to x and y axis
#' @param cond row numbers of points used for conditioning
#' @return ggplot
#' @export
plotChi2 <- function(wc, chi2, x, y, cond){
  dplyr::mutate(wc[cond,], chi2 = chi2[cond]) %>%
    ggplot2::ggplot(ggplot2::aes_string(x, y, fill="chi2")) +
    ggplot2::geom_tile() +
    ggplot2::guides(fill = FALSE) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Chi2 values") +
    ggplot2::theme(aspect.ratio = 1)
}

#' Plot sigma bins in parameter space
#'
#' Parameter values need to be on a regular grid for correct appearance of
#' this plot.
#'
#' @param wc parameter values as matrix
#' @param sm,bf index values for the SM and BF point
#' @param bmID index values for the benchmark points
#' @param sigmabins binning in sigma
#' @param x,y variables names (as string) to map to x and y axis
#' @param cond row numbers of points used for conditioning
#' @return ggplot
#' @export
plotSigBin <- function(wc, sm, bf, bmID, sigmabins, x, y, cond){
  palSig <- RColorBrewer::brewer.pal(length(sigmabins), "Set2")
  colSig <- palSig[sigmabins]
  ggplot2::ggplot(wc[cond,], ggplot2::aes_string(x, y)) +
    ggplot2::geom_tile(fill= colSig[cond]) +
    ggplot2::geom_point(data=wc[sm,],  shape=1, size=3) +
    ggplot2::geom_point(data=wc[bf,], shape=8, size=3) +
    ggplot2::geom_point(data = wc[bmID,], shape=5, size=3) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Sigma bins in parameter space") +
    ggplot2::theme(aspect.ratio = 1, legend.position = "none")
}

#' Make GIF of tour animation.
#'
#' Renders the tour animation, stored in file tour_animation.gif
#'
#' @param coord coordinate representation of the points
#' @param col color vector according to group assignment
#' @param pch vector of plotting symbols
#' @export
tourGif <- function(coord, col, pch){
  set.seed(2021)
  colnames(coord) <- paste0("O",1:ncol(coord))
  tourr::render_gif(coord, tourr::grand_tour(),
                tourr::display_xy(col = col, pch = pch),
                gif_file = "tour_animation.gif",
                frames = 100, rescale = FALSE)
}

#' Make plot as specified in settings
#'
#' Used as interface to generate a specific graph seen when using the GUI.
#' Settings include: coord, useCov, metric, linkage, k, plotType
#'
#' @param pred prediction matrix
#' @param covInv inverse covariance matrix
#' @param wc matrix specifying the model parameters (on a grid)
#' @param exp observable reference value (e.g. experimental measurement)
#' @param settings list specifying parameters usually selected in the app
#' @param user_coord input coordinate matrix (optional)
#' @param user_dist input distance matrix (optional)
#' @param c specification for conditioning
#' @return ggplot
#' @export
makePlots <- function(pred, covInv, wc, exp, settings,
                      user_coord=NULL, user_dist=NULL, c=NULL){
  n <- nrow(pred)
  chi2 <- computeChi2(pred, covInv, exp)
  sig <- computeSigma(chi2, 2)
  bf <- which.min(chi2)
  sm <- which.min(rowSums(abs(wc)))
  if(!is.null(c)){
    x <- c$x
    y <- c$y
    cond <- which(wc[[c$pz1]]==as.numeric(c$vz1))
    if(! is.null(c$pz2)){
      cond2 <- which(wc[[c$pz2]]==as.numeric(c$vz2))
      cond <- intersect(cond, cond2)
    }
  }
  else{
    x <- colnames(wc)[1]
    y <- colnames(wc)[2]
    cond <- 1:nrow(wc)
  }
  coord <- getCoords(settings$coord, settings$useCov, pred, covInv, exp, user_coord)
  dists <- getDists(coord, settings$metric, user_dist)
  fit <- stats::hclust(dists, settings$linkage)
  groups <- stats::cutree(fit, k=settings$k)
  lvl <- unique(groups[stats::order.dendrogram(stats::as.dendrogram(fit))])
  groups <- as.numeric(factor(groups, levels= lvl))
  pal <- RColorBrewer::brewer.pal(settings$k, "Dark2")
  col <- pal[groups]
  pch <- rep(20, n)
  pch[sm] <- 8
  pch[bf] <- 15
  sigmabins <- chi2bins(chi2, 2, settings$k)
  benchmarks <- getBenchmarkInformation(as.matrix(dists), groups)
  if(settings$plotType == "PC"){ return( plotPC(coord, groups, benchmarks$id))}
  else if(settings$plotType == "PCscaled"){ return(plotPC(coord, groups, benchmarks$id, TRUE))}
  else if(settings$plotType == "WC") {return(plotWC(wc, x, y, sm, bf, benchmarks$id, col, cond))}
  else if(settings$plotType == "chi2") {return(plotChi2(wc, chi2, x, y, cond))}
  else if(settings$plotType == "sigBins"){return(plotSigBin(wc, sm, bf, benchmarks$id,
                                                            sigmabins, x, y, cond))}
  else if(settings$plotType == "heatmap"){return(plotHeatmap(as.matrix(dists), fit, settings$k, pal))}
  else if(settings$plotType %in% names(cstat_names)){return(plotCstat(dists, fit, chi2,
                                                                      settings$plotType))}
  else if(startsWith(settings$plotType, "O")){return(plotObs(coord, x, y, wc, settings$plotType, cond))}
  else if(settings$plotType== "tour"){return(tourGif(coord, col, pch))}

  "plotType unknown"
}
