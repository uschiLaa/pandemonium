# dat is coordinate representation in observable space
# gr is grouping returned from clustering
# benchmarkIds gives index values of benchmarks
# if s=TRUE rescale each variable separately
# a is alpha level to draw non-benchmark points with
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

# wc is the full WC data
# x, y are the variables to plot (as string), in the app these are the first two columns of wc
# sm, bf are the index values for the (approximate) SM and BF point
# benchmarkIds is vector of index values for the benchmark points
# col is vector of colors according to cluster assignment
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

plotHeatmap <- function(dat, fit, k, pal){
  dendo <- as.dendrogram(fit) %>%
    dendextend::set("branches_lwd", 3) %>%
    dendextend::color_branches(k = k, col = pal)

  heatmap(dat, Rowv = dendo, Colv = rev(dendo), scale = "none")
}

plotCstat <- function(dist, fit, chivals, stat, kmax=8){
  cstats <- getClusterStats(dist, fit, chivals, kmax)
  ggplot2::ggplot(cstats, ggplot2::aes_string("k", stat)) +
    ggplot2::geom_line() +
    ggplot2::xlab("# clusters") +
    ggplot2::ylab(cstat_names[[stat]]) +
    ggplot2::theme_bw()
}

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

tourGif <- function(coord, col, pch){
  set.seed(2021)
  colnames(coord) <- paste0("O",1:ncol(coord))
  tourr::render_gif(coord, tourr::grand_tour(),
                tourr::display_xy(col = col, pch = pch),
                gif_file = "tour_animation.gif",
                frames = 100, rescale = FALSE)
}

# settings is list with: coord (Pull, pVal), useCov (T/F),
# metric (euclidean, manhatten, ...), linkage (single, complete, ...),
# k (number of clusters), plotType (PC, PCscaled, WC)
makePlots <- function(pred, covInv, wc, exp, settings, user_coord=NULL, c=NULL){
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
  fit <- hclust(dists, settings$linkage)
  groups <- cutree(fit, k=settings$k)
  lvl <- unique(groups[order.dendrogram(as.dendrogram(fit))])
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
