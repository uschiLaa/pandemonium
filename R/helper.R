getCoords <- function(coord, useCov, pred, covInv, exp, user_coord){

  if(coord=="User") return(as.matrix(user_coord))

  n <- nrow(pred)
  pred <- as.matrix(pred)
  nc <- ncol(pred)
  coord_mat <- matrix(nrow = n, ncol = nc)

  if (coord=="normalised"){
    for (i in 1:n){
      for (j in 1:nc){
        if(useCov) coord_mat[i, j] <- sum(covInv[j,] * pred[i,]) / sqrt(covInv[j, j])
        else coord_mat[i, j] <- as.numeric(pred[i, j]) / sqrt(solve(covInv)[j, j])
      }
    }
    return(coord_mat)
  }
  if (coord=="Pull"){
    for (i in 1:n){
      for (j in 1:nc){
        if(useCov) coord_mat[i, j] <- sum(covInv[j,] * (pred[i,] - exp$value)) / sqrt(covInv[j, j])
        else coord_mat[i, j] <- as.numeric((pred[i, j] - exp$value[j]) / sqrt(solve(covInv)[j, j]))
      }
    }
    return(coord_mat)
  }
  if (coord=="p-val"){
    for (i in 1:n){
      for (j in 1:nc){
        if(useCov) coord_mat[i, j] <- pchisq((sum(covInv[j,] * (pred[i,] - exp$value)) / sqrt(covInv[j, j]))^2, df=1)
        else coord_mat[i, j] <- pchisq((as.numeric((pred[i, j] - exp$value[j]) / sqrt(solve(covInv)[j, j])))^2, df=1)
      }
    }
    return(coord_mat)
  }
}

getDists <- function(coord, metric, user_dist){
  if(metric == 'euclidean2') dists <- dist(coord)^2
  else if(metric == "user") dists <- as.dist(user_dist)
  else dists <- dist(coord, method = metric)
  return(dists)
}

computeChi2 <- function(pred, covInv, exp){
  chi2 <- double(nrow(pred))
  for (i in 1:nrow(pred)){
    chi2[i] <- as.matrix(exp$value - pred[i,]) %*% covInv %*% t(as.matrix(exp$value - pred[i,]))
  }
  return(chi2)
}

getBenchmarkInformation <- function(dmat, groups){
  k <- length(unique(groups))
  ret <- tibble::tibble(id = numeric(length = k), r = numeric(length = k), d = numeric(length = k),
                        group = numeric(length = k))
  i <- 1
  for (gr in unique(groups)){
    idx <- which(groups==gr)
    d_gr <- dmat[idx, idx]
    if(length(idx) == 1) d_vec <- 0
    else d_vec <- colSums(d_gr^2)
    id <- idx[which.min(d_vec)]
    d <- max(d_gr)
    r <- max(d_gr[which.min(d_vec)])
    ret[i,] <- t(c(id, r, d, gr))
    i <- i+1
  }
  ret
}

getClusterDists <- function(dmat, groups, benchmarks){
  k <- length(unique(groups))
  n <- choose(k, 2)
  ret <- tibble::tibble(grA = numeric(length = n),  grB = numeric(length = n),
                        d1 = numeric(length = n), d2 = numeric(length = n),
                        d3 = numeric(length = n))

  ni <- 1
  for (i in 1:(k-1)){
    for (j in (i+1):k){

      id1 <- dplyr::filter(benchmarks, group == i)$id
      id2 <- dplyr::filter(benchmarks, group == j)$id

      d1 <- dmat[id1, id2]

      idx1 <- which(groups==i)
      idx2 <- which(groups==j)


      d2 <- min(dmat[idx1, idx2])
      d3 <- max(dmat[idx1, idx2])


      ret[ni,] <- t(c(i, j, d1, d2, d3))
      ni <- ni + 1
    }
  }
  ret
}

getClusterStats <- function(dist, fit, chivals, kmax=10){
  ret <- tibble::tibble(k = numeric(length = kmax-1),
                        wb.ratio = numeric(length = kmax-1),
                        ch = numeric(length = kmax-1),
                        pearsongamma = numeric(length = kmax-1),
                        dunn = numeric(length = kmax-1),
                        dchi2rand = numeric(length = kmax-1),
                        rmin = numeric(length = kmax-1),
                        rmax = numeric(length = kmax-1),
                        dmax = numeric(length = kmax-1),
                        dmin = numeric(length = kmax-1))
  for(k in 2:kmax){
    gr <- cutree(fit, k)
    chibins <- chi2bins(chivals, 2, k)

    x <- fpc::cluster.stats(dist, gr, alt.clustering = chibins)
    bmInfo <- getBenchmarkInformation(as.matrix(dist), gr)
    bmDists <- getClusterDists(as.matrix(dist), gr, bmInfo)
    bmMinDist <- min(bmDists$d1)

    ret[k-1,] <- t(c(k, x$wb.ratio,
                   x$ch, x$pearsongamma, x$dunn,
                   x$corrected.rand,
                   min(bmInfo$r), max(bmInfo$r), max(bmInfo$d),
                   bmMinDist))
  }
  ret
}

# arguments are list of chi2 values (chivals),
# number of parameters (ndf) and number of clusters (k)
# returns list of integer bin codes from lowest to highest values of sigma
chi2bins <- function(chivals, ndf, k){
  chimin <- min(chivals)
  # map chivals to sigmas
  sigvals <- sqrt(qchisq(pchisq(chivals-chimin, ndf), 1))
  sigvals <- pmin(sigvals, 5)
  # get bins in sigma
  sigmabins <- seq(0, max(sigvals), length.out = k+1)
  sigbinned <- cut(sigvals, sigmabins,
                   include.lowest = TRUE,
                   labels = FALSE)
  sigbinned
}

computeSigma <- function(chivals, ndf){
  chimin <- min(chivals)
  # map chivals to sigmas, cutoff at 5
  pmin(sqrt(qchisq(pchisq(chivals-chimin, ndf), 1)), 5)
}

cstat_names <- list(
  "wb.ratio" = "WB ratio",
  "pearsongamma" = "Normalized gamma",
  "dunn" = "Dunn index",
  "ch" = "Calinski and Harabasz index",
  "rmin" = "Minimum radius",
  "rmax" = "Maximum radius",
  "dmax" = "Maximum diameter",
  "dmin" = "Minimum benchmark distance",
  "dchi2rand" = "ARI with CI binning"
)

cstat_labeller <- function(variable,value){
  return(cstat_names[value])
}
