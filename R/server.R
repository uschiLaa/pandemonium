#' Shiny app for exploring clustering solutions
#'
#' Opening the GUI to cluster the data points based on the predictions.
#' Coordinates and distances are computed on the fly, or can be entered
#' in the function call.
#'
#' @param pred prediction matrix
#' @param covInv inverse covariance matrix
#' @param wc matrix specifying the model parameters (on a grid)
#' @param exp observable reference value (e.g. experimental measurement)
#' @param user_coord input coordinate matrix (optional)
#' @param user_dist input distance matrix (optional)
#' @export
#' @examples \dontrun{
#' pandemonium(b_anomaly$pred, b_anomaly$covInv, b_anomaly$wc, b_anomaly$exp)
#' }
pandemonium <- function(pred, covInv, wc, exp, user_coord = NULL, user_dist = NULL){

  rv <- shiny::reactiveValues()
  n <- nrow(pred)
  ndf <- ncol(wc)
  chi2 <- computeChi2(pred, covInv, exp)
  sig <- computeSigma(chi2, ndf)
  bf <- which.min(chi2)
  sm <- which.min(rowSums(abs(wc)))
  rv$x <- colnames(wc)[1]
  rv$y <- colnames(wc)[2]
  rv$vz1 <- NULL
  rv$vz2 <- NULL
  rv$cond <- 1:nrow(wc)
  tmp <- tempdir()

  server <- function(input, output, session) {

    if(!is.null(user_coord)){
      shiny::updateSelectInput(session, "coord",
                               choices = c("Pull", "p-val", "User"))
      shiny::updateSelectInput(session, "coordA",
                               choices = c("Pull", "p-val", "User"))
      shiny::updateSelectInput(session, "coordB",
                               choices = c("Pull", "p-val", "User"))
    }

    if(!is.null(user_dist)){
      shiny::updateSelectInput(session, "metric",
                               choices = c("euclidean", "maximum", "manhattan", "canberra",
                                           "binary", "minkowski", "user"))
      shiny::updateSelectInput(session, "metricA",
                               choices = c("euclidean", "maximum", "manhattan", "canberra",
                                           "binary", "minkowski", "user"))
      shiny::updateSelectInput(session, "metricB",
                               choices = c("euclidean", "maximum", "manhattan", "canberra",
                                           "binary", "minkowski", "user"))
    }

    shiny::observeEvent(c(input$metric, input$coord, input$useCov, input$linkage), {
      n_points <- nrow(pred)
      rv$coord <- getCoords(input$coord, input$useCov, pred, covInv, exp, user_coord)
      dists <- getDists(rv$coord, input$metric, user_dist)
      rv$d_mat <- as.matrix(dists)
      rv$fit <- stats::hclust(dists, input$linkage)
      rv$cstats <- getClusterStats(stats::as.dist(rv$d_mat), rv$fit, chi2, 8)
    }, priority = 99
    )

    shiny::observeEvent(c(input$kC, rv$fit), {
      rv$kC <- as.numeric(input$kC)
      groups <- stats::cutree(rv$fit, k=rv$kC)
      rv$lvl <- unique(groups[stats::order.dendrogram(stats::as.dendrogram(rv$fit))])
      rv$groups <- as.numeric(factor(groups, levels= rv$lvl))
      rv$pal <- RColorBrewer::brewer.pal(rv$kC, "Dark2")
      rv$col <- rv$pal[rv$groups]
      rv$sigmabins <- chi2bins(chi2, ndf, rv$kC)
      rv$benchmarks <- getBenchmarkInformation(rv$d_mat, rv$groups)
      rv$bpDists <- getClusterDists(rv$d_mat, rv$groups, rv$benchmarks)
      rv$pointcol <- rv$col
      rv$pointcol[c(sm,bf)] <- "black"
      rv$pch <- rep(20, n)
      rv$pch[sm] <- 8
      rv$pch[bf] <- 15
      rv$alpha <- rep(0.6, n)
      rv$alpha[c(sm,bf)] <- 1
    }, priority = 80)

    shiny::observeEvent(c(input$px, input$py),{
      rv$x <- input$px
      rv$y <- input$py
      add_params <- colnames(wc)[!colnames(wc) %in% c(rv$x, rv$y)]
      if(length(add_params) > 0) {
        rv$pz1 <- add_params[1]
        output$conditions1 <- shiny::renderUI(
          shiny::selectInput("vz1", label = add_params[1],
                             choices = unique(wc[[add_params[1]]]))
        )
      }
      if(length(add_params) > 1) {
        rv$pz2 <- add_params[2]
        output$conditions2 <- shiny::renderUI(
          shiny::selectInput("vz2", label = add_params[2],
                             choices = unique(wc[[add_params[2]]]))
        )
      }
      else rv$pz2 <- NULL
    })


    shiny::observeEvent(c(input$vz1, input$vz2), {
      rv$cond1 <- which(wc[[rv$pz1]]==as.numeric(input$vz1))
      if(! is.null(rv$pz2)){
        rv$cond2 <- which(wc[[rv$pz2]]==as.numeric(input$vz2))
        rv$cond <- intersect(rv$cond1, rv$cond2)
      }
      else rv$cond <- rv$cond1
    })

    output$heatmap <- shiny::renderPlot({
      plotHeatmap(rv$d_mat, rv$fit, rv$kC, rv$pal)
      }, height = 325)

    output$chi2 <- shiny::renderPlot({
      plotChi2(wc, chi2, rv$x, rv$y, rv$cond)
    }, height = 325)

    output$clusterstats <- shiny::renderPlot({
      rv$cstats %>%
        tidyr::pivot_longer(cols=within.cluster.ss:dmin, names_to="stat") %>%
        ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x=k, y=value)) +
        ggplot2::xlab("# clusters") +
        ggplot2::ylab("") +
        ggplot2::facet_wrap(~stat, ncol=3, scales = "free_y", labeller=cstat_labeller()) +
        ggplot2::theme_bw()
    }, height = 600)

    output$scalar <- shiny::renderPlot({
      dat <- rv$coord[rv$cond,]
      colnames(dat) <- paste0("O",1:ncol(pred))
      dat %>%
        tourr::rescale() %>%
        tibble::as_tibble() %>%
        cbind(wc[rv$cond,]) %>%
        tidyr::pivot_longer(cols = tidyselect::starts_with("O"), names_to="observable") %>%
        dplyr::mutate(observable = factor(observable, levels = paste0("O",1:ncol(pred)))) %>%
        ggplot2::ggplot(ggplot2::aes_string(rv$x, rv$y, fill="value")) +
        ggplot2::geom_tile() +
        ggplot2::facet_wrap(~observable, scales = "free", ncol=7) +
        ggplot2::guides(fill = FALSE) +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::theme_bw() +
        ggplot2::ggtitle("Centered coordinate values for all observables") +
        ggplot2::theme(aspect.ratio = 1)
    })

    output$pc <- shiny::renderPlot({
      plotPC(rv$coord, rv$groups, rv$benchmarks$id)
    }, height = 300)

    output$pc2 <- shiny::renderPlot({
      plotPC(rv$coord, rv$groups, rv$benchmarks$id, TRUE)
    }, height = 300)

    output$wc <- shiny::renderPlot({
      plotWC(wc, rv$x, rv$y, sm, bf, rv$benchmarks$id, rv$col, rv$cond)
      }, height = 325)


    output$sigbins <- shiny::renderPlot({
      plotSigBin(wc, sm, bf, rv$benchmarks$id, rv$sigmabins,
                 rv$x, rv$y, rv$cond)

    }, height = 325)

    output$distText <- shiny::renderText(paste0("Average distance: ",
                                         round(mean(as.vector(rv$d_mat)),1),
                                         ", Maximum distance: ",
                                         round(max(as.vector(rv$d_mat)), 1)))


    output$hist <- shiny::renderPlot({

      # to get vector of distances without double-counting we
      # transform to distance keeping only entries below the diagonal
      dist_vec <- as.vector(rv$d_mat)
      gr1 <- double(length(dist_vec))
      gr2 <- double(length(dist_vec))
      diff_m <- matrix(nrow = length(dist_vec), ncol = ncol(pred))
      ctr <- 1
      #FIXME faster implementation of this?
      for(i in 1:n){
        for(j in 1:n){
          gr1[ctr] <- rv$groups[i]
          gr2[ctr] <- rv$groups[j]
          diff_m[ctr,] <- rv$coord[i,] - rv$coord[j,]
          ctr <- ctr + 1
        }
      }

      colnames(diff_m) <- paste0("diff_",1:ncol(pred))

      dist_tib <- tibble::as_tibble(diff_m) %>%
        dplyr::mutate(dist = dist_vec, gr1 = as.factor(gr1), gr2 = as.factor(gr2)) %>%
        dplyr::mutate(match = dplyr::if_else(gr1==gr2, "within", "between"))

      p1 <- ggplot2::ggplot(dist_tib,
                            ggplot2::aes(x = dist, y = stat(count / sum(count)))) +
        ggplot2::geom_histogram(data=dplyr::select(dist_tib, -gr1, -gr2, -match),
                                fill=NA, color="grey", position="identity") +
        ggplot2::geom_histogram(mapping = ggplot2::aes(color= gr1),
                                fill=NA, position="identity") +
        ggplot2::scale_color_brewer(palette="Dark2") +
        ggplot2::facet_grid(gr1~match) +
        ggplot2::theme_bw() +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::ggtitle(
          "Distribution of distances within and between clusters"
          ) +
        ggplot2::theme(legend.position = "none",
                       strip.text.y = ggplot2::element_blank())

      dist_long <- dist_tib %>%
        dplyr::select(-dist) %>%
        tidyr::pivot_longer(cols = tidyselect::starts_with("diff"),
                            names_to="observable", values_to="difference",
                            names_prefix = "diff_")
      dist_long$observable <- factor(as.integer(dist_long$observable),
                                     levels=1:ncol(pred))

      p2 <-  ggplot2::ggplot(dist_long,
                             ggplot2::aes( x = gr1, y = difference,
                                           color = gr1, fill = match)) +
        ggplot2::geom_boxplot() +
        ggplot2::scale_color_brewer(palette="Dark2") +
        ggplot2::scale_fill_manual(values = c("grey", "white")) +
        ggplot2::facet_wrap(.~observable, ncol = 3) +
        ggplot2::theme_bw() +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::ggtitle("Distribution of coordinate differences between (left)
                         and within (right) clusters") +
        ggplot2::theme(legend.position = "none",
                       strip.text.y = ggplot2::element_blank())

     gridExtra::grid.arrange(p1, p2, ncol=2)


    },height = 600, width = 1000)

    output$benchmarks = DT::renderDT({

      dt <- wc[rv$benchmarks$id,] %>%
        tibble::add_column(cluster = rv$benchmarks$group) %>%
        tibble::add_column(r = rv$benchmarks$r) %>%
        tibble::add_column(d = rv$benchmarks$d) %>%
        tibble::add_column(sigma = sig[rv$benchmarks$id])

      DT::formatRound(DT::datatable(
        dt[,c(ncol(wc)+1, 1:ncol(wc), (ncol(wc)+2):ncol(dt))],
        rownames = FALSE),
        2:ncol(dt), digits=2) %>%
        DT::formatStyle(
          "cluster",
          backgroundColor = DT::styleEqual(rv$benchmarks$group, unique(rv$col)))
    }
    )


    output$tsne <- shiny::renderPlot({
      dist_tsne <- Rtsne::Rtsne(rv$d_mat, is_distance=TRUE)$Y
      colnames(dist_tsne) <- c("tsne1", "tsne2")
      ggplot2::ggplot(tibble::as_tibble(dist_tsne),
                      ggplot2::aes(tsne1, tsne2)) +
        ggplot2::geom_point(size=3, color= rv$pointcol,
                   shape=rv$pch, alpha=rv$alpha) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle("t-SNE embedding") +
        ggplot2::theme(aspect.ratio = 1, legend.position = "none")
    })

    output$umap <- shiny::renderPlot({
      dist_umap <- as.matrix(uwot::umap(stats::as.dist(rv$d_mat)))
      colnames(dist_umap) <- c("umap1", "umap2")
      ggplot2::ggplot(tibble::as_tibble(dist_umap),
                      ggplot2::aes(umap1, umap2)) +
        ggplot2::geom_point(size=3, color= rv$pointcol,
                   shape=rv$pch, alpha=rv$alpha) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle("UMAP embedding") +
        ggplot2::theme(aspect.ratio = 1, legend.position = "none")
    })


    output$lle <- shiny::renderPlot({
      dist_lle <- lle::lle(rv$coord, 2, 10)$Y
      colnames(dist_lle) <- c("lle1", "lle2")
      ggplot2::ggplot(tibble::as_tibble(dist_lle),
                      ggplot2::aes(lle1, lle2)) +
        ggplot2::geom_point(size=3, color= rv$pointcol,
                   shape=rv$pch, alpha=rv$alpha) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle("Local linear embedding") +
        ggplot2::theme(aspect.ratio = 1, legend.position = "none")
    })

    shiny::observeEvent(rv$groups,{
      tour_data <- tour_coord(rv$coord, rv$pointcol, rv$pch)
      tourr::render(tour_data$coord, tourr::grand_tour(),
                    tourr::display_xy(col = tour_data$col, pch = tour_data$pch),
                    'png', paste0(tmp, "/tour-%03d.png"),
                    frames = 100, rescale = F)
      dist_pca <- stats::prcomp(rv$coord)$x[,1:min(5, ncol(pred))]
      colnames(dist_pca) <- paste0("pc",1:min(5, ncol(pred)))
      tourr::render(dist_pca, tourr::grand_tour(),
                    tourr::display_xy(col = rv$pointcol, pch = rv$pch),
                    'png', paste0(tmp, "/tour-pca-%03d.png"),
                    frames = 100, rescale = F)
    })


    output$tourImg <- shiny::renderImage({

      filename <- sprintf(paste0(tmp, '/tour-%03d.png'), input$n)
      # Return a list containing the filename and alt text
      list(src = filename,
           width = 400,
           height = 400,
           alt = paste("Parameter space", input$n))},
      deleteFile = FALSE)

    output$tourImgPca <- shiny::renderImage({

      filename <- sprintf(paste0(tmp, '/tour-pca-%03d.png'), input$npca)
      # Return a list containing the filename and alt text
      list(src = filename,
           width = 400,
           height = 400,
           alt = paste("Parameter space", input$n))},
      deleteFile = FALSE)


    shiny::observeEvent(c(input$metricA, input$linkageA, input$linkageB, input$metricB,
                          input$coordA, input$coordB, input$kA, input$kB, input$useCovA,
                          input$useCovB), {

      kA <- as.numeric(input$kA)
      kB <- as.numeric(input$kB)


      coordA <- getCoords(input$coordA, input$useCovA, pred, covInv, exp, user_coord)
      coordB <- getCoords(input$coordB, input$useCovB, pred, covInv, exp, user_coord)

      distsA <- getDists(coordA, input$metricA, user_dist)
      d_matA <- as.matrix(distsA)

      distsB <- getDists(coordB, input$metricB, user_dist)
      d_matB <- as.matrix(distsB)

      fitA <- stats::hclust(distsA, input$linkageA)
      fitB <- stats::hclust(distsB, input$linkageB)

      groupsA <- stats::cutree(fitA, kA)
      groupsB <- stats::cutree(fitB, kB)

      palA <- RColorBrewer::brewer.pal(kA, "Dark2")
      groupsA <- as.numeric(factor(groupsA,
                                   levels=unique(groupsA[stats::order.dendrogram(stats::as.dendrogram(fitA))])))
      colA <- palA[groupsA]

      palB <- RColorBrewer::brewer.pal(kB, "Set2")
      groupsB <- as.numeric(factor(groupsB,
                                   levels=unique(groupsB[stats::order.dendrogram(stats::as.dendrogram(fitB))])))
      colB <- palB[groupsB]

      output$tableAB <- shiny::renderPlot({
        table(groupsA, groupsB) %>%
          as.data.frame() %>%
          ggplot2::ggplot(ggplot2::aes(groupsA, groupsB)) +
          ggplot2::geom_tile(ggplot2::aes(fill = Freq)) +
          ggplot2::geom_text(ggplot2::aes(label=Freq)) +
          ochRe::scale_fill_ochre(palette="galah", discrete=FALSE) +
          ggplot2::theme_minimal() +
          ggplot2::xlab("") +
          ggplot2::ylab("") +
          ggplot2::theme(legend.position = "none",
                axis.text.x =
                  ggplot2::element_text(colour = palA, face = "bold", size=15),
                axis.text.y =
                  ggplot2::element_text(colour = palB, face = "bold", size=15))
        }, height = 325)


      output$wcA <- shiny::renderPlot({
        plotWC(wc, rv$x, rv$y, sm, bf, c(), colA, rv$cond)
      }, height = 325)

      output$wcB <- shiny::renderPlot({
        plotWC(wc, rv$x, rv$y, sm, bf, c(), colB, rv$cond)
      }, height = 325)

      output$heatmapA <- shiny::renderPlot({
        plotHeatmap(d_matA, fitA, kA, palA)
      }, height = 325)
      output$heatmapB <- shiny::renderPlot({
        plotHeatmap(d_matB, fitB, kB, palB)
      }, height = 325)

    }
    )

    shiny::observeEvent(rv$groups,{
      detour_data <- as.data.frame(pred)
      detour_data$colour <- rv$groups
      rv$shared_detour <- crosstalk::SharedData$new(detour_data)
      rv$depal <- if (length(unique(rv$groups)) == 2) rv$pal[-3] else rv$pal
    })
    
    output$detour1   <- detourr::shinyRenderDisplayScatter2d(detourr::detour(rv$shared_detour, 
                                                                             detourr::tour_aes(projection = colnames(pred),colour = colour)) |>
                                                               detourr::tour_path(detourr::grand_tour(2), fps = 60, 
                                                                                  max_bases=20) |>
                                                               detourr::show_scatter(alpha = 0.7, 
                                                                                     axes = TRUE, scale_factor = 0.4, palette = rv$depal), quoted = FALSE)
    output$detour2 <- detourr::shinyRenderDisplayScatter2d(detourr::detour(rv$shared_detour, 
                                                                           detourr::tour_aes(projection = colnames(pred),colour = colour)) |>
                                                             detourr::tour_path(detourr::local_tour(tourr::basis_random(n=ncol(pred))), fps = 60, 
                                                                                max_bases=20) |>
                                                             detourr::show_scatter(alpha = 0.7, 
                                                                                   axes = TRUE, scale_factor = 0.4, palette = rv$depal), quoted = FALSE)
    

  }

  shiny::shinyApp(ui(colnames(wc)), server)

}

