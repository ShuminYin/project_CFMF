#' Run the complete CFMF workflow
#'
#' This is the main function. It calculates gene percentages, runs purity analysis,
#' and then calculates silhouette scores.
#'
#' @param seu A Seurat object (must have PCA run).
#' @param neighbor Number of neighbors for KNN (default 50).
#' @param ncores Number of cores for parallel processing (default 10).
#' @return A data frame containing gene, purity stats, and silhouette scores.
#' @importFrom dplyr filter mutate
#' @importFrom Seurat GetAssayData
#' @importFrom Matrix rowMeans rowSums
#' @export
RunCFMF <- function(seu, neighbor = 50, ncores = 10) {

  # 1. check PCA
  if (is.null(seu@reductions$pca)) {
    stop("The Seurat object must have PCA calculated (RunPCA) before running CFMF.")
  }

  message("Step 1: Calculating gene percentages...")


  mat <- Seurat::GetAssayData(seu, slot = "data")
  cutoff <- 0


  gene.pct <- data.frame(
    gene = rownames(seu),
    pct = Matrix::rowMeans(mat > cutoff),
    count = Matrix::rowSums(mat > cutoff)
  )


  gene.pct <- gene.pct |> dplyr::filter(count > 3, pct < 1)

  message(paste("      Found", nrow(gene.pct), "genes passing criteria."))

  # 2. run Purity
  message("Step 2: Calculating Purity (this may take a while)...")
  purity_res <- calculate.purity(seu, gene.pct, neighbor, ncores)

  # 3. run Silhouette
  message("Step 3: Calculating Silhouette scores...")
  final_res <- calculate.silhouette(seu, purity_res, ncores)

  message("Done!")
  return(final_res)
}


#' Calculate Gene Purity (Internal)
#'
#' @param seu Seurat object
#' @param gene.pct Data frame with gene filters
#' @param neighbor K for KNN
#' @param ncores Cores
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom BiocNeighbors findKNN
#' @importFrom dplyr mutate
#' @importFrom stats fisher.test pchisq p.adjust
#' @export
calculate.purity <- function(seu, gene.pct, neighbor, ncores = 10) {

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)


  data_mat <- as.matrix(seu@reductions$pca@cell.embeddings)


  knn_res <- BiocNeighbors::findKNN(data_mat, k = neighbor, get.distance = TRUE)

  cutoff <- 0
  genes <- gene.pct$gene

  gene.purity <- foreach::foreach(i = genes, .combine = "rbind", .packages = c('Seurat', 'doParallel', 'stats')) %dopar% {


    gene.exp <- Seurat::GetAssayData(seu, slot = "data")[i, ]
    clusters <- gene.exp > cutoff

    cells <- 1:ncol(seu)
    pos.cells <- which(clusters)


    df <- foreach::foreach(j = pos.cells, .combine = 'rbind') %do% {
      other.cells <- cells[-j]
      is.neighbor <- other.cells %in% knn_res$index[j, ]
      is.pos      <- other.cells %in% pos.cells

      fit   <- stats::fisher.test(x = is.neighbor, y = is.pos)
      p.val <- fit$p.value

      data.frame(cell = j, p.val, row.names = NULL)
    }

    if (is.null(df) || nrow(df) == 0) {
      data.frame(gene = i, meta.p = 1)
    } else {
      chisq_stat <- -2 * sum(log(df$p.val + 1e-300))
      v          <- 2 * length(df$p.val)
      meta.p     <- stats::pchisq(chisq_stat, v, lower.tail = FALSE)
      data.frame(gene = i, meta.p)
    }
  }

  parallel::stopCluster(cl)


  gene.purity <- gene.purity |> dplyr::mutate(FDR = stats::p.adjust(meta.p, method = "BH"))
  gene.purity <- merge(gene.purity, gene.pct, by = "gene")

  return(gene.purity)
}


#' Calculate Silhouette Scores (Internal)
#'
#' @param seu Seurat object
#' @param purity Output from calculate.purity
#' @param ncores Cores
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom cluster silhouette
#' @importFrom stats dist
#' @importFrom dplyr filter
#' @importFrom Seurat GetAssayData
#' @export
calculate.silhouette <- function(seu, purity, ncores = 10) {

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)


  data_dist <- stats::dist(as.matrix(seu@reductions$pca@cell.embeddings))

  genes <- purity$gene

  df <- foreach::foreach(i = genes, .combine = "rbind", .packages = c("Seurat", "cluster", "dplyr")) %dopar% {


    expr     <- Seurat::GetAssayData(seu, slot = "data")[i, ]
    clusters <- ifelse(expr > 0, 1, 2)

    if (length(unique(clusters)) < 2) {
      data.frame(gene = i, silhouette = -1)
    } else {
      sil_scores <- cluster::silhouette(x = as.numeric(clusters), dist = data_dist)
      sil_df <- as.data.frame(sil_scores)

      sil_pos <- sil_df |> dplyr::filter(cluster == 1)

      pos.silhouette <- median(sil_pos$sil_width)

      data.frame(gene = i, silhouette = pos.silhouette)
    }
  }

  parallel::stopCluster(cl)

  silhouette.df <- merge(df, purity, by = "gene")
  return(silhouette.df)
}
