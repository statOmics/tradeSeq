#' @title Extract info from Monocle models
#'
#' @description This function extracts info that will be used downstream to make
#' \code{CellDataSet} objects compatible with a \code{tradeSeq}
#' analysis
#'
#' @rdname extract_monocle_info
#' @param cds A \code{CellDataSet} object.
#' @details For now, this only works for the DDRTree dimentionality reduction.
#' It is the one recommanded by the Monocle developers.
#' @return
#' A list with four objects. A \code{pseudotime} matrix and a \code{cellWeights}
#' matrix that can be used as input to \code{\link{fitGAM}} or
#' \code{\link{evaluateK}}, the reduced dimension matrix for the cells, and a
#' list of length the number of lineages, containing the reduced dimension of
#' each lineage.
#' @importFrom magrittr %>%
#' @importFrom Biobase exprs
#' @import monocle
#' @importFrom igraph degree shortest_paths
#' @export
extract_monocle_info <- function(cds) {
  if (cds@dim_reduce_type != "DDRTree") {
    stop(paste0("For now tradeSeq only support Monocle with DDRTree",
                "reduction. If you want to use another type",
                "please use another format for tradeSeq inputs."))
  }
  # Get the reduced dimension of DDRT
  rd <- t(monocle::reducedDimS(cds)) %>% as.data.frame()
  
  # Get the various lineages info for weights and pseudotime
  y_to_cells <- cds@auxOrderingData[["DDRTree"]]
  y_to_cells <- y_to_cells$pr_graph_cell_proj_closest_vertex %>%
    as.data.frame()
  y_to_cells$cells <- rownames(y_to_cells)
  y_to_cells$Y <- y_to_cells$V1
  root <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
  root <- y_to_cells$Y[y_to_cells$cells == root]
  mst <- monocle::minSpanningTree(cds)
  endpoints <- names(which(igraph::degree(mst) == 1))
  endpoints <- endpoints[endpoints != paste0("Y_", root)]
  cellWeights <- lapply(endpoints, function(endpoint) {
    path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
    path <- as.character(path)
    df <- y_to_cells[y_to_cells$Y %in% path, ]
    df <- data.frame(weights = as.numeric(colnames(cds) %in% df$cells))
    colnames(df) <- endpoint
    return(df)
  }) %>% do.call(what = 'cbind', args = .)
  pseudotime <- sapply(cellWeights, function(w) cds$Pseudotime)
  rownames(cellWeights) <- rownames(pseudotime) <- colnames(cds)
  # Get the lineages representation
  edges_rd <- t(monocle::reducedDimK(cds)) %>% as.data.frame()
  rd_lineages <- lapply(endpoints, function(endpoint){
    path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
    path <- as.character(path)
    path <- paste("Y", path, sep = "_")
    return(edges_rd[path, ])
  })
  return(list("pseudotime" = pseudotime,
              "cellWeights" = as.matrix(cellWeights)))
}
