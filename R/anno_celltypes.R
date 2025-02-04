#' Annotate each cell with the corresponding cell type
#'
#' @param object Seurat object
#'
#' @return A Seurat object
#' @export
#'
#' @examples
#' \dontrun{
#' sobj <- anno_celltypes(object = sobj, anno_level = 2, species = "Hs")
#' }
anno_celltypes <- function(object, anno_level = 2, selfClusters = NULL, species = "Hs", seed = 42, ...) {
  default_assay <- Seurat::DefaultAssay(object)
  Seurat::DefaultAssay(object) <- "RNA"

  ## load Human Panglao database
  base::load(system.file("data", "Human_PanglaoDB.Rdata", package = "scMRMA"))


  anno_res <- scMRMA::scMRMA(
    input = object,
    species = species,
    db = "panglaodb",
    selfClusters = selfClusters,
    ...
  )
  if (length(anno_level) == 1) {
    object@meta.data$cell_type <-
      anno_res$multiR$annotationResult[[paste0("Level", anno_level)]]
  } else {
    for (i in anno_level) {
      object@meta.data[[paste0("cell_type_", i)]] <-
        anno_res$multiR$annotationResult[[paste0("Level", i)]]
    }
    ## additionaly add as cell type (for highest (lowest number) provided level)
    object@meta.data[[paste0("cell_type")]] <-
      anno_res$multiR$annotationResult[[paste0("Level", min(anno_level))]]
  }

  Seurat::DefaultAssay(object) <- default_assay

  return(object)
}



#' Create A UMAP from the annotated Seurat object
#'
#' @param object A Seurat object
#' @param group.by Paraneter by witch to group visualization by
#' @param ndims Number of dimensions for dimensional reduction
#'
#' @return p A Dim plot
#' @export
#'
#' @examples
#' \dontrun{
#' p <- visualize_data(sobj, group.by = "cell_type", ndims = 30)
#' }
visualize_data <- function(object, group.by = "cell_type", ndims = NULL, ...) {
  object <- object %>% Seurat::RunPCA(npcs = max(ndims, 100))

  ## estimate dimensionality of data if no ndims is supplied
  if (is.null(ndims)) {
    ndims <- ceiling(intrinsicDimension::maxLikGlobalDimEst(object@reductions[[paste0("pca")]]@cell.embeddings, k = 20)[["dim.est"]])
  }

  object <- object %>% Seurat::RunUMAP(dims = 1:ndims)

  p <- Seurat::DimPlot(object, group.by = group.by, ...)
  base::print(p)

  return(p)
}




#' Title
#'
#' @param object A
#' @param method A
#' @param samples A
#' @param resolution A
#'
#' @return object
#' @export
#'
#' @examples
#' \dontrun{
#' #' integrate_samples()
#' }
integrate_samples <- function(object, method = "rpca", samples = "samples", seed = 42, npcs = 100, k.weight = 100, verbose = FALSE) {
  set.seed(42)
  if (typeof(object) == "list") {
    object_list <- object
  } else {
    object_list <- Seurat::SplitObject(object, split.by = samples)
  }

  ## normalize data and find variable features for each sample
  object_list <- lapply(object_list, function(ob) {
    ob <- ob %>%
      Seurat::NormalizeData(verbose = verbose) %>%
      Seurat::FindVariableFeatures(verbose = verbose)
    ob
  })

  ## select features that are repeatedly variable across samples
  features <- Seurat::SelectIntegrationFeatures(object.list = object_list, verbose = verbose)

  ## scale and calculate PCs for each sample
  object_list <- lapply(object_list, function(ob) {
    ob <- ob %>%
      Seurat::ScaleData(verbose = verbose, features = features) %>%
      Seurat::RunPCA(verbose = verbose, npcs = npcs, features = features)
    ob
  })


  ## find integration anchors
  integration.anchors <- Seurat::FindIntegrationAnchors(object.list = object_list, anchor.features = features, reduction = "rpca", verbose = verbose)

  object_integrated <- Seurat::IntegrateData(anchorset = integration.anchors, k.weight = k.weight, verbose = verbose)

  Seurat::DefaultAssay(object_integrated) <- "integrated"

  # ## run default processing steps also with rna assay if later used
  Seurat::DefaultAssay(object_integrated) <- "RNA"
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    object_integrated <- object_integrated %>% SeuratObject::JoinLayers()
  }
  object_integrated <- object_integrated %>%
    Seurat::NormalizeData(verbose = verbose) %>%
    Seurat::FindVariableFeatures(verbose = verbose) %>%
    Seurat::ScaleData(verbose = verbose)

  Seurat::DefaultAssay(object_integrated) <- "integrated"

  return(object_integrated)
}


#' Title
#'
#' @param object A
#' @param resolution A
#'
#' @return object
#' @export
#'
#' @examples
#' \dontrun{
#' cluster_data(object, resolution = 0.8)
#' }
cluster_data <- function(object, resolution = 0.8, npcs_calculate = 100, npcs = NULL, seed = 42, verbose = FALSE) {
  set.seed(seed)
  default_assay <- Seurat::DefaultAssay(object)

  if (default_assay == "integrated") {
    object <- object %>%
      Seurat::ScaleData(verbose = verbose) %>%
      Seurat::RunPCA(assay = default_assay, npcs = npcs_calculate, verbose = verbose)
  } else {
    object <- object %>%
      Seurat::NormalizeData(verbose = verbose) %>%
      Seurat::ScaleData(verbose = verbose) %>%
      Seurat::FindVariableFeatures(verbose = verbose) %>%
      Seurat::RunPCA(assay = default_assay, npcs = npcs_calculate, verbose = verbose)
  }

  if (is.null(npcs)) {
    ndims <- ceiling(intrinsicDimension::maxLikGlobalDimEst(object@reductions[[paste0("pca")]]@cell.embeddings,
      k = 20
    )[["dim.est"]])
  } else {
    ndims <- npcs
  }

  print(paste0("Number of used dimensions for clustering: ", ndims))
  object <- object %>%
    Seurat::RunUMAP(dims = 1:ndims, verbose = verbose) %>%
    Seurat::FindNeighbors(dims = 1:ndims, reduction = "pca", verbose = verbose) %>%
    Seurat::FindClusters(resolution = resolution, verbose = verbose)

  return(object)
}
