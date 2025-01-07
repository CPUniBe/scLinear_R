

#' Prepare data for modality prediction
#'
#' @param object A Seurat object
#' @param remove_doublets Should doublets be removed? (TRUE,FALSE)
#' @param low_qc_cell_removal Should low quality cells be removed using median absolute deviation?
#' @param anno_level Level of annotation used for cell type annotation. Either a single number or a vector with multiple numbers c(1,2,3,4).
#' @param samples Variable that contains the sample names. If set, will use rpca to integrate the samples.
#' @param integrate_data Should the sample be integrated through the rpca method?
#' @param remove_empty_droplets If the raw unfiltered matrix was used to create the Seurat object, empty droplets should be filtered out. (TRUE, FALSE)
#' @param lower Lower boundary to define empty droplets. All cells with a lower amount of nFeatures are assumed to be empty.
#' @param FDR FDR threshold to define non empty cells.
#' @param annotation_selfCluster Should the clusters determined by Seurat be used for cell type annotation.
#' @param resolution Resolution for louvain clustering.
#' @param seed Used seed.
#' @param return_plots Should plots be returned from function
#' @param print_plots Print plots. (TRUE, FALSE)
#' @param species Species. Relevant for cell type annotation. ("Hs", "Mm")
#' @param min.features Minimum ammount of features per cell. Replaces automatically determined threshold if bigger.
#'
#' @return object A pre-processed Seurat object  with annotated cell types
#' @export
#'
#' @examples
#' \dontrun{
#' sobj <- scLinear(object = sobj, remove_doublets = TRUE, low_qc_cell_removal = TRUE, anno_level = 2, samples = NULL, integrate_data = FALSE, resolution = 0.8)
#' }

prepare_data <- function(object, remove_doublets = TRUE, low_qc_cell_removal = TRUE, anno_level = 2, samples = NULL, integrate_data = FALSE,remove_empty_droplets = FALSE, lower = 100, FDR = 0.01, annotation_selfCluster = TRUE, resolution = 0.8, seed = 42, return_plots = FALSE, print_plots = TRUE, species = "Hs", min.features = NULL, verbose = FALSE){
  set.seed(seed)

  plot_list <- list()

  Seurat::DefaultAssay(object) <- "RNA"

  if(remove_empty_droplets){
    object <- empty_drops(object = object, lower = lower, FDR = FDR, samples = samples, seed = seed)
    #plot_list[["empty_dropts"]] <- object[[2]]
    #object <- object[[1]]
  }

  if(!("mito_percent" %in% names(object@meta.data))){
    if(species == "Hs"){
      object$mito_percent <- Seurat::PercentageFeatureSet(object, pattern = "^MT-")
    }else{
      object$mito_percent <- Seurat::PercentageFeatureSet(object, pattern = "^mt-")
    }

  }

  if(remove_doublets){
    print("Start remove doublets")
    object <- object %>% remove_doublets(samples = samples, print_plots = print_plots, seed = seed, verbose = verbose)
    plot_list[["doublets"]] <- object[[2]]
    object <- object[[1]]
  }

  if(low_qc_cell_removal){
    print("Start low quality cell removal")
    object <- object %>% mad_filtering(samples = samples, print_plots = print_plots, seed = seed, min.features = min.features, verbose = verbose)
    plot_list[["low_qc_cells"]] <- object[[2]]
    object <- object[[1]]
  }

  if(integrate_data){
    print("Start integrate data")
    object <- integrate_samples(object, samples = samples, seed = seed, verbose = verbose)
  }

  print("Start clustering data")
  object <- cluster_data(object, resolution = resolution, seed = seed)
  Seurat::Idents(object) <- object@meta.data[["seurat_clusters"]]

  print("Start cell type annotation")
  if(annotation_selfCluster){
    object <- object %>% anno_celltypes(anno_level = anno_level, selfClusters = Seurat::Idents(.), species = species, seed = seed)
  }else{
    object <- object %>% anno_celltypes(anno_level = anno_level, selfClusters = NULL, species = species, seed = seed)
  }

  p1 <- Seurat::DimPlot(object, group.by = "cell_type", label = TRUE, repel = TRUE) + ggplot2::theme(legend.position = "null")
  if(print_plots){base::print(p1)}

  if(return_plots){
    return_object <- list(object = object, plots = plot_list)
  }else{
    return_object <- object
  }


  return(return_object)

}


#' Predict modalities based on gene expression data
#'
#' @param object A Seurat object
#'
#' @return object A Seurat object containing additional single cell modalities
#' @export
#'
#' @examples
#' \dontrun{
#' sobj <- scLinear(object = sobj)
#' }
scLinear <- function(object, remove_doublets = TRUE, low_qc_cell_removal = TRUE, anno_level = 2, samples = NULL, integrate_data = FALSE, remove_empty_droplets = FALSE, lower = 100, FDR = 0.01, annotation_selfCluster = TRUE, resolution = 0.8, seed = 42, return_plots = FALSE, model = "all", assay_name = "RNA", print_plots = FALSE, species = "Hs", min.features = NULL, verbose = FALSE){
  set.seed(seed)
  object <- prepare_data(object,
                         remove_doublets = remove_doublets,
                         low_qc_cell_removal = low_qc_cell_removal,
                         anno_level = anno_level,
                         samples = samples,
                         integrate_data = integrate_data,
                         remove_empty_droplets = remove_empty_droplets,
                         lower = lower,
                         FDR = FDR,
                         annotation_selfCluster = annotation_selfCluster,
                         resolution = resolution,
                         seed = seed,
                         return_plots = FALSE,
                         print_plots = print_plots,
                         species = species,
                         min.features = min.features,
                         verbose = verbose)

  pipe <- create_adt_predictor()
  pipe <- load_pretrained_model(pipe, model = model)

  object[["predicted_ADT"]] <-  adt_predict(pipe = pipe,
                                                  gexp = Seurat::GetAssay(object, assay = assay_name),
                                                  normalize = TRUE)

  return(object)

}



fit_predictor <- function(gexp_train,adt_train, gexp_test = NULL,
                            layer_gex = "counts", layer_adt = "counts",
                            normalize_gex = TRUE,normalize_adt = TRUE, margin = 2,
                            n_components = 300, zscore_relative_to_tsvd = 'after'){
  
  # If objects are passed as matrices, leave as is. If passed as Seurat objects --> extract count data for assays
  if(class(gexp_train)[1] == "Assay" | class(gexp_train)[1] == "Assay5"){ gexp_train <- Seurat::GetAssayData(gexp_train, layer = layer_gex) }
  if(class(adt_train)[1] == "Assay" | class(adt_train)[1] == "Assay5"){ adt_train <- Seurat::GetAssayData(adt_train, layer = layer_adt) }
  if(!is.null(gexp_test)){
    if(class(gexp_test)[1] == "Assay" |class(gexp_test)[1] == "Assay5"){ gexp_test <- Seurat::GetAssayData(gexp_test, layer = layer_gex) }
  }
}

# Why set prediction to 0 if negative? We're trying to predict the CLR of ADT not ADT itself --> can be negative
# For some reason this never seems to be the case
predict <- function(predictor,gexp,layer="counts",normalize_gex=TRUE){
  
  if(any(class(gexp) %in% c("Seurat", "Assay", "Assay5"))){
    gexp <- Seurat::GetAssayData(gexp, layer = layer)
  }else{ # assume it is a matrix type
    gexp <- Matrix::Matrix(gexp, sparse = TRUE)
  }
  
  if(normalize_gex){
    gexp <- gexp_normalize(gexp)
  }
  
  # Bring into shape cells x genes (done during training in the step which combines train & test)
  gexp <- filter_input_genes(gexp,predictor)
  gexp <- Matrix::t(gexp)
  
  if(predictor$zscore_relative_to_tsvd == 'before'){
    gexp <- Matrix::t(apply(gexp, 1, function(x) {
      (x - mean(x)) / sd(x)
    }
    ))  
  }
  
  gexp_projected <- gexp %*% predictor$tsvd$v
  
  if(predictor$zscore_relative_to_tsvd == 'after'){
    gexp_projected <- Matrix::t(apply(gexp_projected, 1, function(x) {
      (x - mean(x)) / sd(x)
    }
    ))
  }
  
  coeff_matrix <- matrix(nrow = ncol(gexp_projected)+1)
  # Last coefficient usually 0 --> reason not quite clear but possibly because there is internal linear dependence
  for(model in predictor$lms){
    coeff_matrix <- cbind(coeff_matrix,model$coefficients)
  }
  coeff_matrix[is.na(coeff_matrix)] <- 0
  # Drop the empty column from the coeff matrix (used only to initialize the object)
  coeff_matrix <- coeff_matrix[,2:ncol(coeff_matrix)]
  
  # Add column of 1s 'to the left' of the tSVD projection of the test data --> adds intercept of each model (intercept is the first coefficient of the lm)
  lm_input <- cbind(rep(1,nrow(gexp_projected)),gexp_projected)
  # Multiply tSVD projected & normed input data with LM coefficients
  res <- lm_input %*% coeff_matrix
  colnames(res) <- names(predictor$lms)
  return(res)
}

evaluate_predictor <- function(predictor,gexp_test,adt_test, normalize_gex = TRUE, normalize_adt = TRUE, margin = 2){
  predicted_adt <- predict(predictor,gexp_test,normalize_gex = normalize_gex)
  if(class(adt_test)[1] == "Assay" |class(adt_test)[1] == "Assay5"){ adt_test <- Seurat::GetAssayData(adt_test, layer = 'counts') }
  real_adt_clr <- t(Seurat::NormalizeData(adt_test, normalization.method = "CLR", margin = margin))
  p_adt <- subset(predicted_adt,features = which(rownames(predicted_adt) %in% rownames(real_adt_clr)) )
  t_adt <- subset(real_adt_clr,features = which(rownames(real_adt_clr) %in% rownames(predicted_adt)) )
  err_sq <- (p_adt-t_adt)^2
  # Not sure how to interpret the means of the single model metrics but for now just replicating the behaviour of python code
  rmse <- mean(sqrt(colSums(err_sq)/nrow(err_sq)))
  # Pearson calculated in a really strange way in python --> why would we average across cells and not across models (if at all?)
  # Here deviating from python version and averaging across models
  # Better --> keep model specific metrics. Some perform really well, others really poorly
  # pearson <- mean(diag(cor(p_adt,t_adt,method = "pearson")))
  pearson <- diag(cor(p_adt,t_adt,method = "pearson"))
  # spearman <- mean(diag(cor(p_adt,t_adt,method = "spearman")))
  spearman <- diag(cor(p_adt,t_adt,method = "spearman"))
  return(list(rmse = rmse,pearson = pearson, spearman = spearman))
}

feature_importance <- function(predictor,gexp,layer_gexp,normalize_gex = TRUE){
  if(!(predictor$zscore_relative_to_tsvd == 'after')){
    print('Feature importance calculation not yet implemented for model with z-score normalization BEFORE tSVD')
    return(NULL)
  }
  
  if(any(class(gexp) %in% c("Seurat", "Assay", "Assay5"))){
    gexp <- Seurat::GetAssayData(gexp, layer = layer_gexp)
  }else{ # assume it is a matrix type
    gexp <- Matrix::Matrix(gexp, sparse = TRUE)
  }
  
  # Preprocessing equivalent to how it is done for fit & predict
  if(normalize_gex){
    gexp <- gexp_normalize(gexp)
  }
  
  gexp <- filter_input_genes(gexp,predictor)
  gexp <- Matrix::t(gexp)
  
  gexp_projected <- gexp %*% predictor$tsvd$v
  
  gexp_projected <- Matrix::t(apply(gexp_projected, 1, function(x) {
    (x - mean(x)) / sd(x)
  }
  ))
  
  # Total gradient = W x J x t(v) where W = LM weights, J = Jacobian of z-score transformation and t(v) = Transpose of 'V' from tSVD
  
  # Get W
  coeff_matrix <- matrix(nrow = ncol(gexp_projected)+1)
  # Last coefficient usually 0 --> reason not quite clear but possibly because there is internal linear dependence (perhaps from z-score since all components should add up to mean 1?)
  for(model in predictor$lms){
    coeff_matrix <- cbind(coeff_matrix,model$coefficients)
  }
  coeff_matrix[is.na(coeff_matrix)] <- 0
  # Drop the empty column from the coeff matrix (used only to initialize the object)
  coeff_matrix <- coeff_matrix[,2:ncol(coeff_matrix)]
  
  # Drop the intercept --> not used for derivative d_prediction/d_gex --> constants drop from derivative
  coeff_matrix <- t(coeff_matrix[2:nrow(coeff_matrix),])
  
  # Auxilliary function --> calculate jacobian of one cell so we can then apply this function across all cells
  cellwise_jacobian <- function(cell_projection){
    jacobian <- numDeriv::jacobian(func = function(x) {(x - mean(x)) / sd(x)}, cell_projection)
  }
  # Slow but that's a looooot of matrix multiplications to run so probably to be expected
  Js <- apply(gexp_projected,1,cellwise_jacobian,simplify = FALSE)
  
  WJ <- lapply(Js,function(J){return(coeff_matrix %*% J)})
  v_t <- Matrix::t(predictor$tsvd$v)
  
  WJV <- abind::abind(lapply(WJ,function(WJ){return(WJ %*% v_t)}), along = 3)
  # Axis 1 = model, Axis 2 = Gene, Axis 3 = Cell --> taking mean 'across cells' = mean over margin of axis 1&2
  feature_importance_means <- apply(WJV, c(1, 2), mean)
  colnames(feature_importance_means) <- colnames(gexp)
  rownames(feature_importance_means) <- names(predictor$lms)
  
  return(feature_importance_means)
}


#' Normalize gene expression matrix with scran and scuttle
#'
#' @param gexp_matrix A gene expression matrix
#' @param center.size.factors A
#' @param log A
#' @param ... For the method, additional arguments passed to logNormCounts.
#'
#' @return Normalized expression matrix
#' @export
#'
#' @examples
#' \dontrun{
#' # Normalize expression matirx
#' normalized_matrix <- gexp_normalize(sobj\@assays[["RNA"]]\@counts)
#' # Add normalized matrix back to RNA assay in Seurat.
#' sobj\@assays[["RNA"]]\@data <- normalized_matrix
#' }
gexp_normalize <- function(gexp_matrix, center.size.factors = FALSE, log = FALSE, ...){
  ## normalize data GEX
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = gexp_matrix))
  clusters <- scran::quickCluster(sce)
  sce <- scran::computeSumFactors(sce, clusters=clusters)
  sce <- scuttle::logNormCounts(sce, center.size.factors = center.size.factors, log = log, ...)
  gexp_matrix <- sce@assays@data@listData[["normcounts"]]
  gexp_matrix <- base::log1p(gexp_matrix)
  return(gexp_matrix)
}

