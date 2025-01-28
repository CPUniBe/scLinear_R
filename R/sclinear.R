

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


#' Forces a gene expression input matrix to match the gene expression input used in training of a given predictor.
#' Expression of genes which were not used in training the predictor is discarded, expression of genes used in training but not present in the input
#' are set to 0. Also reorders columns so they match the expected order for the tSVD projection
#' @param gexp gene expression matrix
#' @param predictor predictor object created by fit_predictor function
#'
#' @return gene expression matrix reshaped to match the set of genes used during training of the predictor (genes not present during training discarded, missing genes set to 0)
filter_input_genes <- function(gexp,predictor){
  to_add <- setdiff(predictor$genes_considered,row.names(gexp))
  artificial_gene_counts <- Matrix::Matrix(0,nrow = length(to_add),ncol = ncol(gexp),
                                           dimnames = list(to_add,colnames(gexp)))
  Matrix::Matrix(as.array(rbind(gexp,artificial_gene_counts))[predictor$genes_considered,],sparse = TRUE)
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
scLinear <- function(object, remove_doublets = TRUE, low_qc_cell_removal = TRUE, anno_level = 2, samples = NULL, integrate_data = FALSE, remove_empty_droplets = FALSE, lower = 100, FDR = 0.01, annotation_selfCluster = TRUE, resolution = 0.8, seed = 42, return_plots = FALSE, model = "model_full_neurips", assay_name = "RNA", print_plots = FALSE, species = "Hs", min.features = NULL, verbose = FALSE){
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

  pipe <- readRDS(model)

  object[["predicted_ADT"]] <-  predict(pipe,gexp = Seurat::GetAssay(object, assay = assay_name))

  return(object)

}

#' Trains a prediction model to predict (CLR-transformed & Seurat normed) ADT expression levels based on gene expression data
#'
#' @param gexp_train Seurat Object containing gene expression data in the layer specified by argument layer_gex
#' @param adt_test Seurat Object containing ADT expression data in the layer specified by argument layer_adt (for the same cells as gexp_train)
#' @param gexp_test Seurat Object of same form as gexp_train. If provided, the 'test-data' is used in the truncated singular value decomposition
#' (affects how input gene expression data is projected into lower dimensional space) but is not used for training the actual ADT-predictor itself.
#' Not yet implemented properly so gexp_test should be left NULL (default value)
#' @param layer_gex From which layer of the Seurat Object should the gene expression data be extracted
#' @param layer_adt From which layer of the Seurat Object should the ADT expression data be extracted
#' @param normalize_gex Should gene expression levels be normalized across cells (see function gexp_normalize for details)
#' @param normalize_adt Should ADT levels be normalized (using the NormalizeData function from the Seurat package with method CLR)
#' @param margin Margin to apply CLR normalization over for ADT levels
#' @param n_components Rank for the truncated singular value decomposition (= desired dimensionality of gene expression data after dimensionality reduction)
#' @param zscore_relative_to_tsvd Specifies whether to apply z-score normalization to the gene expression data before or after tSVD projection.
#' ('before'/'after' otherwise skips z-score normalization entirely)
#' @param n_cores Number of cores made available for parallelization of lm-fitting. Defaults to all 'available' cores minus 4.
#' @return A predictor in the form of a list containing the following components:
#' 1. 'tsvd_v': Right singular vector of tsvd used for dimensionality reduction
#' 2. 'lm_coefficients': List of the actual 'prediction models' with weights for each 'component' of the tSVD (one model per ADT)
#' 3. 'zscore_relative_to_tsvd': See explanation of input parameter with the same name. Saved in the output so that the same setting
#' will also be used during the preprocessing of data during any subsequent uses of the predictor
#' 4. 'genes_considered': Genes used during training of the predictor (genes with any non-zero counts in the training data)
#' @export
fit_predictor <- function(gexp_train,adt_train, gexp_test = NULL,
                            layer_gex = "counts", layer_adt = "counts",
                            normalize_gex = TRUE,normalize_adt = TRUE, margin = 2,
                            n_components = 300, zscore_relative_to_tsvd = 'after',n_cores = NULL){
  # If no number is specified, get number of available cores and omit 4 to avoid using up all system resources
  if(is.null(n_cores)){n_cores <- parallelly::availableCores(omit = 4)}
  # If objects are passed as matrices, leave as is. If passed as Seurat objects --> extract count data for assays
  if(class(gexp_train)[1] == "Assay" | class(gexp_train)[1] == "Assay5"){ gexp_train <- Seurat::GetAssayData(gexp_train, layer = layer_gex) }
  if(class(adt_train)[1] == "Assay" | class(adt_train)[1] == "Assay5"){ adt_train <- Seurat::GetAssayData(adt_train, layer = layer_adt) }
  if(!is.null(gexp_test)){
    if(class(gexp_test)[1] == "Assay" |class(gexp_test)[1] == "Assay5"){ gexp_test <- Seurat::GetAssayData(gexp_test, layer = layer_gex) }
  }

  #TODO: Change it so normalization happens for both train and test at the same time
  if(normalize_gex){
    gexp_train <- gexp_normalize(gexp_train)
    if( !is.null(gexp_test)){gexp_test <- gexp_normalize(gexp_test)}
  }
  if(normalize_adt){
    adt_train <- Seurat::NormalizeData(adt_train, normalization.method = "CLR", margin = margin)
  }

  # The way the sets are generated, same gene names should be a given (stem from same dataset just split)
  if(!is.null(gexp_test)){
    training_set <- Matrix::t(cbind(gexp_train,gexp_test))
  }else{training_set <- Matrix::t(gexp_train)}

  # Keep all genes expressed in at least one cell --> might work poorly if it needs to filter out columns (filter out cells instead of genes)
  training_set <- training_set[unique(Matrix::summary(training_set)$i),unique(Matrix::summary(training_set)$j)]
  keep_genes <- colnames(training_set)
  # z-score normalization --> transpose after applying z-score is necessary to get THE SAME dimension as the input
  if(zscore_relative_to_tsvd == 'before'){
    training_set <- Matrix::t(apply(training_set, 1, function(x) {
      (x - mean(x)) / sd(x)
    }
    ))
  }
  # Create tSVD decomposition
  print('Calculating truncated singular value decomposition - for large input matrices this may take several minutes')
  trained_tsvd <- sparsesvd::sparsesvd(training_set,rank = n_components)
  print('tSVD done')
  # apply tSVD projection on input data
  training_set <- training_set %*% trained_tsvd$v
  # z-score normalizaton

  if(zscore_relative_to_tsvd == 'after'){
    training_set <- Matrix::t(apply(training_set, 1, function(x) {
      (x - mean(x)) / sd(x)
    }
    ))
  }

  # If lm-test-data was used for tsvd training step, reduce lm input back to only train data
  training_set <- training_set[intersect(colnames(gexp_train),rownames(training_set)),]

  # Filter out adt data for cells dropped due to RNA 0-counts
  # transposed to match expected input format for lm
  # here the indexing works rather quickly --> no conversion to array necessary
  adt_train_modelling <- Matrix::t(adt_train)[intersect(rownames(training_set),colnames(adt_train)),]


  rm(gexp_train,gexp_test)
  print('Fitting linear models')
  cl <- parallel::makeCluster(n_cores,outfile = 'lm_log.txt')
  parallel::clusterExport(cl,list("training_set","adt_train_modelling"),envir = environment())
  parallel::clusterEvalQ(cl,library(Matrix))
  results <- pbapply::pblapply(cl = cl, X = 1:ncol(adt_train_modelling), FUN = function(i) {
    tryCatch({
      return(lm(adt_train_modelling[, i] ~ training_set)$coefficients)  # Fit a linear model
    }, error = function(e) {
      message("Error fitting model for response ", i, ": ", e$message)
      return(NULL)  # Return NULL if there's an error
    })
  })
  parallel::stopCluster(cl)
  file.remove('lm_log.txt')
  names(results) <- colnames(adt_train_modelling)
  return(list(tsvd_v=trained_tsvd$v,lm_coefficients=results,zscore_relative_to_tsvd = zscore_relative_to_tsvd,
              genes_considered = keep_genes))
}

# Why set prediction to 0 if negative? We're trying to predict the CLR of ADT not ADT itself --> can be negative
# For some reason this never seems to be the case
#' Predicts ADT levels based on the gene expression levels.
#' @param predictor predictor trained with 'fit_predictor' function
#' @param gexp Seurat Object containing gene expression data in the layer specified by argument layer
#' @param layer From which layer of the Seurat Object should the gene expression data be extracted
#' @param normalize_gex Should gene expression levels be normalized across cells (see function gexp_normalize for details)
#'
#' @return Predicted ADT levels. If ADT levels have been normalized during predictor training (default behaviour), the predicted ADT levels should also be considered 'normalized' in the same way.
#'
#' @export
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

  gexp_projected <- gexp %*% predictor$tsvd_v

  if(predictor$zscore_relative_to_tsvd == 'after'){
    gexp_projected <- Matrix::t(apply(gexp_projected, 1, function(x) {
      (x - mean(x)) / sd(x)
    }
    ))
  }

  coeff_matrix <- matrix(nrow = ncol(gexp_projected)+1)
  # Last coefficient usually 0 --> reason not quite clear. Does z-score normalization introduce linear dependence since mean must be 0
  for(model in predictor$lm_coefficients){
    coeff_matrix <- cbind(coeff_matrix,model)
  }
  coeff_matrix[is.na(coeff_matrix)] <- 0
  # Drop the empty column from the coeff matrix (used only to initialize the object)
  coeff_matrix <- coeff_matrix[,2:ncol(coeff_matrix)]

  # Add column of 1s 'to the left' of the tSVD projection of the test data --> adds intercept of each model (intercept is the first coefficient of the lm)
  lm_input <- cbind(rep(1,nrow(gexp_projected)),gexp_projected)
  # Multiply tSVD projected & normed input data with LM coefficients
  res <- lm_input %*% coeff_matrix
  colnames(res) <- names(predictor$lm_coefficients)
  return(res)
}

#' Reports the mean RMSE of all ADT-specific regression models as well as the Pearson-/ & Spearman correlation coefficients
#' between predicted and measured ADT for each model separately
#' Unlike other predictor related functions, layers
#' @param predictor ADT-predictor (of form as produced by fit_predictor function) whose performance should be evaluated
#' @param gexp_test Seurat object containing gene expression data to use for ADT level prediction
#' @param adt_test Measured ADT levels to which the predicted ADT levels are to be compared
#' @param gexp_layer From which layer of the Seurat Object should the gene expression data be extracted
#' @param adt_layer From which layer of the Seurat Object should the gene expression data be extracted
#' @param normalize_gex Should gene expression levels be normalized across cells (see function gexp_normalize for details)
#' @param normalize_adt Should ADT levels be normalized (using the NormalizeData function from the Seurat package with method CLR) --> this needs
#' to be the same as the setting used during training of the prediction model.
#' @param margin Margin to apply CLR normalization over for ADT levels
#'
#' @return List of metrics for model performance (Mean RMSE of all models and ADT-specific correlation coefficients between predicted and measured ADT Values)
#' @export
evaluate_predictor <- function(predictor,gexp_test,adt_test,gexp_layer = 'counts', adt_layer = 'counts', normalize_gex = TRUE, normalize_adt = TRUE, margin = 2){

  # ToFix --> can't compare CLR if different subsets were considered since it's compositional data
  predicted_adt <- predict(predictor,gexp_test,layer = gexp_layer, normalize_gex = normalize_gex)
  if(class(adt_test)[1] == "Assay" |class(adt_test)[1] == "Assay5"){ adt_test <- Seurat::GetAssayData(adt_test, layer = adt_layer) }
  if(normalize_adt){
    adt_test <- t(Seurat::NormalizeData(adt_test, normalization.method = "CLR", margin = margin))
  }else{
      # Converting from Matrix::Matrix to base Matrix (and transposing)
      adt_test <- as.matrix(Matrix::t(adt_test))
      }

  p_adt <- subset(predicted_adt,features = which(colnames(predicted_adt) %in% colnames(adt_test)) )
  t_adt <- subset(adt_test,features = which(colnames(adt_test) %in% colnames(predicted_adt)) )

  t_adt <- t_adt[match(rownames(p_adt),rownames(t_adt)),match(colnames(p_adt),colnames(t_adt))]
  err_sq <- (p_adt-t_adt)^2
  # Not sure how to interpret the means of the single model metrics but for now just replicating the behaviour of python code
  rmse <- mean(sqrt(colSums(err_sq)/nrow(err_sq)))
  # Pearson calculated in a really strange way in python --> why would we average across cells and not across models (if at all?)
  # Here reporting the correlation coefficients of each ADT prediction separately
  pearson <- diag(cor(t_adt,p_adt,method = "pearson"))
  spearman <- diag(cor(t_adt,p_adt,method = "spearman"))
  # Original behaviour:
  mean_pearson <- mean(unlist(lapply(1:nrow(t_adt),function(i) cor(t_adt[i,],p_adt[i,],method = 'pearson'))))
  mean_spearman <- mean(unlist(lapply(1:nrow(t_adt),function(i) cor(t_adt[i,],p_adt[i,],method = 'spearman'))))
  return(list(rmse = rmse,pearson = pearson, spearman = spearman, mean_pearson = mean_pearson, mean_spearman = mean_spearman))
}

# Auxilliary function --> calculate jacobian of one cell so we can then apply this function across all cells
# Pass N as input so it doesn't have to be calculated on each call and all projections have the same number of components
# Unclear if better to use N or N-1 for standard deviation --> in that case only the N*sd would change to (N-1) * sd --> the other Ns come from mean not sd
cellwise_jacobian <- function(cell_projection){
  N <- length(cell_projection)
  sd <- sd(cell_projection)
  m <- mean(cell_projection)
  non_constant <- Matrix::tcrossprod(cell_projection-m)
  # Derivative of standard score of input x_i by x_j --> diagonal is correction term for when i = j
  return((((-sd/N) - (non_constant/(N*sd)))/(sd^2)) + diag(sd, nrow = N))
}

# defining matrix product as a function instead of operator just because operator %*% can't be passed easily to apply function
# without slowing down parallelization (no FUN = function()...)
# Need f(y,x) = x %*% y --> order swapped
# Function will be called by apply and m2 is the big object to iterate over so it will be the first argument
matprod_rl <- function(m2,m1){
  return (m1 %*% m2)
}

matprod_lr <- function(m1,m2){
  return (m1 %*% m2)
}

# Helper function for feature importance
cross_cell_average_fi_c <- function(WJ_single_model,v){
  gc()
  WJV <- matrix_product(t(WJ_single_model),v)
  return(colMeans(WJV))
}

#' Feature Importance
#'
#' Only implemented for the case where z-score normalization was applied AFTER dimensionality reduction. \cr
#' Calculates the derivative of the predicted level of a given ADT as a function of the expression levels of a single gene 'dADT/dGene'\cr
#' Obtained by decomposing d(ADT)/d(Gene) = WJV where\cr
#' W = d(ADT)/d(zscores)\cr
#' J = d(zscores)/d(tSVD components)\cr
#' V = d(tSVD components)/d(gene expression)\cr
#' W and V directly correspond to the weights of the linear regression model & the right singular vectors of the tSVD.\cr
#' However, the effect of the z-score normalization had to be numerically approximated using the function 'jacobian' from the package 'numDeriv'.
#' @param predictor Predictor created by fit_predictor function
#' @param gexp Seurat object containing gene expression data
#' @param layer_gexp From which layer of the gexp Seurat Object should the gene expression data be extracted
#' @param normalize_gex Should gene expression levels be normalized across cells (see function gexp_normalize for details)
#'
#' @return Effect of the expression level of each gene on the prediction of the ADT values (averaged across all cells present in the input).
#' @export
feature_importance <- function(predictor,gexp,layer_gexp,normalize_gex = TRUE,n_cores = NULL){
  if(is.null(n_cores)){n_cores <- parallelly::availableCores(omit = 4)}
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
  gexp_names <- colnames(gexp)

  gexp_projected <- gexp %*% predictor$tsvd_v

  # No normalization here --> we want dz/dprojection evaluated AT THE CURRENT PROJECTION
  # gexp_projected <- Matrix::t(apply(gexp_projected, 1, function(x) {
  #   (x - mean(x)) / sd(x)
  # }
  # ))

  # Total gradient = W x J x t(v) where W = LM weights, J = Jacobian of z-score transformation and t(v) = Transpose of 'V' from tSVD

  # Get W
  coeff_matrix <- matrix(nrow = ncol(gexp_projected)+1)
  # Last coefficient usually 0 --> reason not quite clear but possibly because there is internal linear dependence (perhaps from z-score since all components should add up to mean 1?)
  for(model in predictor$lm_coefficients){
    coeff_matrix <- cbind(coeff_matrix,model)
  }
  coeff_matrix[is.na(coeff_matrix)] <- 0
  # Drop the empty column from the coeff matrix (used only to initialize the object)
  coeff_matrix <- coeff_matrix[,2:ncol(coeff_matrix)]

  # Drop the intercept --> not used for derivative d_prediction/d_gex --> constants drop from derivative
  coeff_matrix <- t(coeff_matrix[2:nrow(coeff_matrix),])


  # environment(cellwise_jacobian) <- environment(feature_importance)
  on.exit(parallel::stopCluster(cl))
  # Slow but that's a looooot of matrix multiplications to run so probably to be expected

  N <- ncol(gexp_projected)
  cl <- parallel::makeCluster(n_cores,outfile = 'feature_importance_log.txt')
  parallel::clusterExport(cl,list('cellwise_jacobian'))
  print('Calculating Jacobian of z-score normalization step')

  Js <- pbapply::pbapply(cl=cl,X = gexp_projected, MARGIN = 1,FUN = cellwise_jacobian,simplify = FALSE)
  parallel::stopCluster(cl)

  #transposing so we can use crossprod for apply function --> circumvents having to do FUN = function(X)
  # ... which slows down parallelization when called from environment other than global
  cl <- parallel::makeCluster(n_cores,outfile = 'feature_importance_log.txt')
  parallel::clusterExport(cl,list('coeff_matrix'),envir = environment())
  parallel::clusterExport(cl,list('matprod_rl'))
  print('Calculating Matrix product WJ')
  WJ <- abind::abind(pbapply::pblapply(cl = cl, X= Js,FUN = matprod_rl,m1 = coeff_matrix ),along = 3)
  parallel::stopCluster(cl)

  # v_t <- Matrix::t(predictor$tsvd_v)
  v <- Matrix::t(predictor$tsvd_v)

  browser()
  print('Calculating Matrix product WJV')
  f <- function(WJ_element){cross_cell_average_fi_c(WJ_element,v)}
  env <- new.env(parent = environment(feature_importance))
  env$v <- v
  environment(f) <- env
  cl <- parallel::makeCluster(n_cores,outfile = 'feature_importance_log.txt')
  parallel::clusterEvalQ(cl,Rcpp::cppFunction('
  NumericMatrix matrix_product(NumericMatrix tm, NumericMatrix tm2) {
    // Convert R matrices to Eigen matrices
    const Eigen::Map<Eigen::MatrixXd> ttm(as<Eigen::Map<Eigen::MatrixXd>>(tm));
    const Eigen::Map<Eigen::MatrixXd> ttm2(as<Eigen::Map<Eigen::MatrixXd>>(tm2));

    // Perform matrix multiplication
    Eigen::MatrixXd prod = ttm * ttm2;

    // Return the result as a NumericMatrix
    return wrap(prod);
  }
',depends="RcppEigen")
  )
  parallel::clusterExport(cl,list('v'),envir = env)
  parallel::clusterExport(cl,list('cross_cell_average_fi_c'))
  WJV <- pbapply::pbapply(cl=cl, X= WJ, MARGIN = 1, FUN= f)
  parallel::stopCluster(cl)
  file.remove('feature_importance_log.txt')
  # Axis 1 = model, Axis 2 = Gene, Axis 3 = Cell --> taking mean 'across cells' = mean over margin of axis 1&2
  rownames(WJV) <- gexp_names
  colnames(WJV) <- names(predictor$lm_coefficients)

  return(WJV)
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
