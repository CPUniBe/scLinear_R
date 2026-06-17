

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
  return(rbind(gexp,artificial_gene_counts)[predictor$genes_considered,])
}

# Auxiliary function
# Get ensg to hgnc mapping
get_gene_map <- function() {
  #Where to store the mapping
  dir <- tools::R_user_dir("scLineaR", "cache")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  path <- file.path(dir, "ensg_hgnc.rds")
  #Read cached mapping table if it exists
  if (file.exists(path)) {
    return(readRDS(path))
  }
  # Otherwise - try to download
  res <- try({
    library(biomaRt)
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    df <- getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol"),
      mart = mart
    )
    saveRDS(df, path)
    df
  }, silent = TRUE)

  if (inherits(res, "try-error")) {
    stop("No cached gene map available and biomaRt failed.")
  }
 return(res)
}




# Similar to filter input genes but with mapping between ENSG IDs and HGNC gene symbols built in
match_genesets <- function(reference_names,newdat){
  # How to import !! for sym without loading the whole library?
  library(dplyr)
  gene_identifiers <- get_gene_map()

  # Drop empty mappings
  gene_identifiers <- gene_identifiers[gene_identifiers$ensembl_gene_id != '' & gene_identifiers$hgnc_symbol != '',]
  # Determine which id type was used for both datasets
  matches_gex_ensembl <- length(intersect(reference_names,gene_identifiers$ensembl_gene_id))
  matches_gex_hgnc <- length(intersect(reference_names,gene_identifiers$hgnc_symbol))
  matches_predictor_input_ensembl <- length(intersect(rownames(newdat),gene_identifiers$ensembl_gene_id))
  matches_predictor_input_hgnc <- length(intersect(rownames(newdat),gene_identifiers$hgnc_symbol))

  id_type_gex <- 'Undetermined'
  id_type_prediction_input <- 'Undetermined'

  if(matches_gex_ensembl > matches_gex_hgnc & matches_gex_ensembl > 0.3 * length(unique(reference_names))){
    id_type_gex <- 'ensembl_gene_id'
  }else if(matches_gex_ensembl < matches_gex_hgnc & matches_gex_hgnc > 0.3 * length(unique(reference_names))){
    id_type_gex <- 'hgnc_symbol'
  }
  if(matches_predictor_input_ensembl > matches_predictor_input_hgnc & matches_predictor_input_ensembl > 0.3 * length(unique(rownames(newdat)))){
    id_type_prediction_input <- 'ensembl_gene_id'
  }else if(matches_predictor_input_ensembl < matches_predictor_input_hgnc & matches_predictor_input_hgnc > 0.3 * length(unique(rownames(newdat)))){
    id_type_prediction_input <- 'hgnc_symbol'
  }

  if('Undetermined' %in% c(id_type_gex,id_type_prediction_input)){
    stop("Unable to determine ID type (ENSG / HGNC / other?) for genes used in gene counts either for the dataset used in model training
          or the one it is being applied to now. Manual translation between gene identifiers required")
  }
  if(id_type_gex != id_type_prediction_input){
    sym_id_type_gex <- sym(id_type_gex)
    sym_id_type_prediction_input <- sym(id_type_prediction_input)
    mapping_options <- gene_identifiers %>% filter(!!sym_id_type_gex %in% reference_names &
                                                   !!sym_id_type_prediction_input %in% rownames(newdat)
    )
    ambiguous <- which(table(mapping_options$hgnc_symbol) > 1)
    if(length(ambiguous) > 0){
      # Avoiding paste so print keeps clean data frame structure
      print('Ambiguous mappings between training data and new input for the following genes: ')
      print(mapping_options %>% filter(hgnc_symbol %in% names(ambiguous)))

      # ENSG IDs seem to be unique , HGNC ids can have different ENSG IDs
      # If training data used HGNC --> aggregate all ENSG ids from new count mat
      # If training data used ENSG ID --> unclear how to split counts --> drop
      # Haven't figured out how to aggregate rows without converting to non-sparse matrix or very time consuming functions
      # For now just drop all ambiguous mappings
      # if(id_type_prediction_input == 'hgnc_symbol'){
      mapping_options <- mapping_options %>% filter(!(hgnc_symbol %in% names(ambiguous)))
      # }
    }

    # Can't left join Csparse matrices and converting them would require tons of space --> handle names separately

    rownames_prediction_input <- data.frame(rownames(newdat))
    colnames(rownames_prediction_input)[1] <- id_type_prediction_input
    rownames_prediction_input <- rownames_prediction_input %>% left_join(mapping_options, by = join_by(!!sym_id_type_prediction_input == !!sym_id_type_prediction_input))
    rownames(newdat) <- rownames_prediction_input %>% pull(!!sym_id_type_gex)
    newdat <- newdat[!is.na(rownames(newdat)),]

  }
  to_add <- setdiff(reference_names,rownames(newdat))
  cat(paste0('Missing counts for ',length(to_add),' genes which were used for training of the ADT prediction model'))
  # Buffer missing gene counts with 0
  artificial_gene_counts <- Matrix::Matrix(0,nrow = length(to_add),ncol = ncol(newdat),
                                           dimnames = list(to_add,colnames(newdat)))
  # Add together original data and buffer 0 and arrange the rows in the same order as the genes were ordered in the training data
  newdat <- rbind(newdat,artificial_gene_counts)[reference_names,]
  return(newdat)
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
                            normalize_gex = TRUE,normalize_adt = TRUE, adt_norm_method = 'CLR', margin = 2,
                            n_components = 300, zscore_relative_to_tsvd = 'after',n_cores = NULL, reused_tsvd = NULL, isos = NULL){
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
    #adding a pseudocount of 0.99 so CLR transform of 0 values doesn't throw an error --> ok if counts overall are not super low
    #Using 0.99 instead of 1 to make it more easily distinguishable from true 1s during debugging
    #Better option would likely to ignore these ADT values entirely for training the LM
    if(adt_norm_method == 'CLR'){
    adt_train[adt_train == 0] <- 0.99
    # drop_mask = adt_train == 0
    # gmeans <- apply(adt_train,2,function(x){exp(mean(log(x[x>0])))})
    # Treat each 0 value as NA for  CLR BUT replace it afterwards
    # by the value it would have gotten if the count had been 1 (ignoring effect of the additional 1 on geometric mean)
    # replace <- log(1/gmeans)
    # replace <- matrix(rep(replace,nrow(adt_train)),nrow = nrow(adt_train),byrow = TRUE)
    adt_train <- Matrix::Matrix((t(compositions::clr(t(as.matrix(adt_train))))),sparse = TRUE)
    # adt_train <- Matrix::Matrix((t(compositions::ilr(t(as.matrix(adt_train))))),sparse = TRUE)
    # adt_train[adt_train == 0 & drop_mask] <- replace[adt_train == 0 & drop_mask]
    }else if(adt_norm_method == 'legacy'){
    adt_train <- Seurat::NormalizeData(adt_train, normalization.method = "CLR", margin = margin)
    }else if(adt_norm_method == 'dsb_no_empty'){
      adt_train <- dsb::ModelNegativeADTnorm(adt_train,isotype.control.name.vec = isos)
    }else if(adt_norm_method == 'log'){
      adt_train[adt_train == 0] <- 0.99
      adt_train <- log(adt_train)
    }
    else{
      Print('If normalize adt is set to TRUE (default), adt_norm_method has to be set to either CLR or legacy')
    }
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

  # When doing train test split to evaluate model performance and then retraining the model on full data
  # --> If gex_test was used in the original tsvd, there is no reason to recalculate it
  # Since tSVD calculation is one of the most time, consuming steps, include option to pass previously calculated svd projector
  if(!is.null(reused_tsvd)){
  print('Reusing tSVD projector matrix')
    projector <- reused_tsvd
  }else{
  # Create tSVD decomposition
  print('Calculating truncated singular value decomposition - for large input matrices this may take several minutes')
  trained_tsvd <- sparsesvd::sparsesvd(training_set,rank = n_components)
  print('tSVD done')
  projector <- trained_tsvd$v
  }
  # apply tSVD projection on input data
  training_set <- training_set %*% projector
  # z-score normalizaton

  if(zscore_relative_to_tsvd == 'after'){
    training_set <- Matrix::t(apply(training_set, 1, function(x) {
      (x - mean(x)) / sd(x)
    }
    ))
  }

  # If not zscore normalized, training set is not a matrix at this point --> take transformation out of the zscore block
  # For not explicit conversion here
  training_set <- as.matrix(training_set)


  # If lm-test-data was used for tsvd training step, reduce lm input back to only train data
  training_set <- training_set[intersect(colnames(gexp_train),rownames(training_set)),]

  # Filter out adt data for cells dropped due to RNA 0-counts
  # transposed to match expected input format for lm
  # here the indexing works rather quickly --> no conversion to array necessary
  adt_train_modelling <- Matrix::t(adt_train)[intersect(rownames(training_set),colnames(adt_train)),]

  rm(gexp_train,gexp_test)
  print('Fitting linear models')
  cl <- parallel::makeCluster(n_cores)
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
  names(results) <- colnames(adt_train_modelling)
  return(list(tsvd_v=projector,lm_coefficients=results,zscore_relative_to_tsvd = zscore_relative_to_tsvd,
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
predict <- function(predictor,gexp,layer="counts",normalize_gex=TRUE,convert = 'None'){
  # TODO: Make it properly part of package requirements
  library(biomaRt)
  library(dplyr)
  if(any(class(gexp) %in% c("Seurat", "Assay", "Assay5"))){
    gexp <- Seurat::GetAssayData(gexp, layer = layer)
  }else{ # assume it is a matrix type
    gexp <- Matrix::Matrix(gexp, sparse = TRUE)
  }

  #Normalization first based on all genes then filter might be more accurate but takes much longer

  # gexp <- filter_input_genes(gexp,predictor)
  gexp <- match_genesets(predictor$genes_considered,gexp)

  if(normalize_gex){
    gexp <- gexp_normalize(gexp)
  }

  # Bring into shape cells x genes (done during training in the step which combines train & test)
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
evaluate_predictor <- function(predictor,gexp_test,adt_test,gexp_layer = 'counts', adt_layer = 'counts', normalize_gex = TRUE,
                               normalize_adt = TRUE, adt_norm_method = 'CLR', margin = 2, convert = 'None', isos = NULL){
  # ToFix --> can't compare CLR if different subsets were considered since it's compositional data
  predicted_adt <- predict(predictor,gexp_test,layer = gexp_layer, normalize_gex = normalize_gex, convert = convert)
  # predicted_adt <- predict_alternative(predictor,gexp_test)
  if(class(adt_test)[1] == "Assay" |class(adt_test)[1] == "Assay5"){ adt_test <- Seurat::GetAssayData(adt_test, layer = adt_layer) }
  if(normalize_adt){
    #adding a pseudocount of 0.99 so CLR transform of 0 values doesn't throw an error --> ok if counts overall are not super low
    #Using 0.99 instead of 1 to make it more easily distinguishable from true 1s during debugging
    #Better option would likely to ignore these ADT values entirely for training the LM
    if(adt_norm_method == 'CLR'){
      adt_test[adt_test == 0] <- 0.99
      adt_test <- Matrix::Matrix(((compositions::clr(t(as.matrix(adt_test))))),sparse = TRUE)
    }else if(adt_norm_method == 'legacy'){
      adt_test <- Matrix::t(Seurat::NormalizeData(adt_test, normalization.method = "CLR", margin = margin))
    }else if(adt_norm_method == 'dsb_no_empty'){
      adt_test <- Matrix::t(dsb::ModelNegativeADTnorm(adt_test,isotype.control.name.vec = isos))
    }else if(adt_norm_method == 'log'){
      adt_test[adt_test == 0] <- 0.99
      adt_test <- Matrix::t(log(adt_test))
    }
    else{
      Print('If normalize adt is set to TRUE (default), adt_norm_method has to be set to either CLR or legacy')
    }
  }else{
      # Converting from Matrix::Matrix to base Matrix (and transposing)
      adt_test <- as.matrix(Matrix::t(adt_test))
  }
  #What's the subsetting for? Protein filter or cell filter?
  # p_adt <- subset(predicted_adt,subset = rownames(predicted_adt) %in% rownames(adt_test))
  # t_adt <- (subset(as.matrix(adt_test),subset = rownames(adt_test) %in% rownames(predicted_adt)))
  p_adt <- predicted_adt
  t_adt <- as.matrix(adt_test)
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
  # Slow but that's a looooot of matrix multiplications to run so probably to be expected

  N <- ncol(gexp_projected)
  # cl <- parallel::makeCluster(n_cores,outfile = 'feature_importance_log.txt')
  # parallel::clusterExport(cl,list('cellwise_jacobian'))
  # print('Calculating Jacobian of z-score normalization step')
  #
  # Js <- pbapply::pbapply(cl=cl,X = gexp_projected, MARGIN = 1,FUN = cellwise_jacobian,simplify = FALSE)
  # parallel::stopCluster(cl)
  Js <- cluster_call_cellwise_jacobian(gexp_projected,n_cores)

  #transposing so we can use crossprod for apply function --> circumvents having to do FUN = function(X)
  # ... which slows down parallelization when called from environment other than global
  # cl <- parallel::makeCluster(n_cores,outfile = 'feature_importance_log.txt')
  # parallel::clusterExport(cl,list('coeff_matrix'),envir = environment())
  # parallel::clusterExport(cl,list('matprod_rl'))
  # print('Calculating Matrix product WJ')
  # WJ <- abind::abind(pbapply::pblapply(cl = cl, X= Js,FUN = matprod_rl,m1 = coeff_matrix ),along = 3)
  # parallel::stopCluster(cl)
  WJ <- cluster_call_WJ(Js,coeff_matrix, n_cores)

  print('Calculating Matrix product WJV')
  v <- Matrix::t(predictor$tsvd_v)
  # f <- function(WJ_element){cross_cell_average_fi_c(WJ_element,v)}
  # env <- new.env(parent = environment(feature_importance))
  # env$v <- v
  # environment(f) <- env
  # cl <- parallel::makeCluster(n_cores)
  # # The parallel workers need to know where to look for the package as they don't automatically inherit the
  # # library of the calling environment
  # libpaths <- .libPaths()
  # parallel::clusterExport(cl, "libpaths")
  # parallel::clusterEvalQ(cl,{
  #   .libPaths(libpaths)
  #   library(scLineaR)
  # }
  #   )
  # # parallel::clusterEvalQ(cl,Rcpp::sourceCpp('C:/Users/mz24b548/Documents/GitRepos_local/scLinear_R/src/matrix_product.cpp'))
  # parallel::clusterExport(cl,list('v'),envir = env)
  # parallel::clusterExport(cl,list('cross_cell_average_fi_c'))
  # WJV <- pbapply::pbapply(cl=cl, X= WJ, MARGIN = 1, FUN= f)
  # parallel::stopCluster(cl)
  WJV <- cluster_call_WJV(WJ,v,n_cores)
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
  #For some reason the computeSumFactors function throws an error if input id dgRMatrix
  if(class(gexp_matrix) == 'dgRMatrix'){gexp_matrix <- as(gexp_matrix, "CsparseMatrix")}

  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = gexp_matrix))
  clusters <- scran::quickCluster(sce)
  sce <- scran::computeSumFactors(sce, clusters=clusters)
  sce <- scuttle::logNormCounts(sce, center.size.factors = center.size.factors, log = log, ...)
  gexp_matrix <- sce@assays@data@listData[["normcounts"]]
  gexp_matrix <- base::log1p(gexp_matrix)
  return(gexp_matrix)
}




# fit_predictor_alternative <- function(gexp_train,adt_train, gexp_test = NULL,
#                           layer_gex = "counts", layer_adt = "counts",
#                           normalize_gex = TRUE,normalize_adt = TRUE, adt_norm_method = 'CLR', margin = 2,
#                           n_components = 300, zscore_relative_to_tsvd = 'after',n_cores = NULL, reused_tsvd = NULL){
#   # If no number is specified, get number of available cores and omit 4 to avoid using up all system resources
#   if(is.null(n_cores)){n_cores <- parallelly::availableCores(omit = 4)}
#   # If objects are passed as matrices, leave as is. If passed as Seurat objects --> extract count data for assays
#   if(class(gexp_train)[1] == "Assay" | class(gexp_train)[1] == "Assay5"){ gexp_train <- Seurat::GetAssayData(gexp_train, layer = layer_gex) }
#   if(class(adt_train)[1] == "Assay" | class(adt_train)[1] == "Assay5"){ adt_train <- Seurat::GetAssayData(adt_train, layer = layer_adt) }
#   if(!is.null(gexp_test)){
#     if(class(gexp_test)[1] == "Assay" |class(gexp_test)[1] == "Assay5"){ gexp_test <- Seurat::GetAssayData(gexp_test, layer = layer_gex) }
#   }
#
#   #TODO: Change it so normalization happens for both train and test at the same time
#   # if(normalize_gex){
#   #   gexp_train <- gexp_normalize(gexp_train)
#   #   if( !is.null(gexp_test)){gexp_test <- gexp_normalize(gexp_test)}
#   # }
#   if(normalize_adt){
#     #adding a pseudocount of 0.99 so CLR transform of 0 values doesn't throw an error --> ok if counts overall are not super low
#     #Using 0.99 instead of 1 to make it more easily distinguishable from true 1s during debugging
#     #Better option would likely to ignore these ADT values entirely for training the LM
#     if(adt_norm_method == 'CLR'){
#       adt_train[adt_train == 0] <- 0.99
#       adt_train <- Matrix::Matrix((t(compositions::clr(t(as.matrix(adt_train))))),sparse = TRUE)
#     }else if(adt_norm_method == 'legacy'){
#       adt_train <- Seurat::NormalizeData(adt_train, normalization.method = "CLR", margin = margin)
#     }else{
#       Print('If normalize adt is set to TRUE (default), adt_norm_method has to be set to either CLR or legacy')
#     }
#   }
#
#   # The way the sets are generated, same gene names should be a given (stem from same dataset just split)
#   if(!is.null(gexp_test)){
#     training_set <- cbind(gexp_train,gexp_test)
#   }else{training_set <- gexp_train}
#
#   gex_trainset_seurat <- Seurat::CreateSeuratObject(counts = training_set, assay = 'RNA', min.cells = 1, min.features = 1)
#   gex_trainset_seurat <- NormalizeData(gex_trainset_seurat)
#   gex_trainset_seurat <- FindVariableFeatures(gex_trainset_seurat)
#   gex_trainset_seurat <- ScaleData(gex_trainset_seurat)
#   gex_trainset_seurat <- RunPCA(gex_trainset_seurat, features = VariableFeatures(object = gex_trainset_seurat),
#                                 return.model = TRUE, npcs = n_components)
#
#   # z-score normalization --> transpose after applying z-score is necessary to get THE SAME dimension as the input
#   # if(zscore_relative_to_tsvd == 'before'){
#   #   training_set <- Matrix::t(apply(training_set, 1, function(x) {
#   #     (x - mean(x)) / sd(x)
#   #   }
#   #   ))
#   # }
#
#   # When doing train test split to evaluate model performance and then retraining the model on full data
#   # --> If gex_test was used in the original tsvd, there is no reason to recalculate it
#   # Since tSVD calculation is one of the most time, consuming steps, include option to pass previously calculated svd projector
#   # if(!is.null(reused_tsvd)){
#   #   print('Reusing tSVD projector matrix')
#   #   projector <- reused_tsvd
#   # }else{
#   #   # Create tSVD decomposition
#   #   print('Calculating truncated singular value decomposition - for large input matrices this may take several minutes')
#   #   trained_tsvd <- sparsesvd::sparsesvd(training_set,rank = n_components)
#   #   print('tSVD done')
#   #   projector <- trained_tsvd$v
#   # }
#   # apply tSVD projection on input data
#   training_set <- gex_trainset_seurat@reductions$pca@cell.embeddings
#   # z-score normalizaton
#
#   # if(zscore_relative_to_tsvd == 'after'){
#   #   training_set <- Matrix::t(apply(training_set, 1, function(x) {
#   #     (x - mean(x)) / sd(x)
#   #   }
#   #   ))
#   # }
#
#   # If not zscore normalized, training set is not a matrix at this point --> take transformation out of the zscore block
#   # For not explicit conversion here
#   # training_set <- as.matrix(training_set)
#
#
#   # If lm-test-data was used for tsvd training step, reduce lm input back to only train data
#   training_set <- training_set[intersect(colnames(gexp_train),rownames(training_set)),]
#
#   # Filter out adt data for cells dropped due to RNA 0-counts
#   # transposed to match expected input format for lm
#   # here the indexing works rather quickly --> no conversion to array necessary
#   adt_train_modelling <- Matrix::t(adt_train)[intersect(Cells(gex_trainset_seurat), colnames(adt_train)),]
#
#
#
#   rm(gexp_train,gexp_test,adt_train)
#   print('Fitting linear models')
#
#   # cl <- parallel::makeCluster(n_cores,outfile = 'predictor_fit_alt_log.txt')
#   # parallel::clusterExport(cl,list("training_set","adt_train_modelling"),envir = environment())
#   # parallel::clusterEvalQ(cl,library(Matrix))
#   # results <- pbapply::pblapply(cl = cl, X = 1:ncol(adt_train_modelling), FUN = function(i) {
#   #   tryCatch({
#   #     return(lm(adt_train_modelling[, i] ~ training_set)$coefficients)  # Fit a linear model
#   #   }, error = function(e) {
#   #     message("Error fitting model for response ", i, ": ", e$message)
#   #     return(NULL)  # Return NULL if there's an error
#   #   })
#   # })
#   # parallel::stopCluster(cl)
#   results <- cluster_call_model_training(adt_train_modelling,training_set,n_cores)
#   names(results) <- colnames(adt_train_modelling)
#   return(list(traindat = gex_trainset_seurat, lm_coefficients = results,
#               zscore_relative_to_tsvd = zscore_relative_to_tsvd))
# }
#
# predict_alternative <- function(train_model,gexp,layer="counts"){
#   # TODO: Make it properly part of package requirements
#   library(dplyr)
#   if(any(class(gexp) %in% c("Seurat", "Assay", "Assay5"))){
#     gexp <- Seurat::GetAssayData(gexp, layer = layer)
#   }else{ # assume it is a matrix type
#     gexp <- Matrix::Matrix(gexp, sparse = TRUE)
#   }
#
#   gex <- match_genesets(Seurat::GetAssayData(train_model$traindat, layer = 'counts'),gexp)
#   gex <- Seurat::CreateSeuratObject(counts = gex, assay = 'RNA')
#   gex <- NormalizeData(gex)
#   gex <- FindVariableFeatures(gex)
#   gex <- ScaleData(gex)
#   IAnchors <- Seurat::FindTransferAnchors(reference = train_model$traindat,query = gex,
#                                           reduction = 'pcaproject', reference.reduction = 'pca',
#                                           dims = 1:ncol(train_model$traindat@reductions$pca@cell.embeddings))
#   gex <- Seurat::IntegrateEmbeddings(IAnchors, reference =  train_model$traindat, query =  gex)
#   coeff_matrix <- matrix(nrow = ncol(gex@reductions$integrated_dr)+1)
#   # Last coefficient usually 0 --> reason not quite clear. Does z-score normalization introduce linear dependence since mean must be 0
#   for(model in train_model$lm_coefficients){
#     coeff_matrix <- cbind(coeff_matrix,model)
#   }
#   coeff_matrix[is.na(coeff_matrix)] <- 0
#   # Drop the empty column from the coeff matrix (used only to initialize the object)
#   coeff_matrix <- coeff_matrix[,2:ncol(coeff_matrix)]
#   # Add column of 1s 'to the left' of the tSVD projection of the test data --> adds intercept of each model (intercept is the first coefficient of the lm)
#   lm_input <- cbind(rep(1,nrow(gex@reductions$integrated_dr@cell.embeddings)),gex@reductions$integrated_dr@cell.embeddings)
#   # Multiply tSVD projected & normed input data with LM coefficients
#   res <- lm_input %*% coeff_matrix
#   colnames(res) <- names(train_model$lm_coefficients)
#   return(res)
# }

#Separate function in the hopes that this will lead to only the variable in THIS cluster's environment being passed to each node
cluster_call_model_training <- function(adt_train_modelling,training_set,n_cores){
  cl <- parallel::makeCluster(n_cores,outfile = 'predictor_fit_alt_log.txt')
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
  return(results)
}

cluster_call_cellwise_jacobian <- function(gexp_projected,n_cores){
  cl <- parallel::makeCluster(n_cores,outfile = 'feature_importance_log.txt')
  parallel::clusterExport(cl,list('cellwise_jacobian'))
  print('Calculating Jacobian of z-score normalization step')
  Js <- pbapply::pbapply(cl=cl,X = gexp_projected, MARGIN = 1,FUN = cellwise_jacobian,simplify = FALSE)
  parallel::stopCluster(cl)
  return(Js)
}

cluster_call_WJ <- function(Js, coeff_matrix, n_cores){
cl <- parallel::makeCluster(n_cores,outfile = 'feature_importance_log.txt')
parallel::clusterExport(cl,list('coeff_matrix'),envir = environment())
parallel::clusterExport(cl,list('matprod_rl'))
print('Calculating Matrix product WJ')
WJ <- abind::abind(pbapply::pblapply(cl = cl, X= Js,FUN = matprod_rl,m1 = coeff_matrix ),along = 3)
parallel::stopCluster(cl)
return(WJ)
}

cluster_call_WJV <- function(WJ,v,n_cores){
  f <- function(WJ_element){cross_cell_average_fi_c(WJ_element,v)}
  # env <- new.env(parent = environment(cluster_call_WJV))
  # env$v <- v
  # environment(f) <- env
  cl <- parallel::makeCluster(n_cores,outfile = 'feature_importance_log.txt')

  # parallel::clusterEvalQ(cl,library(scLineaR))
  libpaths <- .libPaths()
  parallel::clusterExport(cl,list('v','cross_cell_average_fi_c','libpaths'))
  # parallel::clusterExport(cl,list('cross_cell_average_fi_c'))
  # parallel::clusterExport(cl, "libpaths")
  parallel::clusterEvalQ(cl,{
    .libPaths(libpaths)
    library(scLineaR)
  }
    )
  WJV <- pbapply::pbapply(cl=cl, X= WJ, MARGIN = 1, FUN= f)
  parallel::stopCluster(cl)
  return(WJV)
}
