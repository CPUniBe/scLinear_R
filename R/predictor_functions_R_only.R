
##TO CHANGE --> make normalization simulataneous for train and test set --> very different size factors calculated so counts 1 lead to very different normcounts
gexp_normalize_modified <- function(gexp_matrix, center.size.factors = FALSE, log = FALSE, ...){
  ## normalize data GEX
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = gexp_matrix))
  clusters <- scran::quickCluster(sce)
  sce <- scran::computeSumFactors(sce, clusters=clusters)
  sce <- scuttle::logNormCounts(sce, center.size.factors = center.size.factors, log = log, ...)
  gexp_matrix <- sce@assays@data@listData[["normcounts"]]
  # Comment the log1p since the logNormCounts already does a log-transformation WITH adding a pseudocount of 1
  # Just for testing put it back to see if rest of pipeline works
  gexp_matrix <- base::log1p(gexp_matrix)
  return(gexp_matrix)
}




fit_predictor_R <- function(gexp_train,adt_train, gexp_test = NULL,
                            slot_gex = "counts", slot_adt = "counts",
                            normalize_gex = TRUE,normalize_adt = TRUE, margin = 2,
                            n_components = 300, zscore_relative_to_tsvd = 'skip'){
  
  # If objects are passed as matrices, leave as is. If passed as Seurat objects --> extract count data for assays
  if(class(gexp_train)[1] == "Assay" | class(gexp_train)[1] == "Assay5"){ gexp_train <- Seurat::GetAssayData(gexp_train, slot = slot_gex) }
  if(class(adt_train)[1] == "Assay" | class(adt_train)[1] == "Assay5"){ adt_train <- Seurat::GetAssayData(adt_train, slot = slot_adt) }
  if(!is.null(gexp_test)){
    if(class(gexp_test)[1] == "Assay" |class(gexp_test)[1] == "Assay5"){ gexp_test <- Seurat::GetAssayData(gexp_test, slot = slot_gex) }
  }

  #TODO: Change it so normalization happens for both train and test at the same time 
  if(normalize_gex){
    gexp_train <- gexp_normalize_modified(gexp_train)
    if( !is.null(gexp_test)){gexp_test <- gexp_normalize_modified(gexp_test)}
  }
  if(normalize_adt){
    adt_train <- Seurat::NormalizeData(adt_train, normalization.method = "CLR", margin = margin)
  }
  
  # The way the sets are generated, same gene names should be a given (stem from same dataset just split)
  if(!is.null(gexp_test)){
  training_set <- Matrix::t(cbind(gexp_train,gexp_test))
  }else{training_set <- Matrix::t(gexp_train)}
   
  # Drop cells with total count (<)=0
  training_set <- training_set[Matrix::rowSums(training_set)>0,]
  
  if(zscore_relative_to_tsvd == 'before'){
    training_set <- Matrix::t(apply(training_set, 1, function(x) {
      (x - mean(x)) / sd(x)
    }
    ))
  }
  
  # Create tSVD decomposition
  trained_tsvd <- sparsesvd::sparsesvd(training_set,rank = 300)
  
  # apply tSVD projection on input data
  training_set_redux <- training_set %*% trained_tsvd$v
  # z-score normalizaton
  
  if(zscore_relative_to_tsvd == 'after'){
  training_set_redux <- Matrix::t(apply(training_set_redux, 1, function(x) {
         (x - mean(x)) / sd(x)
     }
     ))
  }
  
  #transposed to match expected input format for lm
  redux_train_modelling <- as.array(training_set_redux)
  adt_train_modelling <- Matrix::t(adt_train)
  
  # Fit linear models for all ADTs
  lm_results <- apply(adt_train_modelling, 2, function(y) lm(y ~ redux_train_modelling))
  return(list(tsvd=trained_tsvd,lms=lm_results))
}

# Why set prediction to 0 if negative? We're trying to predict the CLR of ADT not ADT itself --> can be negative
predict_R <- function(predictor,gexp,slot="counts",normalize_gex=TRUE,zscore_relative_to_tsvd = 'skip'){
  
  if(any(class(gexp) %in% c("Seurat", "Assay", "Assay5"))){
    gexp <- Seurat::GetAssayData(gexp, slot = slot)
  }else{ # assume it is a matrix type
    gexp <- Matrix::Matrix(gexp, sparse = TRUE)
  }
  
  if(normalize_gex){
    gexp <- Matrix::t(gexp_normalize_modified(gexp))
  }else {gexp <- Matrix::t(gexp)}
  
  if(zscore_relative_to_tsvd == 'before'){
    gexp <- Matrix::t(apply(gexp, 1, function(x) {
      (x - mean(x)) / sd(x)
    }
    ))  
  }
  
  gexp_projected <- gexp %*% predictor$tsvd$v
  
  if(zscore_relative_to_tsvd == 'after'){
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

#Discuss this function --> what metrics do we want and why
evaluate_predictor <- function(predictor,gexp_test,adt_test, normalize_gex = TRUE, normalize_adt = TRUE, margin = 2){
predicted_adt <- predict_R(predictor,gexp_test,normalize_gex = normalize_gex,zscore_relative_to_tsvd = 'after')
if(class(adt_test)[1] == "Assay" |class(adt_test)[1] == "Assay5"){ adt_test <- Seurat::GetAssayData(adt_test, slot = 'counts') }
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