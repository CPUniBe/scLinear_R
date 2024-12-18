
#TO CHANGE --> make normalization simulataneous for train and test set --> very different size factors calculated so counts 1 lead to very different normcounts
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



fit_predictor_R <- function(gexp_train,adt_train, gexp_test = NULL,
                            layer_gex = "counts", layer_adt = "counts",
                            normalize_gex = TRUE,normalize_adt = TRUE, margin = 2,
                            n_components = 300, zscore_relative_to_tsvd = 'after'){
  
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
   
  # Drop cells with total count (<)=0
  training_set <- training_set[Matrix::rowSums(training_set)>0,]
  
  # z-score normalization --> transpose after applying z-score is necessary to get THE SAME dimension as the input
  if(zscore_relative_to_tsvd == 'before'){
    training_set <- Matrix::t(apply(training_set, 1, function(x) {
      (x - mean(x)) / sd(x)
    }
    ))
  }
  
  # Create tSVD decomposition
  trained_tsvd <- sparsesvd::sparsesvd(training_set,rank = n_components)
  
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
  return(list(tsvd=trained_tsvd,lms=lm_results,zscore_relative_to_tsvd = zscore_relative_to_tsvd))
}

# Why set prediction to 0 if negative? We're trying to predict the CLR of ADT not ADT itself --> can be negative
predict_R <- function(predictor,gexp,layer="counts",normalize_gex=TRUE){
  
  if(any(class(gexp) %in% c("Seurat", "Assay", "Assay5"))){
    gexp <- Seurat::GetAssayData(gexp, layer = layer)
  }else{ # assume it is a matrix type
    gexp <- Matrix::Matrix(gexp, sparse = TRUE)
  }
  
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

#Discuss this function --> what metrics do we want and why
evaluate_predictor <- function(predictor,gexp_test,adt_test, normalize_gex = TRUE, normalize_adt = TRUE, margin = 2){
predicted_adt <- predict_R(predictor,gexp_test,normalize_gex = normalize_gex)
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