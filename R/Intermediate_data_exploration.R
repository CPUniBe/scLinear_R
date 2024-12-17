pbmc10k <- readRDS("pbmc10k_prepared.rds")
## Create a training and a test set
set.seed(42)
# After running the new prediction function once, this block doesn't work anymore until session reload --> somehow overwrites some behaviour from before?
# --> the trigger that causes this seems to be running any function from Matrix or ExtraMatrix 
indx <- sample(1:length(colnames(pbmc10k)), size = length(colnames(pbmc10k)), replace = FALSE)
pbmc10k_train <- pbmc10k[,indx[1:5000]]
pbmc10k_test <- pbmc10k[,indx[5001:length(colnames(pbmc10k))]]



rpredictor_normed_z <- fit_predictor_R(gexp_train = pbmc10k_train@assays[["RNA"]],
                                   adt_train = pbmc10k_train@assays[["ADT"]],
                                   zscore_relative_to_tsvd = 'after'
                              )
r_predicted_adt_normed_z <- predict_R(rpredictor_normed_z,pbmc10k_test@assays[["RNA"]],zscore_relative_to_tsvd = 'after')

reticulate::use_condaenv("C:/Users/mz24b548/AppData/Local/miniconda3/envs/scLinearDependencies")
devtools::load_all("scLinear/")
pipe <- create_adt_predictor()
pipe <- fit_predictor(pipe = pipe,
                      gexp_train = pbmc10k_train@assays[["RNA"]],
                      adt_train = pbmc10k_train@assays[["ADT"]],
                      normalize_gex = TRUE,
                      normalize_adt = TRUE)

py_predicted_adt <- t(adt_predict(pipe = pipe,
            gexp = pbmc10k_test@assays[["RNA"]],
            normalize = TRUE)$data)


real_adt_normed <- t(Seurat::NormalizeData(pbmc10k_test@assays[["ADT"]], normalization.method = "CLR", margin = 2)$data)
real_adt_unnormed <- pbmc10k_test@assays[["ADT"]]['counts']

library(ggplot2)
plot_comparison <- function(adt1,adt2,xaxis = 'ADT Values 1', yaxis = 'ADT Values 2',title = ''){
  plot_dat <- cbind(adt1,adt2)
  colnames(plot_dat) <- c('ADT1','ADT2')
  comp_plot <- ggplot(plot_dat,aes(x=ADT1,y=ADT2))+
    geom_point()+
    geom_abline(intercept=0,slope=1,color = "red")+
    xlab(xaxis) + ylab(yaxis) +
    ggtitle(title)
  return(comp_plot)
}

r_vs_py <- list()
r_vs_real <- list()
py_vs_real <- list()
for(model in 1:ncol(r_predicted_adt_normed_z)){
  comparison_r_py <- plot_comparison(py_predicted_adt[,model],r_predicted_adt_normed_z[,model],xaxis = 'Prediction of Python model', yaxis = 'Prediction of R model',
                                title = names(rpredictor_normed_z$lms)[model])
  comparison_r_real <- plot_comparison(real_adt_normed[,model],r_predicted_adt_normed_z[,model],xaxis = 'Real (normed) ADT value', yaxis = 'Prediction of R model',
                                title = names(rpredictor_normed_z$lms)[model])
  comparison_py_real <- plot_comparison(real_adt_normed[,model],py_predicted_adt[,model],xaxis = 'Real (normed) ADT value', yaxis = 'Prediction of Python model',
                                        title = names(rpredictor_normed_z$lms)[model])
  
  r_vs_py <- append(r_vs_py,list(comparison_r_py))
  r_vs_real <- append(r_vs_real, list(comparison_r_real))
  py_vs_real <- append(py_vs_real, list(comparison_py_real))
}

# ggpubr::ggarrange(plotlist = r_vs_py, nrow = ceiling(sqrt(length(plotlist))), ncol = floor(sqrt(length(plotlist))))
# ggpubr::ggarrange(plotlist = r_vs_real, nrow = ceiling(sqrt(length(plotlist))), ncol = floor(sqrt(length(plotlist))))
# ggpubr::ggarrange(plotlist = py_vs_real, nrow = ceiling(sqrt(length(plotlist))), ncol = floor(sqrt(length(plotlist))))


gexp_fi <- pbmc10k_test@assays[["RNA"]]
slot_gex <- 'counts'
normalize_gex <- TRUE
predictor <- rpredictor_normed_z
# Feature importance --> first part the same as for predict_R --> project data and get model coefficients
if(any(class(gexp_fi) %in% c("Seurat", "Assay", "Assay5"))){
  gexp_fi <- Seurat::GetAssayData(gexp_fi, slot = slot_gex)
}else{ # assume it is a matrix type
  gexp_fi <- Matrix::Matrix(gexp_fi, sparse = TRUE)
}

if(normalize_gex){
  gexp_fi <- Matrix::t(gexp_normalize_modified(gexp_fi))
}else {gexp_fi <- Matrix::t(gexp_fi)}

if(zscore_relative_to_tsvd == 'before'){
  gexp_fi <- Matrix::t(apply(gexp_fi, 1, function(x) {
    (x - mean(x)) / sd(x)
  }
  ))  
}

gexp_fi_projected <- gexp_fi %*% predictor$tsvd$v

if(zscore_relative_to_tsvd == 'after'){
  gexp_fi_projected <- Matrix::t(apply(gexp_fi_projected, 1, function(x) {
    (x - mean(x)) / sd(x)
  }
  ))
}

coeff_matrix <- matrix(nrow = ncol(gexp_fi_projected)+1)
# Last coefficient usually 0 --> reason not quite clear but possibly because there is internal linear dependence
for(model in predictor$lms){
  coeff_matrix <- cbind(coeff_matrix,model$coefficients)
}
coeff_matrix[is.na(coeff_matrix)] <- 0
# Drop the empty column from the coeff matrix (used only to initialize the object)
coeff_matrix <- coeff_matrix[,2:ncol(coeff_matrix)]

# Drop the intercept --> not used for derivative d_prediction/d_gex --> constants drop from derivative
coeff_matrix_fi <- t(coeff_matrix[2:nrow(coeff_matrix),])

# input x = rowvector (corresponding to one cell) of tSVD projection matrix
cellwise_jacobian <- function(cell_projection){
jacobian <- numDeriv::jacobian(func = function(x) {(x - mean(x)) / sd(x)}, cell_projection)
}
# Slow but that's a looooot of matrix multiplications to run so probably to be expected
Js <- apply(gexp_fi_projected,1,cellwise_jacobian,simplify = FALSE)
# JV calculation too slow and object becomes so big that the space for it cannot be allocated
# JV <- lapply(Js_test,function(J){return(J %*% Matrix::t(predictor$tsvd$v))})
# Trying WJ first (resulting in 1741 17x300 matrices instead of 300x300 matrices)
WJ <- lapply(Js_test,function(J){return(coeff_matrix_fi %*% J)})
v_t <- Matrix::t(predictor$tsvd$v)
WJV <- abind::abind(lapply(WJ,function(WJ){return(WJ %*% v_t)}), along = 3)
feature_importance_means <- apply(WJV, c(1, 2), mean)
colnames(feature_importance_means) <- colnames(gexp_fi)
rownames(feature_importance_means) <- names(predictor$lms)

" Getting quite different results compared to python
There is very rough match in magnitudes i.e. the features which have < e-10 importance are generally set to 0 here
and the ranks of feature importance seem superficially similar but the actual values often differ by a factor of 10 or potentially much more
However, it is reassuring to note that the RNAs corresponding to each ADT protein are typically amongst the most important features
to predict the measured ADT level --> increases trust in the calculated feature importance for the other RNAs
"

# Depending on z-score after or skip the output is a simple array or a dgCMatrix


# I don't get the CLR back from the Seurat Normalization of ADT --> not quite sure WHAT it returns
# CLR not of counts but of scaled counts?
# Check if ADTs somewhat proportional to total RNA counts --> would allow us to skip a lot of scaling

evaluate_predictor(rpredictor_normed_z,pbmc10k_test@assays[["RNA"]],pbmc10k_test@assays[["ADT"]])


# evaluate_predictor(pipe = pipe,
#                    gexp_test = pbmc10k_test@assays[["RNA"]],
#                    adt_test = pbmc10k_test@assays[["ADT"]],
#                    normalize_gex = TRUE,
#                    normalize_adt = TRUE)