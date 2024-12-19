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
r_predicted_adt_normed_z <- predict_R(rpredictor_normed_z,pbmc10k_test@assays[["RNA"]])

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

# ggpubr::ggarrange(plotlist = r_vs_py, nrow = ceiling(sqrt(length(r_vs_py))), ncol = floor(sqrt(length(r_vs_py))))
# ggpubr::ggarrange(plotlist = r_vs_real, nrow = ceiling(sqrt(length(r_vs_real))), ncol = floor(sqrt(length(r_vs_real))))
# ggpubr::ggarrange(plotlist = py_vs_real, nrow = ceiling(sqrt(length(py_vs_real))), ncol = floor(sqrt(length(py_vs_real))))

feature_imp <- feature_importance(rpredictor_normed_z,pbmc10k_test@assays[["RNA"]],layer_gexp = 'counts')
# intersect(row.names(fi),colnames(fi)) --> only CD4, CD14, CD19 & TIGIT have a directly corresponding RNA measured in the dataset
# The RNA is generally within the top 10-20 most important genes for predicting protein abundance (as expected) --> for CD19 it's only at rank 59th
# Still seems to indicate that the process 'works'


evaluate_predictor(rpredictor_normed_z,pbmc10k_test@assays[["RNA"]],pbmc10k_test@assays[["ADT"]])


# evaluate_predictor(pipe = pipe,
#                    gexp_test = pbmc10k_test@assays[["RNA"]],
#                    adt_test = pbmc10k_test@assays[["ADT"]],
#                    normalize_gex = TRUE,
#                    normalize_adt = TRUE)