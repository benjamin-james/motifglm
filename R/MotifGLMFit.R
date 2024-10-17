
#' Fit Motifs
#' @importFrom Matrix rowMeans
#' @importFrom S4Vectors metadata
#' @export
motif_fit <- function(motif_se, feat_mat=NULL, covariates=NULL, use_expected=TRUE, size_factors="poscounts", verbose=TRUE, model="glmGamPoi", BPPARAM=bpparam(), ...) {
  require(Matrix)
  require(sparseMatrixStats)
  require(SummarizedExperiment)
  comm_feat <- intersect(rownames(feat_mat), rownames(motif_se))
  stopifnot(length(comm_feat) > 1)
  M <- t(as.matrix(assays(motif_se)$counts[comm_feat,]))
  nzero <- sum(colSums(M) < 1)
  if (nzero > 0) {
    warning(paste0(nzero, " genes have zero total motif counts, discarding..."))
    comm_feat <- intersect(colnames(M)[colSums(M) > 0], comm_feat) ### Throw out zero-sum genes
    M <- M[,comm_feat]
  }
  if (is.null(covariates)) {
      covariates <- metadata(motif_se)$covariates
  }
  col_data <- as.data.frame(rowData(motif_se))[comm_feat, covariates, drop=FALSE]
  if (use_expected && (ncol(feat_mat) > 1)) {
    col_data$expected <- rowMeans(feat_mat[comm_feat,])
    covariates <- c(covariates, "expected")
  }
  if (model == "glmGamPoi") {
    model <- glmGamPoi::glm_gp(M, as.formula(paste0(c("~1", covariates), collapse="+")), col_data=col_data, size_factors=size_factors, verbose=verbose, ...)
  } else if (model == "DESeq2") {
    dds <- DESeq2::DESeqFromMatrix(M, design=as.formula(paste0(c("~1", covariates), collapse="+")), colData=col_data)
    model <- DESeq2::DESeq(dds, fitType="glmGamPoi", sfType=size_factors, BPPARAM=BPPARAM, parallel=TRUE, ...)
  } else if (model == "edgeR") {
    design <- model.matrix(as.formula(paste0(c("~1", covariates), collapse="+")), col_data)
    dgel <- edgeR::calcNormFactors(edgeR::DGEList(counts=M))
    dgel <- edgeR::estimateDisp(dgel, design)
    model <- edgeR::glmQLFit(dgel, model$model_matrix, ...)
  } else {
    stop("Unknown model name provided")
  }
  return(model)
}
                      
