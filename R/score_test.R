
#' @importFrom Matrix Diagonal colSums
score_test_inner <- function(X, Z, res, mu) {
  W <- Diagonal(x=mu)
  XtW <- as.matrix(t(X) %*% W)
  cholNegH <- chol(solve(XtW %*% X))
  lhs <- colSums((Diagonal(x=sqrt(mu)) %*% Z)**2)
  rhs <- colSums((cholNegH %*% XtW %*% Z)**2)
  stat <- as.numeric(res %*% Z)**2 / (lhs - rhs)
  return(stat)
}
#' Score test for a fitted object.
#' Passing a BPPARAM and Z matrix (named) are required arguments.
#' The first argument can be a glmGamPoi object, DGEGLM object, or DESeqDataSet object.
#' 
#' @export
score_test <- function(X, ...) {
  UseMethod("score_test")
}

#' Score test for default (numeric)
#' 
#' @method score_test default
#' @export
score_test.default <- function(X, Z, Res, Mu, BPPARAM=bpparam()) {
  comm <- Reduce(intersect, list(rownames(X), rownames(Z), colnames(Res), colnames(Mu)))
  X <- X[comm,]
  Z <- Z[comm,]
  Res <- Res[,comm]
  Mu <- Mu[,comm]
  res <- bplapply(setNames(seq_len(nrow(Mu)), rownames(Mu)), function(i) {
    stat <- score_test_inner(X, Z, Res[i, ], Mu[i, ])
    return(stat)
  }, BPPARAM=BPPARAM)
  return(t(simplify2array(res)))
}

#' Score test for DESeqDataSet
#' 
#' @method score_test DESeqDataSet
#' @export
score_test.DESeqDataSet <- function(dds, Z, BPPARAM=bpparam()) {
### dds <- DESeq(DESeqFromMatrix(M, design=formula, colData=rowData(prom)), sfType="poscounts", fitType="glmGamPoi")
  score_test.default(X=DESeq2::design(dds),
                     Z=Z,
                     Res=SummarizedExperiment::assays(dds)$counts - SummarizedExperiment::assays(dds)$mu,
                     Mu=SummarizedExperiment::assays(dds)$mu,
                     BPPARAM=BPPARAM)
}

#' Score test for glmGamPoi
#' 
#' @method score_test glmGamPoi
#' @export
score_test.glmGamPoi <- function(fit, Z, BPPARAM=bpparam()) {
  return(score_test.default(X=fit$model_matrix,
                            Z=Z,
                            Res=residuals(fit, type="response"),
                            Mu=fit$Mu,
                            BPPARAM=BPPARAM))
}

#' Score test for DGEGLM
#' 
#' @method score_test DGEGLM
#' @export
score_test.DGEGLM <- function(fit, Z, BPPARAM=bpparam()) {
  return(score_test.default(X=fit$design,
                            Z=Z,
                            Res=fit$counts - fit$fitted.values,
                            Mu=fit$fitted.values,
                            BPPARAM=BPPARAM))
}

