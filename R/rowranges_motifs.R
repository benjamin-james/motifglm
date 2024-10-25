
#' Match motifs against a generic RangedSummarizedExperiment
#'
#' @param se A RangedSummarizedExperiment object with a rowRanges to be counted against
#' @param JASPAR A JASPAR SQLite file such as JASPAR2024::JASPAR2024()@db
#' @param genome Genome string e.g. BSgenome.Hsapiens.UCSC.hg38
#' @param counts Whether to use counts or binary matches
#' @param species JASPAR species parameter. By default will be set to the BSgenome species
#' @param collection JASPAR collection
#' @param width Width for motifmatchr
#' @param cutoff P-value cutoff for motifmatchr
#' @param bg Background frequency for motifmatchr
#' @param kmer K-mer length to use for covariate adjustment
#' @param percent_var Cumulative percent variance captured of PCA of K-mers to use as a covariate
#' @return SummarizedExperiment object of genes by TF matrix of counts and associated metadata
#'
#' @importFrom BSgenome getBSgenome seqinfo
#' @importFrom S4Vectors metadata
#' @importFrom GenomicRanges GRanges seqinfo
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicFeatures makeTxDbFromGFF intronicParts
#' @importFrom IRanges IRanges
#' @importFrom Biostrings alphabetFrequency getSeq oligonucleotideFrequency
#' @importFrom SummarizedExperiment rowRanges rowData colData assays SummarizedExperiment
#' @export
rowRanges_motifs <- function(se, JASPAR, genome=NULL,
                             counts=TRUE,
                             species=NULL,
                             collection = "CORE",
                             width = 7L, cutoff = 5e-05, bg="even",
                             kmer = 2L, percent_var=99.) {
    if (is.null(genome)) {
        warning("BSgenome unspecified. Using BSgenome.Hsapiens.UCSC.hg38")
        genome <- "BSgenome.Hsapiens.UCSC.hg38"
    }
    bsg <- getBSgenome(genome)
    gr <- rowRanges(se)[as.character(seqnames(rowRanges(se))) %in% seqlevels(bsg)]
    seqlevels(gr) <- seqlevels(bsg)
    seqinfo(gr) <- seqinfo(bsg)
    mo <- motif_overlap(trim(gr), JASPAR, counts=counts, BSgenome=bsg,
                        species=species, collection=collection, width=width, cutoff=cutoff, bg=bg)
    af <- alphabetFrequency(getSeq(bsg, rowRanges(mo)), as.prob=TRUE)
    rowData(mo)$GC <- rowSums(af[,c("C", "G")])
    metadata(mo) <- list(genome=genome,
                                    species=species,
                                    collection=collection,
                                    width=width, cutoff=cutoff,
                                    bg=bg, 
                                    kmer=kmer, 
                                    covariates="GC")
  if (kmer > 0) {
    dn <- oligonucleotideFrequency(getSeq(bsg, rowRanges(mo)), kmer, as.prob=TRUE)
    pca <- prcomp(dn)
    frac_var <- cumsum(pca$sdev**2)/sum(colVars(dn))
    n_pcs <- sum(frac_var <= percent_var / 100.)
    for (i in seq_len(ncol(pca$x))) {
      rowData(mo)[[paste0("PC", i)]] <- pca$x[,i]
    }
    metadata(mo)$percent_var <- frac_var*100
    metadata(mo)$covariates <- paste0("PC", seq_len(n_pcs))
  }
  return(mo)

}
