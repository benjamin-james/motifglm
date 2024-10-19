
#' Match motifs against an annotation
#'
#' @param JASPAR A JASPAR SQLite file such as JASPAR2024::JASPAR2024()@db
#' @param gff GFF or GTF file to be read in by rtracklayer::readGFF
#' @param genome Genome string e.g. BSgenome.Hsapiens.UCSC.hg38
#' @param upstream Integer upstream
#' @param downstream Downstream base pairs of gene
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
#' @importFrom Matrix Diagonal sparseMatrix
#' @importFrom BSgenome getBSgenome seqinfo
#' @importFrom rtracklayer readGFF
#' @importFrom S4Vectors metadata
#' @importFrom GenomicRanges seqnames GRanges seqinfo trim promoters
#' @importFrom GenomicFeatures makeTxDbFromGFF intronicParts
#' @importFrom IRanges IRanges
#' @importFrom Biostrings alphabetFrequency getSeq oligonucleotideFrequency
#' @importFrom SummarizedExperiment rowRanges rowData colData assays
#' @export
intron_motifs <- function(JASPAR, gff,
                          genome=NULL,
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
  if (is.null(species)) {
    species <- S4Vectors::metadata(bsg)$organism
  }
  txdb <- makeTxDbFromGFF(gff, organism=species)
  introns <- intronicParts(txdb, linked.to.single.gene.only=TRUE)
  introns <- introns[as.character(seqnames(introns)) %in% seqlevels(bsg)]
  names(introns) <- make.unique(paste0(introns$gene_id, "-", introns$intronic_part), sep="#")
  seqlevels(introns) <- seqlevels(bsg)
  seqinfo(introns) <- seqinfo(bsg)
  mo <- motif_overlap(introns, JASPAR, counts=counts, BSgenome=bsg,
                      species=species, collection=collection, width=width, cutoff=cutoff, bg=bg)
  gff <- as.data.frame(rtracklayer::readGFF(gff))
  gff <- gff[gff$seqid %in% seqnames(bsg),]
  gene_gr <- with(gff[gff$type == "gene",], GRanges(seqnames=seqid,
                                                    ranges=IRanges(start, end),
                                                    gene_name=make.unique(gene_name, sep="-"),
                                                    gene_id=gene_id,
                                                    strand=strand,
                                                    gene_type=gene_type,
                                                    gene_length=end-start,
                                                    seqinfo=seqinfo(bsg)))
  AF <- alphabetFrequency(getSeq(bsg, rowRanges(mo)),
                          as.prob=FALSE)
  M <- sparseMatrix(i=match(rownames(mo), gene_gr$gene_id),
                    j=seq_len(nrow(mo)),
                    x=rep(1L, nrow(mo)),
                    dimnames=list(gene_gr$gene_name, rownames(mo)),
                    dims=c(length(gene_gr), nrow(mo)))
  AFM <- as.matrix(M %*% AF)
  gene_gr$GC <- rowSums(AFM[,c("C", "G")])/(1e-300 + rowSums(AFM))
  gene_gr$introns <- as.numeric(M %*% rep(1L, nrow(mo)))
  names(gene_gr) <- gene_gr$gene_name
  go <- SummarizedExperiment(list(counts=M %*% assays(mo)$counts),
                             metadata=list(genome=genome, species=species,
                                           collection=collection,
                                           width=width, cutoff=cutoff, bg=bg,
                                           gff=setNames(metadata(txdb)$value, metadata(txdb)$name)[["Data source"]],
                                           kmer=kmer, covariates=c("CG")),
                             colData=colData(mo),
                             rowData=gene_gr)
  if (kmer > 0) {
    dn_intron <- oligonucleotideFrequency(getSeq(bsg, rowRanges(mo)), kmer, as.prob=FALSE)
    dn_pergene <- as.matrix(M %*% dn_intron)
    dn <- as.matrix(Diagonal(x=1/(1e-300 + rowSums(dn_pergene))) %*% dn_pergene)
    pca <- prcomp(dn)
    for (i in seq_len(ncol(pca$x))) {
      rowData(go)[[paste0("PC", i)]] <- pca$x[,i]
    }
    frac_var <- cumsum(pca$sdev**2)/sum(colVars(dn))
    n_pcs <- sum(frac_var <= percent_var / 100.)
    metadata(go)$percent_var <- frac_var*100
    metadata(go)$covariates <- paste0("PC", seq_len(n_pcs))
  }
  return(go)
}

