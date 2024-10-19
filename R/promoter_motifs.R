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
#' @importFrom BSgenome getBSgenome
#' @importFrom rtracklayer readGFF
#' @importFrom S4Vectors metadata
#' @importFrom GenomicRanges seqnames GRanges seqinfo trim promoters
#' @importFrom IRanges IRanges
#' @importFrom Biostrings alphabetFrequency getSeq oligonucleotideFrequency
#' @importFrom SummarizedExperiment rowRanges rowData
#' @export
promoter_motifs <- function(JASPAR, gff,
                            genome=NULL,
                            upstream=2000L, downstream=200L, counts=TRUE,
                            species=NULL,
                            collection = "CORE",
                            width = 7L, cutoff = 5e-05, bg="even",
                            kmer = 2L, percent_var=99.) {
  if (is.null(genome)) {
    warning("BSgenome unspecified. Using BSgenome.Hsapiens.UCSC.hg38")
    genome <- "BSgenome.Hsapiens.UCSC.hg38"
  }
  bsg <- BSgenome::getBSgenome(genome)
  if (is.null(species)) {
      species <- S4Vectors::metadata(bsg)$organism
  }
  gff_file <- normalizePath(gff)
  gff <- as.data.frame(rtracklayer::readGFF(gff))
  gff <- gff[gff$type == "gene",]
  gff$gene_name <- make.unique(gff$gene_name, sep="-")
  gr <- with(gff[(gff$seqid %in% GenomicRanges::seqnames(bsg)),],
             GenomicRanges::GRanges(seqnames=seqid,
                                    ranges=IRanges::IRanges(start, end),
                                    gene_name=gene_name,
                                    gene_id=gene_id,
                                    strand=strand,
                                    gene_type=gene_type,
                                    gene_length=end-start,
                                    seqinfo=GenomicRanges::seqinfo(bsg)))
  names(gr) <- gr$gene_name
  prom.gr <- GenomicRanges::trim(GenomicRanges::promoters(gr, upstream=upstream, downstream=downstream))
  mo <- motif_overlap(prom.gr, JASPAR, counts=counts, BSgenome=bsg,
                      species=species, collection=collection, width=width, cutoff=cutoff, bg=bg)
  af <- Biostrings::alphabetFrequency(Biostrings::getSeq(bsg, SummarizedExperiment::rowRanges(mo)), as.prob=TRUE)
  SummarizedExperiment::rowData(mo)$GC <- rowSums(af[,c("C", "G")])
  S4Vectors::metadata(mo) <- list(genome=genome,
                                  species=species,
                                  collection=collection,
                                  width=width, cutoff=cutoff,
                                  bg=bg, gff=attr(gff, "pragmas"),
                                  upstream=upstream, downstream=downstream,
                                  kmer=kmer, 
                                  covariates="GC")
  if (kmer > 0) {
    dn <- Biostrings::oligonucleotideFrequency(Biostrings::getSeq(bsg, SummarizedExperiment::rowRanges(mo)), kmer, as.prob=TRUE)
    pca <- prcomp(dn)
    frac_var <- cumsum(pca$sdev**2)/sum(colVars(dn))
    n_pcs <- sum(frac_var <= percent_var / 100.)
    for (i in seq_len(ncol(pca$x))) {
      SummarizedExperiment::rowData(mo)[[paste0("PC", i)]] <- pca$x[,i]
    }
    S4Vectors::metadata(mo)$percent_var <- frac_var*100
    S4Vectors::metadata(mo)$covariates <- paste0("PC", seq_len(n_pcs))
  }
  return(mo)
}

