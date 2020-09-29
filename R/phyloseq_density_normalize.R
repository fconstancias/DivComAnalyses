#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'library(phyloseq)
#'data("GlobalPatterns")
#'sample_data(GlobalPatterns)$norm = sample(5000000:1000000, size=nsamples(GlobalPatterns), replace=TRUE)
#'phyloseq_density_normalize(GlobalPatterns, "norm") -> out
#'HTSSIP::OTU_qPCR_trans(GlobalPatterns, GlobalPatterns %>% sample_data() %>% data.frame() , "X.SampleID", value_idx = "norm") -> out2
#'otu_table(GlobalPatterns)["540305","TRRsed3"] #2305
#'sample_data(GlobalPatterns)["TRRsed3",] #2664122
#'sample_sums(GlobalPatterns)["TRRsed3"] #279704
#'(2305/279704)  * 2664122
#'otu_table(out)["540305","TRRsed3"] #
#'otu_table(out2)["540305","TRRsed3"] #


phyloseq_density_normalize <-  function(physeq = physeq,
                                        qPCR = physeq %>% sample_data() %>% data.frame() %>% rownames_to_column("tmp"), 
                                        sample_idx = "tmp", 
                                        value_idx = "ReadNorm")   


{
  require(tidyverse)
  require(phyloseq)
  require(HTSSIP)
  
    stopifnot(class(qPCR) == "data.frame" | class(qPCR) == "matrix")
    # stopifnot(!is.null(qPCR$Sample))
    df_OTU_col = colnames(phyloseq::otu_table(physeq))
    df_OTU = phyloseq2df(physeq, phyloseq::otu_table)
    df_OTU_rn = rownames(df_OTU)
    df_OTU = as.data.frame(apply(df_OTU, 2, HTSSIP::as.Num))
    rownames(df_OTU) = df_OTU_rn
    df_OTU = tss(df_OTU)
    rownames(qPCR) = make.names(qPCR[, sample_idx])
    qPCR = qPCR[colnames(df_OTU), ]
    qPCR_vals = qPCR[, value_idx]
    if (length(qPCR_vals) != ncol(df_OTU)) {
      stop("length qPCR_vals (", length(qPCR_vals), ") != ncol df_OTU (", 
           ncol(df_OTU), ")")
    }
    df_OTU = sweep(df_OTU %>% as.data.frame, 2, qPCR_vals, "*")
    df_OTU = apply(df_OTU, 2, function(x) round(x, 0))
    colnames(df_OTU) = df_OTU_col
    df_OTU[is.na(df_OTU)] = 0
    tree = phyloseq::phy_tree(physeq, errorIfNULL = FALSE)
    tax = phyloseq::tax_table(physeq, errorIfNULL = FALSE)
    sam = phyloseq::sample_data(physeq, errorIfNULL = FALSE)
    physeq2 = phyloseq::phyloseq(phyloseq::otu_table(df_OTU, 
                                                     taxa_are_rows = TRUE), phyloseq::phy_tree(tree, errorIfNULL = FALSE), 
                                 phyloseq::tax_table(tax, errorIfNULL = FALSE), phyloseq::sample_data(sam, 
                                                                                                      errorIfNULL = FALSE))
    return(physeq2)
  }
  
# https://bioconductor.org/packages/release/data/experiment/vignettes/curatedMetagenomicData/inst/doc/curatedMetagenomicData.html
  # 
  # sweep(physeq %>%
  #         otu_table() %>%
  #         magrittr::divide_by(sample_sums(physeq)), 2,
  #       sample_data(physeq) %>% data.frame() %>% pull(ReadNorm), "*") %>%
  #   # sample_data(physeq)$norm, "*") %>%
  #   round() %>%
  #   as.data.frame() %>%
  #   # replace_na(list(0))
  #   select_if(~ !any(is.na(.))) -> counts
  # 
  # counts %>%
  #   as_tibble()
  # 
  # # plot(counts %>%
  # #        colSums(),
  # #      sample_data(physeq)$ReadNorm) -> p
  # 
  # physeq@otu_table = NULL
  # 
  # merge_phyloseq(
  #   otu_table(counts, taxa_are_rows = TRUE),
  #   physeq) -> physeq


  # out <- list("physeq" = physeq,
  #             "plot" = p

#   return(out2)
#   
#   
# }
