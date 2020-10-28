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
#'data("enterotype")
#'
#'phyloseq_check_lib_size(enterotype,"SeqTech","Project", 1000, 10) -> out
#'
#'

phyloseq_check_lib_size <- function(physeq, data_color, data_facet, nreads_display, first_n)
{
  require(tidyverse)
  require(speedyseq)

  if ('SampleID' %in% (physeq %>% sample_data() %>% colnames()))
      {
        physeq %>%
          sample_data() %>%
          as.matrix() %>%
          as.data.frame() -> df
  }else
  {
    physeq %>%
      sample_data() %>%
      as.matrix() %>%
      as.data.frame() %>%
      rownames_to_column('SampleID') -> df
  }
  df$LibrarySize <- sample_sums(physeq)
  df <- df[order(df$LibrarySize),]
  df$Index <- seq(nrow(df))

  p <- ggplot(data=df,  #ggplot(data=subset(df,Index<400),
              aes_string(x="Index", y="LibrarySize", color = data_color , label = "SampleID")) +
    geom_point(size=0.5) +
    ggrepel::geom_text_repel(
      data = df %>% filter(LibrarySize < nreads_display)  #sample_type == "NC" | sample_type ==  "MOCK")
    ) 
  if(is.null(data_facet)==FALSE){
    p <- p+ facet_wrap(~ get(data_facet) , ncol = 1)
  }
  

  p + coord_cartesian(xlim=c(0,first_n)) +
    geom_hline(yintercept = nreads_display, color="red" , size = 0.5) -> p

  out <- list("plot" = p,
              "df" = df)

  return(out)

}


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
#'data("enterotype")
#'
#'phyloseq_rarefaction_curves(GlobalPatterns, 10000,  "Primer", "SampleType") -> out
#'
#'

phyloseq_rarefaction_curves<- function (physeq, stepsize, color_data, facet_data)
{
  require(ampvis2)
  #devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6") # to source phyloseq_to_ampvis2()

  tax_table(physeq) <- tax_table(physeq)[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")] # because amp_rarecurve() doesn't want to deal with our "Strain" level taxonomy

  physeq %>%
    # prune_samples(samples  =  sample_sums(physeq_tmp) > 1000) %>%
    # prune_taxa(taxa =  taxa_sums(physeq_tmp) > 1) %>%
    # rarefy_even_depth(rngseed = seed) %>%
    phyloseq_to_ampvis2() %>%
    amp_rarecurve(stepsize = stepsize,
                  color_by = color_data,
                  facet_by = facet_data,
                  facet_scales = "fixed") +
    ylab("Number of Observed ASVs") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) -> p

  return(p)
}

phyloseq_to_ampvis2 <- function(physeq) {
  #check object for class
  if(!any(class(physeq) %in% "phyloseq"))
    stop("physeq object must be of class \"phyloseq\"", call. = FALSE)
  
  #ampvis2 requires taxonomy and abundance table, phyloseq checks for the latter
  if(is.null(physeq@tax_table))
    stop("No taxonomy found in the phyloseq object and is required for ampvis2", call. = FALSE)
  
  #OTUs must be in rows, not columns
  if(phyloseq::taxa_are_rows(physeq))
    abund <- as.data.frame(phyloseq::otu_table(physeq)@.Data)
  else
    abund <- as.data.frame(t(phyloseq::otu_table(physeq)@.Data))
  
  #tax_table is assumed to have OTUs in rows too
  tax <- phyloseq::tax_table(physeq)@.Data
  
  #merge by rownames (OTUs)
  otutable <- merge(
    abund,
    tax,
    by = 0,
    all.x = TRUE,
    all.y = FALSE,
    sort = FALSE
  )
  colnames(otutable)[1] <- "OTU"
  
  #extract sample_data (metadata)
  if(!is.null(physeq@sam_data)) {
    metadata <- data.frame(
      phyloseq::sample_data(physeq),
      row.names = phyloseq::sample_names(physeq), 
      stringsAsFactors = FALSE, 
      check.names = FALSE
    )
    
    #check if any columns match exactly with rownames
    #if none matched assume row names are sample identifiers
    samplesCol <- unlist(lapply(metadata, function(x) {
      identical(x, rownames(metadata))}))
    
    if(any(samplesCol)) {
      #error if a column matched and it's not the first
      if(!samplesCol[[1]])
        stop("Sample ID's must be in the first column in the sample metadata, please reorder", call. = FALSE)
    } else {
      #assume rownames are sample identifiers, merge at the end with name "SampleID"
      if(any(colnames(metadata) %in% "SampleID"))
        stop("A column in the sample metadata is already named \"SampleID\" but does not seem to contain sample ID's", call. = FALSE)
      metadata$SampleID <- rownames(metadata)
      
      #reorder columns so SampleID is the first
      metadata <- metadata[, c(which(colnames(metadata) %in% "SampleID"), 1:(ncol(metadata)-1L)), drop = FALSE]
    }
  } else
    metadata <- NULL
  
  #extract phylogenetic tree, assumed to be of class "phylo"
  if(!is.null(physeq@phy_tree)) {
    tree <- phyloseq::phy_tree(physeq)
  } else
    tree <- NULL
  
  #extract OTU DNA sequences, assumed to be of class "XStringSet"
  if(!is.null(physeq@refseq)) {
    #convert XStringSet to DNAbin using a temporary file (easiest)
    fastaTempFile <- tempfile(pattern = "ampvis2_", fileext = ".fa")
    Biostrings::writeXStringSet(physeq@refseq, filepath = fastaTempFile)
  } else
    fastaTempFile <- NULL
  
  #load as normally with amp_load
  ampvis2::amp_load(
    otutable = otutable,
    metadata = metadata,
    tree = tree,
    fasta = fastaTempFile
  )
}

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
#'data("enterotype")
#'
#'phyloseq_rarefied_unrerefied_richness(GlobalPatterns, 123,  "Primer", "Primer", "SampleType") -> out
#'
#'
#'
#'
#'
#'

phyloseq_rarefied_unrarefied_richness <- function(physeq, sample_size, seed, color_data, fill_data, shape_data)
{
  require(phyloseq)
  physeq %>%
    rarefy_even_depth(rngseed = seed,
                      sample.size = sample_size) -> physeq_rarefied

  prune_samples(sample_sums(physeq)>= sample_size, physeq) -> physeq


  data.frame(a =  estimate_richness(physeq_rarefied, measures = "Observed")[, 1],
             b = estimate_richness(physeq, measures = "Observed")[, 1],
             sample_data(physeq)) %>%
    ggplot(aes(x = a, y = b)) +
    geom_point(aes_string(color = color_data, fill = fill_data, shape = shape_data)) +
    geom_smooth(method="lm", level=0.95) +
    labs(x = "\nRarefied Richness", y = "UN-Rarefied Richness\n") +
    theme_minimal() -> p

  return(p)
}


#' @title Identify and Filter ASV from phyloseq object / dataframe with samples as rows.
#'
#' @param physeq phyloseq object / dataframe with samples as rows
#' @param thrs threshold in %
#' @author Florentin Constancias
#' @note For each sample, determine cumulated sum of percentage of sorted ASV.
#' @note Then, define Common ASV for wich cumulated sum of percentag is less or equal to percentage threshold.
#' @note As ASV can be rare / common in different samples, ASV are defined as common if they are at least common in on samples. The rares ASV are the remaining ones.
#' @return list of 2 data frames / phyloseq objects of rare - common ASV
#' @export
#' @examples
#' library(phyloseq)
#' library(tidyverse)
#' data("esophagus")
#'
#' # Data table with taxa samples as rows
#' esophagus %>%
#' otu_table() %>%
#' as.data.frame() -> data
#'
#' phyloseq_filter_samples(data, 50) -> results
#'

phyloseq_filter_samples <- function(physeq, thrs)
{
  require(phyloseq)
  require(tidyverse)

  if (class(physeq) == "phyloseq")
  {
    physeq %>%
      otu_table() %>%
      data.frame() %>%
      t() %>%
      as.data.frame() -> df
    tax <- taxa_names(physeq)
    ntax <- ntaxa(physeq)
    samples <- sample_names(physeq)
  }
  else{
    df = physeq
    samples <- colnames(df)
    ntax <- nrow(df)
    tax <- rownames(df)
  }

  df %>%
    rownames_to_column('ASV')%>%
    select(ASV) -> results_all
  results_all_2 <- results_all
  results_all_3 <- results_all
  for (c in samples)
  {
    # print(c)
    df %>%
      rownames_to_column('ASV') %>%
      select(ASV, c) %>%
      arrange(-get(c)) %>%
      mutate(percent = get(c)/sum(get(c))) %>%
      mutate(Accumulated_pc = cumsum(percent)) %>%
      mutate(Common_Rare = ifelse(Accumulated_pc <= thrs/100, "Common", "Rare")) -> tpm

    tpm %>%
      mutate(new_count = ifelse(Common_Rare == "Common", get(c), 0)) %>%
      select(- percent, - Accumulated_pc, -!!c, -Common_Rare) %>%
      dplyr::rename(!!c := new_count) -> filtter_sample_common

    tpm %>%
      mutate(new_count = ifelse(Common_Rare == "Rare", get(c), 0)) %>%
      select(- percent, - Accumulated_pc, -!!c, -Common_Rare) %>%
      dplyr::rename(!!c := new_count) -> filtter_sample_rare

    tpm %>%
      select(- percent, - Accumulated_pc, -!!c) %>%
      dplyr::rename(!!c := Common_Rare) -> results

    full_join(results_all,
              results) -> results_all

    full_join(results_all_2,
              filtter_sample_common) -> results_all_2

    full_join(results_all_3,
              filtter_sample_rare) -> results_all_3
  }

  results_all %>%
    pivot_longer(samples) %>%
    # group_by(ASV) %>%
    filter(value == "Common") %>%
    pull(ASV) %>%
    unique() -> commons

  # results_all %>%
  #   pivot_longer(samples) %>% group_by(ASV) %>%
  #   filter(value == "Rare") %>% distinct(ASV) %>% pull() -> rares

  setdiff(tax,
          commons) -> rare

  # union(commons,
  #         rare) %>% setdiff(tax)

  # print(paste0("Results : ", length(commons)," common OTUs ", length(rare), " rare OTUs - over a total of ", ntax, " OTUs"))

  # common defined as never rare:

  if (class(physeq) == "phyloseq")
  {
    prune_taxa(commons,
               physeq) -> physeq_common

    prune_taxa(rare,
               physeq) -> physeq_rare

    merge_phyloseq(
      results_all_2 %>% column_to_rownames('ASV') %>% as.matrix() %>% otu_table(taxa_are_rows = TRUE),
      physeq %>% tax_table(),
      physeq %>% refseq(),
      physeq %>% phy_tree(),
      physeq %>% sample_data()
    ) %>%
      filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> physeq_persample_common

    merge_phyloseq(
      results_all_3 %>% column_to_rownames('ASV') %>% as.matrix() %>% otu_table(taxa_are_rows = TRUE),
      physeq %>% tax_table(),
      physeq %>% refseq(),
      physeq %>% phy_tree(),
      physeq %>% sample_data()
    ) %>%
      filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> physeq_persample_rare


    out <- list("global_common" = physeq_common,
                "global_rare" = physeq_rare,
                "per_sample_common" = physeq_persample_common,
                "per_sample_rare" = physeq_persample_rare)

    return(out)
  }else
  {
    df %>%
      rownames_to_column('ASV') %>%
      filter(ASV %in% commons) -> df_common

    df %>%
      rownames_to_column('ASV') %>%
      filter(ASV %in% rare) -> df_rare

    out <- list("global_common" = df_common,
                "global_rare" = df_rare,
                "per_sample_common" = results_all_2,
                "per_sample_rare" = results_all_3
    )

    return(out)
  }

}

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
#'
#'phyloseq_get_strains(GlobalPatterns) -> out
#'
#'

phyloseq_get_strains <- function(physeq)
{

  require(tidyverse)
  physeq_tmp = physeq

  as(tax_table(physeq), "matrix") %>%
    data.frame() %>%
    rownames_to_column(var = "ASV") %>%
    mutate_at(vars(everything()), na_if, "unknown") -> tmp1
  
  tmp1 %>%
   column_to_rownames("ASV")  %>%
   as.matrix() -> tax_table(physeq_tmp)

  physeq_tmp %>%
    get_strains(label = "unknown",
              species = TRUE) -> tmp2
  
  as(tax_table(tmp2), "matrix") %>%
    data.frame() %>%
    rownames_to_column(var = "ASV") %>%
    mutate_if(is.factor, as.character) %>%
    unite(Strain, Species, ASV,
          sep = " ", remove = FALSE, na.rm = TRUE) %>%
    column_to_rownames("ASV") %>%
    select(-Strain, Strain) -> tax_tbl_tmp

  # full_join(tmp1,
  #           tax_tbl_tmp %>%
  #             select(Strain) %>%
  #             rownames_to_column("ASV")) %>%
  #   column_to_rownames("ASV") %>%
  #   # replace(is.na(.), "unknown") %>%
  #   as.matrix() -> tax_table(physeq)
  tmp1 %>%
    as.matrix() -> tax_table(physeq) 

  return(physeq)
}

get_strains <- function(physeq_obj, label = "Unannotated", other_label = NULL, 
          species = FALSE, unique_rank = NULL, unique_sep = " ") 
{
  tax_tbl <- phyloseq::tax_table(physeq_obj)
  tax_names <- colnames(tax_tbl)
  tax_tbl <- t(apply(tax_tbl, 1, function(x) {
    n <- length(x)
    if (sum(is.na(x)) == n) {
      tax_ranks <- rep("Unknown", n)
    }
    else {
      if (sum(is.na(x)) != 0) {
        i <- max(which(!is.na(x)))
        rank <- x[i]
        x[which(is.na(x))] <- sprintf("%s %s (%s)", label, 
                                      rank, names(x)[i])
      }
      else {
        if (!is.null(other_label)) {
          if (sum(other_label %in% x) > 0) {
            tax_ranks <- x
          }
          else {
            if (species) {
              x[n] <- sprintf("%s %s", x[n - 1], x[n])
            }
          }
        }
        else {
          if (species) {
            x[n] <- sprintf("%s %s", x[n - 1], x[n])
          }
        }
      }
      tax_ranks <- x
      return(tax_ranks)
    }
  }))
  if (!is.null(unique_rank)) {
    ind <- which(colnames(tax_tbl) == unique_rank)
    tax_tbl[, ind] <- as.character(tax_tbl[, ind])
    tax_tbl[, ind] <- as.character(gen_uniq_lbls(tax_tbl[, 
                                                         ind], sep_char = unique_sep))
  }
  phyloseq::tax_table(physeq_obj) <- tax_tbl
  colnames(phyloseq::tax_table(physeq_obj)) <- tax_names
  return(physeq_obj)
}

#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .remove ASV assigned to Chloroplast and Mitochondria
#' @note .
#' @note .
#' @return .
#' @export
#' @examples
#'
#'library(phyloseq)
#'data("enterotype")
#'
#'
#'
#'

phyloseq_remove_chloro_mitho <- function(physeq)
{
  physeq %>%
    subset_taxa(Order != "Chloroplast" |
                  Family != "Mitochondria") %>%
    filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> your_phyloseq_clean
  return(your_phyloseq_clean)
}

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
  df_OTU = HTSSIP::phyloseq2df(physeq, phyloseq::otu_table)
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
#'
#' #' @title ...
#' #' @param .
#' #' @param ..
#' #' @author Florentin Constancias
#' #' @note .
#' #' @note .
#' #' @note .
#' #' @return .
#' #' @export
#' #' @examples
#' #'
#' #'
#' phyloseq_remove_contaminants <- function(physeq, sample_type, NTC_label, batch = "PCR_plate", batch.combine = "minimum", normalize = TRUE, threshold = 0.1, facet_plot = NULL, taxa_plot="Family")
#' 
#' {
#' # Remove ASV as long as they occur once in the NTC samples                
#' ## get ASV id found in the Nuclease_free_H2O = NTC samples
#' physeq %>%
#'   subset_samples(origin %in% c("Nuclease_free_H2O")) %>%
#'   filter_taxa(function(x) sum(x) > 0, TRUE) %>%
#'   taxa_names() -> ASV_NTC
#'
#' physeq %>%
#'   transform_sample_counts(function(x) x/sum(x) *100)  -> physeq_tmp
#'
#' ## export p1 wich is a diagnostic plot                            
#' prune_taxa(ASV_NTC, physeq_tmp) %>%
#'   plot_bar(fill= taxa_plot) +
#'   facet_grid(~ facet_plot ,scales = "free_x", space = "free") +
#'   theme(plot.title = element_text(hjust = 0.5)) +
#'   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 3)) -> p1
#'
#'
#' # Use the decontam appraoch # https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
#'
#' sample_data(physeq)$is.neg <- sample_data(physeq)$origin == "Nuclease_free_H2O"
#'
#' contamdf.prev <- decontam::isContaminant(physeq,
#'                                          method = "prevalence",
#'                                          neg = "is.neg",
#'                                          batch = batch,
#'                                          batch.combine = batch.combine,
#'                                          normalize = normalize,
#'                                          threshold = threshold,
#'                                          detailed = T)
#'
#'
#'
#' # Make phyloseq object of presence-absence in negative controls and true samples
#' ps.pa <- transform_sample_counts(physeq, function(abund) 1*(abund>0))
#' ps.pa.neg <- prune_samples(sample_data(physeq)$origin == "Nuclease_free_H2O", physeq)
#' ps.pa.pos <- prune_samples(sample_data(physeq)$origin != "Nuclease_free_H2O", physeq)
#'
#' # Make data.frame of prevalence in positive and negative samples
#' df.pa <- data.frame(pa.pos = taxa_sums(ps.pa.pos),
#'                     pa.neg = taxa_sums(ps.pa.neg),
#'                     contaminant = contamdf.prev$contaminant)#,
#' # tax = tax_table(ps.pa.pos))
#'
#' ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color = contaminant)) + geom_point() +
#'   # scale_y_continuous(trans = 'log2') +
#'   # scale_x_continuous(trans = 'log2') +
#'   # xlim(0,1000) + ylim(0,1000) + coord_equal() +
#'   xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
#'
#'
#' prune_taxa(subset(df.pa , contaminant == TRUE) %>%
#'              rownames(), physeq %>%
#'              transform_sample_counts(function(x) x/sum(x) * 100) ) %>%
#'   plot_bar(fill="Strain") +
#'   facet_grid(as.formula(paste0(taxa_plot, "~ ", facet_plot)),scales = "free_x", space = "free") +
#'   theme(plot.title = element_text(hjust = 0.5)) +
#'   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 1)) -> p2
#'
#' physeq_save <- physeq
#'
#' physeq %>%
#'   prune_taxa(taxa =! contamdf.prev$contaminant) %>%
#'   #subset_taxa(! Class == "unknown") %>%
#'   #subset_samples(! origin %in% c("GDC-MOCK", "Nuclease_free_H2O")) %>%
#'   filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> physeq # more than 0 in at least two samples (more than 1) # here we should check based on prevalence plot
#'
#' # export all outputs: p1, p2, phyloseq-object filterd with p1, with p2
#' out <- list("plot_ASV" = p1,
#'             "plot_decontam" = p2,
#'              "phyloseq_ASV" = ,
#'              "phyloseq_decontam" = )
#' return(out)                
#'}
