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

## Set of graphical methods for phyloseq objects (mainly related to ordination,
## and library size after normalisation)
require(ggplot2)
require(scales)
require(reshape2)

## Rarefaction curve, ggplot style
ggrare <- function(physeq, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  ## Args:
  ## - physeq: phyloseq class object, from which abundance data are extracted
  ## - step: Step size for sample size in rarefaction curves
  ## - label: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - color: (Optional). Default ‘NULL’. Character string. The name of the
  ##          variable to map to colors in the plot. This can be a sample
  ##          variable (among the set returned by
  ##          ‘sample_variables(physeq)’ ) or taxonomic rank (among the set
  ##          returned by ‘rank_names(physeq)’).
  ##
  ##          Finally, The color scheme is chosen automatically by
  ##          ‘link{ggplot}’, but it can be modified afterward with an
  ##          additional layer using ‘scale_color_manual’.
  ## - color: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - plot:  Logical, should the graphic be plotted.
  ## - parallel: should rarefaction be parallelized (using parallel framework)
  ## - se:    Default TRUE. Logical. Should standard errors be computed.
  # require(vegan)
  x <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  ## Get sample data
  if (!is.null(sample_data(physeq, FALSE))) {
    sdf <- as(sample_data(physeq), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  ## Add, any custom-supplied plot-mapped variables
  if( length(color) > 1 ){
    data$color <- color
    names(data)[names(data)=="color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  if( length(label) > 1 ){
    labels$label <- label
    names(labels)[names(labels)=="label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot(data = data, aes_string(x = "Size", y = ".S", group = "Sample", color = color))
  p <- p + labs(x = "Sample Size", y = "Species Richness")
  if (!is.null(label)) {
    p <- p + geom_text(data = labels, aes_string(x = "x", y = "y", label = label, color = color),
                       size = 4, hjust = 0)
  }
  p <- p + geom_line()
  if (se) { ## add standard error if available
    p <- p + geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
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
  
  require(phyloseq); require(tidyverse)
  physeq_tmp = physeq
  
  physeq %>%
    tax_table() %>%
    as.data.frame() %>%
    rownames_to_column(var = "ASV") %>%
    mutate_at(vars(everything()), na_if, "unknown") %>%
    column_to_rownames("ASV") %>%
    as.matrix() -> tax_table(physeq_tmp)
  
  physeq_tmp %>%
    get_strains(label = "unknown",
                species = TRUE) %>%
    tax_table() %>%
    as.data.frame() %>%
    rownames_to_column("ASV") %>%
    mutate_if(is.factor, as.character) %>%
    unite(Strain, Species, ASV,
          sep = " ", remove = FALSE, na.rm = TRUE) %>%
    column_to_rownames("ASV") %>%
    dplyr::select(-Strain, Strain) -> tax_tbl_tmp
  
  full_join(physeq %>% tax_table() %>%
              as.data.frame() %>%
              rownames_to_column("ASV"),
            tax_tbl_tmp %>%
              select(Strain) %>%
              rownames_to_column("ASV")) %>%
    column_to_rownames("ASV") %>%
    # replace(is.na(.), "unknown") %>%
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

phyloseq_get_strains_fast <- function(physeq)
{
  
  require(phyloseq);require(tidyverse)
  physeq_tmp = physeq
  
  as(tax_table(physeq), "matrix") %>%
    data.frame() %>%
    rownames_to_column('ASV') %>% 
    mutate_at(vars(everything()), na_if, "unknown") -> tmp1
  
  
  tmp1 %>%
    column_to_rownames("ASV") %>%
    as.matrix() -> tax_table(physeq_tmp)
  
  physeq_tmp %>%
    get_strains(label = "unknown",
                species = TRUE) -> tmp2
  
  as(tax_table(tmp2), "matrix") %>%
    data.frame() %>%
    rownames_to_column('ASV')  %>%
    mutate_if(is.factor, as.character) %>%
    unite(Strain, Species, ASV,
          sep = " ", remove = FALSE, na.rm = TRUE) %>%
    column_to_rownames("ASV") %>%
    dplyr::select(-Strain, Strain) -> tax_tbl_tmp
  
  full_join(tmp1,
            tax_tbl_tmp %>%
              dplyr::select(Strain) %>%
              rownames_to_column("ASV")) %>%
    column_to_rownames("ASV") %>%
    # replace(is.na(.), "unknown") %>%
    as.matrix() -> tax_table(physeq)
  
  return(physeq)
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
    subset_taxa(Order != "Chloroplast" &
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
#'sample_data(GlobalPatterns)$toto = sample(5000000:1000000, size=nsamples(GlobalPatterns), replace=TRUE)
#'phyloseq_density_normalize(GlobalPatterns, "toto") -> out
#'HTSSIP::OTU_qPCR_trans(GlobalPatterns, GlobalPatterns %>% sample_data() %>% data.frame() , "X.SampleID", value_idx = "norm") -> out2
#'otu_table(GlobalPatterns)["540305","TRRsed3"] #2305
#'sample_data(GlobalPatterns)["TRRsed3",] #2664122
#'sample_sums(GlobalPatterns)["TRRsed3"] #279704
#'(2305/279704)  * 2664122
#'otu_table(out)["540305","TRRsed3"] #
#'otu_table(out2)["540305","TRRsed3"] #

phyloseq_density_normalize <-  function(physeq = physeq,
                                        value_idx = "norm",
                                        remove.na=TRUE,
                                        set.na.0=FALSE)
  
  
{
  require(tidyverse)
  require(phyloseq)
  
  # stopifnot(!is.null(qPCR$Sample))
  
  
  #if we want to remove the samples for which we have no qPCR data to avoid future NAs in the phyloseq object
  if(remove.na==TRUE){

    prune_samples(!is.na(physeq %>% get_variable(value_idx)),physeq) -> physeq
  }
  
  df_OTU_col = colnames(phyloseq::otu_table(physeq))
  df_OTU = phyloseq2df(physeq, phyloseq::otu_table)
  df_OTU_rn = rownames(df_OTU)
  df_OTU = as.data.frame(apply(df_OTU, 2, as.Num))
  rownames(df_OTU) = df_OTU_rn
  df_OTU = tss(df_OTU)
  
  sample_idx="tmp"
  physeq %>% sample_data() %>% data.frame() %>% rownames_to_column(sample_idx) -> qPCR
  
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
  
  #if we want to replace NAs with 0 
  if(set.na.0==TRUE){
    df_OTU[is.na(df_OTU)] = 0
  }
  
  tree = phyloseq::phy_tree(physeq, errorIfNULL = FALSE)
  tax = phyloseq::tax_table(physeq, errorIfNULL = FALSE)
  sam = phyloseq::sample_data(physeq, errorIfNULL = FALSE)
  physeq2 = phyloseq::phyloseq(phyloseq::otu_table(df_OTU,
                                                   taxa_are_rows = TRUE), phyloseq::phy_tree(tree, errorIfNULL = FALSE),
                               phyloseq::tax_table(tax, errorIfNULL = FALSE), phyloseq::sample_data(sam,
                                                                                                    errorIfNULL = FALSE))
  return(physeq2)
}


# phyloseq_density_normalize_2 <-  function(physeq = physeq,
#                                         qPCR = physeq %>% sample_data() %>% data.frame() %>% rownames_to_column("tmp"),
#                                         sample_idx = "tmp",
#                                         value_idx = "ReadNorm")
# 
# 
# {
#   require(tidyverse)
#   require(phyloseq)
#   require(HTSSIP)
# 
#   stopifnot(class(qPCR) == "data.frame" | class(qPCR) == "matrix")
#   # stopifnot(!is.null(qPCR$Sample))
#   df_OTU_col = colnames(phyloseq::otu_table(physeq))
#   df_OTU = HTSSIP::phyloseq2df(physeq, phyloseq::otu_table)
#   df_OTU_rn = rownames(df_OTU)
#   df_OTU = as.data.frame(apply(df_OTU, 2, HTSSIP::as.Num))
#   rownames(df_OTU) = df_OTU_rn
#   df_OTU = tss(df_OTU)
#   rownames(qPCR) = make.names(qPCR[, sample_idx])
#   qPCR = qPCR[colnames(df_OTU), ]
#   qPCR_vals = qPCR[, value_idx]
#   if (length(qPCR_vals) != ncol(df_OTU)) {
#     stop("length qPCR_vals (", length(qPCR_vals), ") != ncol df_OTU (",
#          ncol(df_OTU), ")")
#   }
#   df_OTU = sweep(df_OTU %>% as.data.frame, 2, qPCR_vals, "*")
#   df_OTU = apply(df_OTU, 2, function(x) round(x, 0))
#   colnames(df_OTU) = df_OTU_col
#   df_OTU[is.na(df_OTU)] = 0
#   tree = phyloseq::phy_tree(physeq, errorIfNULL = FALSE)
#   tax = phyloseq::tax_table(physeq, errorIfNULL = FALSE)
#   sam = phyloseq::sample_data(physeq, errorIfNULL = FALSE)
#   physeq2 = phyloseq::phyloseq(phyloseq::otu_table(df_OTU,
#                                                    taxa_are_rows = TRUE), phyloseq::phy_tree(tree, errorIfNULL = FALSE),
#                                phyloseq::tax_table(tax, errorIfNULL = FALSE), phyloseq::sample_data(sam,
#                                                                                                     errorIfNULL = FALSE))
#   return(physeq2)
# }
# 

phyloseq2df <- function (physeq, table_func) 
{
  physeq.md = table_func(physeq)
  physeq.md = suppressWarnings(as.data.frame(as.matrix(physeq.md)))
  physeq.md = as.matrix(data.frame(lapply(physeq.md, as.character)))
  physeq.md = as.data.frame(apply(physeq.md, 2, trimws))
  rownames(physeq.md) = rownames(table_func(physeq))
  return(physeq.md)
}

as.Num <- function (x) 
{
  as.numeric(as.character(x))
}

tss <- function (x, MARGIN = 2, na.rm = FALSE) 
{
  k = min(x, na.rm = na.rm)
  tmp = pmax(k, apply(x, MARGIN, sum, na.rm = na.rm))
  x = sweep(x, MARGIN, tmp, "/")
  return(x)
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



#' @title Decontaminate the phyloseq object by removing all the ASVs present in the negative controls
#' @author Sneha Sundar
#' @param physeq The phyloseq object
#' @param sample_type The metadata column that indicates whether a sample is a negative control or not
#' @param NTC_label The label given to the negative control samples in the \code{sample_type} column
#' @param facet_plot The name of the column in the metadata that should be used for faceting in the diagnostic bar plot
#' @param taxa_plot The taxa level to plot in the diagnostic bar plot
#' @param Strain logical. TRUE gives the strain-level information of contaminant ASV rather than just the ASV IDs (like ASV01,ASV04) of contaminants. If TRUE, can also use 'Strain' in taxa_plot argument.
#' @return 
#' A list containing the contaminant ASVs (\code{contaminant.ASV}),
#' the decontaminated phyloseq object (\code{physeq.decontaminated}) 
#' and a diagnostic bar plot of the relative abundances of ASVs present 
#' in the negative controls for all the samples (\code{diagnostic.plot})
#' @export
#' @examples 
#' physeq_remove_contaminants_crude(physeq=physeq,
#'                                  sample_type = "Experiment",
#'                                  ,NTC_label = "NTC",
#'                                  facet_plot = 'Reactor', 
#'                                  taxa_plot="Strain"
#'                                  Strain = TRUE)

physeq_remove_contaminants_crude <- function(physeq, 
                                             sample_type = 'Type', 
                                             NTC_label,
                                             facet_plot, 
                                             taxa_plot,
                                             Strain){
  
  
  
  
  # Crude Approach: Remove ASV as long as they occur once in the NTC samples
  
  require(phyloseq)
  require(tidyverse)

  
  #adding the strain level annotation to the asvs (ie. instead of ASV1,ASV2..,etc, we now have strain level annotation as taxa_names())
  
  
  if("Strain" %in% rank_names(physeq)  && Strain == TRUE)
  {
    taxa_names(physeq)  <- tax_table(physeq)[,"Strain"]
    
    
  }
  
  #defining an operator 'not in'
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  if("Strain" %!in% rank_names(physeq)  && Strain == TRUE)
  {
    physeq %>%
      phyloseq_get_strains_fast() -> physeq
    
    taxa_names(physeq)  <- tax_table(physeq)[,"Strain"]
    
  }
  

  #create logical vector indicating whether sample is negative or not
  is.neg <- sample_data(physeq)[[sample_type]]==NTC_label
  
  ##Get the ASVs that are present in the negative controls  
  prune_samples(is.neg,physeq) %>% 
    filter_taxa(function(x) sum(x) > 0, TRUE) %>%
    taxa_names() -> ASV_NTC
  
  ## creating a temporary physeq object for use in plotting that is normalized to relative abundaces 
  physeq %>%
    transform_sample_counts(function(x) x/sum(x) *100)  -> physeq_tmp
  
  
  ##Creating a diagnostic bar plot of 
  # relative abundances of the ASVs (at the taxonomic level indicated in `taxa_plot`) 
  # present in the negative controls for all the samples. Facet wrap according to `facet_plot`.
  prune_taxa(ASV_NTC, physeq_tmp) %>%
    plot_bar(fill= taxa_plot) +
    facet_wrap(.~get(facet_plot),scales="free_x") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 3)) -> diagnostic_plot
  
  ASV_to_keep <- taxa_names(physeq)[!(taxa_names(physeq) %in% ASV_NTC)]
  physeq.decontaminated <- prune_taxa(ASV_to_keep,physeq)
  
  out <- list("contaminant.ASV" = ASV_NTC,
              "physeq.decontaminated" = physeq.decontaminated,
              "diagnostic.plot" = diagnostic_plot)
  
  return(out)
}


#' @title Decontaminate the phyloseq object by using the prevalence approach in the decontam package: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
#' @author Sneha Sundar
#' @param physeq The phyloseq object
#' @param sample_type The metadata column that indicates whether a sample is a negative control or not
#' @param NTC_label The label given to the negative control samples in the \code{sample_type} column
#' @param facet_plot The name of the column in the metadata that should be used for faceting in the diagnostic bar plot
#' @param taxa_plot The taxa level to plot in the diagnostic bar plot
#' @param Strain logical. TRUE gives the strain-level information of contaminant ASV rather than just the ASV IDs (like ASV01,ASV04) of contaminants. If TRUE, can also use 'Strain' in taxa_plot argument.#' @param batch The name of the metadata column indicating batch information (like sequencing run or plate number)
#' @param batch.combine Default "minimum". For each input sequence variant (or OTU) the probabilities calculated in each batch are combined into a single probability that is compared to \code{threshold} to classify contaminants. Valid values: "minimum", "product", "fisher".
#' @param threshold Default 0.1. The probability threshold below which (strictly less than) the null-hypothesis (not a contaminant) should be rejected in favor of the alternate hypothesis (contaminant).
#' @param normalize Default TRUE. If TRUE, the input seqtab is normalized so that each row sums to 1 (converted to frequency). If FALSE, no normalization is performed (the data should already be frequencies or counts from equal-depth samples).
#' @param detailed Default TRUE. If TRUE, the return value is a data.frame containing diagnostic information on the contaminant decision. If FALSE, the return value is a logical vector containing the binary contaminant classifications.
#' @return 
#' A list containing the contaminant ASVs (\code{contaminant.ASV}),
#' the decontaminated phyloseq object (\code{physeq.decontaminated}) 
#' and two diagnostic plots.\code{diagnostic.plot1} is a a scatter plot with the prevalence of ASVs using the  Negative controls on the x axis and the prevalence of ASVs using the TRUE samples
#' on the y axis. Each point is coloured according to whether it is a contaminant or not. \code{diagnostic.plot2} is a bar plot of the relative abundances of contaminant ASVs in all samples.
#' @export
#' @examples 
#' phyloseq_remove_contaminants_decontam(physeq, 
#'                                       sample_type="Experiment", 
#'                                       NTC_label="NTC",
#'                                       facet_plot = 'Reactor', 
#'                                       taxa_plot="Strain", 
#'                                       Strain = TRUE,
#'                                       batch = 'plate', 
#'                                       batch.combine = "minimum", 
#'                                      normalize = TRUE, 
#'                                      threshold = 0.1,
#'                                       detailed=TRUE)
#'@note No need to normalize to frequencies before running function


phyloseq_remove_contaminants_decontam <- function(physeq, 
                                                  sample_type, 
                                                  NTC_label,
                                                  facet_plot = NULL, 
                                                  taxa_plot,
                                                  Strain,
                                                  batch = NULL, 
                                                  batch.combine = "minimum", 
                                                  normalize = TRUE, 
                                                  threshold = 0.1,
                                                  detailed=TRUE){
  
  require(phyloseq)
  require(tidyverse)
  require(decontam)
  
  
  # decontam prevalence approach # https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
  
  
  #adding the strain level annotation to the asvs (ie. instead of ASV1,ASV2..,etc, we now have strain level annotation as taxa_names())
  
  if("Strain" %in% rank_names(physeq)  && Strain == TRUE)
  {
    taxa_names(physeq)  <- tax_table(physeq)[,"Strain"]
    
    
  }
  
  #defining an operator 'not in'
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  if("Strain" %!in% rank_names(physeq)  && Strain == TRUE)
  {
    physeq %>%
      phyloseq_get_strains_fast() -> physeq
    
    taxa_names(physeq)  <- tax_table(physeq)[,"Strain"]
    
  }
  
  #Create logical vector indicating whether a sample is a negative control or not: TRUE for Negative control
  sample_data(physeq)$is.neg <- sample_data(physeq)[[sample_type]]==NTC_label
  
  #decontam function
  contamdf.prev <- decontam::isContaminant(physeq,
                                           method = "prevalence",
                                           neg = "is.neg",
                                           batch = batch,
                                           batch.combine = batch.combine,
                                           normalize = normalize,
                                           threshold = threshold,
                                           detailed = detailed)
  
  #get the ASVs identified as contaminants by this approach
  ASV_decontam_prev<-rownames(contamdf.prev)[contamdf.prev$contaminant]
  
  # Make phyloseq object of presence-absence in negative controls and true samples (otu table is now a binary matrix of 0s and 1s)
  
  ps.pa <- transform_sample_counts(physeq, function(abund) 1*(abund>0))
  

  ps.pa.neg <- prune_samples(sample_data(physeq)[[sample_type]]==NTC_label,physeq)
  ps.pa.pos <- prune_samples(sample_data(physeq)[[sample_type]]!=NTC_label, physeq)
  
  # Make data.frame of prevalence in positive and negative samples
  df.pa <- data.frame(Prevalence_true_samples = taxa_sums(ps.pa.pos),
                      Prevalence_negative_controls = taxa_sums(ps.pa.neg),
                      contaminant = contamdf.prev$contaminant)
  
  df.pa %>% rownames_to_column("ASV") -> df.pa
  
  #diagnostic plot: a scatter plot with the prevalence of ASVs using the  Negative controls on the x axis and the prevalence of ASVs using the TRUE samples
  #on the y axis. Each point is coloured according to whether it is a contaminant or not.
  ggplot(data=df.pa, aes(x=Prevalence_negative_controls, y=Prevalence_true_samples, color = contaminant,label=ASV)) + geom_point() +
    # scale_y_continuous(trans = 'log2') +
    # scale_x_continuous(trans = 'log2') +
    # xlim(0,1000) + ylim(0,1000) + coord_equal() +
    xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") -> p1.decontam
  
  
  ## creating a temporary physeq object for use in plotting that is normalized to relative abundaces 
  physeq %>%
    transform_sample_counts(function(x) x/sum(x) *100)  -> physeq_tmp
  
  
  prune_taxa(ASV_decontam_prev, physeq_tmp) %>%
    #subset_samples(Experiment == "Continuous") %>% 
    plot_bar(fill= taxa_plot) +
    facet_wrap(.~get(facet_plot),scales="free_x") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 3)) -> p2.decontam
  
  ASV_to_keep <- taxa_names(physeq)[!(taxa_names(physeq) %in% ASV_decontam_prev)]
  physeq.decontaminated <- prune_taxa(ASV_to_keep,physeq)
  
  out <- list("contaminant.ASV" = ASV_decontam_prev,
              "physeq.decontaminated" = physeq.decontaminated,
              "diagnostic.plot1" = p1.decontam,
              "diagnostic.plot2" = p2.decontam)
  
  return(out)
  
  
}



#' @title Decontaminate the phyloseq object using the microdecon package (https://github.com/donaldtmcknight/microDecon )
#' @author Sneha Sundar
#' @param physeq The phyloseq object
#' @param grouping_var The name of the metadata column you want to group your samples by: e.g sequencing run or plate .
#' @param sample_type_var The metadata column that indicates whether a sample is a negative control or not 
#' @param NTC_label The label given to the negative control samples in the \code{sample_type} column
#' @param taxa_info Default TRUE. logical vector indicating that you want taxonomic information of the otus as well
#' @param taxa_level What level of taxonomic information is needed ? e.g "Strain","Genus"
#' @param facet_plot The name of the column in the metadata that should be used for faceting in the diagnostic bar plot
#' @param taxa_plot The taxa level to plot in the diagnostic bar plot
#' @param run_groupwise logical. If TRUE , will run the decon function separately for every group. See note for why this might be helpful. 
#' @return p
#' A list containing the result of the decon function (\code{result.microDecon}),
#' and a diagnostic bar plot of the relative abundances of ASVs that are totally removed from all the groups
#' for all the samples (\code{diagnostic.plot} and the decontaminated phyloseq object5)
#' If \code{run_groupwise} is TRUE, it will return a list and each element of the list is the output (as defined above with the microDecon function result, the diagnostic plot and the decontaminated phyloseq object) of doing the decontamination for each group separately . 
#' @export
#' @examples 
#' phyloseq_remove_contaminants_microDecon(physeq=physeq, 
#'                                        grouping_var=plate,
#'                                        sample_type_var=Experiment,
#'                                        NTC_label="NTC",
#'                                        taxa_info=TRUE,
#'                                        taxa_level='Strain',
#'                                        facet_plot='Reactor',
#'                                        taxa_plot='Order')
#' @note 
#' For the microdecontam approach we need to provide the data to the function in a specifically formatted dataframe:
#' rows are OTUs each column is individual sample . 
#' First column: OTU ID  
#' last column: taxonomic info (optional)
#' From column 2 : negative control samples (if there is more than one negative control , they will be averaged to produce a single mean negative control)* 
#'   
#'   The rest of the samples should be grouped by population ID, species, or some other sensible a priori grouping criteria (i.e., treat these groups as experimental blocks). For example PCR_plate can be the grouping criterion because that is where we are anticipating batch effects to play out. 
#' 
#' It is a good idea to see if the negative controls are similar or heterogenous. If they are heterogenous it is a good idea to run the function separately on the different groups (provided of course you have negative controls specifically for each group). You can do this by setting argument \code{run_groupwise} to TRUE

#' DON"T NORMALIZE READS BEFORE RUNNING FUNCTION

#' 

phyloseq_remove_contaminants_microDecon <- function(physeq,
                                                    grouping_var,
                                                    sample_type_var,
                                                    NTC_label,
                                                    taxa_info=TRUE,
                                                    taxa_level,
                                                    facet_plot = NULL , 
                                                    taxa_plot,
                                                    run_groupwise=FALSE){
  
  require(phyloseq)
  require(tidyverse)
  require(microDecon)
  
  #adding the strain level annotation to the asvs (ie. instead of ASV1,ASV2..,etc, we now have strain level annotation as taxa_names())
  
  if("Strain" %in% rank_names(physeq)  && taxa_level == "Strain")
  {
    taxa_names(physeq)  <- tax_table(physeq)[,"Strain"]
    
    
  }
  
  #defining an operator 'not in'
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  if("Strain" %!in% rank_names(physeq)  && taxa_level == "Strain")
  {
    physeq %>%
      phyloseq_get_strains_fast() -> physeq
    
    taxa_names(physeq)  <- tax_table(physeq)[,"Strain"]
    
  }
  
  
  if(run_groupwise==TRUE){
    
    #run the decontamination separately for each group
    
    result.run_groupwise <- vector("list",length(unique(physeq %>% get_variable(grouping_var))))
    
    groups <- unique(physeq %>% get_variable(grouping_var))
    for(i in 1:length(groups)){
      physeq.subset<-prune_samples(sample_data(physeq)[[grouping_var]]==groups[i],physeq)
      
      #sample_data(physeq) as a dataframe
      as(sample_data(physeq.subset),"matrix") %>%
        data.frame(check.names=FALSE) %>%
        rownames_to_column("Sample_ID") -> sample.data
      
      #group and arrange `sample.data`according to the grouping variable (e.g Plate or Batch number or sequencing run)
      sample.data<-sample.data %>% group_by(.data[[grouping_var]]) %>% arrange(.data[[grouping_var]])
      
      #get the grouped sample ids: we will use this to format the otu table
      grouped_sample_ids<-sample.data$Sample_ID
      
      #get otu table as dataframe
      as(otu_table(physeq.subset),"matrix") %>% data.frame(check.names=FALSE) %>% rownames_to_column("OTU_ID") -> otu.df
      
      #checking if we need taxonomic information of the ASVs and adding them if needed
      if(taxa_info == TRUE){
        tax.table <- as(tax_table(physeq.subset),"matrix") %>% data.frame(check.names=FALSE ) %>% rownames_to_column('ASV')
        taxa <- tax.table[,taxa_level]
        otu.df<-otu.df %>% mutate(Taxa = taxa)
      }
      
      #formatting the otu table in the way required by microDecon
      otu.df<-otu.df %>% relocate(all_of(grouped_sample_ids),.after=OTU_ID)
      
      NTC.df <- sample.data %>% filter(.data[[sample_type_var]] == "NTC")
      
      NTC_samples<- NTC.df$Sample_ID
      
      otu.df<-otu.df %>% relocate(all_of(NTC_samples),.after=OTU_ID)  #formatted table
      
      #number of negative controls
      n.blanks = length(NTC_samples)
      
      #vector of numbers listing the number of individuals in each user-specified group
      numb.ind<-sample.data %>% filter(.data[[sample_type_var]]!='NTC') %>% group_size()
      
      #performing decontamination:
      
      result<-decon(data=otu.df,
                    numb.blanks = n.blanks,
                    numb.ind = numb.ind,
                    taxa = T,
                    runs=2,
                    thresh = 0.7,
                    prop.thresh = 0.00005,
                    regression = 0,
                    low.threshold = 40,
                    up.threshold = 400)
      
      #type ?decon() in the console for more info about other parameters and how the result is formatted
      
      
      otu_removed_df<-result$OTUs.removed
      
      #otus removed from all groups
      otus_totally_removed <- otu_removed_df[otu_removed_df$All.groups=="Totally.removed","OTU_ID"]
      
      
      physeq.subset %>%
        transform_sample_counts(function(x) x/sum(x) *100)  -> physeq.subset_tmp
      
      #Creating a diagnostic bar plot of relative abundances of the ASVs (at the taxonomic level indicated in `taxa_plot`) removed from all groups. Facet wrap according to `facet_plot`.
      
      prune_taxa(otus_totally_removed, physeq.subset_tmp) %>%
        plot_bar(fill= taxa_plot) +
        facet_wrap(.~get(facet_plot),scales="free_x") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 3)) -> p1.microDecon
      
      #updating the otu table of physeq.subset object with decontaminated otu table
      physeq.decon <- physeq.subset
      otu_table(physeq.decon) <- otu_table(result$decon.table %>% select(-any_of(c("Mean.blank","Taxa"))) %>% remove_rownames() %>% column_to_rownames("OTU_ID") %>% data.matrix(),taxa_are_rows = TRUE)
      
      result.run_groupwise[[i]]<- list("result.microDecon" = result,
                                       "diagnostic.plot" = p1.microDecon,
                                       "physeq.decontaminated" = physeq.decon)
      
    }
    
    #name the result
    naming_fun <- function(x){return(paste(grouping_var,x,sep = "_"))}
    
    names(result.run_groupwise) <-c(sapply(as.character(groups),naming_fun,USE.NAMES = FALSE))
    
    
    return(result.run_groupwise) 
    
    
  }else{
    #sample_data(physeq) as a dataframe
    as(sample_data(physeq),"matrix") %>%
      data.frame(check.names=FALSE) %>%
      rownames_to_column("Sample_ID") -> sample.data
    
    #group and arrange `sample.data`according to the grouping variable (e.g Plate or Batch number or sequencing run)
    sample.data<-sample.data %>% group_by(.data[[grouping_var]]) %>% arrange(.data[[grouping_var]])
    
    #get the grouped sample ids: we will use this to format the otu table
    grouped_sample_ids<-sample.data$Sample_ID
    
    #get otu table as dataframe
    as(otu_table(physeq),"matrix") %>% data.frame(check.names=FALSE) %>% rownames_to_column("OTU_ID") -> otu.df
    
    #checking if we need taxonomic information of the ASVs and adding them if needed
    if(taxa_info == TRUE){
      tax.table <- as(tax_table(physeq),"matrix") %>% data.frame(check.names=FALSE ) %>% rownames_to_column('ASV')
      taxa <- tax.table[,taxa_level]
      otu.df<-otu.df %>% mutate(Taxa = taxa)
    }
    
    #formatting the otu table in the way required by microDecon
    otu.df<-otu.df %>% relocate(all_of(grouped_sample_ids),.after=OTU_ID)
    
    NTC.df <- sample.data %>% filter(.data[[sample_type_var]] == "NTC")
    
    NTC_samples<- NTC.df$Sample_ID
    
    otu.df<-otu.df %>% relocate(all_of(NTC_samples),.after=OTU_ID)  #formatted table
    
    #number of negative controls
    n.blanks = length(NTC_samples)
    
    #vector of numbers listing the number of individuals in each user-specified group
    numb.ind<-sample.data %>% filter(.data[[sample_type_var]]!='NTC') %>% group_size()
    
    #performing decontamination:
    
    result<-decon(data=otu.df,
                  numb.blanks = n.blanks,
                  numb.ind = numb.ind,
                  taxa = T,
                  runs=2,
                  thresh = 0.7,
                  prop.thresh = 0.00005,
                  regression = 0,
                  low.threshold = 40,
                  up.threshold = 400)
    
    #type ?decon() in the console for more info about other parameters and how the result is formatted
    
    
    otu_removed_df<-result$OTUs.removed
    
    #otus removed from all groups
    otus_totally_removed <- otu_removed_df[otu_removed_df$All.groups=="Totally.removed","OTU_ID"]
    
    
    physeq %>%
      transform_sample_counts(function(x) x/sum(x) *100)  -> physeq_tmp
    
    #Creating a diagnostic bar plot of relative abundances of the ASVs (at the taxonomic level indicated in `taxa_plot`) removed from all groups. Facet wrap according to `facet_plot`.
    
    prune_taxa(otus_totally_removed, physeq_tmp) %>%
      plot_bar(fill= taxa_plot) +
      facet_wrap(.~get(facet_plot),scales="free_x") +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 3)) -> p1.microDecon
    
    #updating the otu table of physeq object with decontaminated otu table
    physeq.decon <- physeq
    otu_table(physeq.decon) <- otu_table(result$decon.table %>% select(-any_of(c("Mean.blank","Taxa"))) %>% remove_rownames() %>% column_to_rownames("OTU_ID") %>% data.matrix(),taxa_are_rows = TRUE)
    
    out <- list("result.microDecon" = result,
                "diagnostic.plot" = p1.microDecon,
                "physeq.decontaminated" = physeq.decon)
    return(out)
  }
  

}

# to add:

# https://github.com/MBARI-BOG/BOG-Banzai-Dada2-Pipeline/blob/e40953dcb4980792d0320d6d1f4c815bfaa7484c/Pipeline_scripts/decon_std_outputs_v1.0.R
# https://github.com/Mettetron/3Species/blob/c36ed383fa0d81aeeec76b8a04bba3c8b588f7c5/DADA2_filterAndNorm.R
# https://github.com/zjgold/gruinard_decon/blob/3b1b2d076f9e74ed9c9c298f17f409621e62f44f/decontamination_utilities.R




