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
#'
#'source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_heatmap.R")
#'
#'library(phyloseq)
#'data("GlobalPatterns")
#'
#'GlobalPatterns %>% 
#'  phyloseq_ampvis_heatmap(physeq = ., 
#'                          transform = "compositional", 
#'                          facet_by = "SampleType" , 
#'                          group_by = "Primer",  
#'                          tax_aggregate = "Genus", 
#'                          tax_add = NULL, 
#'                          ntax =  5) -> heat_overall
#'                          
#'source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_varia.R")
#'GlobalPatterns %>% 
#' physeq_most_abundant(group_var = "SampleType",
#'                     ntax = 5,
#'                     tax_level = "Genus") -> top_taxa_per_group
#'
#'GlobalPatterns %>% 
#'  transform_sample_counts(function(x) x/sum(x) * 100) %>%
#'  subset_taxa(Genus %in% top_taxa_per_group) %>% 
#'  phyloseq_ampvis_heatmap(physeq = ., 
#'                          transform = FALSE, 
#'                          facet_by = "SampleType" , 
#'                          group_by = "Primer",  
#'                          tax_aggregate = "Genus", 
#'                          tax_add = NULL, 
#'                          ntax =  Inf) -> heat_top_taxa_per_group
#'
#'
#'
#'
phyloseq_ampvis_heatmap <- function(physeq,
                                    transform = 'compositional', 
                                    group_by, 
                                    facet_by, 
                                    tax_aggregate = FALSE, 
                                    tax_add = NULL, 
                                    ntax = 10,
                                    plot_values = TRUE,
                                    order_y_by = NULL,
                                    order_x_by = NULL)
{
  require(tidyverse)
  require(ampvis2)
  require(phyloseq)
  
  #devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6") # to source phyloseq_to_ampvis2()
  
  
  # if ('SampleID' !%in% (physeq %>% sample_data() %>% colnames()))
  # {
  physeq %>%
    sample_data() %>%
    data.frame() %>%
    rownames_to_column('temp') %>%
    mutate("SampleID" = temp) %>%
    dplyr::select(SampleID, everything()) %>%
    column_to_rownames("temp") %>%
    sample_data() -> df
  
  physeq@sam_data = NULL
  
  physeq <- merge_phyloseq(physeq,
                           df)
  # }else 
  # {
  #   physeq %>%
  #     sample_data() %>%
  #     as.matrix() %>%
  #     as.data.frame() %>%
  #     rownames_to_column('SampleID') -> df
  # }
  # if(any(rank_names(physeq) == "Strain")){
  # tax_table(physeq) <- tax_table(physeq)[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Strain")]
  # colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  # }
  if (transform != FALSE)
  {
    physeq %>%
      microbiome::transform(transform = transform) -> physeq
    # }else{
    #   physeq %>%
    #     filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> physeq
  }
  
  if ('Strain' %in% (physeq %>% tax_table() %>% colnames()))
  {
    tax_table(physeq) <- tax_table(physeq)[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Strain")]
    colnames(tax_table(physeq)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  }
  physeq %>%
    phyloseq_to_ampvis2() %>% 
    amp_heatmap(group_by = group_by,#treatment
                facet_by = facet_by,
                normalise = FALSE,
                plot_values = plot_values,
                plot_values_size = 2,
                tax_show = ntax,
                order_x_by = order_x_by,
                order_y_by = order_y_by,
                min_abundance = 0,
                tax_aggregate = tax_aggregate,
                tax_add = tax_add,
                plot_na = FALSE,
                color_vector = c("white", "red"),
                plot_colorscale = "sqrt",
                plot_legendbreaks = c(1, 10, 20)
    ) -> p 
  
  
  p + theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) + 
    theme(axis.text.y = element_text(angle = 0,  size = 8)) -> p
  
  return(p)
  
  detach("package:ampvis2", unload=TRUE)
  
  
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
