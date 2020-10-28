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
#'phyloseq_ampvis_heatmap(GlobalPatterns,"percent", "Primer","SampleType",  "Species", "Phylum", 10) -> out
#'
#'
#'
#'
#'
#'
phyloseq_ampvis_heatmap <- function(physeq,transform, group_by, facet_by, tax_aggregate, tax_add, ntax)
{
  require(tidyverse)
  require(ampvis2)
  require(phyloseq)
  
  devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6") # to source phyloseq_to_ampvis2()
  
  
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
    tax_table(physeq) <- tax_table(physeq)[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Strain")]
    colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    if (transform == "percent")
    {
  physeq %>%
    transform_sample_counts(function(x) x/sum(x) * 100) %>%
    filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> physeq
    }else{
      physeq %>%
        filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> physeq
    }
  
    if ('Strain' %in% (physeq %>% tax_table() %>% colnames()))
    {
  tax_table(physeq) <- tax_table(physeq)[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Strain")]
  colnames(tax_table(physeq)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    }else{
      
    }
  physeq %>%
    phyloseq_to_ampvis2() %>% 
    amp_heatmap(group_by = group_by,#treatment
                facet_by = facet_by,
                normalise = FALSE,
                plot_values = TRUE,
                plot_values_size = 2,
                tax_show = ntax,
                # order_x_by = "cluster",
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
  
}
