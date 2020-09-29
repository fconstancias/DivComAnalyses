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
#'
#'
#'
#'
phyloseq_rarefaction_curves<- function (physeq, stepsize, color_data, facet_data)
{
  require(ampvis2)
  devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6") # to source phyloseq_to_ampvis2()
  
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