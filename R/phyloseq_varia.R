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
#'
#'
#'
#'
# library(phyloseq)
# data("soilrep")
# 
# enterotype %>% 
#   physeq_most_abundant(group_var = "warmed",
#                      ntax = 2,
#                      tax_level = "Genus") -> tmp


physeq_most_abundant <- function(physeq,
                                 group_var,
                                 ntax = 10,
                                 tax_level = "Species"){
  
  require(tidyverse); require(phyloseq)
  
  taxa_top_all = NULL  
  
  for(tp in physeq %>% 
      get_variable(group_var) %>%
      unique()){
    # print(tp)
    prune_samples(get_variable(physeq, group_var) == tp,
                  physeq) %>%
      fantaxtic::get_top_taxa(n = ntax, 
                              relative = FALSE, 
                              discard_other = TRUE) -> tmp2
    
    as(tax_table(tmp2), "matrix") %>%
      data.frame() %>%
      # dplyr::filter(grepl("ose",Gene_name)) %>%
      pull(tax_level) -> spc
    
    c(spc, taxa_top_all) %>%
      # discard(is.na)  %>%
      unique() %>%
      sort() -> taxa_top_all
    
  }
  return(taxa_top_all)
}
