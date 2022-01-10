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



'%!in%' <- function(x,y)!('%in%'(x,y))

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

physeq_get_unique <- function(ps, var){
  
  ps %>% 
    speedyseq::psmelt() %>% 
    distinct( get(var)) %>%  # could be anything: taxa, metadata, ...
    pull() -> out
  
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
#'


generate_color_palette <- function(ps, var, seed = 123456, pal = "randomcoloR", runTsne = FALSE, altCol = FALSE, print = TRUE){
  
  ps %>% 
    physeq_get_unique(var) -> var_uniq
  
  if(pal == "randomcoloR"){
    set.seed(seed)
    var_uniq %>% 
      length() %>% 
      randomcoloR::distinctColorPalette(k = ., altCol = altCol, runTsne = runTsne) -> col
  }else{
    var_uniq %>% 
      length() %>% 
      ggpubr::get_palette(k = ., palette = pal) -> col
  }
  
  if(print == TRUE){
    pie(rep(1, length(col)),                                
        col = col)
  }
  
  names(col) <- var_uniq
  
  return(col)
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
#'
#'
#'
#'

physeq_simplify_tax <- function(ps, tax_sel, round_otu = FALSE){
  
  ps %>% 
    tax_table() %>% 
    data.frame() -> tax_mapping
  
  
  if(round_otu == TRUE){
    
    round(otu_table(ps)) -> otu_table(ps)
    
  }
  
  ps %>% 
    tax_table() %>% 
    data.frame() %>% 
    select(!!tax_sel) %>% 
    as.matrix() -> tax_table(ps)
  
  
  ps %>% 
    speedyseq::tax_glom(tax_sel[1]) -> ps_glom
  
  
  
  taxa_names(ps_glom) <- tax_table(ps_glom)[,tax_sel[1]]
  
  tax_sel[1] -> sel
  
  
  ps_glom %>% 
    tax_table() %>% 
    data.frame() %>% 
    left_join(tax_mapping %>% distinct(!!sel, .keep_all = TRUE),
              by = setNames(tax_sel[1], tax_sel[1]),
              suffix = c("_x", "")) %>% 
    # mutate(str_remove_all(:=) %>% 
    column_to_rownames(tax_sel[1])  %>% 
    as.matrix() -> tax_table(ps_glom)
  
  return(ps_glom)
}
