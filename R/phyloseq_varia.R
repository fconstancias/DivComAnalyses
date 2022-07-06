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
#'source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_heatmap.R")
#'
#'library(phyloseq); library(tidyverse)
#'
#'data("GlobalPatterns")
#'
#'GlobalPatterns %>%
#'  phyloseq_ampvis_heatmap(physeq = .,
#'                          transform = "compositional", # = percentage transformation
#'                          facet_by = "SampleType" ,
#'                          group_by = "Primer",
#'                          tax_aggregate = "Species",
#'                          tax_add = NULL,
#'                          ntax =  5) -> heat_overall
#'
#'heat_overall
#'
#'
#'# How many sample types?
#'GlobalPatterns %>%
#'  physeq_get_unique("SampleType")
#'
#'
#'# 5 Most abundant Species based on ASV classification - i.e., not agglomerated at the Species level - per sample types:
#'GlobalPatterns %>%
#'  physeq_most_abundant(group_var = "SampleType",
#'                       ntax = 5,
#'                       tax_level = "Species") -> top_taxa_per_group
#'
#'# looking at the five most abondant per SampleType we obtain 8 Species
#'
#'top_taxa_per_group
#'
#'
#'GlobalPatterns %>%
#'  transform_sample_counts(function(x) x/sum(x) * 100) %>% # transform as percentage before filtering
#'  subset_taxa(Species %in% top_taxa_per_group) %>% # extract only the taxa to display - after percentage normalisation
#'  phyloseq_ampvis_heatmap(physeq = .,
#'                          transform = FALSE, # extract only the taxa to display - after percentage normalisation
#'                          facet_by = "SampleType" ,
#'                          group_by = "Primer",
#'                          tax_aggregate = "Species",
#'                          tax_add = NULL,
#'                          ntax =  Inf) -> heat_top_taxa_per_group
#'
#'heat_top_taxa_per_group

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
#'library(phyloseq)
#'
#'data("GlobalPatterns")
#'
#'# Applied on tax_table information
#'GlobalPatterns %>%
#'  physeq_get_unique("Kingdom")
#'
#'GlobalPatterns %>%
#'  physeq_get_unique("Phylum") %>%
#'  length()
#'
#'# Applied on sample_metadata
#'# First check the variable names
#'GlobalPatterns %>%
#'  sample_variables()
#'
#'# apply
#'GlobalPatterns %>%
#'  physeq_get_unique("Primer")
#'

physeq_get_unique <- function(ps, var){

  ps %>%
    speedyseq::psmelt() %>%
    distinct(get(var)) %>%  # could be anything: taxa, metadata, ...
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
#'library(phyloseq)
#'data("GlobalPatterns")
#'
#'# First check the variable names
#'GlobalPatterns %>%
#'  sample_variables()
#'
#'# apply
#'GlobalPatterns %>%
#'  generate_color_palette(var = "SampleType",
#'                         pal = "npg") -> my_pal
#'my_pal

generate_color_palette <- function(ps,
                                   var,
                                   seed = 123456,
                                   pal = "randomcoloR",
                                   runTsne = FALSE,
                                   altCol = FALSE,
                                   print = TRUE){

  ps %>%
    physeq_get_unique(var) -> var_uniq

  if(pal == "randomcoloR"){
    require(randomcoloR)
    set.seed(seed)
    var_uniq %>%
      length() %>%
      randomcoloR::distinctColorPalette(k = ., altCol = altCol, runTsne = runTsne) -> col

    detach("package:randomcoloR", unload=TRUE)

  }else{
    require(ggpubr)

    var_uniq %>%
      length() %>%
      ggpubr::get_palette(k = ., palette = pal) -> col

    detach("package:ggpubr", unload=TRUE)

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
#'ps_CARD %>%
#'physeq_simplify_tax(round_otu = TRUE, tax_sel = c("Best_Hit_ARO")) -> ps_AMRgn
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


#' @title Perform agglomeration at a particular level and rename OTU based on that level or rename at highest level
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
#'#'library(phyloseq)
#'data("GlobalPatterns")
#'GlobalPatterns %>%  physeq_glom_rename(taxrank = "Family") %>%  taxa_names()


physeq_glom_rename <- function(phyloseq,
                               speedyseq = FALSE,
                               taxrank = FALSE,
                               rename_ASV = taxrank,
                               taxnames_rm = c("unknown", "Incertae Sedis")){

  ##---------------------------------------------
  require(tidyverse); require(phyloseq)
  if(speedyseq == TRUE){require(speedyseq)}
  ##---------------------------------------------

  if (taxrank %in% rank_names(phyloseq)){
    phyloseq %>%
      tax_glom(taxrank = taxrank) -> phyloseq

    prune_taxa(data.frame(tax_table(phyloseq)[,taxrank])  %>%
                 dplyr::filter(!get(taxrank) %in% taxnames_rm) %>% rownames(),
               phyloseq) -> phyloseq

    # taxa_names(phyloseq) <-  tax_table(phyloseq)[,taxrank]

  }

  if (rename_ASV != FALSE){

    taxa_names(phyloseq) <-  tax_table(phyloseq)[,rename_ASV]
  }

  ##---------------------------------------------

  return(phyloseq)
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
#'ps_CARD %>%
#'physeq_simplify_tax(round_otu = TRUE, tax_sel = c("Best_Hit_ARO")) -> ps_AMRgn
#'
#'

physeq_sel_tax_table <- function(ps, tax_sel){

  # tax_table(ps)[,tax_sel] -> tax_table(ps)

  ps %>%
    tax_table() %>%
    data.frame() %>%
    rownames_to_column('tmp_id') %>%
    dplyr::filter(all_of(tmp_id,tax_sel)) %>%
    mutate_if(is.character,as.factor) -> tax_table(ps)

  return(ps)
}
