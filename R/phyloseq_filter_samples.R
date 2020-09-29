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