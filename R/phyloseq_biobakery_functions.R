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
#'metaphlan_2phyloseq(merged_metaphlan = here::here("data/processed/metaphlan/default_merged_metagenome.txt"),
#'metadata = here::here("data/metadata_all_DNA_RNA.xlsx") %>%
#'readxl::read_xlsx() %>%
#'  rownames_to_column('tmp') %>%
#'  filter(Type == "DNA") %>%
#'  column_to_rownames('Sample')) -> ps

#'metaphlan_2phyloseq(merged_metaphlan = here::here("data/processed/humann/merged_metaphlan3_humann3.tsv"),
#'                    metadata = here::here("data/metadata_all_DNA_RNA.xlsx") %>%
#'                      readxl::read_xlsx() %>%
#'                      rownames_to_column('tmp') %>%
#'                      filter(Type == "DNA") %>%
#'                      column_to_rownames('Sample'),
#'                    skip_col = 1,
#'                    id = "clade_name") -> ps

metaphlan_2phyloseq <- function(merged_metaphlan,
                                metadata,
                                tree = "https://raw.githubusercontent.com/biobakery/MetaPhlAn/master/metaphlan/utils/mpa_v30_CHOCOPhlAn_201901_species_tree.nwk",
                                skip_col = 0,
                                id = "#SampleID",
                                tax_label = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                                tax_sep = "\\|"){
  
  ## ------------------------------------------------------------------------
  require(tidyverse); require(speedyseq)
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  cat(paste0('\n##',"You are using speedyseq version ", packageVersion('speedyseq'),'\n\n'))
  `%!in%` = Negate(`%in%`)
  
  cat('################################\n\n')
  ## ------------------------------------------------------------------------
  merged_metaphlan %>% 
    read_tsv(col_names = TRUE,
             skip = skip_col) %>%
    dplyr::select_if(names(.) %!in% c('NCBI_tax_id')) %>%
    tidyr::separate(id,tax_label ,sep = tax_sep) %>%
    dplyr::filter(!is.na(Species)) -> df
  ## ------------------------------------------------------------------------
  
  df %>%
    dplyr::select_if(is_character) %>%
    as.matrix() -> tax
  
  df %>%
    dplyr::select_if(is.double) %>% 
    as.matrix() -> count
  ## ------------------------------------------------------------------------
  
  merge_phyloseq(otu_table(count, taxa_are_rows = TRUE),
                 tax_table(tax)) -> physeq
  
  tax_table(physeq) <- tax_table(physeq) %>% gsub(pattern="[a-s]__",replacement="")
  
  taxa_names(physeq) <- tax_table(physeq)[,"Species"]
  ## ------------------------------------------------------------------------
  if (file.exists(metadata) == TRUE){
    merge_phyloseq(physeq,
                   metadata %>% phyloseq::sample_data()) -> physeq
  }
  if (tree != FALSE){
    tree %>%
      ape::read.tree() -> tree_file
    
    tree_file$tip.label <- gsub(".+\\|s__", "", tree_file$tip.label)
    
    filt_tree <- ape::keep.tip(tree_file, intersect(taxa_names(physeq),tree_file$tip.label))
    
    merge_phyloseq(physeq,
                   filt_tree %>% phy_tree()) -> physeq
  }
  return(physeq) 
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
#' here::here("data/processed/humann/DNA/genefamilies_joined_tables_uniref90_ko_renamed_kegg-orthology.tsv")
#' 
#'
#'
#'

humann_2df <- function(humann_renamed = here::here("data/processed/humann/DNA/genefamilies_joined_tables_uniref90_ko_renamed_kegg-orthology.tsv"),
                       type = '# Gene Family',#`# Gene Family` | `# Pathway`)
                       n_rows = Inf) # for testing purpose
{
  
  ## ------------------------------------------------------------------------
  require(tidyverse)
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  
  cat('################################\n\n')
  ## ------------------------------------------------------------------------
  humann_renamed %>% 
    readr::read_tsv(col_names = TRUE, n_max = n_rows) %>%
    dplyr::rename(Feature = !!type) %>%
    dplyr::mutate(id = Feature) %>%
    dplyr::select(id, everything(.)) %>%
    tidyr::separate(Feature, c("Feature", "organism"), sep = "\\|", fill = "right") %>%
    tidyr::separate(organism, c("Genus", "Species"), sep = "\\.", fill = "right", remove = FALSE) %>%
    dplyr::mutate(organism = str_replace(organism, ".s__", "_")) %>% 
    dplyr::mutate(organism = str_replace(organism, "g__|", "")) %>% 
    dplyr::mutate(Genus = str_replace(Genus, "g__", "")) %>% 
    dplyr::mutate(Species = str_replace(Species, "s__", "")) -> df
  
  ## ------------------------------------------------------------------------
  
  return(df)
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
#' here::here("data/processed/humann/DNA/genefamilies_joined_tables_uniref90_ko_renamed_kegg-orthology.tsv") %>% humann_2df() %>% clean_humann_df() -> dff
#' 
#'
#'
#'

clean_humann_df <- function(humann_df){
  ## ------------------------------------------------------------------------
  
humann_df %>%
    dplyr::mutate(Feature = str_replace(Feature, "biosynthesis", "bios.")) %>%
    dplyr::mutate(Feature = str_replace(Feature, "degradation", "deg.")) %>%
    dplyr::mutate(Feature = str_replace(Feature, "superpathway", "spw")) %>%
    dplyr::mutate(Feature = gsub("\\s*\\([^\\)]+\\)","",as.character(Feature))) -> df
  ## ------------------------------------------------------------------------
  
  return(df)
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
#'here::here("data/processed/humann/DNA/pathabundance_cpm_joined_tables.tsv") %>% humann_2df(type = '# Pathway') -> DNA
#'here::here("data/processed/humann/RNA/pathabundance_cpm_joined_tables.tsv") %>% humann_2df(type = '# Pathway') -> RNA
#'
#'humann_DNA_RNA_2phyloseq(DNA,RNA) -> ps

humann_DNA_RNA_2phyloseq <- function(DNA_humann_2df,
                                     RNA_humann_2df)
{
  ## ------------------------------------------------------------------------
  require(tidyverse); require(speedyseq)
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  cat(paste0('\n##',"You are using speedyseq version ", packageVersion('speedyseq'),'\n\n'))
  
  cat('################################\n\n')
  ## ------------------------------------------------------------------------
  DNA_humann_2df %>%
    dplyr::rename_with( ~ paste( .x, "DNA", sep = "_")) -> DNA
  RNA_humann_2df %>%
    dplyr::rename_with( ~ paste( .x, "RNA", sep = "_")) -> RNA
  
  dplyr::left_join(DNA,
                   RNA, 
                   by= c("id_DNA"="id_RNA"), 
                   suffix = c("_DNA", "_RNA"), 
                   keep = FALSE) %>% #-> DNA_RNA
    dplyr::mutate_if(is.numeric, ~ replace(., is.na(.), 0)) %>%
    dplyr::mutate(Species_DNA = if_else(is.na(Species_DNA), "unclassified", Species_DNA)) -> DNA_RNA
  
  DNA_RNA %>%
    select(c("id_DNA", "Feature_DNA","organism_DNA", "Genus_DNA", "Species_DNA")) %>%
    # separate(col=Description_DNA,
    #      into=c("Description","EC_number"),
    #      sep = "_EC_", remove = TRUE) %>%
    dplyr::rename(id = id_DNA,
                  Feature = Feature_DNA, 
                  organism = organism_DNA, 
                  Genus = Genus_DNA, 
                  Species = Species_DNA) -> tax
  
  DNA_RNA %>%
    dplyr::select_if(is.numeric) -> count
  # colnames(df) <- gsub("[^[:alnum:] ]", "",colnames(df))
  
  if(dim(tax)[1] != dim(count)[1]) stop ("Something went wrong...")
  
  
  merge_phyloseq(otu_table(count %>% as.matrix(), 
                           taxa_are_rows = TRUE),
                 tax_table(tax %>% as.matrix())) -> physeq
  ## ------------------------------------------------------------------------
  
  return(physeq)
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
#'here::here("data/processed/humann/DNA/genefamilies_joined_tables_uniref90_ko_renamed_kegg-orthology.tsv") %>% humann_2df() %>% clean_humann_df() %>% humann_2phyloseq() -> ps
#'

humann_2phyloseq <- function(humann_2df)
{
  ## ------------------------------------------------------------------------
  require(tidyverse); require(speedyseq)
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  cat(paste0('\n##',"You are using speedyseq version ", packageVersion('speedyseq'),'\n\n'))
  
  cat('################################\n\n')
  ## ------------------------------------------------------------------------
  humann_2df %>%
    dplyr::select_if(is_character) -> tax
  
  humann_2df %>%
    dplyr::select_if(is.numeric) -> count
  # colnames(df) <- gsub("[^[:alnum:] ]", "",colnames(df))
  
  if(dim(tax)[1] != dim(count)[1]) stop ("Something went wrong...")
  
  
  merge_phyloseq(otu_table(count %>% as.matrix(), 
                           taxa_are_rows = TRUE),
                 tax_table(tax %>% as.matrix())) -> physeq
  ## ------------------------------------------------------------------------
  
  return(physeq)
}

#' @title ...
#' @param .
#' @param ..
#' @author Florentin Constancias
#' @note .
#' @note .
#' @note .
#' #'As seen here <https://github.com/fconstancias/omnibus-and-maaslin2-rscripts-and-hmp2-data/blob/patch-1/Maaslin2.R>
#'  - No unmapped and unintegrated pathway 
#' - ~ 1 normalised.
#' --> we have to renormalise
#'--> get rid of stratification
#' <https://forum.biobakery.org/t/compatibility-of-humann2-output-files-for-maaslin2/412/9>
#'  **However, they (unmapped/unintegrated) can be useful to retain for linear modeling as they limit the potential for housekeeping functions to artificially inflate in less-well-characterized communities.**
#' @return .
#' @export
#' @examples
#'


phyloseq_get_humann_strat_un_output <- function(physeq,
                                                output = "stratified", # stratified / unstratified
                                                remove_unmapped_unintegrated = FALSE, 
                                                transform = "compositional",# from microbiome:: 'compositional' (ie relative abundance), 'Z', 'log10', 'log10p', 'hellinger', 'identity', 'clr', or any method from the vegan::decostand function.
                                                export_long_df = TRUE)
  ## ------------------------------------------------------------------------

{
  if (remove_unmapped_unintegrated == TRUE){
    physeq %>%
      subset_taxa(!(grepl("^UN", Gene))) -> physeq
  }
  if (output == "stratified"){
    physeq %>%
      subset_taxa(!is.na(organism)) -> physeq
  }
  if (output == "unstratified"){
    physeq %>%
      subset_taxa(is.na(organism)) -> physeq
  }
  
  if (transform != FALSE){
    physeq %>%
      microbiome::transform(transform) -> physeq
  }
  
  ## ------------------------------------------------------------------------
  
  return(physeq)
  
  if (export_long_df == TRUE){
    out <- list("physeq" = physeq,
                "df" = physeq %>% 
                  speedyseq::psmelt())
    return(out)
    
  }
}

