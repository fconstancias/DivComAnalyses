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
#'here::here("data/processed/humann/DNA/pathabundance_cpm_joined_tables.tsv") %>% humann_2df(type = '# Pathway') -> DNA_humann_2df
#'here::here("data/processed/humann/RNA/pathabundance_cpm_joined_tables.tsv") %>% humann_2df(type = '# Pathway') -> RNA_humann_2df
#'
#'humann_DNA_RNA_2phyloseq(DNA_humann_2df,RNA_humann_2df) -> ps

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
    dplyr::rename_with( ~ paste( .x, "DNA", sep = "_")) %>%
    dplyr::rename(id = id_DNA) -> DNA
  
  RNA_humann_2df %>%
    dplyr::rename_with( ~ paste( .x, "RNA", sep = "_")) %>%
    dplyr::rename(id = id_RNA) %>%
    dplyr::mutate(id_RNA = id) -> RNA
  
  dplyr::full_join(DNA,
                   RNA, 
                   by= c("id"="id"), 
                   suffix = c("", "_RNA"), 
                   keep = FALSE) %>% #-> DNA_RNA
    dplyr::mutate_if(is.numeric, ~ replace(., is.na(.), 0)) %>% # **tight here, we need to combine the column from DNA and RNA since not all DNA features are in RNA. sth like mutate if**
    dplyr::mutate(Species_DNA = if_else(is.na(Species_DNA), "unclassified", Species_DNA)) %>%
    dplyr::rename(#id = id_DNA,
                  Feature = Feature_DNA, 
                  organism = organism_DNA, 
                  Genus = Genus_DNA, 
                  Species = Species_DNA) -> DNA_RNA
  
  cat(paste0('##',"DNA table: ",ncol(DNA %>% dplyr::select_if(is.numeric)),' samples and ',nrow(DNA),' features'))
  cat(paste0('##',"RNA table: ",ncol(RNA %>% dplyr::select_if(is.numeric)),' samples and ',nrow(RNA),' features'))
  cat(paste0('##',"Merged table: ",ncol(DNA_RNA %>% dplyr::select_if(is.numeric)),' samples and ',nrow(DNA_RNA),' features \n'))
  
  DNA_RNA %>%
    select(c("id", "Feature","organism", "Genus", "Species")) %>%
    # separate(col=Description_DNA,
    #      into=c("Description","EC_number"),
    #      sep = "_EC_", remove = TRUE) %>%
    # dplyr::rename(id = id_DNA,
    #               Feature = Feature_DNA, 
    #               organism = organism_DNA, 
    #               Genus = Genus_DNA, 
    #               Species = Species_DNA) %>%
    column_to_rownames('id') -> tax
  
  DNA_RNA %>%
    dplyr::select_if(is.numeric) %>%
    colnames() -> samples
  
  DNA_RNA %>%
    dplyr::select(c(id, samples)) %>%
    column_to_rownames('id') -> count
  # colnames(df) <- gsub("[^[:alnum:] ]", "",colnames(df))
  
  if(dim(tax)[1] != dim(count)[1]) stop ("Something went wrong...")
  
  
  merge_phyloseq(otu_table(count %>% as.matrix(), 
                           taxa_are_rows = TRUE),
                 tax_table(tax %>% as.matrix())) -> physeq
  
  
  cat(paste0('##',"Created phyloseq object ",nsamples(physeq),' samples and ', ntaxa(physeq),' features'))


  if (visualize == TRUE){
    
    # ggVennDiagram::ggVennDiagram(list("RNA" = RNA$id,
    #                                   "DNA" = DNA$id)) -> p
    
    ggVennDiagram::ggVennDiagram(list("DNA" = DNA_RNA$id,
                                      "RNA" = DNA_RNA$id_RNA_RNA)) -> p
    
    out <- list("physeq" = physeq,
                "venn" = p)
    return(out)
    
  }else{
    return(physeq)
  }

  
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
    dplyr::select_if(is_character) %>%
    column_to_rownames('id') -> tax
  
  humann_2df %>%
    dplyr::select_if(is.numeric) %>%
    colnames() -> samples
  
  humann_2df %>%
    dplyr::select(c(id, samples)) %>%
    column_to_rownames('id') -> count
  
  # colnames(df) <- gsub("[^[:alnum:] ]", "",colnames(df))
  
  if(dim(tax)[1] != dim(count)[1]) stop ("Something went wrong...")
  
  
  merge_phyloseq(otu_table(count %>% as.matrix(), 
                           taxa_are_rows = TRUE),
                 tax_table(tax %>% as.matrix())) -> physeq
  
  cat(paste0('##',"Created phyloseq object ",nsamples(physeq),' samples and ', ntaxa(physeq),' features'))
  
  # taxa_names(physeq) <- speedyseq::tax_table(physeq)[,"id"]
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
#'here::here("data/processed/humann/DNA/genefamilies_joined_tables_uniref90_ko_renamed_kegg-orthology.tsv") %>% humann_2df() %>% clean_humann_df() %>% humann_2phyloseq() -> tmp
#'tmp %>% phyloseq_get_humann_strat_un_output(output = "unstratified" ,transform = "clr",  export_long_df = TRUE, remove_unmapped_unintegrated = TRUE) -> ps



phyloseq_get_humann_strat_un_output <- function(physeq,
                                                output = "stratified", # stratified / unstratified
                                                remove_unmapped_unintegrated = FALSE, 
                                                transform = "compositional",# from microbiome:: 'compositional' (ie relative abundance), 'Z', 'log10', 'log10p', 'hellinger', 'identity', 'clr', or any method from the vegan::decostand function.
                                                export_long_df = TRUE){
  ## ------------------------------------------------------------------------
  require(tidyverse); require(phyloseq)
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  cat(paste0('\n##',"You are using phyloseq version ", packageVersion('phyloseq'),'\n\n'))
  
  cat('################################\n\n')
  ## ------------------------------------------------------------------------  
  if (remove_unmapped_unintegrated == TRUE){
    physeq %>%
      subset_taxa(!(grepl("^UN", Feature))) -> physeq
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
  
  
  
  if (export_long_df == TRUE){
    
    physeq %>% 
      speedyseq::psmelt() -> df
    out <- list("physeq" = physeq,
                "df" = df)
    return(out)
    
  }else{
    return(physeq)
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
#' here::here("data/processed/humann/DNA/genefamilies_joined_tables_uniref90_ko_renamed_kegg-orthology.tsv") %>% humann_2df() %>% clean_humann_df() %>% humann_2phyloseq() %>% phyloseq_get_humann_strat_un_output(output = "unstratified" ,transform = "clr",  export_long_df = TRUE, remove_unmapped_unintegrated = TRUE) -> ps
#' 
#' sample_names(ps$physeq) <- str_replace(sample_names(ps$physeq), "_DNA_cat_Abundance-RPKs", "")
#' 
#' physeq_add_metadata(ps$physeq, here::here("data/metadata_all_DNA_RNA.xlsx") %>%  readxl::read_xlsx() %>% filter(Type == "DNA"), sample_column = "Sample") -> test

physeq_add_metadata <- function(physeq,
                                metadata,
                                sample_column = "Sample"){
  
  ## ------------------------------------------------------------------------
  require(tidyverse); require(phyloseq)
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  cat(paste0('\n##',"You are using phyloseq version ", packageVersion('phyloseq'),'\n\n'))
  
  cat('################################\n\n')
  ## ------------------------------------------------------------------------  
  
  phyloseq::merge_phyloseq(physeq,
                           metadata %>%
                             column_to_rownames( sample_column ) %>% 
                             phyloseq::sample_data()) -> physeq
  
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
#' here::here("data/processed/humann/DNA/genefamilies_joined_tables_uniref90_ko_renamed_kegg-orthology.tsv") %>% humann_2df() %>% clean_humann_df() %>% humann_2phyloseq() %>% phyloseq_get_humann_strat_un_output(output = "stratified" ,transform = "compositional",  export_long_df = FALSE, remove_unmapped_unintegrated = TRUE) -> ps
#'
#' sample_names(ps) <- str_replace(sample_names(ps), "_DNA_cat_Abundance-RPKs", "")
#' 
#' physeq_add_metadata(ps, here::here("data/metadata_all_DNA_RNA.xlsx") %>%  readxl::read_xlsx() %>% filter(Type == "DNA"), sample_column = "Sample") -> physeq
#' 
here::here("data/processed/humann/DNA/pathabundance_cpm_joined_tables.tsv") %>% humann_2df(type = '# Pathway') -> DNA; here::here("data/processed/humann/RNA/pathabundance_cpm_joined_tables.tsv") %>% humann_2df(type = '# Pathway') -> RNA
humann_DNA_RNA_2phyloseq(DNA,RNA) -> physeq
sample_names(physeq) <- str_replace(sample_names(physeq), "_DNA_cat_Abundance-CPM_DNA", "_DNA"); sample_names(physeq) <- str_replace(sample_names(physeq), "_cat_Abundance-CPM_RNA", "_RNA")
physeq_add_metadata(physeq, here::here("data/metadata_all_DNA_RNA.xlsx") %>%  readxl::read_xlsx(), sample_column = "Sample_ID") -> physeq


humann2_species_contribution <- function(physeq,
                                         meta_data_var = c("Sample", "Subject", "Type", "Oral_Site", "Health_status"),
                                         filter_cut = 0,
                                         type = "RNA_DNA")
  {
  
  
}

physeq %>%
  subset_taxa(!(grepl("^UN", Feature))) %>% # to make sure we are working with stratefied data
  subset_taxa(!is.na(organism)) %>% # to make sure we are working with stratefied data
  speedyseq::psmelt() %>%
  dplyr::select(-organism) %>%
  dplyr::select(id, Feature, 
                Genus, Species, meta_data_var,
                Abundance) -> long_strat

  long_strat %>%
    dplyr::filter(Abundance > filter_cut) %>%
    tidyr::pivot_wider(names_from  = Type,
                       values_from = Abundance,
                       values_fill = list(Abundance = 0)) -> tmp
  
if (type == "RNA_DNA"){
  tmp %>%
    dplyr::mutate(RNA_DNA = RNA/DNA)  -> strat_DNA_RNA
}


  strat_DNA_RNA %>%
    filter(Feature %in% "HISDEG-PWY: L-histidine degradation I")
  filter(Species %in% c("Streptococcus_parasanguinis")) %>%
  filter(DNA  > 0 , RNA > 0) -> toto

ggpubr::compare_means(RNA_DNA ~ Health_status,
                      group.by = c("Oral_Site","Gene"),
                      data = toto,
                      method = "wilcox.test",
                      p.adjust.method = "fdr") %>%
  filter(p.adj < 0.05) %>%
  select(Gene, group1, group2, p.adj) -> RNA_DAN_signif

RNA_DAN_signif %>% pull(Gene) %>% unique() %>% as.vector() -> RNA_DAN_signif_gn

RNA_DAN_signif %>%
  DT::datatable()



toto %>%
  filter(Gene %in% RNA_DAN_signif_gn) %>%
  ggplot(aes(x = log10(DNA), y = log10(RNA),
             fill = Gene, color = Gene)) + 
  geom_point(alpha = 0.8, show.legend = FALSE, size = 1)  +
  facet_grid(Oral_Site ~ Health_status, drop=TRUE,
             scale="fixed",space="free_x") +
  ggConvexHull::geom_convexhull(aes(group = Gene), alpha = 0.1,
                                size = 0.08, linetype = "dotted",
                                show.legend = FALSE) +
  # geom_polygon(aes(group = Species), alpha= 0.11,
  #              size = 0.08, linetype = "dotted",
  #              show.legend = FALSE, rule = "winding") +
  geom_abline(slope=1, linetype = "dashed", size = 0.08, intercept=0) +
  theme_light() +
  # coord_equal() +
  ggrepel::geom_text_repel(cex = 1.5,
                           force = 2,
                           aes(label=Gene),
                           show.legend = FALSE,
                           segment.size = 0.05, segment.alpha = 0.5,
                           direction = "both",
                           # nudge_x = 1,
                           nudge_y = 0,
                           toto %>%
                             filter(Subject == "K8")) +
  guides(fill=guide_legend(ncol= 1)) +
  theme(legend.key.size = unit(0.2,"cm")) -> p

p