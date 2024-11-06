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

# metaphlan_2phyloseq(merged_metaphlan = here::here("data/processed/humann/merged_metaphlan3_humann3.tsv"),
#                    metadata = here::here("data/metadata_all_DNA_RNA.xlsx") %>%
#                      readxl::read_xlsx() %>%
#                      rownames_to_column('tmp') %>%
#                      filter(Type == "DNA") %>%
#                      column_to_rownames('Sample'),
#                    skip_col = 1,
#                    id = "clade_name") -> ps

metaphlan_2phyloseq <- function(merged_metaphlan="~/metaphlan_analysis/merged_abundance_table.txt",
                                metadata="no",
                                rm_unclassified = FALSE,
                                tree = "https://raw.githubusercontent.com/biobakery/MetaPhlAn/refs/heads/master/metaphlan/utils/mpa_vJun23_CHOCOPhlAnSGB_202403.nwk",
                                skip_col = 0,
                                metaphlan_sample_names_to_rm = "",
                                tax_label = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"),
                                tax_sep = "\\|"){
  
  ## ------------------------------------------------------------------------
  require(tidyverse);  require(phyloseq)
  
  `%!in%` = Negate(`%in%`)
  ## ------------------------------------------------------------------------
  merged_metaphlan %>%
    read.table(sep ="\t", header  = TRUE) %>%
    # read_tsv(col_names = TRUE,
    # skip = skip_col) %>%
    # dplyr::select_if(names(.) %!in% c('NCBI_tax_id')) %>%
    dplyr::filter(., grepl('t__|UNCLASSIFIED', clade_name)) %>%
    # dplyr::filter(clade_name != NA) %>% 
    # dplyr::filter(., !grepl('t__', clade_name)) %>%
    tidyr::separate(clade_name,tax_label ,sep = tax_sep)  -> df
  # dplyr::filter(!is.na(Strain) | !Kingdom == "UNCLASSIFIED") -> df
  
  if(rm_unclassified == TRUE){
    df %>% 
      dplyr::filter(!is.na(Strain) | !Kingdom == "UNCLASSIFIED") -> df
  }
  
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
  
  tax_table(physeq) <- tax_table(physeq) %>% gsub(pattern="[a-t]__",replacement="") %>%  data.frame() %>%  replace(is.na(.), "UNCLASSIFIED") %>%  as.matrix() %>%  tax_table()
  
  ## ------------------------------------------------------------------------
  
  taxa_names(physeq) <- paste0(tax_table(physeq)[,"Species"], "_", tax_table(physeq)[,"Strain"])
  
  ## ------------------------------------------------------------------------
  
  # physeq %>%
  #   clean_phyloseq_sample_names(sub_pat = metaphlan_sample_names_to_rm) -> physeq
  
  ## ------------------------------------------------------------------------
  if (file.exists(metadata) == TRUE){
    merge_phyloseq(physeq,
                   metadata %>% phyloseq::sample_data()) -> physeq
  }
  # if (tree != FALSE){
  #   tree %>%
  #     ape::read.tree() -> tree_file
  # 
  #   tree_file$tip.label <- gsub(".+\\|s__", "", tree_file$tip.label)
  # 
  #   filt_tree <- ape::keep.tip(tree_file, intersect(taxa_names(physeq),tree_file$tip.label))
  # 
  #   merge_phyloseq(physeq,
  #                  filt_tree %>% phy_tree()) -> physeq
  # }
  return(physeq)
}
## ------------------------------------------------------------------------------------------------------------------------------------------------

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

clean_phyloseq_sample_names <- function(physeq,
                                        sub_pat =  "_DNA_cat_Abundance-CPM_DNA",
                                        str_replace = ""){

  ## ------------------------------------------------------------------------
  require(tidyverse); require(phyloseq)

  physeq %>%
    sample_names() %>%
    # gsub("[^-]*-(.*)", "\\1", .) %>%
    sub(sub_pat, str_replace, .) -> new_names

  sample_names(physeq) <- new_names

  #sample_names(physeq) <- str_replace(sample_names(physeq), str_rm, str_replace)
  ## ------------------------------------------------------------------------
  return(physeq)
}


## ------------------------------------------------------------------------------------------------------------------------------------------------

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

## ------------------------------------------------------------------------------------------------------------------------------------------------

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

## ------------------------------------------------------------------------------------------------------------------------------------------------

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
                                     RNA_humann_2df,
                                     visualize = FALSE)
{
  ## ------------------------------------------------------------------------
  require(tidyverse); require(speedyseq)

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

  cat(paste0('##',"DNA table: ",ncol(DNA %>% dplyr::select_if(is.numeric)),' samples and ',nrow(DNA),' features \n'))
  cat(paste0('##',"RNA table: ",ncol(RNA %>% dplyr::select_if(is.numeric)),' samples and ',nrow(RNA),' features \n'))
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

## ------------------------------------------------------------------------------------------------------------------------------------------------

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

## ------------------------------------------------------------------------------------------------------------------------------------------------

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
                                                export_long_df = TRUE,
                                                rm_un_sp = TRUE){
  ## ------------------------------------------------------------------------
  require(tidyverse); require(phyloseq)

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

  if(rm_un_sp == TRUE){
    physeq %>%
      subset_taxa(Species != "unclassified") -> physeq
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
## ------------------------------------------------------------------------------------------------------------------------------------------------


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

# here::here("data/processed/humann/DNA/pathabundance_cpm_joined_tables.tsv") %>% humann_2df(type = '# Pathway') -> DNA; here::here("data/processed/humann/RNA/pathabundance_cpm_joined_tables.tsv") %>% humann_2df(type = '# Pathway') -> RNA
# humann_DNA_RNA_2phyloseq(DNA,RNA) -> physeq
# sample_names(physeq) <- str_replace(sample_names(physeq), "_DNA_cat_Abundance-CPM_DNA", "_DNA"); sample_names(physeq) <- str_replace(sample_names(physeq), "_cat_Abundance-CPM_RNA", "_RNA")
#
# physeq_add_metadata(physeq,
#                     here::here("data/metadata_all_DNA_RNA.xlsx") %>% readxl::read_xlsx() %>% dplyr::rename(Sample_ID2 = Sample),
#                     sample_column = "Sample_ID") %>%
#   humann2_species_contribution(meta_data_var = c("Subject", "Sample_ID2" ,"Type", "Oral_Site", "Health_status")) -> df


humann2_species_contribution <- function(physeq,
                                         meta_data_var,
                                         DNA_RNA_meta = "Type",
                                         filter_cut = 0,
                                         transform = 'identity')#transform after removeinf ^un and na(organisms)
{
  physeq %>%
    subset_taxa(!(grepl("^UN", Feature))) %>% # to make sure we are working with stratefied data
    subset_taxa(!is.na(organism)) %>% # to make sure we are working with stratefied data
    microbiome::transform(transform = transform) %>%
    speedyseq::psmelt() %>%
    dplyr::select(-organism) %>%
    dplyr::rename(id = OTU) %>%
    dplyr::select(id, Feature,
                  Genus, Species, meta_data_var,
                  Abundance) %>%
    dplyr::filter(Abundance > filter_cut) %>%
    tidyr::pivot_wider(names_from  = DNA_RNA_meta,
                       values_from = Abundance,
                       values_fill = list(Abundance = 0)) %>%
    dplyr::mutate(RNA_DNA = RNA/DNA) -> strat_DNA_RNA

  return(strat_DNA_RNA)
}


## ------------------------------------------------------------------------------------------------------------------------------------------------

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
# here::here("data/processed/humann/DNA/pathabundance_cpm_joined_tables.tsv") %>% humann_2df(type = '# Pathway') -> DNA; here::here("data/processed/humann/RNA/pathabundance_cpm_joined_tables.tsv") %>% humann_2df(type = '# Pathway') -> RNA
# humann_DNA_RNA_2phyloseq(DNA,RNA) -> physeq
# sample_names(physeq) <- str_replace(sample_names(physeq), "_DNA_cat_Abundance-CPM_DNA", "_DNA"); sample_names(physeq) <- str_replace(sample_names(physeq), "_cat_Abundance-CPM_RNA", "_RNA")
# physeq_add_metadata(physeq,
#                     here::here("data/metadata_all_DNA_RNA.xlsx") %>% readxl::read_xlsx() %>% dplyr::rename(Sample_ID2 = Sample),
#                     sample_column = "Sample_ID") %>%
#   humann2_species_contribution(meta_data_var = c("Subject", "Sample_ID2" ,"Type", "Oral_Site", "Health_status")) -> df
# df %>%
#   humann2_RNA_DNA_plot(facet_formula = "Oral_Site ~ Health_status",
#                        filter_feature = c("PWY-6737: starch degradation V",
#                                           "PWY-7111: pyruvate fermentation to isobutanol (engineered)"),
#                        filter_species = c("Streptococcus_oralis", )) -> p

humann2_RNA_DNA_plot <- function(df,
                                 rm_un_sp = TRUE,
                                 filter_feature = FALSE,
                                 x_plot = "log10(DNA)",
                                 y_plot = "log10(RNA)",
                                 color = "Feature",
                                 fill = "Feature",
                                 group = "Feature",
                                 shape = NULL,
                                 filter_genus = FALSE,
                                 filter_species = FALSE,
                                 only_pos = TRUE, # dplyr::filter(DNA  > 0 , RNA > 0)
                                 facet_formula = FALSE){

  require(tidyverse); require(ggConvexHull)

  if(rm_un_sp == TRUE){
    df %>%
      dplyr::filter(Species != "unclassified") -> df
  }
  if(filter_feature != FALSE){
    df %>%
      dplyr::filter(Feature %in% filter_feature) -> df
  }

  if(filter_genus != FALSE){
    df %>%
      dplyr::filter(Genus %in% filter_genus) -> df
  }

  if(filter_species != FALSE){
    df %>%
      dplyr::filter(Species %in% filter_species) -> df
  }

  if(only_pos == TRUE){
    df %>%
      dplyr::filter(DNA  > 0 , RNA > 0) -> df
  }

  df %>%
    ggplot(aes_string(x = x_plot, y = y_plot, color = color, fill = fill, shape = shape)) +
    geom_point(alpha = 0.8, show.legend = TRUE, size = 1)  +
    ggConvexHull::geom_convexhull(aes_string(group = group), alpha = 0.1,
                                  size = 0.08, linetype = "dotted",
                                  show.legend = FALSE) +
    geom_abline(slope=1, linetype = "dashed", size = 0.08, intercept = 0) +
    theme_light() +
    # coord_equal() +
    # ggrepel::geom_text_repel(cex = 1.5,
    #                          force = 2,
    #                          aes(label=Feature),
    #                          show.legend = FALSE,
    #                          segment.size = 0.05, segment.alpha = 0.5,
    #                          direction = "both",
    #                          # nudge_x = 1,
    #                          nudge_y = 0,
    #                          toto %>%
    #                            filter(Subject == "K8")) +
  # guides(fill=guide_legend(ncol= 1)) +
  theme(legend.key.size = unit(0.2,"cm")) -> p

  if(facet_formula != FALSE){
    p + facet_grid(as.formula(facet_formula), drop=TRUE,
                   scale="fixed",space="free_x") -> p
  }

  out <- list("legend" = p %>% ggpubr::get_legend() %>% ggpubr::as_ggplot(),
              "plot" = p + theme(legend.position = "none"))

  return(out)
}

## ------------------------------------------------------------------------------------------------------------------------------------------------

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

# here::here("data/processed/humann/DNA/pathabundance_cpm_joined_tables.tsv") %>% humann_2df(type = '# Pathway') -> DNA; here::here("data/processed/humann/RNA/pathabundance_cpm_joined_tables.tsv") %>% humann_2df(type = '# Pathway') -> RNA
# humann_DNA_RNA_2phyloseq(DNA,RNA) -> physeq
# sample_names(physeq) <- str_replace(sample_names(physeq), "_DNA_cat_Abundance-CPM_DNA", "_DNA"); sample_names(physeq) <- str_replace(sample_names(physeq), "_cat_Abundance-CPM_RNA", "_RNA")
# physeq_add_metadata(physeq,
#                     here::here("data/metadata_all_DNA_RNA.xlsx") %>% readxl::read_xlsx() %>% dplyr::rename(Sample_ID2 = Sample),
#                     sample_column = "Sample_ID") %>%
#   humann2_species_contribution(meta_data_var = c("Subject", "Sample_ID2" ,"Type", "Oral_Site", "Health_status")) %>%
#   humann2_RNA_DNA_ratio_plot(x_plot = "Oral_Site",
#                              y_plot = "log10(RNA_DNA)",
#                              color = "Oral_Site",
#                              fill = "Oral_Site",
#                              filter_feature = c("PWY-6737: starch degradation V",                                                                                                           "PWY-7111: pyruvate fermentation to isobutanol (engineered)"),
#                              filter_genus = c("Rothia"),
#                              facet_formula = ". ~ Health_status") -> plot


humann2_RNA_DNA_ratio_plot <- function(df,
                                       rm_un_sp = TRUE,
                                       x_plot,
                                       y_plot,
                                       color,
                                       fill,
                                       shape = NULL,
                                       filter_feature = FALSE,
                                       filter_genus = FALSE,
                                       filter_species = FALSE,
                                       only_pos = TRUE, # dplyr::filter(DNA  > 0 , RNA > 0)
                                       facet_formula = FALSE,#". ~ Health_status"
                                       export_legend = FALSE,
                                       box_width = 0.8){ #

  require(tidyverse); require(ggConvexHull)

  if(rm_un_sp == TRUE){
    df %>%
      dplyr::filter(Species != "unclassified") -> df
  }

  if(filter_feature != FALSE){
    df %>%
      dplyr::filter(Feature %in% filter_feature) -> df
  }

  if(filter_genus != FALSE){
    df %>%
      dplyr::filter(Genus %in% filter_genus) -> df
  }

  if(filter_species != FALSE){
    df %>%
      dplyr::filter(Species %in% filter_species) -> df
  }

  if(only_pos == TRUE){
    df %>%
      dplyr::filter(DNA  > 0 , RNA > 0) -> df
  }

  df %>%
    ggplot(aes_string(x = x_plot, y = y_plot, color = color, fill = fill, shape = shape)) +
    geom_boxplot(outlier.colour = NA, alpha = 0.2,
                 position = position_dodge(width=box_width)) +
    # geom_jitter(size=1, alpha=0.2) +
    # ggbeeswarm::geom_beeswarm(size=1, alpha=0.2) +
    # geom_violin(size=1, alpha=0.2) +
    geom_hline(yintercept = 0,
               col = "red",
               linetype = "dotted",
               size = 0.5) +
    theme_light() +
    guides(fill=guide_legend(ncol=1)) +
    theme(axis.title.x = element_blank()) +
    # ylab(paste0("Relative Abundance (Top ",n," (RNA)) \n")) + #scale_y_continuous(trans='sqrt') +
    # theme(legend.text = element_text(size= 6)) +
    # theme(axis.text.x = element_text(size = 4)) +
    theme(legend.key.size = unit(0.2,"cm")) -> p

  if(facet_formula != FALSE){
    p + facet_grid(as.formula(facet_formula), drop=TRUE,
                   scale="fixed",space="free_x") -> p
  }

  if(export_legend == TRUE){

    out <- list("legend" = p %>% ggpubr::get_legend() %>% ggpubr::as_ggplot(),
                "plot" = p  + ggpubr::rotate_x_text(90) +
                  ggpubr::rotate() + theme(legend.position = "none"))

  }else{
    out <- p  + ggpubr::rotate_x_text(90) +
      ggpubr::rotate() + theme(legend.position = "none")
  }

  return(out)
}

## ------------------------------------------------------------------------------------------------------------------------------------------------

# here::here("data/processed/humann/DNA/pathabundance_cpm_joined_tables.tsv") %>% humann_2df(type = '# Pathway') -> DNA; here::here("data/processed/humann/RNA/pathabundance_cpm_joined_tables.tsv") %>% humann_2df(type = '# Pathway') -> RNA
# humann_DNA_RNA_2phyloseq(DNA,RNA) -> physeq
# sample_names(physeq) <- str_replace(sample_names(physeq), "_DNA_cat_Abundance-CPM_DNA", "_DNA"); sample_names(physeq) <- str_replace(sample_names(physeq), "_cat_Abundance-CPM_RNA", "_RNA")
# physeq_add_metadata(physeq,
#                     here::here("data/metadata_all_DNA_RNA.xlsx") %>% readxl::read_xlsx() %>% dplyr::rename(Sample_ID2 = Sample),
#                     sample_column = "Sample_ID") %>%
#   humann2_species_contribution(meta_data_var = c("Subject", "Sample_ID2" ,"Type", "Oral_Site", "Health_status")) -> df

## ------------------------------------------------------------------------------------------------------------------------------------------------

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

humann2_RNA_DNA_ratio_plot_lapply <- function(filter_species,
                                              df,
                                              rm_un_sp = TRUE,
                                              x_plot,
                                              y_plot,
                                              color,
                                              fill,
                                              shape = NULL,
                                              filter_feature = FALSE,
                                              filter_genus = FALSE,
                                              only_pos = TRUE, # dplyr::filter(DNA  > 0 , RNA > 0)
                                              facet_formula = FALSE,
                                              export_legend = FALSE){ #  ". ~ Health_status"

  require(tidyverse); require(ggConvexHull)

  if(rm_un_sp == TRUE){
    df %>%
      dplyr::filter(Species != "unclassified") -> df
  }

  if(filter_feature != FALSE){
    df %>%
      dplyr::filter(Feature %in% filter_feature) -> df
  }

  if(filter_genus != FALSE){
    df %>%
      dplyr::filter(Genus %in% filter_genus) -> df
  }
  df %>%
    dplyr::filter(Species %in% filter_species) -> df
  if(only_pos == TRUE){
    df %>%
      dplyr::filter(DNA  > 0 , RNA > 0) -> df
  }

  df %>%
    # dplyr::mutate(Feature = fct_reorder(Feature, RNA, .desc = TRUE)) %>%
    ggplot(aes_string(x = x_plot, y = y_plot, color = color, fill = fill, shape = shape)) +
    geom_boxplot(outlier.colour = NA, alpha = 0.2,
                 position = position_dodge(width=0.7)) +
    # geom_jitter(size=1, alpha=0.2) +
    ggbeeswarm::geom_beeswarm(size=1, alpha=0.2) +
    # geom_violin(size=1, alpha=0.2) +
    geom_hline(yintercept = 0,
               col = "red",
               linetype = "dotted",
               size = 0.5) +
    theme_light() +
    guides(fill=guide_legend(ncol=1)) +
    theme(axis.title.x = element_blank()) +
    # ylab(paste0("Relative Abundance (Top ",n," (RNA)) \n")) + #scale_y_continuous(trans='sqrt') +
    # theme(legend.text = element_text(size= 6)) +
    # theme(axis.text.x = element_text(size = 4)) +
    theme(legend.key.size = unit(0.2,"cm")) -> p

  if(facet_formula != FALSE){
    p + facet_grid(as.formula(facet_formula), drop=TRUE,
                   scale="fixed",space="free_x") -> p
  }
  if(export_legend == TRUE){

    out <- list("legend" = p %>% ggpubr::get_legend() %>% ggpubr::as_ggplot(),
                "plot" = p  + ggpubr::rotate_x_text(90) +
                  ggpubr::rotate() + theme(legend.position = "none") + ggtitle(paste0(filter_species)))

  }else{
    out <- p  + ggpubr::rotate_x_text(90) +
      ggpubr::rotate() + ggtitle(paste0(filter_species))
  }
  return(out)
}


# all_long_strat %>%
#   select(Full, Gene,Genus, Species, Subject, Abundance, Type, Oral_Site, Health_status) %>%
#   rename(cpm = Abundance) %>%
#   filter(cpm > 0) %>%
#   group_by(Species, Genus, Oral_Site, Health_status, Subject, Type) %>%
#   summarize(count = n(),
#             Mean = mean(cpm, na.rm = T),
#             Sum = sum(cpm>0, na.rm = T)) %>%
#   ungroup() %>%
#   select(-count, -Sum) %>%
#   group_by(Species, Genus, Oral_Site, Health_status, Subject, Type) %>%
#   pivot_wider(names_from  = Type,
#               values_from = Mean,
#               values_fill = list(mean = 0)) %>%
#   mutate(RNA_DNA = RNA/DNA) %>%
#   ungroup() -> species_DNA_RNA
#
# all_long_strat %>%
#   select(Full, Gene,Genus, Species, Subject, Abundance, Type, Oral_Site, Health_status) %>%
#   rename(cpm = Abundance) %>%
#   filter(cpm > 0) %>%
#   group_by(Species, Genus, Oral_Site, Health_status, Subject, Type) %>%
#   summarize(count = n(),
#             Mean = mean(cpm, na.rm = T),
#             Sum = sum(cpm>0, na.rm = T)) %>%
#   ungroup() %>%
#   select(-count, -Sum) %>%
#   group_by(Species, Genus, Oral_Site, Health_status, Subject, Type) %>%
#   pivot_wider(names_from  = Type,
#               values_from = Mean,
#               values_fill = list(mean = 0)) %>%
#   mutate(RNA_DNA = RNA/DNA) %>%
#   ungroup() -> species_DNA_RNA




# ## ------------------------------------------------------------------------------------------------------------------------------------------------
#
# # contributional_diversity
# humann2 <- contributional_diversity()
#
# ## ------------------------------------------------------------------------------------------------------------------------------------------------
#
# # also correlation
# all_long_strat_DNA_RNA %>%
#   # head(1000) %>%
#   group_by(Species, Gene) %>%
#   summarize(Cor_species = cor(DNA,RNA, method = "spearman"),
#             Mn_DNA = mean(DNA),
#             Mn_RNA = mean(RNA),
#             Sm_DNA = sum(DNA),
#             Sm_RNA = sum(DNA)) %>%
#   select(Species, Gene, Cor_species) -> Species_PW_corDNA_RNA
#
# all_long_strat_DNA_RNA %>%
#   # head(1000) %>%
#   group_by(Subject, Gene) %>%
#   summarize(Cor_subject = cor(DNA,RNA, method = "spearman"),
#             Mn_DNA = mean(DNA),
#             Mn_RNA = mean(RNA),
#             Sm_DNA = sum(DNA),
#             Sm_RNA = sum(DNA)) %>%
#   select(Subject, Gene, Cor_subject) -> Subject_PW_corDNA_RNA
#
# all_long_strat_DNA_RNA %>%
#   # tail(10000) %>%
#   group_by(Gene, Oral_Site, Health_status) %>%
#   summarize(Cor = cor(DNA,RNA, method = "spearman"),
#             Mn_DNA = mean(DNA),
#             Mn_RNA = mean(RNA),
#             Sm_DNA = sum(DNA),
#             Sm_RNA = sum(DNA)) %>%
#   select(Gene, Oral_Site, Health_status) -> Gene_Site_health_corDNA_RNA
#
#
# full_join(Gene_Site_health_corDNA_RNA,
#           Subject_PW_corDNA_RNA) %>%
#   full_join(Species_PW_corDNA_RNA) -> pw_all_corr
#
# pw_all_corr %>%
#   arrange(Cor_subject)
#
# ## ------------------------------------------------------------------------------------------------------------------------------------------------
#
# # boxplot per species average per patway transcriptional activity (figure c preprint)
# species_PWY_DNA_RNA %>%
#   dplyr::filter(Species %in% c("Actinobaculum_sp_oral_taxon_183",
#                                "Actinomyces_graevenitzii",
#                                "Actinomyces_massiliensis",
#                                "Actinomyces_naeslundii",
#                                "Actinomyces_oris",
#                                "Actinomyces_sp_ICM47",
#                                "Actinomyces_sp_oral_taxon_448",
#                                "Corynebacterium_matruchotii",
#                                "Fusobacterium_nucleatum",
#                                "Gemella_sanguinis",
#                                "Neisseria_flavescens",
#                                "Parvimonas_micra",
#                                "Prevotella_histicola",
#                                "Prevotella_melaninogenica",
#                                "Prevotella_nigrescens",
#                                "Prevotella_oris",
#                                "Rothia_dentocariosa",
#                                "Rothia_mucilaginosa",
#                                "Streptococcus_infantis",
#                                "Streptococcus_oralis",
#                                "Streptococcus_parasanguinis",
#                                "Streptococcus_salivarius",
#                                "Streptococcus_sanguinis",
#                                "Tannerella_forsythia",
#                                "Veillonella_atypica",
#                                "Veillonella_dispar",
#                                "Veillonella_parvula")) %>%
#   mutate(RNA_DNA = RNA/DNA) -> tmp
#
# tmp %>%
#   ggplot(aes(x = reorder(Species, -log10(RNA_DNA), FUN=median, na.rm = TRUE), y = log10(RNA_DNA), color = Genus, fill = Genus)) +
#   geom_boxplot(outlier.colour = NA, alpha = 0.2,
#                position = position_dodge(width=0.7)) +
#   # geom_jitter(size=1, alpha=0.2) +
#   ggbeeswarm::geom_beeswarm(size=1, aes(alpha = abs(Cor_species))) +
#   # geom_violin(size=1, alpha=0.2) +
#   geom_hline(yintercept = 0,
#              col = "red",
#              linetype = "dotted",
#              size = 0.5) +
#   theme_light() +
#   guides(fill=guide_legend(ncol=1)) +
#   theme(axis.title.x = element_blank()) +
#   # theme(axis.text.x=element_text(angle=75,hjust=0,vjust=0)) +
#   # ylab(paste0("Relative Abundance (Top ",n," (RNA)) \n")) + #scale_y_continuous(trans='sqrt') +
#   theme(legend.text = element_text(size= 6)) +
#   theme(axis.text.x = element_text(size = 4)) +
#   # scale_y_log10() +
#   theme(legend.key.size = unit(0.2,"cm")) + coord_cartesian(ylim = c(-1 ,1 )) -> p
#
# p  +  ggpubr::rotate_x_text(75) + xlab(NULL) +
#   ggrepel::geom_text_repel(cex = 1.5,
#                            force = 1.5,
#                            aes(label=PW),
#                            fill = "black",
#                            color = "black",
#                            show.legend = FALSE,
#                            segment.size = 0.05, segment.alpha = 0.5,
#                            direction = "both",
#                            # nudge_x = 1,
#                            nudge_y = 0,
#                            tmp %>%
#                              separate(Gene,
#                                       sep = ": ",
#                                       into = "PW",
#                                       remove = FALSE) %>%
#                              dplyr::filter(abs(log10(RNA_DNA)) > 0.6))
#
# ## ------------------------------------------------------------------------------------------------------------------------------------------------
