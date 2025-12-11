#' Extract Plot and Legend as Separate ggplot Objects
#'
#' This function extracts the legend from a ggplot object and prepares it as a separate 
#' ggplot object. It also modifies the original plot to remove its legend.
#'
#' @param p A ggplot object from which the legend will be extracted.
#' @return A list with two elements:
#'   \itemize{
#'     \item `p`: The original ggplot object without a legend.
#'     \item `p_leg`: The legend of the original plot as a separate ggplot object.
#'   }
#' @examples
#' library(ggplot2)
#' library(ggpubr)
#' p <- ggplot(mtcars, aes(x = wt, y = mpg, color = factor(cyl))) +
#'   geom_point() +
#'   theme_minimal()
#' result <- get_plotandlegend(p)
#' result$p # The plot without a legend
#' result$p_leg # The legend as a separate ggplot object
#' @export
get_plotandlegend <- function(p){
  # Suppress warnings and messages during legend extraction
  suppressWarnings(
    suppressMessages(
      {
        # Extract the legend from the ggplot object and convert it to a ggplot object
        p_leg <- p %>% 
          ggpubr::get_legend() %>% 
          ggpubr::as_ggplot()
      }
    )
  )
  
  # Suppress warnings and messages during plot modification
  suppressWarnings(
    suppressMessages(
      {
        # Remove the legend from the original ggplot object
        p <- p + theme(legend.position = "none")
      }
    )
  )
  
  # Return the modified plot and the extracted legend as a list
  list(p_leg = p_leg, p = p)
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

#' @title Extract Unique Values from a Phyloseq Object
#' @description Extracts the unique values of a specified variable from a phyloseq object. 
#'   Works on taxonomic table information or sample metadata.
#' @param ps A `phyloseq` object.
#' @param var Name of the variable/column to extract unique values from. This can be a taxonomic rank
#'   (e.g., "Kingdom", "Phylum") or a sample metadata variable (e.g., "Primer").
#' @author Florentin Constancias
#' @note Uses `speedyseq::psmelt()` to melt the phyloseq object into a long format tibble.
#' @note The function is general and works with any column present in the melted phyloseq object.
#' @note Returns a vector of unique values, preserving the original class (character, factor, etc.).
#' @return A vector containing the unique values of the specified variable.
#' @export
#' @examples
#' library(phyloseq)
#' data("GlobalPatterns")
#'
#' # Applied on tax_table information
#' physeq_get_unique(GlobalPatterns, "Kingdom")
#'
#' physeq_get_unique(GlobalPatterns, "Phylum") %>% length()
#'
#' # Applied on sample_metadata
#' # First check the variable names
#' sample_variables(GlobalPatterns)
#'
#' # Apply
#' physeq_get_unique(GlobalPatterns, "Primer")

physeq_get_unique <- function(ps, var){
  
  ps %>%
    speedyseq::psmelt() %>%
    distinct(get(var)) %>%  # could be anything: taxa, metadata, ...
    pull() -> out
  
  return(out)
  
}



#' @title Generate a Color Palette for a Variable
#' @description Generate a distinct color palette for a variable from a phyloseq object, dataframe, or vector. 
#'   Can use `randomcoloR` or `ggpubr` palettes. Optionally displays a pie chart of the colors.
#' @param input A `phyloseq` object, `data.frame`, or vector containing the variable of interest.
#' @param var Name of the variable/column to generate colors for (ignored if `input` is a vector).
#' @param seed Numeric seed for reproducibility (default = 123456).
#' @param pal Character string: "randomcoloR" (default) or any palette name supported by `ggpubr::get_palette`.
#' @param runTsne Logical; only used if `pal = "randomcoloR"` (default = FALSE).
#' @param altCol Logical; only used if `pal = "randomcoloR"` (default = FALSE).
#' @param print Logical; if TRUE, displays a pie chart of the colors (default = TRUE).
#' @author Florentin Constancias
#' @note Works with phyloseq objects via `physeq_get_unique()`.
#' @note For dataframes, the variable must exist as a column.
#' @note For vectors, `var` is ignored and the vector unique values are used.
#' @return A named vector of colors, with names corresponding to the unique values of the variable.
#' @export
#' @examples
#' library(phyloseq)
#' data("GlobalPatterns")
#' 
#' # Check the sample variables
#' sample_variables(GlobalPatterns)
#' 
#' # Apply to phyloseq object
#' my_pal <- generate_color_palette(GlobalPatterns, var = "SampleType", pal = "npg")
#' my_pal
#' 
#' # Apply to dataframe
#' df <- data.frame(group = c("A", "B", "C", "A"))
#' generate_color_palette(df, var = "group", pal = "npg")
#' 
#' # Apply to vector
#' generate_color_palette(c("X", "Y", "Z", "X"), pal = "npg")

generate_color_palette <- function(input,
                                   var = NULL,
                                   seed = 123456,
                                   pal = "randomcoloR",
                                   runTsne = FALSE,
                                   altCol = FALSE,
                                   print = TRUE){
  
  # Determine unique values
  unique_vals <- NULL
  if ("phyloseq" %in% class(input)) {
    if (is.null(var)) stop("For phyloseq objects, 'var' must be specified")
    unique_vals <- physeq_get_unique(input, var)
  } else if (is.data.frame(input)) {
    if (is.null(var)) stop("For data frames, 'var' must be specified")
    if (!var %in% colnames(input)) stop("Variable not found in dataframe")
    unique_vals <- unique(input[[var]])
  } else if (is.vector(input)) {
    unique_vals <- unique(input)
  } else {
    stop("Input must be a phyloseq object, dataframe, or vector")
  }
  
  # Generate color palette
  if (pal == "randomcoloR") {
    set.seed(seed)
    col <- randomcoloR::distinctColorPalette(k = length(unique_vals), altCol = altCol, runTsne = runTsne)
  } else {
    col <- ggpubr::get_palette(k = length(unique_vals), palette = pal)
  }
  
  # Optionally display pie chart
  if (print) pie(rep(1, length(col)), col = col)
  
  # Name colors
  names(col) <- unique_vals
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
    dplyr::select(c("tmp_id",tax_sel)) %>%
    mutate_if(is.character,as.factor) %>%
    column_to_rownames('tmp_id') %>%
    as.matrix() -> tax_table(ps)

  return(ps)
}
