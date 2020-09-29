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
#'phyloseq_check_lib_size(enterotype,"SeqTech","Project", 1000, 10) -> out
#'
#'
#'
#'
#'
#'
phyloseq_check_lib_size <- function(physeq, data_color, data_facet, nreads_display, first_n)
{
  require(tidyverse)
  require(speedyseq)
  
  if ('SampleID' %in% (physeq %>% sample_data() %>% colnames()))
      {
        physeq %>%
          sample_data() %>%
          as.matrix() %>%
          as.data.frame() -> df 
  }else 
  {
    physeq %>%
      sample_data() %>%
      as.matrix() %>%
      as.data.frame() %>%
      rownames_to_column('SampleID') -> df
  }
  df$LibrarySize <- sample_sums(physeq)
  df <- df[order(df$LibrarySize),]
  df$Index <- seq(nrow(df))
  
  p <- ggplot(data=df,  #ggplot(data=subset(df,Index<400), 
              aes_string(x="Index", y="LibrarySize", color = data_color , label = "SampleID")) + 
    geom_point(size=0.5) +
    ggrepel::geom_text_repel(
      data = df %>% filter(LibrarySize < nreads_display)  #sample_type == "NC" | sample_type ==  "MOCK")
    ) + facet_wrap(~ get(data_facet) , ncol = 1)
  
  p + coord_cartesian(xlim=c(0,first_n)) +
    geom_hline(yintercept = nreads_display, color="red" , size = 0.5) -> p
  
  out <- list("plot" = p, 
              "df" = df)
  
  return(out)
  
}
