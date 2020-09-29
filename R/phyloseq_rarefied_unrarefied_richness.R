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
#'phyloseq_rarefied_unrerefied_richness(GlobalPatterns, 123,  "Primer", "Primer", "SampleType") -> out
#'
#'
#'
#'
#'
#'
phyloseq_rarefied_unrarefied_richness <- function (physeq, sample_size, seed, color_data, fill_data, shape_data)
{
  require(phyloseq)
  physeq %>%
    rarefy_even_depth(rngseed = seed,
                      sample.size = sample_size) -> physeq_rarefied

  prune_samples(sample_sums(physeq)>= sample_size, physeq) -> physeq
  
  
  data.frame(a =  estimate_richness(physeq_rarefied, measures = "Observed")[, 1],
             b = estimate_richness(physeq, measures = "Observed")[, 1],
             sample_data(physeq)) %>%
    ggplot(aes(x = a, y = b)) +
    geom_point(aes_string(color = color_data, fill = fill_data, shape = shape_data)) +
    geom_smooth(method="lm", level=0.95) +
    labs(x = "\nRarefied Richness", y = "UN-Rarefied Richness\n") + 
    theme_minimal() -> p 
  
  return(p)
}