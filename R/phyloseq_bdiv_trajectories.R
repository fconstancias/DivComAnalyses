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
#'data(esophagus)
#'phyloseq_compute_bdiv(esophagus, 100) -> dist
#'
#'
#'
#'
#'
#'
phyloseq_bdiv_trajectories <- function(bdist, phylo_tmp, mouse_label)
{
  require(tidyverse)
  require(speedyseq)
  
  # from distance matrix (e.g., bray, select sample belonging to same mice)
  df_all_t0 = NULL
  df_all_tprev = NULL
  
  bdist %>%
    as.matrix() -> mtest
  
  for(mice in phylo_tmp %>% 
      get_variable(var) %>%
      levels())
  {
    print(mice)
    
    tmp <- prune_samples(get_variable(phylo_tmp, get(var)) == mice,
                         phylo_tmp)
    
    
    # condition : we need to have at least 3 samples, and/or day_-1.
    if (nsamples(tmp) >= 2)
    {
      # paste0(nsamples(tmp)," samples for mice : ", mice) %>%
      #     print()
      
      tmp %>%
        sample_data() %>%
        rownames_to_column("SampleID") %>%
        as_tibble() %>%
        select(SampleID, day, get(var), treatment) %>%
        drop_na() -> tmp2
      
      # a bit ugly bot does the work
      
      left_join(
        reshape2::melt(mtest[sample_names(tmp),sample_names(tmp)], 
                       varnames = c("Sample1", "Sample2"), 
                       value.name = "distance") %>% 
          filter(distance != 0 ),
        tmp2,
        by = c("Sample1" = "SampleID")
      ) %>%
        as_tibble() %>%
        dplyr::rename(Day1 = day) %>%
        left_join(
          tmp2,
          by = c("Sample2" = "SampleID")
        )  %>%
        as_tibble() %>%
        dplyr::rename(Day2 = day) %>%
        select(starts_with("Samp"), starts_with("Day"), distance, treatment.x, mouse_label.x) %>%
        arrange(Day1, Day2) %>%
        # mutate(Days = paste0(Day1,"_", Day2)) %>%
        dplyr::rename(treatment = treatment.x, mouse_label = mouse_label.x) -> tmp3
      
      # select values of interest :
      
      tmp3 %>%
        filter(Day1 %in% c("Day_-1")) %>%
        mutate(metric = i) -> tmp4
      
      bind_rows(tmp4, 
                df_all_t0) -> df_all_t0
      
      bind_rows(
        tmp3 %>%
          filter(Day1 %in% c("Day_-1") & Day2 %in% c("Day_6")),
        tmp3 %>%
          filter(Day1 %in% c("Day_6") & Day2 %in% c("Day_12")),
        tmp3 %>%
          filter(Day1 %in% c("Day_12") & Day2 %in% c("Day_14"))
      ) %>%
        mutate(metric = i) -> tmp5
      
      bind_rows(tmp5,
                df_all_tprev) -> df_all_tprev
      
    }else{
      # paste0(nsamples(tmp)," samples for mice : ", mice, " : not enough") %>%
      #     print()
    }
  }
  
  out <- list("dist_to_prev" = df_all_tprev, 
              "dist_to_t0" = df_all_t0)
  
  return(out)
  
}
