##Code to summarize qpAdm results from 2 or more way 

#function to read from qpAdm
read_qpAdm <- function(log_lines) {
  # parse the admixture proportions and standard errors
  stats <- stringr::str_subset(log_lines, "(best coefficients|std. errors):") %>%
    stringr::str_replace("(best coefficients|std. errors): +", "") %>%
    stringr::str_replace_all(" +", " ") %>%
    stringr::str_replace_all("^ | $", "") %>%
    stringr::str_split(" ") %>%
    lapply(as.numeric) %>%
    stats::setNames(c("proportion", "stderr"))
  
  #parse the nested value
  #nest <- stringr::str_subset(log_lines, "p-value for nested model:") %>%
    #grep(pattern = "infeasible",invert = T)
  nest <- stringr::str_subset(log_lines, "p-value for nested model:") %>%
    stringr::str_replace("best pat: +", "") %>%
    stringr::str_replace("p-value for nested model: +", "") %>%
    stringr::str_replace("chi\\(nested\\): +", "") %>%
    stringr::str_replace_all(" +", " ") %>%
    stringr::str_replace_all("^ | $", "")
    #stringr::str_split(" ") %>% 
    #lapply('[',c(1,4,5)) %>%
    #stats::setNames(c("nested", "P","infeasible"))
  nest <- sapply(nest, USE.NAMES = FALSE, function(l)
    if (stringr::str_detect(l, "infeasible")) l else paste0(l, " -"))
  
  if (length(nest)==1){
    nest_df <- nest %>% 
      c(.,"\n") %>%
      paste0(collapse = "\n") %>%
      readr::read_delim(delim = " ", col_names = FALSE) %>%
      stats::setNames(c("pat","p","chi","nestp","comment"))
  }else{
    nest_df <- nest %>%
      paste0(collapse = "\n") %>%
      readr::read_delim(delim = " ", col_names = FALSE) %>%
      stats::setNames(c("pat","p","chi","nestp","comment"))
  }
  
     
  
    
  # parse the population names
  leftpops <- stringr::str_locate(log_lines, "(left|right) pops:") %>%
    .[, 1] %>%  { !is.na(.) } %>% which
  target <- log_lines[leftpops[1] + 1]
  source <- log_lines[(leftpops[1] + 2) : (leftpops[2] - 2)]
  
  # parse the SNP count
  snp_count <- stringr::str_subset(log_lines, paste0("coverage: +", target)) %>%
    stringr::str_replace(paste0("coverage: +", target, " +"), "") %>%
    as.numeric
  
#  proportions <- rbind(c(target, stats$proportion, stats$stderr, snp_count)) %>% #this line do not work well in results.R
#    tibble::as_tibble(.name_repair = "minimal") %>%
#    subset(select=-1) %>%
#    as.numeric %>%
#    stats::setNames(c(source, paste0("stderr_", source), "nsnps"))
  
#  proportions <- rbind(c(target, stats$proportion, stats$stderr, snp_count)) %>%
#    tibble::as_tibble(.name_repair = "minimal") %>%
#    stats::setNames(c("target", source, paste0("stderr_", source), "nsnps")) %>%
#    dplyr::mutate_at(dplyr::vars(-target), as.numeric)
  
  proportions <- rbind(c(target, stats$proportion, stats$stderr, snp_count)) %>%
    tibble::as_tibble(.name_repair = "minimal") %>%
    stats::setNames(c("target", source, paste0("stderr_", source), "nsnps")) 
#    dplyr::mutate_at(dplyr::vars(-target), as.numeric)  #do not mutate to numerics, just keep them as characters
  
  # parse the population combination patterns into a data.frame
  pat_start <- stringr::str_detect(log_lines, "fixed pat") %>% which
  pat_end <- stringr::str_detect(log_lines, "best pat") %>% which
  patterns <- log_lines[pat_start : (pat_end[1] - 1)] %>%
    stringr::str_replace(" fixed", "") %>%
    stringr::str_replace(" prob", "") %>%
    stringr::str_replace(" pattern", "") %>%
    stringr::str_replace_all(" +", " ") %>%
    stringr::str_replace_all("^ | $", "")
  pat_header <- c(strsplit(patterns[1], " ")[[1]], source) %>%
    stringr::str_replace("pat", "pattern")
  pat_header <- c(pat_header, "comment")
  #if (any(stringr::str_detect(patterns, "infeasible"))) {
    patterns[-1] <- sapply(patterns[-1], USE.NAMES = FALSE, function(l)
      if (stringr::str_detect(l, "infeasible")) l else paste0(l, " -"))
  #}
  pat_df <- patterns[-1] %>%
    paste0(collapse = "\n") %>%
    readr::read_delim(delim = " ", col_names = FALSE) %>%
    stats::setNames(pat_header)
  
  # parse the rank test results
  ranks <- read_qpWave(log_lines)
  
  list(proportions = proportions, ranks = ranks, subsets = pat_df,nest=nest_df)
}


