
# #function to read from qpWave
read_qpWave <- function(log_lines, details = FALSE) {
  test_pos  <- which(stringr::str_detect(log_lines, "f4info:"))
  b_pos <- which(stringr::str_detect(log_lines, "B:"))
  a_pos <- which(stringr::str_detect(log_lines, "A:"))
  a_end <- c(test_pos[-c(1, 2)], which(stringr::str_detect(log_lines, "## end of run")))
  
  test_df <- log_lines[test_pos + 1] %>%
    stringr::str_replace_all(" *[a-z0-9]+: ", "") %>%
    stringr::str_replace_all(" +", "\t") %>%
    paste0(collapse = "\n") %>%
    readr::read_tsv(col_names = c("rank", "df", "chisq", "tail", "dfdiff",
                                  "chisqdiff", "taildiff"))
  
  if (details) {
    B_matrix <- lapply(seq_along(b_pos), function(i) {
      parse_matrix(log_lines[(b_pos[i] + 1) : (a_pos[i] - 1)])
    }) %>% stats::setNames(paste0(seq_along(.)))
    
    A_matrix <- lapply(seq_along(a_pos), function(i) {
      parse_matrix(log_lines[(a_pos[i] + 1) : (a_end[i] - 2)])
    }) %>% stats::setNames(paste0(seq_along(.)))
    
    matrices <- lapply(seq_along(B_matrix), function(rank) {
      list(A = A_matrix[[rank]], B = B_matrix[[rank]])
    })
    
    return(list(ranks = test_df, matrices = matrices))
  } else {
    return(test_df)
  }
  
}

