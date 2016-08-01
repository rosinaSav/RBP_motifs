my_chisq = function(observed, expected) {
  chi_sq = sum(((observed - expected)^2)/expected)
  df = length(observed) - 1
  chisq_table = read.csv("general/chisq_table.csv", header = TRUE, sep = "\t")
  current_line = chisq_table[df, ]
  if (chi_sq < current_line[1]) {
    p = "p > 0.95"
  }
  else {
    for (chi_pos in 1:length(current_line)) {
      if (chi_sq > current_line[chi_pos]) {
        p_string = colnames(current_line)[chi_pos]
        p = paste("p < ", substr(p_string, start = 2, stop = nchar(p_string)), sep = "")
      }
    }
  }
  result = c(chi_sq, p)
  return(result)
}

input_data = commandArgs(TRUE)
input_data = strsplit(input_data, split = "|", fixed = TRUE)
observed = as.numeric(unlist(strsplit(input_data[[1]][1], split = ",")))
expected = as.numeric(unlist(strsplit(input_data[[1]][2], split = ",")))
result = my_chisq(observed, expected)
sprintf("%s", result)
