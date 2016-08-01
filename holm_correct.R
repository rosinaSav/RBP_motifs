p_values = commandArgs(TRUE)
p_values = unlist(strsplit(p_values, split = ","))
method = p_values[length(p_values)]
p_values = p_values[-length(p_values)]
p_values = as.numeric(p_values)
adjusted_p_values = p.adjust(p_values, method = method)
sprintf("%.5f",adjusted_p_values)