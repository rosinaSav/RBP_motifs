input_data = commandArgs(TRUE)
input_data = strsplit(input_data, split = "_", fixed = TRUE)
alt = input_data[[1]][2]
input_data = input_data[[1]][1]
input_data = strsplit(input_data, split = "|", fixed = TRUE)
vector1 = as.numeric(unlist(strsplit(input_data[[1]][1], split = ",")))
vector2 = as.numeric(unlist(strsplit(input_data[[1]][2], split = ",")))
result = wilcox.test(vector1, vector2, alternative = alt, paired = TRUE)$p.value
sprintf("%s", result)