cont_table = commandArgs(TRUE)
cont_table = unlist(strsplit(cont_table, split = ","))
alt = cont_table[5]
cont_table = as.numeric(cont_table[1:4])
test = fisher.test(matrix(cont_table,  ncol=2), alternative = alt)
results = c(test$estimate, test$p.value)
sprintf("%.5f", results)
