library('meta')

# load all data
all_data <- read.csv('testsets/test_input.csv')

# get all outcomes
ocns = unique(all_data$outcome)

# the final dataframe for output
df_rsts = data.frame()

# start count time
start_time <- Sys.time()

for (i in 1:length(ocns)) {
  ocn = ocns[i]
  df <- all_data[all_data$outcome == ocn, c('study', 'Et', 'Nt', 'Ec', 'Nc')]
  
  # get the result of this study
  r <- metabin(
    Et, 
    Nt, 
    Ec, 
    Nc, 
    data = df, 
    studlab = study, 
    sm = "OR", 
    method = "Inverse", 
    method.tau = "DL", 
    prediction = FALSE, 
    hakn = FALSE
  )
  
  df_rst = data.frame(
    "Outcome" = c(paste(ocn)), 
    "sm" = c('OR'),
    "TE.fixed" = c(r$TE.fixed), 
    "lower.fixed" = c(r$lower.fixed), 
    "upper.fixed" = c(r$upper.fixed),
    "TE.random" = c(r$TE.random), 
    "lower.random" = c(r$lower.random), 
    "upper.random" = c(r$upper.random),
    "I2" = c(r$I2), 
    "tau2" = c(r$tau2), 
    "pval.Q" = c(r$pval.Q)
  )

  # merge this record
  df_rsts = rbind(df_rsts, df_rst)

}
end_time <- Sys.time()
print(paste0('Spent: ', end_time - start_time))
print(df_rsts)

write.csv(
  df_rsts, 
  file='testsets/test_result_OR.csv', 
  row.names = FALSE
)

forest(r, backtransf = FALSE)